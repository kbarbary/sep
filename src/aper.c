/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
* This file is part of SEP
*
* Copyright 1993-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
* Copyright 2014 SEP developers
*
* SEP is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* SEP is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with SEP.  If not, see <http://www.gnu.org/licenses/>.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/* Note: was photom.c in SExtractor. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sep.h"
#include "sepcore.h"

/* determine the extent of the box that just contains the circle with
 * parameters x, y, r. xmin is inclusive and xmax is exclusive.
 * Ensures that box is within image bound and sets a flag if it is not.
 */


inline void boxextent(double x, double y, double rx, double ry, int w, int h,
		      int *xmin, int *xmax, int *ymin, int *ymax,
		      short *flag)
{
  *xmin = (int)(x - rx + 0.5);
  *xmax = (int)(x + rx + 1.499999);
  *ymin = (int)(y - ry + 0.5);
  *ymax = (int)(y + ry + 1.499999);
  if (*xmin < 0)
    {
      *xmin = 0;
      *flag |= SEP_APER_TRUNC;
    }
  if (*xmax > w)
    {
      *xmax = w;
      *flag |= SEP_APER_TRUNC;
    }
  if (*ymin < 0)
    {
      *ymin = 0;
      *flag |= SEP_APER_TRUNC;
    }
  if (*ymax > h)
    {
      *ymax = h;
      *flag |= SEP_APER_TRUNC;
    }
}


inline void boxextent_ellipse(double x, double y,
			      double cxx, double cyy, double cxy, double r,
			      int w, int h,
			      int *xmin, int *xmax, int *ymin, int *ymax,
			      short *flag)
{
  double dxlim, dylim;

  dxlim = cxx - cxy*cxy/(4.0*cyy);
  dxlim = dxlim>0.0 ? r/sqrt(dxlim) : 0.0;
  dylim = cyy - cxy*cxy/(4.0*cxx);
  dylim = dylim > 0.0 ? r/sqrt(dylim) : 0.0;
  boxextent(x, y, dxlim, dylim, w, h, xmin, xmax, ymin, ymax, flag);
}


int sep_apercirc(void *data, void *error, void *mask,
		 int dtype, int edtype, int mdtype, int w, int h,
		 double maskthresh, double gain, short inflag,
		 double x, double y, double r, int subpix,
		 double *sum, double *sumerr, double *area, short *flag)
{
  float dx, dy, dx1, dy2, r2, rpix2, overlap, offset, rin, rout, rin2, rout2;
  float scale, scale2, pix, varpix, tmp;
  double tv, sigtv, okarea, totarea;
  int ix, iy, xmin, xmax, ymin, ymax, sx, sy, status, size, esize, msize;
  long pos;
  short errisarray, errisstd, maskignore;
  BYTE *datat, *errort, *maskt;
  converter convert, econvert, mconvert;

  /* get data converter(s) for input array(s) */
  status = get_converter(dtype, &convert, &size);
  if (status)
    return status;
  if (error)
    {
      status = get_converter(edtype, &econvert, &esize);
      if (status)
	return status;
    }
  if (mask)
    {
      status = get_converter(mdtype, &mconvert, &msize);
      if (status)
	return status;
    }

  /* initializations */
  size = esize = msize = 0;
  tv = sigtv = 0.0;
  overlap = totarea = okarea = 0.0;
  datat = maskt = NULL;
  errort = error;
  *flag = 0;
  r2 = r*r;
  varpix = 0.0;
  scale = 1.0/subpix;
  scale2 = scale*scale;
  offset = 0.5*(scale-1.0);
  rin = r - 0.7072;  /* Internal radius of oversampled annulus */
  rout = r + 0.7072;  /* external radius of oversampled annulus */
  rin2 = (rin>0.0)? rin*rin: 0.0;
  rout2 = rout*rout;

  /* get data converter(s) for input array(s) */
  if ((status = get_converter(dtype, &convert, &size)))
    return status;
  if (error && (status = get_converter(edtype, &econvert, &esize)))
    return status;
  if (mask && (status = get_converter(mdtype, &mconvert, &msize)))
    return status;

  /* get options for how to interpret error array and
   * how to treat masked pixels.
   */
  errisarray = inflag & SEP_ERROR_IS_ARRAY;
  if (!error)
    errisarray = 0;  /* in case user set flag but error is NULL */
  errisstd = !(inflag & SEP_ERROR_IS_VAR);
  maskignore = inflag & SEP_MASK_IGNORE;

  /* If error exists and is scalar, set the pixel variance now */
  if (error && !errisarray)
    {
      varpix = econvert(errort);
      if (errisstd)
	varpix *= varpix;
    }

  /* get extent of box */
  boxextent(x, y, r, r, w, h, &xmin, &xmax, &ymin, &ymax, flag);

  /* loop over rows in the box */
  for (iy=ymin; iy<ymax; iy++)
    {
      /* set pointers to the start of this row */
      pos = (iy%h) * w + xmin;
      datat = data + pos*size;
      if (errisarray)
	errort = error + pos*esize;
      if (mask)
	maskt = mask + pos*msize;

      /* loop over pixels in this row */
      for (ix=xmin; ix<xmax; ix++)
	{
	  dx = ix - x;
	  dy = iy - y;
	  if ((rpix2=dx*dx+dy*dy) < rout2)
	    {
	      if (rpix2 > rin2)
		/* might be partially in the aperture; calculate overlap */
		{
		  dx += offset;
		  dy += offset;
		  overlap = 0.0;
		  for (sy=subpix; sy--; dy+=scale)
		    {
		      dx1 = dx;
		      dy2 = dy*dy;
		      for (sx=subpix; sx--; dx1+=scale)
			if (dx1*dx1 + dy2 < r2)
			  overlap += scale2;
		    }
		}
	      else
		/* definitely fully in aperture */
		overlap = 1.0;
	      
	      /* get pixel value */
	      pix = convert(datat);

	      /* only update varpix if error is an array */
	      if (errisarray)
		{
		  varpix = econvert(errort);
		  if (errisstd)
		    varpix *= varpix;
		}

	      /* if mask is given and mask value is above thresh, enter
	         masking procedure */
	      if (mask)
		{
		  if (mconvert(maskt) > maskthresh)
		    { 
		      *flag |= SEP_APER_HASMASKED;
		    }
		  else
		    {
		      tv += pix*overlap;
		      sigtv += varpix*overlap;
		      okarea += overlap;
		    }
		}
	      else
		{
		  tv += pix*overlap;
		  sigtv += varpix*overlap;
		}

	      totarea += overlap;

	    } /* closes "if pixel within rout" */
	  
	  /* increment pointers by one element */
	  datat += size;
	  if (errisarray)
	    errort += esize;
	  maskt += msize;
	}
    }

  /* correct for masked values */
  if (mask)
    {
      if (maskignore)
	totarea = okarea;
      else
	{
	  tv *= (tmp = totarea/okarea);
	  sigtv *= tmp;
	}
    }

  /* add poisson noise, only if gain > 0 */
  if (gain > 0.0 && tv>0.0)
    sigtv += tv/gain;

  *sum = tv;
  *sumerr = sqrt(sigtv);
  *area = totarea;

  return status;
}

int sep_apercircann(void *data, void *error, void *mask,
		    int dtype, int edtype, int mdtype, int w, int h,
		    double maskthresh, double gain, short inflag,
		    double x, double y, double rin, double rout, int subpix,
		    double *sum, double *sumerr, double *area, short *flag)
{
  float dx, dy, dx1, dy2, rpix2, overlap, offset, rin2, rout2;
  float rinin, rinout, routin, routout, rinin2, rinout2, routin2, routout2;
  float scale, scale2, pix, varpix, tmp;
  double tv, sigtv, okarea, totarea;
  int ix, iy, xmin, xmax, ymin, ymax, sx, sy, status, size, esize, msize;
  long pos;
  short errisarray, errisstd, maskignore;
  BYTE *datat, *errort, *maskt;
  converter convert, econvert, mconvert;

  /* initializations */
  size = esize = msize = 0;
  tv = sigtv = 0.0;
  overlap = totarea = okarea = 0.0;
  datat = maskt = NULL;
  errort = error;
  *flag = 0;
  varpix = 0.0;
  scale = 1.0/subpix;
  scale2 = scale*scale;
  offset = 0.5*(scale-1.0);
  rin2 = rin*rin;
  rout2 = rout*rout;

  /* inner oversampled annulus */
  rinin = rin - 0.7072;
  rinout = rin + 0.7072;
  rinin2 = (rinin>0.0)? rinin*rinin: 0.0;
  rinout2 = rinout*rinout;

  /* outer oversampled annulus */
  routin = rout - 0.7072;
  routout = rout + 0.7072;
  routin2 = (routin>0.0)? routin*routin: 0.0;
  routout2 = routout*routout;

  /* get data converter(s) for input array(s) */
  if ((status = get_converter(dtype, &convert, &size)))
    return status;
  if (error && (status = get_converter(edtype, &econvert, &esize)))
    return status;
  if (mask && (status = get_converter(mdtype, &mconvert, &msize)))
    return status;

  /* get options */
  errisarray = inflag & SEP_ERROR_IS_ARRAY;
  if (!error)
    errisarray = 0; /* in case user set flag but error is NULL */
  errisstd = !(inflag & SEP_ERROR_IS_VAR);
  maskignore = inflag & SEP_MASK_IGNORE;

  /* If error exists and is scalar, set the pixel variance now */
  if (error && !errisarray)
    {
      varpix = econvert(errort);
      if (errisstd)
	varpix *= varpix;
    }

  /* get extent of box */
  boxextent(x, y, rout, rout, w, h, &xmin, &xmax, &ymin, &ymax, flag);

  /* loop over rows in the box */
  for (iy=ymin; iy<ymax; iy++)
    {
      /* set pointers to the start of this row */
      pos = (iy%h) * w + xmin;
      datat = data + pos*size;
      if (errisarray)
	errort = error + pos*esize;
      if (mask)
	maskt = mask + pos*msize;

      /* loop over pixels in this row */
      for (ix=xmin; ix<xmax; ix++)
	{
	  dx = ix - x;
	  dy = iy - y;
	  rpix2 = dx*dx + dy*dy;

	  if ((rpix2 < routout2) && (rpix2 > rinin2))
	    {
	      /* check if pixel is in either the outer or inner oversampled
                 annulus */
	      if ((rpix2 > routin2) || (rpix2 < rinout2))
		{
		  dx += offset;
		  dy += offset;
		  overlap = 0.0;
		  for (sy=subpix; sy--; dy+=scale)
		    {
		      dx1 = dx;
		      dy2 = dy*dy;
		      for (sx=subpix; sx--; dx1+=scale)
			{
			  rpix2 = dx1*dx1 + dy2;
			  if ((rpix2 < rout2) && (rpix2 > rin2))
			    overlap += scale2;
			}
		    }
		}
	      else
		/* definitely fully in aperture */
		overlap = 1.0;
	      
	      /* get pixel value */
	      pix = convert(datat);

	      /* only update varpix if error is an array */
	      if (errisarray)
		{
		  varpix = econvert(errort);
		  if (errisstd)
		    varpix *= varpix;
		}

	      /* if mask is given and mask value is above thresh, enter
	         masking procedure */
	      if (mask)
		{
		  if (mconvert(maskt) > maskthresh)
		    { 
		      *flag |= SEP_APER_HASMASKED;
		    }
		  else
		    {
		      tv += pix*overlap;
		      sigtv += varpix*overlap;
		      okarea += overlap;
		    }
		}
	      else
		{
		  tv += pix*overlap;
		  sigtv += varpix*overlap;
		}

	      totarea += overlap;

	    } /* closes "if pixel within rout" */
	  
	  /* increment pointers by one element */
	  datat += size;
	  if (errisarray)
	    errort += esize;
	  maskt += msize;
	}
    }

  /* correct for masked values */
  if (mask)
    {
      if (maskignore)
	totarea = okarea;
      else
	{
	  tv *= (tmp = totarea/okarea);
	  sigtv *= tmp;
	}
    }

  /* add poisson noise, only if gain > 0 */
  if (gain > 0.0 && tv>0.0)
    sigtv += tv/gain;

  *sum = tv;
  *sumerr = sqrt(sigtv);
  *area = totarea;

  return status;
}

int sep_aperellip(void *data, void *error, void *mask,
		  int dtype, int edtype, int mdtype, int w, int h,
		  double maskthresh, double gain, short inflag,
		  double x, double y, double cxx, double cyy, double cxy,
		  double r, int subpix,
		  double *sum, double *sumerr, double *area, short *flag)
{
  float dx, dy, dx1, dy2, r2, rpix2, overlap, offset, rin, rout, rin2, rout2;
  float scale, scale2, pix, varpix, tmp;
  double tv, sigtv, okarea, totarea;
  int ix, iy, xmin, xmax, ymin, ymax, sx, sy, status, size, esize, msize;
  long pos;
  short errisarray, errisstd, maskignore;
  BYTE *datat, *errort, *maskt;
  converter convert, econvert, mconvert;

  /* initializations */
  size = esize = msize = 0;
  tv = sigtv = 0.0;
  overlap = totarea = okarea = 0.0;
  datat = maskt = NULL;
  errort = error;
  *flag = 0;
  r2 = r*r;
  varpix = 0.0;
  scale = 1.0/subpix;
  scale2 = scale*scale;
  offset = 0.5*(scale-1.0);
  rin = r - 0.7072;  /* Internal radius of oversampled annulus */
  rout = r + 0.7072;  /* external radius of oversampled annulus */
  rin2 = (rin>0.0)? rin*rin: 0.0;
  rout2 = rout*rout;

  /* get converter(s) for input array(s) */
  if ((status = get_converter(dtype, &convert, &size)))
    return status;
  if (error && (status = get_converter(edtype, &econvert, &esize)))
    return status;
  if (mask && (status = get_converter(mdtype, &mconvert, &msize)))
    return status;

  /* options for interpreting error array and treating masked pixels. */
  errisarray = inflag & SEP_ERROR_IS_ARRAY;
  if (!error)
    errisarray = 0;  /* in case user set flag but error is NULL */
  errisstd = !(inflag & SEP_ERROR_IS_VAR);
  maskignore = inflag & SEP_MASK_IGNORE;

  /* If error exists and is scalar, set the pixel variance now */
  if (error && !errisarray)
    {
      varpix = econvert(errort);
      if (errisstd)
	varpix *= varpix;
    }

  /* get extent of box */
  boxextent_ellipse(x, y, cxx, cyy, cxy, r, w, h,
		    &xmin, &xmax, &ymin, &ymax, flag);

  /* loop over rows in the box */
  for (iy=ymin; iy<ymax; iy++)
    {
      /* set pointers to the start of this row */
      pos = (iy%h) * w + xmin;
      datat = data + pos*size;
      if (errisarray)
	errort = error + pos*esize;
      if (mask)
	maskt = mask + pos*msize;

      /* loop over pixels in this row */
      for (ix=xmin; ix<xmax; ix++)
	{
	  dx = ix - x;
	  dy = iy - y;
	  if ((rpix2 = cxx*dx*dx + cyy*dy*dy + cxy*dx*dy) < rout2)
	    {
	      if (rpix2 > rin2)
		/* might be partially in the aperture; calculate overlap */
		{
		  dx += offset;
		  dy += offset;
		  overlap = 0.0;
		  for (sy=subpix; sy--; dy+=scale)
		    {
		      dx1 = dx;
		      dy2 = dy*dy;
		      for (sx=subpix; sx--; dx1+=scale)
			if (cxx*dx1*dx1 + cyy*dy2 + cxy*dx1*dy < r2)
			  overlap += scale2;
		    }
		}
	      else
		/* definitely fully in aperture */
		overlap = 1.0;
	      
	      /* get pixel value */
	      pix = convert(datat);

	      /* only update varpix if error is an array */
	      if (errisarray)
		{
		  varpix = econvert(errort);
		  if (errisstd)
		    varpix *= varpix;
		}

	      /* if mask is given and mask value is above thresh, enter
	         masking procedure */
	      if (mask)
		{
		  if (mconvert(maskt) > maskthresh)
		    { 
		      *flag |= SEP_APER_HASMASKED;
		    }
		  else
		    {
		      tv += pix*overlap;
		      sigtv += varpix*overlap;
		      okarea += overlap;
		    }
		}
	      else
		{
		  tv += pix*overlap;
		  sigtv += varpix*overlap;
		}

	      totarea += overlap;

	    } /* closes "if pixel within rout" */
	  
	  /* increment pointers by one element */
	  datat += size;
	  if (errisarray)
	    errort += esize;
	  maskt += msize;
	}
    }

  /* correct for masked values */
  if (mask)
    {
      if (maskignore)
	totarea = okarea;
      else
	{
	  tv *= (tmp = totarea/okarea);
	  sigtv *= tmp;
	}
    }

  /* add poisson noise, only if gain > 0 */
  if (gain > 0.0 && tv>0.0)
    sigtv += tv/gain;

  *sum = tv;
  *sumerr = sqrt(sigtv);
  *area = totarea;

  return status;
}

int sep_aperellipann(void *data, void *error, void *mask,
		     int dtype, int edtype, int mdtype, int w, int h,
		     double maskthresh, double gain, short inflag,
		     double x, double y, double cxx, double cyy, double cxy,
		     double rin, double rout, int subpix,
		     double *sum, double *sumerr, double *area, short *flag)
{
  float dx, dy, dx1, dy2, rpix2, overlap, offset, rin2, rout2;
  float rinin, rinout, routin, routout, rinin2, rinout2, routin2, routout2;
  float scale, scale2, pix, varpix, tmp;
  double tv, sigtv, okarea, totarea;
  int ix, iy, xmin, xmax, ymin, ymax, sx, sy, status, size, esize, msize;
  long pos;
  short errisarray, errisstd, maskignore;
  BYTE *datat, *errort, *maskt;
  converter convert, econvert, mconvert;

  /* initializations */
  size = esize = msize = 0;
  tv = sigtv = 0.0;
  overlap = totarea = okarea = 0.0;
  datat = maskt = NULL;
  errort = error;
  *flag = 0;
  varpix = 0.0;
  scale = 1.0/subpix;
  scale2 = scale*scale;
  offset = 0.5*(scale-1.0);
  rin2 = rin*rin;
  rout2 = rout*rout;

  /* inner oversampled annulus */
  rinin = rin - 0.7072;
  rinout = rin + 0.7072;
  rinin2 = (rinin>0.0)? rinin*rinin: 0.0;
  rinout2 = rinout*rinout;

  /* outer oversampled annulus */
  routin = rout - 0.7072;
  routout = rout + 0.7072;
  routin2 = (routin>0.0)? routin*routin: 0.0;
  routout2 = routout*routout;

  /* get data converter(s) for input array(s) */
  if ((status = get_converter(dtype, &convert, &size)))
    return status;
  if (error && (status = get_converter(edtype, &econvert, &esize)))
    return status;
  if (mask && (status = get_converter(mdtype, &mconvert, &msize)))
    return status;

  /* get options */
  errisarray = inflag & SEP_ERROR_IS_ARRAY;
  if (!error)
    errisarray = 0; /* in case user set flag but error is NULL */
  errisstd = !(inflag & SEP_ERROR_IS_VAR);
  maskignore = inflag & SEP_MASK_IGNORE;

  /* If error exists and is scalar, set the pixel variance now */
  if (error && !errisarray)
    {
      varpix = econvert(errort);
      if (errisstd)
	varpix *= varpix;
    }

  /* get extent of box */
  boxextent_ellipse(x, y, cxx, cyy, cxy, rout, w, h,
		    &xmin, &xmax, &ymin, &ymax, flag);

  /* loop over rows in the box */
  for (iy=ymin; iy<ymax; iy++)
    {
      /* set pointers to the start of this row */
      pos = (iy%h) * w + xmin;
      datat = data + pos*size;
      if (errisarray)
	errort = error + pos*esize;
      if (mask)
	maskt = mask + pos*msize;

      /* loop over pixels in this row */
      for (ix=xmin; ix<xmax; ix++)
	{
	  dx = ix - x;
	  dy = iy - y;
	  rpix2 = cxx*dx*dx + cyy*dy*dy + cxy*dx*dy;

	  if ((rpix2 < routout2) && (rpix2 > rinin2))
	    {
	      /* check if pixel is in either the outer or inner oversampled
                 annulus */
	      if ((rpix2 > routin2) || (rpix2 < rinout2))
		{
		  dx += offset;
		  dy += offset;
		  overlap = 0.0;
		  for (sy=subpix; sy--; dy+=scale)
		    {
		      dx1 = dx;
		      dy2 = dy*dy;
		      for (sx=subpix; sx--; dx1+=scale)
			{
			  rpix2 = cxx*dx1*dx1 + cyy*dy2 + cxy*dx1*dy;
			  if ((rpix2 < rout2) && (rpix2 > rin2))
			    overlap += scale2;
			}
		    }
		}
	      else
		/* definitely fully in aperture */
		overlap = 1.0;
	      
	      /* get pixel value */
	      pix = convert(datat);

	      /* only update varpix if error is an array */
	      if (errisarray)
		{
		  varpix = econvert(errort);
		  if (errisstd)
		    varpix *= varpix;
		}

	      /* if mask is given and mask value is above thresh, enter
	         masking procedure */
	      if (mask)
		{
		  if (mconvert(maskt) > maskthresh)
		    { 
		      *flag |= SEP_APER_HASMASKED;
		    }
		  else
		    {
		      tv += pix*overlap;
		      sigtv += varpix*overlap;
		      okarea += overlap;
		    }
		}
	      else
		{
		  tv += pix*overlap;
		  sigtv += varpix*overlap;
		}

	      totarea += overlap;

	    } /* closes "if pixel within rout" */
	  
	  /* increment pointers by one element */
	  datat += size;
	  if (errisarray)
	    errort += esize;
	  maskt += msize;
	}
    }

  /* correct for masked values */
  if (mask)
    {
      if (maskignore)
	totarea = okarea;
      else
	{
	  tv *= (tmp = totarea/okarea);
	  sigtv *= tmp;
	}
    }

  /* add poisson noise, only if gain > 0 */
  if (gain > 0.0 && tv>0.0)
    sigtv += tv/gain;

  *sum = tv;
  *sumerr = sqrt(sigtv);
  *area = totarea;

  return status;
}

/* calculate Kron radius from pixels within an ellipse. */
int sep_kronrad(void *data, void *mask, int dtype, int mdtype, int w, int h,
		double maskthresh,
		double x, double y, double cxx, double cyy, double cxy,
		double r, double *kronrad, short *flag)
{
  float pix;
  double r1, v1, r2, area, rpix2, dx, dy;
  int ix, iy, xmin, xmax, ymin, ymax, status, size, msize;
  long pos;
  BYTE *datat, *maskt;
  converter convert, mconvert;

  r2 = r*r;
  r1 = v1 = 0.0;
  area = 0.0;
  *flag = 0;
  datat = maskt = NULL;
  size = msize = 0;

  /* get data converter(s) for input array(s) */
  if ((status = get_converter(dtype, &convert, &size)))
    return status;
  if (mask && (status = get_converter(mdtype, &mconvert, &msize)))
      return status;


  /* get extent of ellipse in x and y */
  boxextent_ellipse(x, y, cxx, cyy, cxy, r, w, h,
		    &xmin, &xmax, &ymin, &ymax, flag);

  /* loop over rows in the box */
  for (iy=ymin; iy<ymax; iy++)
    {
      /* set pointers to the start of this row */
      pos = (iy%h) * w + xmin;
      datat = data + pos*size;
      if (mask)
	maskt = mask + pos*msize;
      
      /* loop over pixels in this row */
      for (ix=xmin; ix<xmax; ix++)
	{
	  dx = ix - x;
	  dy = iy - y;
	  rpix2 = cxx*dx*dx + cyy*dy*dy + cxy*dx*dy;
	  if (rpix2 <= r2)
	    {
	      pix = convert(datat);
	      if ((pix < -BIG) || (mask && mconvert(maskt) > maskthresh))
		{
		  *flag |= SEP_APER_HASMASKED;
		}
	      else
		{
		  r1 += sqrt(rpix2)*pix;
		  v1 += pix;
		  area++;
		}
	    }

	  /* increment pointers by one element */
	  datat += size;
	  maskt += msize;
	}
    }

  if (area == 0)
    {
      *flag |= SEP_APER_ALLMASKED;
      *kronrad = 0.0;
    }
  else if (r1 <= 0.0 || v1 <= 0.0)
    {
      *flag |= SEP_APER_NONPOSITIVE;
      *kronrad = 0.0;
    }
  else
    {
      *kronrad = r1 / v1;
    }

  return RETURN_OK;
}


/* set array values within an ellipse (uc = unsigned char array) */
void sep_setellip_uc(unsigned char *arr, int w, int h,
		     double x, double y, double cxx, double cyy, double cxy,
		     double r, unsigned char val)
{
  unsigned char *arrt;
  int xmin, xmax, ymin, ymax, xi, yi;
  double r2, dx, dy, dy2;
  short flag; /* not actually used, but input to boxextent */

  flag = 0;
  r2 = r*r;

  boxextent_ellipse(x, y, cxx, cyy, cxy, r, w, h,
		    &xmin, &xmax, &ymin, &ymax, &flag);

  for (yi=ymin; yi<ymax; yi++)
    {
      arrt = arr + (yi*w + xmin);
      dy = yi - y;
      dy2 = dy*dy;
      for (xi=xmin; xi<xmax; xi++, arrt++)
	{
	  dx = xi - x;
	  if ((cxx*dx*dx + cyy*dy2 + cxy*dx*dy) <= r2)
	    *arrt = val;
	}
    }
}
