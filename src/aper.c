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
  size = esize = msize = 0;
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
  rin = r - 0.7072; /* Internal radius of oversampled annulus: r - sqrt(2)/2 */
  rin2 = (rin>0.0)? rin*rin: 0.0;
  rout = r + 0.7072; /* external radius of oversampled annulus */
  rout2 = rout*rout;

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

  /* set extent of box to loop over */
  xmin = (int)(x - r + 0.5);
  xmax = (int)(x + r + 1.499999);
  ymin = (int)(y - r + 0.5);
  ymax = (int)(y + r + 1.499999);
  if (xmin < 0)
    {
      xmin = 0;
      *flag |= SEP_APER_TRUNC;
    }
  if (xmax > w)
    {
      xmax = w;
      *flag |= SEP_APER_TRUNC;
    }
  if (ymin < 0)
    {
      ymin = 0;
      *flag |= SEP_APER_TRUNC;
    }
  if (ymax > h)
    {
      ymax = h;
      *flag |= SEP_APER_TRUNC;
    }

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
		{
		  /* might be partially in aperture; get 'overlap' */
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


/* NOTE: This is mostly a copy of the circular aperture function above.
         In the future, we may wish to refactor these functions to reduce
         code duplication, if possible, as long as performance is not
         affected too much. */
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

  /* get data converter(s) for input array(s) */
  size = esize = msize = 0;
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
  tv = sigtv = 0.0;
  overlap = totarea = okarea = 0.0;
  datat = maskt = NULL;
  errort = error;
  *flag = 0;
  rin2 = rin*rin;
  rout2 = rout*rout;
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

  /* set extent of box to loop over */
  xmin = (int)(x - rout + 0.5);
  xmax = (int)(x + rout + 1.499999);
  ymin = (int)(y - rout + 0.5);
  ymax = (int)(y + rout + 1.499999);
  if (xmin < 0)
    {
      xmin = 0;
      *flag |= SEP_APER_TRUNC;
    }
  if (xmax > w)
    {
      xmax = w;
      *flag |= SEP_APER_TRUNC;
    }
  if (ymin < 0)
    {
      ymin = 0;
      *flag |= SEP_APER_TRUNC;
    }
  if (ymax > h)
    {
      ymax = h;
      *flag |= SEP_APER_TRUNC;
    }

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
