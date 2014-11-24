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

/*****************************************************************************/
/* Convenience functions for aperture functions */

typedef struct
{
  void *data;
  void *error;
  void *mask;
  converter convert;
  converter econvert;
  converter mconvert;
  int size;
  int esize;
  int msize;
  short errisarray;
  short errisstd;
} sepimage;

inline int get_image_info(void *data, void *error, void *mask,
			  int dtype, int edtype, int mdtype,
			  short inflag, sepimage *sepim)
{
  int status = 0;

  sepim->data = data;
  sepim->error = error;
  sepim->mask = mask;
  sepim->size = 0;
  sepim->esize = 0;
  sepim->msize = 0;

  /* get data converter(s) for input array(s) */
  if ((status = get_converter(dtype, &(sepim->convert), &(sepim->size))))
    return status;
  if (error)
    if ((status = get_converter(edtype, &(sepim->econvert), &(sepim->esize))))
      return status;
  if (mask)
    if ((status = get_converter(mdtype, &(sepim->mconvert), &(sepim->msize))))
      return status;

  /* get options for how to interpret error array and
   * how to treat masked pixels.
   */
  sepim->errisarray = inflag & SEP_ERROR_IS_ARRAY;
  if (!error)
    sepim->errisarray = 0;  /* in case user set flag but error is NULL */
  sepim->errisstd = !(inflag & SEP_ERROR_IS_VAR);

  return status;
}


inline void update_var(sepimage *sepim, void *errort, PIXTYPE *var)
{
  *var = sepim->econvert(errort);
  if (sepim->errisstd)
    *var *= *var;
}
  

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

/*****************************************************************************/
/* exact mode subpixel sampling */

/* Return area of a circle arc between (x1, y1) and (x2, y2) with radius r */
/* reference: http://mathworld.wolfram.com/CircularSegment.html */
static inline double area_arc(double x1, double y1, double x2, double y2,
			      double r)
{
  double a, theta;

  a = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
  theta = 2. * asin(0.5 * a / r);
  return 0.5 * r * r * (theta - sin(theta));
}

/* Area of a triangle defined by three verticies */
static inline double area_triangle(double x1, double y1, double x2, double y2,
				   double x3, double y3)
{
  return 0.5 * fabs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));
}

/* Core of circular overlap routine.
 * Assumes that xmax >= xmin >= 0.0, ymax >= ymin >= 0.0.
 * (can always modify input to conform to this).
 */
static inline double circoverlap_core(double xmin, double ymin,
				      double xmax, double ymax, double r)
{
  double a, b, x1, x2, y1, y2, r2, xmin2, ymin2, xmax2, ymax2;

  xmin2 = xmin*xmin;
  ymin2 = ymin*ymin;
  r2 = r*r;
  if (xmin2 + ymin2 > r2)
    return 0.;

  xmax2 = xmax*xmax;
  ymax2 = ymax*ymax;
  if (xmax2 + ymax2 < r2)
    return (xmax-xmin)*(ymax-ymin);

  a = xmax2 + ymin2;  /* (corner 1 distance)^2 */
  b = xmin2 + ymax2;  /* (corner 2 distance)^2 */

  if (a < r2 && b < r2)
    {
      x1 = sqrt(r2 - ymax2);
      y1 = ymax;
      x2 = xmax;
      y2 = sqrt(r2 - xmax2);
      return ((xmax-xmin)*(ymax-ymin) -
	      area_triangle(x1, y1, x2, y2, xmax, ymax) +
	      area_arc(x1, y1, x2, y2, r));
    }

  if (a < r2)
    {
      x1 = xmin;
      y1 = sqrt(r2 - xmin2);
      x2 = xmax;
      y2 = sqrt(r2 - xmax2);
      return (area_arc(x1, y1, x2, y2, r) +
	      area_triangle(x1, y1, x1, ymin, xmax, ymin) +
	      area_triangle(x1, y1, x2, ymin, x2, y2));
    }

  if (b < r2)
    {
      x1 = sqrt(r2 - ymin2);
      y1 = ymin;
      x2 = sqrt(r2 - ymax2);
      y2 = ymax;
      return (area_arc(x1, y1, x2, y2, r) +
	      area_triangle(x1, y1, xmin, y1, xmin, ymax) +
	      area_triangle(x1, y1, xmin, y2, x2, y2));
    }

  /* else */
  x1 = sqrt(r2 - ymin2);
  y1 = ymin;
  x2 = xmin;
  y2 = sqrt(r2 - xmin2);
  return (area_arc(x1, y1, x2, y2, r) +
	  area_triangle(x1, y1, x2, y2, xmin, ymin));
}



/* Area of overlap of a rectangle and a circle */
static double circoverlap(double xmin, double ymin, double xmax, double ymax,
			  double r)
{
  if (0. <= xmin)
    {
      if (0. <= ymin)
	return circoverlap_core(xmin, ymin, xmax, ymax, r);
      else if (0. >= ymax)
	return circoverlap_core(-ymax, xmin, -ymin, xmax, r);
      else
	return (circoverlap(xmin, ymin, xmax, 0., r) +
                circoverlap(xmin, 0., xmax, ymax, r));
    }
  else if (0. >= xmax)
    {
      if (0. <= ymin)
	return circoverlap_core(-xmax, ymin, -xmin, ymax, r);
      else if (0. >= ymax)
	return circoverlap_core(-xmax, -ymax, -xmin, -ymin, r);
      else:
	return (circoverlap(xmin, ymin, xmax, 0., r) +
                circoverlap(xmin, 0., xmax, ymax, r));
    }
  else
    {
      if (0. <= ymin)
	return (circoverlap(xmin, ymin, 0., ymax, r) +
                circoverlap(0., ymin, xmax, ymax, r));
      if (0. >= ymax)
	return (circoverlap(xmin, ymin, 0., ymax, r) +
                circoverlap(0., ymin, xmax, ymax, r));
      else
	return (circoverlap(xmin, ymin, 0., 0., r) +
                circoverlap(0., ymin, xmax, 0., r) +
                circoverlap(xmin, 0., 0., ymax, r) +
                circoverlap(0., 0., xmax, ymax, r));
    }
}

/*
Start of new circular overlap routine that might be faster.

double circoverlap_new(double dx, double dy, double r)
{
  double xmin, xmax, ymin, ymax, xmin2, xmax2, ymin2, ymax2, r2;

  if (dx < 0.)
    dx = -dx;
  if (dy < 0.)
    dy = -dy;
  if (dy > dx)
    {
      r2 = dy;
      dy = dx;
      dx = r2;
    }

  xmax = dx + 0.5;
  ymax = dy + 0.5;
  xmax2 = xmax*xmax;
  ymax2 = ymax*ymax;
  r2 = r*r;

  if (xmax2 + ymax2 < r2)
    return 1.;

  xmin2 = xmin*xmin;
  if (xmin2 +
}
*/

/*****************************************************************************/

int sep_apercirc(void *data, void *error, void *mask,
		 int dtype, int edtype, int mdtype, int w, int h,
		 double maskthresh, double gain, short inflag,
		 double x, double y, double r, int subpix,
		 double *sum, double *sumerr, double *area, short *flag)
{
  float dx, dy, dx1, dy2, r2, rpix2, overlap, offset, rin, rout, rin2, rout2;
  float scale, scale2, pix, varpix, tmp;
  double tv, sigtv, totarea, maskarea;
  int ix, iy, xmin, xmax, ymin, ymax, sx, sy, status;
  long pos;
  BYTE *datat, *errort, *maskt;
  sepimage sepim;

  /* initializations */
  tv = sigtv = 0.0;
  overlap = totarea = maskarea = 0.0;
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

  /* get image info */
  if ((status = get_image_info(data, error, mask, dtype, edtype, mdtype,
			       inflag, &sepim)))
    return status;

  /* set varpix once if error is not an array */
  if (error && !sepim.errisarray)
    update_var(&sepim, error, &varpix);

  /* get extent of box */
  boxextent(x, y, r, r, w, h, &xmin, &xmax, &ymin, &ymax, flag);

  /* loop over rows in the box */
  for (iy=ymin; iy<ymax; iy++)
    {
      /* set pointers to the start of this row */
      pos = (iy%h) * w + xmin;
      datat = data + pos*sepim.size;
      if (sepim.errisarray)
	errort = error + pos*sepim.esize;
      if (mask)
	maskt = mask + pos*sepim.msize;

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
	      
	      pix = sepim.convert(datat);

	      if (sepim.errisarray)
		update_var(&sepim, errort, &varpix);

	      if (mask && (sepim.mconvert(maskt) > maskthresh))
		{ 
		  *flag |= SEP_APER_HASMASKED;
		  maskarea += overlap;
		}
	      else
		{
		  tv += pix*overlap;
		  sigtv += varpix*overlap;
		}

	      totarea += overlap;

	    } /* closes "if pixel within rout" */
	  
	  /* increment pointers by one element */
	  datat += sepim.size;
	  if (sepim.errisarray)
	    errort += sepim.esize;
	  maskt += sepim.msize;
	}
    }

  /* correct for masked values */
  if (mask)
    {
      if (inflag & SEP_MASK_IGNORE)
	totarea -= maskarea;
      else
	{
	  tv *= (tmp = totarea/(totarea-maskarea));
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
