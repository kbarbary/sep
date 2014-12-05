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
#include "overlap.h"

/****************************************************************************/
/* conversions between ellipse representations */

/* return ellipse semi-major and semi-minor axes and position angle, given
   representation: cxx*x^2 + cyy*y^2 + cxy*x*y = 1
   derived from http://mathworld.wolfram.com/Ellipse.html

   Input requirements:
   cxx*cyy - cxy*cxy/4. > 0.
   cxx + cyy > 0.
*/
int sep_ellipse_axes(double cxx, double cyy, double cxy,
                     double *a, double *b, double *theta)
{
  double p, q, t;

  p = cxx + cyy;
  q = cxx - cyy;
  t = sqrt(q*q + cxy*cxy);

  /* Ensure that parameters actually describe an ellipse. */
  if ((cxx*cyy - cxy*cxy/4. <= 0.) || (p <= 0.))
    return NON_ELLIPSE_PARAMS;

  *a = sqrt(2. / (p - t));
  *b = sqrt(2. / (p + t));

  /* theta = 0 if cxy == 0, else (1/2) acot(q/cxy) */
  *theta = (cxy == 0.) ? 0. : (q == 0. ? 0. : atan(cxy/q))/2.;
  if (cxx>cyy)
    *theta += PI/2.;
  if (*theta > PI/2.)
    *theta -= PI;

  return RETURN_OK;
}

void sep_ellipse_coeffs(double a, double b, double theta,
                        double *cxx, double *cyy, double *cxy)
{
  double costheta, sintheta;

  costheta = cos(theta);
  sintheta = sin(theta);

  *cxx = costheta*costheta/(a*a) + sintheta*sintheta/(b*b);
  *cyy = sintheta*sintheta/(a*a) + costheta*costheta/(b*b);
  *cxy = 2.*costheta*sintheta * (1./(a*a) - 1./(b*b));
}

/*****************************************************************************/
/* Helper functions for aperture functions */

/* determine the extent of the box that just contains the circle with
 * parameters x, y, r. xmin is inclusive and xmax is exclusive.
 * Ensures that box is within image bound and sets a flag if it is not.
 */
static void boxextent(double x, double y, double rx, double ry, int w, int h,
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


static void boxextent_ellipse(double x, double y,
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

/* determine oversampled annulus for a circle */
static void oversamp_ann_circle(double r, double *r_in2, double *r_out2)
{
   *r_in2 = r - 0.7072;
   *r_in2 = (*r_in2 > 0.0) ? (*r_in2)*(*r_in2) : 0.0;
   *r_out2 = r + 0.7072;
   *r_out2 = (*r_out2) * (*r_out2);
}

/* determine oversampled "annulus" for an ellipse */
static void oversamp_ann_ellipse(double r, double b, double *r_in2,
				 double *r_out2)
{
   *r_in2 = r - 0.7072/b;
   *r_in2 = (*r_in2 > 0.0) ? (*r_in2)*(*r_in2) : 0.0;
   *r_out2 = r + 0.7072/b;
   *r_out2 = (*r_out2) * (*r_out2);
}

/*****************************************************************************/

#define APER_NAME sep_sum_circle
#define APER_ARGS double r
#define APER_DECL double r2, r_in2, r_out2
#define APER_INIT				\
  r2 = r*r;					\
  oversamp_ann_circle(r, &r_in2, &r_out2)
#define APER_BOXEXTENT boxextent(x, y, r, r, w, h,			\
                                 &xmin, &xmax, &ymin, &ymax, flag)
#define APER_EXACT circoverlap(dx-0.5, dy-0.5, dx+0.5, dy+0.5, r)
#define APER_RPIX2 dx*dx + dy*dy
#define APER_RPIX2_SUBPIX dx1*dx1 + dy2
#define APER_COMPARE1 rpix2 < r_out2
#define APER_COMPARE2 rpix2 > r_in2
#define APER_COMPARE3 rpix2 < r2
#include "aperbody.h"
#undef APER_NAME
#undef APER_ARGS
#undef APER_DECL
#undef APER_INIT
#undef APER_BOXEXTENT
#undef APER_EXACT
#undef APER_RPIX2
#undef APER_RPIX2_SUBPIX
#undef APER_COMPARE1
#undef APER_COMPARE2
#undef APER_COMPARE3

/* TODO: require that a>b, a and b positive, theta in range -PI/2, PI/2 */
#define APER_NAME sep_sum_ellipse
#define APER_ARGS double a, double b, double theta, double r
#define APER_DECL double cxx, cyy, cxy, r2, r_in2, r_out2
#define APER_INIT						\
  r2 = r*r;							\
  oversamp_ann_ellipse(r, b, &r_in2, &r_out2);			\
  sep_ellipse_coeffs(a, b, theta, &cxx, &cyy, &cxy);		\
  a *= r;	       						\
  b *= r
#define APER_BOXEXTENT boxextent_ellipse(x, y, cxx, cyy, cxy, r, w, h, \
		                         &xmin, &xmax, &ymin, &ymax, flag)
#define APER_EXACT ellipoverlap(dx-0.5, dy-0.5, dx+0.5, dy+0.5, a, b, theta)
#define APER_RPIX2 cxx*dx*dx + cyy*dy*dy + cxy*dx*dy
#define APER_RPIX2_SUBPIX cxx*dx1*dx1 + cyy*dy2 + cxy*dx1*dy
#define APER_COMPARE1 rpix2 < r_out2
#define APER_COMPARE2 rpix2 > r_in2
#define APER_COMPARE3 rpix2 < r2
#include "aperbody.h"
#undef APER_NAME
#undef APER_ARGS
#undef APER_DECL
#undef APER_INIT
#undef APER_BOXEXTENT
#undef APER_EXACT
#undef APER_RPIX2
#undef APER_RPIX2_SUBPIX
#undef APER_COMPARE1
#undef APER_COMPARE2
#undef APER_COMPARE3

#define APER_NAME sep_sum_circann
#define APER_ARGS double rin, double rout
#define APER_DECL double rin2, rin_in2, rin_out2, rout2, rout_in2, rout_out2
#define APER_INIT					\
  rin2 = rin*rin;					\
  oversamp_ann_circle(rin, &rin_in2, &rin_out2);	\
  rout2 = rout*rout;					\
  oversamp_ann_circle(rout, &rout_in2, &rout_out2)
#define APER_BOXEXTENT boxextent(x, y, rout, rout, w, h, \
				 &xmin, &xmax, &ymin, &ymax, flag)
#define APER_EXACT (circoverlap(dx-0.5, dy-0.5, dx+0.5, dy+0.5, rout) - \
		    circoverlap(dx-0.5, dy-0.5, dx+0.5, dy+0.5, rin))
#define APER_RPIX2 dx*dx + dy*dy
#define APER_RPIX2_SUBPIX dx1*dx1 + dy2
#define APER_COMPARE1 (rpix2 < rout_out2) && (rpix2 > rin_in2)
#define APER_COMPARE2 (rpix2 > rout_in2) || (rpix2 < rin_out2)
#define APER_COMPARE3 (rpix2 < rout2) && (rpix2 > rin2)
#include "aperbody.h"
#undef APER_NAME
#undef APER_ARGS
#undef APER_DECL
#undef APER_INIT
#undef APER_BOXEXTENT
#undef APER_EXACT
#undef APER_RPIX2
#undef APER_RPIX2_SUBPIX
#undef APER_COMPARE1
#undef APER_COMPARE2
#undef APER_COMPARE3

#define APER_NAME sep_sum_ellipann
#define APER_ARGS double a, double b, double theta, double rin, double rout
#define APER_DECL double cxx, cyy, cxy;				\
  double rin2, rin_in2, rin_out2, rout2, rout_in2, rout_out2
#define APER_INIT						\
  rin2 = rin*rin;						\
  oversamp_ann_ellipse(rin, b, &rin_in2, &rin_out2);		\
  rout2 = rout*rout;						\
  oversamp_ann_ellipse(rout, b, &rout_in2, &rout_out2);		\
  sep_ellipse_coeffs(a, b, theta, &cxx, &cyy, &cxy)
#define APER_BOXEXTENT boxextent_ellipse(x, y, cxx, cyy, cxy, rout, w, h, \
		                         &xmin, &xmax, &ymin, &ymax, flag)
#define APER_EXACT							\
  (ellipoverlap(dx-0.5, dy-0.5, dx+0.5, dy+0.5, a*rout, b*rout, theta) - \
   ellipoverlap(dx-0.5, dy-0.5, dx+0.5, dy+0.5, a*rin, b*rin, theta))
#define APER_RPIX2 cxx*dx*dx + cyy*dy*dy + cxy*dx*dy
#define APER_RPIX2_SUBPIX cxx*dx1*dx1 + cyy*dy2 + cxy*dx1*dy
#define APER_COMPARE1 (rpix2 < rout_out2) && (rpix2 > rin_in2)
#define APER_COMPARE2 (rpix2 > rout_in2) || (rpix2 < rin_out2)
#define APER_COMPARE3 (rpix2 < rout2) && (rpix2 > rin2)
#include "aperbody.h"
#undef APER_NAME
#undef APER_ARGS
#undef APER_DECL
#undef APER_INIT
#undef APER_BOXEXTENT
#undef APER_EXACT
#undef APER_RPIX2
#undef APER_RPIX2_SUBPIX
#undef APER_COMPARE1
#undef APER_COMPARE2
#undef APER_COMPARE3

/*****************************************************************************/
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
