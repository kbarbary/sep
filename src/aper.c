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

void sep_apercirc(PIXTYPE *im, PIXTYPE *var, int w, int h,
		  PIXTYPE gain, PIXTYPE varthresh,
		  double cx, double cy, double r,int subpix,
		  double *flux, double *fluxerr, short *flag)
{

  /*
   float		r2, raper,raper2, rintlim,rintlim2,rextlim2,
			mx,my,dx,dx1,dy,dy2,
			offsetx,offsety,scalex,scaley,scale2, ngamma, locarea;
   double		tv, sigtv, area, pix, var, backnoise2, gain;
   int			x,y, x2,y2, xmin,xmax,ymin,ymax, sx,sy, w,h,
			fymin,fymax, pflag,corrflag, gainflag;
			long			pos;
  */
  float dx, dy, dx1, dy2, r2, rpix2, locarea, offset, rin, rin2, rout2;
  float scale, scale2;
  double pix, varpix, tv, sigtv;
  int x, y, xmin, xmax, ymin, ymax, sx, sy;
  long pos;
  PIXTYPE *imt, *vart;

   /*
  if (wfield)
    wthresh = wfield->weight_thresh;
  wstrip = wstript = NULL;
  mx = obj->mx;
  my = obj->my;
  w = field->width;
  h = field->stripheight;
  fymin = field->ymin;
  fymax = field->ymax;
  ngamma = field->ngamma;
  pflag = (prefs.detect_type==PHOTO)? 1:0;
  corrflag = (prefs.mask_type==MASK_CORRECT);
  gainflag = wfield && prefs.weightgain_flag;
  var = backnoise2 = field->backsig*field->backsig;
  gain = field->gain;
   */

  tv = sigtv = locarea = 0.0;
  vart = NULL;
  *flag = 0;
  r2 = r*r;

  /* Internal & external radii of the oversampled annulus: r +/- sqrt(2)/2 */
  rin = r - 0.75;
  rin2 = (rin>0.0)? rin*rin: 0.0;
  rout2 = (r + 0.75)*(r + 0.75);

  varpix = 0.0;
  scale = 1.0/subpix;
  scale2 = scale*scale;
  offset = 0.5*(scale-1.0);

  xmin = (int)(cx - r + 0.499999);
  xmax = (int)(cx + r + 1.499999);
  ymin = (int)(cy - r + 0.499999);
  ymax = (int)(cy + r + 1.499999);

  if (xmin < 0)
    {
      xmin = 0;
      *flag |= SEP_OBJ_APERT_PB;
    }
  if (xmax > w)
    {
      xmax = w;
      *flag |= SEP_OBJ_APERT_PB;
    }
  if (ymin < 0)
    {
      ymin = 0;
      *flag |= SEP_OBJ_APERT_PB;
    }
  if (ymax > h)
    {
      ymax = h;
      *flag |= SEP_OBJ_APERT_PB;
    }

  for (y=ymin; y<ymax; y++)
    {
      imt = im + (pos = (y%h)*w + xmin);
      if (var)
	vart = var + pos;
      for (x=xmin; x<xmax; x++, imt++, vart++)
	{
	  dx = x - cx;
	  dy = y - cy;
	  if ((rpix2=dx*dx+dy*dy) < rout2)
	    {
	      if (rpix2 > rin2)
		{
		  /* might be partially in aperture; get overlap area */
		  dx += offset;
		  dy += offset;
		  locarea = 0.0;
		  for (sy=subpix; sy--; dy+=scale)
		    {
		      dx1 = dx;
		      dy2 = dy*dy;
		      for (sx=subpix; sx--; dx1+=scale)
			if (dx1*dx1 + dy2 < r2)
			  locarea += scale2;
		    }
		}
	      else
		locarea = 1.0;
	      
	      /* if the values are crazy, set to 0 */
	      if ((pix=*imt)<=-BIG || (var && (varpix=*vart)>=varthresh))
		{
		  pix = 0.0;
		  if (var)
		    varpix = 0.0;
		}
	      
	      tv += pix*locarea;
	      sigtv += varpix*locarea;
	    } /* closes "if pixel within rout" */
	}
    }

    /* add poisson noise, only if gain > 0 */
    if (gain > 0.0 && tv>0.0)
      sigtv += tv/gain;

    *flux = tv;
    *fluxerr = sqrt(sigtv);
    
    return;
}
