/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
* This file is part of SEP and distributed under an MIT license.
* 
* Copyright 2014 SEP developers
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <math.h>
#include "sep.h"
#include "sepcore.h"


void sep_setellipse_ucc(unsigned char *arr, int w, int h,
                        float x, float y, float cxx, float cyy, float cxy,
		        float r, unsigned char val)
{
  unsigned char *arrt;
  int xmin, xmax, ymin, ymax, xi, yi;
  float dxlim, dylim, r2, dx, dy, dy2;

  /* get extent of ellipse in x and y */
  dxlim = cxx - cxy*cxy/(4.0*cyy);
  dxlim = dxlim>0.0 ? r/sqrt(dxlim) : 0.0;
  dylim = cyy - cxy*cxy/(4.0*cxx);
  dylim = dylim > 0.0 ? r/sqrt(dylim) : 0.0;
  r2 = r*r;

  xmin = (int)(x - dxlim + 0.5);
  xmax = (int)(x + dxlim + 1.499999);
  ymin = (int)(y - dylim + 0.5);
  ymax = (int)(y + dylim + 1.499999);
  if (xmin < 0)
    xmin = 0;
  if (xmax > w)
    xmax = w;
  if (ymin < 0)
    ymin = 0;
  if (ymax > h)
    ymax = h;

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
