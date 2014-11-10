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

/* convolve_*() functions and get_convolver()
 *
 * There is a separate convolve function for each image datatype.  In
 * sep_extract(), we select the correct function for the given array
 * type at runtime, using get_convolver(). The convolve functions are
 * identical except for name and type information.
 *
 * NOTE: This file uses macro tricks to avoid duplicating the source
 * code for the convolve function. The body is defined in the first
 * part of this file, using macros CONVOLVE_FN (function name) and
 * CONVOLVE_TYPE (image type).  The second half of the file sets these
 * macros and #includes the first half of the file (multiple times).
 */

/* Function body definition. Upon first reading this file, this block will
 * be skipped, as _CONVOLVE_FUNCTION_BODY is not yet defined. */
#ifdef _CONVOLVE_FUNCTION_BODY

/* Convolve one line of an image with a given kernel.
 *
 * image : full input array
 * w, h : width and height of image
 * y : line to convolve in image
 * conv : convolution kernel
 * convw, convh : width and height of conv
 * buf : output convolved line (`w` elements long)
 */
void CONVOLVE_FN(void *image, int w, int h, int y,
	         float *conv, int convw, int convh, PIXTYPE *buf)
{
  int	  convw2, cx, dcx, y0, dy;
  float	  *convend;
  PIXTYPE *bufend, *buft, *buftend, cval;

  CONVOLVE_TYPE *im, *line, *imt;
  im = (CONVOLVE_TYPE *)image;

  line = NULL;        /* To avoid gcc -Wall warnings */
  convw2 = convw/2;
  bufend = buf+w;     /* pointer to end of output buffer */
  y0 = y - (convh/2); /* starting y line for convolution */

  /* reduce height of convolution kernel if it extends beyond image */
  dy = h - y0;        /* distance from starting line to top of image */
  if (convh > dy)
    convh = dy;
  convend = conv + (convw*convh);

  /* Set start position in convolution kernel (and start line in image) */
  if (y0 < 0)
    {
      conv += convw*(-y0);
      y0 = 0;
    }

  /* initialize output buffer to zero */
  memset(buf, 0, w*sizeof(PIXTYPE));

  /* loop over pixels in the mask */
  for (cx=0; conv<convend; conv++, cx++)
    {
      cval = *conv;

      /* get the x position in the mask */
      if (cx==convw)
	cx = 0;

      /* when cx goes to zero, increment the start line in the image */
      if (!cx)
	line = im + w*((y0++)%h);

      /* get start and end positions in the source and target line */
      if ((dcx = cx-convw2)>=0)
	{
	  imt = line + dcx;
	  buft = buf;
	  buftend = bufend - dcx;
	}
      else
	{
	  imt = line;
	  buft = buf - dcx;
	  buftend = bufend;
	}

      while (buft < buftend)
	*(buft++) += cval * *(imt++);
    }

  return;
}

#else /* _CONVOLVE_FUNCTION_BODY not defined */

/* Here, we just entered this file for the first time.
 * Define the macro so that when this self-same file is #include'd,
 * it will include the only the body of the convolve function and not
 * this part of the file. */
#define _CONVOLVE_FUNCTION_BODY 1

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sep.h"
#include "sepcore.h"
#include "extract.h"

#define CONVOLVE_FN convolve_flt
#define CONVOLVE_TYPE float
#include "convolve.c"
#undef CONVOLVE_FN
#undef CONVOLVE_TYPE

#define CONVOLVE_FN convolve_dbl
#define CONVOLVE_TYPE double
#include "convolve.c"
#undef CONVOLVE_FN
#undef CONVOLVE_TYPE

#define CONVOLVE_FN convolve_int
#define CONVOLVE_TYPE int
#include "convolve.c"
#undef CONVOLVE_FN
#undef CONVOLVE_TYPE

/* return the correct converter depending on the datatype code */
int get_convolver(int dtype, convolver *f)
{
  int status = RETURN_OK;
  char errtext[80];

  if (dtype == SEP_TFLOAT)
    {
      *f = convolve_flt;
    }
  else if (dtype == SEP_TDOUBLE)
    {
      *f = convolve_dbl;
    }
  else if (dtype == SEP_TINT)
    {
      *f = convolve_int;
    }
  else
    {
      *f = NULL;
      status = ILLEGAL_DTYPE;
      sprintf(errtext, "unsupported data type in get_convolver(): %d",
	      dtype);
      put_errdetail(errtext);
    }
  return status;
}

/* cleanup */
#undef _CONVOLVE_FUNCTION_BODY

#endif
