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
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sep.h"
#include "sepcore.h"
#include "extract.h"

#define CONVOLVE_FN convolve_flt
#define CONVOLVE_TYPE float
#include "convbody.h"
#undef CONVOLVE_FN
#undef CONVOLVE_TYPE

#define CONVOLVE_FN convolve_dbl
#define CONVOLVE_TYPE double
#include "convbody.h"
#undef CONVOLVE_FN
#undef CONVOLVE_TYPE

#define CONVOLVE_FN convolve_int
#define CONVOLVE_TYPE int
#include "convbody.h"
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
      sprintf(errtext, "%d (in get_convolver())", dtype);
      put_errdetail(errtext);
    }
  return status;
}
