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

/* matched_filter_*() functions and get_matched_filter()
 *
 * There is a separate matched filter function for each image datatype.  In
 * sep_extract(), we select the correct function for the given array
 * type at runtime, using get_matched_filter(). The matched filter functions
 * are identical except for name and type information.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "sep.h"
#include "sepcore.h"
#include "extract.h"

#define MATCHED_FILTER_FN matched_filter_flt_flt
#define MATCHED_FILTER_IMAGE_TYPE float
#define MATCHED_FILTER_NOISE_TYPE float
#include "matchedfilterbody.c.inc"
#undef MATCHED_FILTER_FN
#undef MATCHED_FILTER_IMAGE_TYPE
#undef MATCHED_FILTER_NOISE_TYPE

#define MATCHED_FILTER_FN matched_filter_flt_dbl
#define MATCHED_FILTER_IMAGE_TYPE float
#define MATCHED_FILTER_NOISE_TYPE double
#include "matchedfilterbody.c.inc"
#undef MATCHED_FILTER_FN
#undef MATCHED_FILTER_IMAGE_TYPE
#undef MATCHED_FILTER_NOISE_TYPE

#define MATCHED_FILTER_FN matched_filter_flt_int
#define MATCHED_FILTER_IMAGE_TYPE float
#define MATCHED_FILTER_NOISE_TYPE int
#include "matchedfilterbody.c.inc"
#undef MATCHED_FILTER_FN
#undef MATCHED_FILTER_IMAGE_TYPE
#undef MATCHED_FILTER_NOISE_TYPE

#define MATCHED_FILTER_FN matched_filter_dbl_flt
#define MATCHED_FILTER_IMAGE_TYPE double
#define MATCHED_FILTER_NOISE_TYPE float
#include "matchedfilterbody.c.inc"
#undef MATCHED_FILTER_FN
#undef MATCHED_FILTER_IMAGE_TYPE
#undef MATCHED_FILTER_NOISE_TYPE

#define MATCHED_FILTER_FN matched_filter_dbl_dbl
#define MATCHED_FILTER_IMAGE_TYPE double
#define MATCHED_FILTER_NOISE_TYPE double
#include "matchedfilterbody.c.inc"
#undef MATCHED_FILTER_FN
#undef MATCHED_FILTER_IMAGE_TYPE
#undef MATCHED_FILTER_NOISE_TYPE

#define MATCHED_FILTER_FN matched_filter_dbl_int
#define MATCHED_FILTER_IMAGE_TYPE double
#define MATCHED_FILTER_NOISE_TYPE int
#include "matchedfilterbody.c.inc"
#undef MATCHED_FILTER_FN
#undef MATCHED_FILTER_IMAGE_TYPE
#undef MATCHED_FILTER_NOISE_TYPE

#define MATCHED_FILTER_FN matched_filter_int_flt
#define MATCHED_FILTER_IMAGE_TYPE int
#define MATCHED_FILTER_NOISE_TYPE float
#include "matchedfilterbody.c.inc"
#undef MATCHED_FILTER_FN
#undef MATCHED_FILTER_IMAGE_TYPE
#undef MATCHED_FILTER_NOISE_TYPE

#define MATCHED_FILTER_FN matched_filter_int_dbl
#define MATCHED_FILTER_IMAGE_TYPE int
#define MATCHED_FILTER_NOISE_TYPE double
#include "matchedfilterbody.c.inc"
#undef MATCHED_FILTER_FN
#undef MATCHED_FILTER_IMAGE_TYPE
#undef MATCHED_FILTER_NOISE_TYPE

#define MATCHED_FILTER_FN matched_filter_int_int
#define MATCHED_FILTER_IMAGE_TYPE int
#define MATCHED_FILTER_NOISE_TYPE int
#include "matchedfilterbody.c.inc"
#undef MATCHED_FILTER_FN
#undef MATCHED_FILTER_IMAGE_TYPE
#undef MATCHED_FILTER_NOISE_TYPE

/* return the correct converter depending on the datatype code */
int get_matched_filter(int dtype, int ndtype, matched_filter *f)
{
  int status = RETURN_OK;
  char errtext[80];

  if (dtype == SEP_TFLOAT && ndtype == SEP_TFLOAT)
    {
      *f = matched_filter_flt_flt;
    }
  else if (dtype == SEP_TFLOAT && ndtype == SEP_TDOUBLE)
    {
      *f = matched_filter_flt_dbl;
    }
  else if (dtype == SEP_TFLOAT && ndtype == SEP_TINT)
    {
      *f = matched_filter_flt_int;
    }
  else if (dtype == SEP_TDOUBLE && ndtype == SEP_TFLOAT)
    {
      *f = matched_filter_dbl_flt;
    }
  else if (dtype == SEP_TDOUBLE && ndtype == SEP_TDOUBLE)
    {
      *f = matched_filter_dbl_dbl;
    }
  else if (dtype == SEP_TDOUBLE && ndtype == SEP_TINT)
    {
      *f = matched_filter_dbl_int;
    }
  else if (dtype == SEP_TINT && ndtype == SEP_TFLOAT)
    {
      *f = matched_filter_int_flt;
    }
  else if (dtype == SEP_TINT && ndtype == SEP_TDOUBLE)
    {
      *f = matched_filter_int_dbl;
    }
  else if (dtype == SEP_TINT && ndtype == SEP_TINT)
    {
      *f = matched_filter_int_int;
    }
  else
    {
      *f = NULL;
      status = ILLEGAL_DTYPE;
      sprintf(errtext, "%d %d (in get_matched_filter())", dtype, ndtype);
      put_errdetail(errtext);
    }
  return status;
}
