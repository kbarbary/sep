/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
* This file is part of SEP
*
* Copyright 2014 Kyle Barbary
*
* SEP is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* SEP is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with SEP.  If not, see <http://www.gnu.org/licenses/>.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
* This file part of: SExtractor
*
* Copyright:         (C) 1993-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
*
* License:           GNU General Public License
*
* SExtractor is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* SExtractor is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* You should have received a copy of the GNU General Public License
* along with SExtractor. If not, see <http://www.gnu.org/licenses/>.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sep.h"


/*i**** fqcmp **************************************************************
PROTO	int	fqcmp(const void *p1, const void *p2)
PURPOSE	Sorting function for floats in qsort().
INPUT	Pointer to first element,
	pointer to second element.
OUTPUT	1 if *p1>*p2, 0 if *p1=*p2, and -1 otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	05/10/2010
 ***/
static int fqcmp(const void *p1, const void *p2)
{
  double f1=*((float *)p1);
  double f2=*((float *)p2);
  return f1>f2? 1 : (f1<f2? -1 : 0);
}


/****** fqmedian *************************************************************
PROTO	float   fqmedian(float *ra, int n)
PURPOSE	Compute the median of an array of floats, using qsort().
INPUT	Pointer to the array,
	Number of array elements.
OUTPUT	Median of the array.
NOTES	Warning: the order of input data is modified!.
*/
float fqmedian(float *ra, int n)
{
  int dqcmp(const void *p1, const void *p2);
  
  qsort(ra, n, sizeof(float), fqcmp);
  if (n<2)
    return *ra;
  else
    return n&1? ra[n/2] : (ra[n/2-1]+ra[n/2])/2.0;
}

/*---------------------------------------------------------------------------*/
/*
  Return a short descriptive error message that corresponds to the input
  error status value.  The message may be up to 60 characters long, plus
  the terminating null character.
*/
void sep_errmsg(int status, char *errtext)
{
  errtext[0] = '\0';
  switch (status)
    {
    case RETURN_OK:
      strcpy(errtext, "OK - no error");
      break;
    case RETURN_ERROR:
      strcpy(errtext, "unspecified error");
      break;
    case MEMORY_PIXSTACK_ERROR:
      strcpy(errtext, "memory pixel stack error");
      break;
    case PIXSTACK_OVERFLOW_ERROR:
      strcpy(errtext, "pixel stack overflow error");
      break;
    case MEMORY_CLEAN_ERROR:
      strcpy(errtext, "memory clean error");
      break;
    case NO_CLEAN_OBJ_ERROR:
      strcpy(errtext, "internal error: no object to remove in subcleanobj");
      break;
    case LUTZ_REALLOC_ERROR:
      strcpy(errtext, "internal error: problem with memory realloc in lutz");
      break;
    case GATHERUP_MEMORY_ERROR:
      strcpy(errtext, "not enough memory to update pixel list in gatherup");
      break;
    case MEMORY_ALLOC_ERROR:
      strcpy(errtext, "memory allocation");
      break;
    case DEBLEND_OVERFLOW_ERROR:
      strcpy(errtext, "deblend overflow");
      break;
    case CLEAN_OVERFLOW_ERROR:
      strcpy(errtext, "clean overflow");
      break;
    default:
       strcpy(errtext, "unknown error status");
       break;
    }
}
