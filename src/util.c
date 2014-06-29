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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "sep.h"
#include "sepcore.h"

#define GETDETAIL 1
#define PUTDETAIL 2
#define DETAILSIZE 512

char *sep_version_string = PACKAGE_VERSION;

/* data type pointer conversion */
float convertd(void *ptr)
{
  return *(double *)ptr;
}

float convertf(void *ptr)
{
  return *(float *)ptr;
}

int sizeof_dtype(int dtype, int *size)
{
  int status = 0;
  if (dtype == TFLOAT)
    {
      *size = sizeof(float);
      return RETURN_OK;
    }
  else if (dtype == TDOUBLE)
    {
      *size = sizeof(double);
      return RETURN_OK;
    }

  *size = 0;
  return ILLEGAL_DTYPE;
}

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
void sep_get_errmsg(int status, char *errtext)
{
  errtext[0] = '\0';
  switch (status)
    {
    case RETURN_OK:
      strcpy(errtext, "OK - no error");
      break;
    case SEP_INTERNAL_ERROR:
      strcpy(errtext, "SEP internal error");
      break;
    case MEMORY_ALLOC_ERROR:
      strcpy(errtext, "memory allocation");
      break;
    default:
       strcpy(errtext, "unknown error status");
       break;
    }
}

void errdetail(int action, char *errtext)
/* Read or write message to a static buffer.

   action == GETDETAIL ==> read
   action == PUTDETAIL ==> write
*/
{
  static char buff[DETAILSIZE];

  if (action == PUTDETAIL)
    strcpy(buff, errtext);
  else if (action == GETDETAIL)
    strcpy(errtext, buff);

  return;
}

void sep_get_errdetail(char *errtext)
{
  errdetail(GETDETAIL, errtext);
  return;
}

void put_errdetail(char *errtext)
{
  errdetail(PUTDETAIL, errtext);
  return;
}
