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

/*
Utilties for dealing with `objstruct`s and `objliststruct`s
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sep.h"
#include "extract.h"

/******************************* addobjshallow *******************************/
/*
Copy an object to an object list, shallowly.
Only top-level information is copied (no pixel list information).
*/
int addobjshallow(objstruct *objin, objliststruct *objlist)
{
  /* If objlist already contains objects, realloc.
     Otherwise, alloc for the first time. */
  if (objlist->nobj)
      objlist->obj = (objstruct *)realloc(objlist->obj,
					  (++objlist->nobj)*sizeof(objstruct));
  else
      objlist->obj = (objstruct *)malloc((++objlist->nobj)*sizeof(objstruct));

  if (!objlist->obj)
    return MEMORY_CLEAN_ERROR;

  /* actually copy the object */
  objlist->obj[objlist->nobj-1] = *objin;

  return RETURN_OK;
}

/******************************* rmobjshallow ********************************/
/* 
remove an object from an objlist
*/
int rmobjshallow(int objnb, objliststruct *objlist)
{
  if (objnb>=objlist->nobj)
    {
      sprintf(errdetail,
	      "tried to remove object index %d from objlist of length %d",
	      objnb, objlist->nobj);
      return SEP_INTERNAL_ERROR;
    }
  
  /* if there are any objects left, need to reallocate memory */
  if (--objlist->nobj)
    {
      /* If we're not removing the last object, copy the last object 
	 to the space currently occupied by the object being removed */
      if (objlist->nobj != objnb)
	objlist->obj[objnb] = objlist->obj[objlist->nobj];

      objlist->obj = (objstruct *)realloc(objlist->obj,
					  objlist->nobj*sizeof(objstruct));
      if (!objlist->obj)
	return MEMORY_CLEAN_ERROR;
    }

  /* otherwise, just free the object list */
  else
    {
      free(objlist->obj);
    }

  return RETURN_OK;
}

/******************************** mergeobject *******************************/
/*
Merge two objects.
*/

void mergeobjshallow(objstruct *objslave, objstruct *objmaster)
{
  objmaster->fdnpix += objslave->fdnpix;
  objmaster->dnpix += objslave->dnpix;
  objmaster->fdflux += objslave->fdflux;
  objmaster->dflux += objslave->dflux;
  objmaster->flux += objslave->flux;
  objmaster->fluxerr += objslave->fluxerr;

  if (objslave->fdpeak>objmaster->fdpeak)
    {
      objmaster->fdpeak = objslave->fdpeak;
      objmaster->peakx = objslave->peakx;
      objmaster->peaky = objslave->peaky;
    }
  if (objslave->dpeak>objmaster->dpeak)
    objmaster->dpeak = objslave->dpeak;
  if (objslave->peak>objmaster->peak)
    objmaster->peak = objslave->peak;

  if (objslave->xmin<objmaster->xmin)
    objmaster->xmin = objslave->xmin;
  if (objslave->xmax>objmaster->xmax)
    objmaster->xmax = objslave->xmax;

  if (objslave->ymin<objmaster->ymin)
    objmaster->ymin = objslave->ymin;
  if (objslave->ymax>objmaster->ymax)
    objmaster->ymax = objslave->ymax;

  objmaster->flag |= (objslave->flag & (~(OBJ_MERGED|OBJ_CROWDED)));

  /* removed next line because it only applies to external flag maps */
  /* mergeflags(objmaster, objslave); */
  
  return;
}



/********** addobjdeep (originally in manobjlist.c) **************************/
/*
Add object number `objnb` from list `objl1` to list `objl2`.
Unlike `addobjshallow` this also copies plist pixels to the second list.
*/

int addobjdeep(int objnb, objliststruct *objl1, objliststruct *objl2)
{
  objstruct	*objl2obj;
  pliststruct	*plist1 = objl1->plist, *plist2 = objl2->plist;
  int		fp, i, j, npx, objnb2;
  
  fp = objl2->npix;      /* 2nd list's plist size in pixels */
  j = fp*plistsize;      /* 2nd list's plist size in bytes */
  objnb2 = objl2->nobj;  /* # of objects currently in 2nd list*/

  /* Allocate space in `objl2` for the new object */
  if (objnb2)
    objl2obj = (objstruct *)realloc(objl2->obj,
				    (++objl2->nobj)*sizeof(objstruct));
  else
    objl2obj = (objstruct *)malloc((++objl2->nobj)*sizeof(objstruct));

  if (!objl2obj)
    goto earlyexit;
  objl2->obj = objl2obj;

  /* Allocate space for the new object's pixels in 2nd list's plist */
  npx = objl1->obj[objnb].fdnpix;
  if (fp)
    plist2 = (pliststruct *)realloc(plist2, (objl2->npix+=npx)*plistsize);
  else
    plist2 = (pliststruct *)malloc((objl2->npix=npx)*plistsize);

  if (!plist2)
    goto earlyexit;
  objl2->plist = plist2;
  
  /* copy the plist */
  plist2 += j;
  for(i=objl1->obj[objnb].firstpix; i!=-1; i=PLIST(plist1+i,nextpix))
    {
      memcpy(plist2, plist1+i, (size_t)plistsize);
      PLIST(plist2,nextpix) = (j+=plistsize);
      plist2 += plistsize;
    }
  PLIST(plist2-=plistsize, nextpix) = -1;
  
  /* copy the object itself */
  objl2->obj[objnb2] = objl1->obj[objnb];
  objl2->obj[objnb2].firstpix = fp*plistsize;
  objl2->obj[objnb2].lastpix = j-plistsize;

  return RETURN_OK;
  
  /* if early exit, reset 2nd list */
 earlyexit:
  objl2->nobj--;
  objl2->npix = fp;
  return RETURN_ERROR;
}
