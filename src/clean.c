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

* SEP is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.

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
#include <stdlib.h>
#include <string.h>
#include "sep.h"
#include "extract.h"

void mergeobject(objstruct *, objstruct *);

/********************************** clean ***********************************
PURPOSE Remove object from frame -buffer and put it in the "CLEANlist".
INPUT   Object number,
        Object list (source).
OUTPUT  0 if the object was CLEANed, 1 otherwise.
 ***/
int clean(int objnb, objliststruct *objlistin, objliststruct *cleanobjlist,
	  LONG *cleanvictim, double clean_param, int *status)
{
  objstruct     *objin, *obj;
  int	        i,j,k;
  double        amp,ampin,alpha,alphain, unitarea,unitareain,beta,val;
  float	       	dx,dy,rlim;

  objin = objlistin->obj+objnb;
  beta = clean_param;
  unitareain = PI*objin->a*objin->b;
  ampin = objin->fdflux/(2*unitareain*objin->abcor);
  alphain = (pow(ampin/objin->dthresh, 1.0/beta)-1)*unitareain/objin->fdnpix;
  j=0;
  obj = cleanobjlist->obj;
  for (i=0; i<cleanobjlist->nobj; i++, obj++)
    {
      dx = objin->mx - obj->mx;
      dy = objin->my - obj->my;
      rlim = objin->a+obj->a;
      rlim *= rlim;
      if (dx*dx+dy*dy<rlim*CLEAN_ZONE*CLEAN_ZONE)
	{
	  if (obj->fdflux < objin->fdflux)
	    {
	      val = 1 + alphain*(objin->cxx*dx*dx + objin->cyy*dy*dy +
				 objin->cxy*dx*dy);
	      if (val>1.0 &&
		  ((float)(val<1e10?ampin*pow(val,-beta) : 0.0) > obj->mthresh))
		/*------- the newcomer puts that object in its menu! */
		cleanvictim[j++] = i;
	    }
	  else
	    {
	      unitarea = PI*obj->a*obj->b;
	      amp = obj->fdflux/(2*unitarea*obj->abcor);
	      alpha = (pow(amp/obj->dthresh, 1.0/beta) - 1) *
		unitarea/obj->fdnpix;
	      val = 1 + alpha*(obj->cxx*dx*dx + obj->cyy*dy*dy +
			       obj->cxy*dx*dy);
	      if (val>1.0 &&
		  ((float)(val<1e10?amp*pow(val,-beta) : 0.0) > objin->mthresh))
		{
		  /*------- the newcomer is eaten!! */
		  mergeobject(objin, obj);
		  return 0;
		}
	    }
	}
    }

  /* the newcomer eats the content of the menu */
  for (i=j; i--;)
    {
      k = cleanvictim[i];
      obj = cleanobjlist->obj + k;
      mergeobject(obj, objin);
      if ((*status = subcleanobj(k, cleanobjlist)) != RETURN_OK);
	return 0;
    }

  return 1;
}

/******************************* addcleanobj ********************************/
/*
Add an object to the "cleanobjlist".
*/
int addcleanobj(objstruct *objin, objliststruct *cleanobjlist)
{
  int	margin, y;
  float	hh1,hh2;

  /*Update the object list */
  if (cleanobjlist->nobj)
    {
      if (!(cleanobjlist->obj =
	    (objstruct *)realloc(cleanobjlist->obj,
				 (++cleanobjlist->nobj) * sizeof(objstruct))))
	return MEMORY_CLEAN_ERROR;
    }
  else
    {
      if (!(cleanobjlist->obj = (objstruct *)malloc((++cleanobjlist->nobj) *
						    sizeof(objstruct))))
	return MEMORY_CLEAN_ERROR;
    }

  /* Compute the max. vertical extent of the object: */
  /* First from 2nd order moments, compute y-limit of the 3-sigma ellips... */
  hh1 = objin->cyy - objin->cxy*objin->cxy/(4.0*objin->cxx);
  hh1 = hh1 > 0.0 ? 1/sqrt(3*hh1) : 0.0;

  /* ... then from the isophotal limit, which should not be TOO different... */
  hh2 = (objin->ymax-objin->ymin+1.0);
  margin = (int)((hh1>hh2?hh1:hh2)*MARGIN_SCALE+MARGIN_OFFSET+0.49999);
  objin->ycmax = objin->ymax+margin;

  /* ... and finally compare with the predefined margin */
  if ((y=(int)(objin->my+0.49999)+CLEAN_MARGIN)>objin->ycmax)
    objin->ycmax = y;
  objin->ycmin = objin->ymin-margin;
  if ((y=(int)(objin->my+0.49999)-CLEAN_MARGIN)<objin->ycmin)
    objin->ycmin = y;

  cleanobjlist->obj[cleanobjlist->nobj-1] = *objin;

  return 0;
  }


/******************************** mergeobject *******************************/
/*
Merge twos objects from "objlist".
*/
/* TODO: are all these fields actually initialized
   (since we're skipping some analysis)?? */

void mergeobject(objstruct *objslave, objstruct *objmaster)
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

/******************************* subcleanobj ********************************/
/* 
remove an object from a "cleanobjlist".
*/
int subcleanobj(int objnb, objliststruct *cleanobjlist)
{
  /* Update the object list */
  if (objnb>=cleanobjlist->nobj)
    return NO_CLEAN_OBJ_ERROR;
  
  if (--cleanobjlist->nobj)
    {
      if (cleanobjlist->nobj != objnb)
	cleanobjlist->obj[objnb] = cleanobjlist->obj[cleanobjlist->nobj];
      if (!(cleanobjlist->obj = (objstruct *)realloc(cleanobjlist->obj,
						     cleanobjlist->nobj * sizeof(objstruct))))
	return MEMORY_CLEAN_ERROR;
    }
  else
    free(cleanobjlist->obj);
  
  return RETURN_OK;
}
