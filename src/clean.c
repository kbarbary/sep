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

/* preallocate this so we don't have to realloacate it in
   every call to clean() */
static LONG cleanvictim[CLEAN_STACKSIZE];

/********************************** clean ***********************************/
/*
Add `objin` to `objlist` shallowly, and with cleaning. "With cleaning" means
that it may be merged into an existing object.
*/
int addobjcleanly(objstruct *objin, objliststruct *objlist,
		  double clean_param)
{
  objstruct     *obj, *objin;
  int	        i,j,k, status;
  double        amp,ampin,alpha,alphain, unitarea,unitareain,beta,val;
  float	       	dx,dy,rlim;

  status = RETURN_OK;

  objin = objl1->obj + objnb;

  beta = clean_param;
  unitareain = PI*objin->a*objin->b;
  ampin = objin->fdflux/(2*unitareain*objin->abcor);
  alphain = (pow(ampin/objin->dthresh, 1.0/beta)-1)*unitareain/objin->fdnpix;
  j=0;
  obj = objl2->obj;
  for (i=0; i<objl2->nobj; i++, obj++)
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
		{
		  /* ensure that cleanvictim doesn't overflow */
		  if (j >= CLEAN_STACKSIZE)
		    {
		      status = CLEAN_OVERFLOW_ERROR;
		      goto exit;
		    }

		  /*------- the newcomer puts that object in its menu! */
		  cleanvictim[j++] = i;
		}
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
		  ((float)(val<1e10?amp*pow(val,-beta):0.0) > objin->mthresh))
		{
		  /*------- the newcomer is eaten!! */
		  mergeobjshallow(objin, obj);
		  return status;
		}
	    }
	}
    }

  /* the newcomer eats the content of the menu */
  for (i=j; i--;)
    {
      k = cleanvictim[i];
      obj = objl2->obj + k;
      mergeobjshallow(obj, objin);
      status = rmobjshallow(k, objl2);
      if (status != RETURN_OK)
	goto exit;
    }

  /* and then gets added to the list! */
  status = addobjshallow(objin, objl2);

 exit:
  return status;
}
