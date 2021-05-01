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

/* Note: was extract.c in SExtractor. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sep.h"
#include "sepcore.h"
#include "extract.h"

#define	NOBJ 256  /* starting number of obj. */

void lutzsort(infostruct *, objliststruct *);

/******************************* lutzalloc ***********************************/
/*
Allocate once for all memory space for buffers used by lutz().
*/
int lutzalloc(int width, int height, lutzbuffers *buffers)
{
  int *discant;
  int stacksize, i, status=RETURN_OK;

  memset(buffers, 0, sizeof(lutzbuffers));

  stacksize = width+1;
  buffers->xmin = buffers->ymin = 0;
  buffers->xmax = width-1;
  buffers->ymax = height-1;
  QMALLOC(buffers->info, infostruct, stacksize, status);
  QMALLOC(buffers->store, infostruct, stacksize, status);
  QMALLOC(buffers->marker, char, stacksize, status);
  QMALLOC(buffers->psstack, pixstatus, stacksize, status);
  QMALLOC(buffers->start, int, stacksize, status);
  QMALLOC(buffers->end, int, stacksize, status);
  QMALLOC(buffers->discan, int, stacksize, status);
  discant = buffers->discan;
  for (i=stacksize; i--;)
    *(discant++) = -1;

  return status;

 exit:
  lutzfree(buffers);

  return status;
}

/******************************* lutzfree ************************************/
/*
Free once for all memory space for buffers used by lutz().
*/
void lutzfree(lutzbuffers *buffers)
{
  free(buffers->discan);
  buffers->discan = NULL;
  free(buffers->info);
  buffers->info = NULL;
  free(buffers->store);
  buffers->store = NULL;
  free(buffers->marker);
  buffers->marker = NULL;
  free(buffers->psstack);
  buffers->psstack = NULL;
  free(buffers->start);
  buffers->start = NULL;
  free(buffers->end);
  buffers->end = NULL;
}

static const infostruct initinfo = {
	.firstpix = -1,
	.lastpix = -1
};

/********************************** lutz *************************************/
/*
C implementation of R.K LUTZ' algorithm for the extraction of 8-connected pi-
xels in an image
*/
int lutz(pliststruct *plistin,
	 int *objrootsubmap, int subx, int suby, int subw,
	 objstruct *objparent, objliststruct *objlist, int minarea,
	 lutzbuffers *buffers)
{
  infostruct	curpixinfo;
  objstruct		*obj;
  pliststruct		*plist,*pixel, *plistint;

  char			newmarker;
  int			cn, co, luflag, pstop, xl,xl2,yl,
                        out, deb_maxarea, stx,sty,enx,eny, step,
                        nobjm = NOBJ,
			inewsymbol, *iscan;
  short		        trunflag;
  PIXTYPE		thresh;
  pixstatus		cs, ps;

  out = RETURN_OK;

  deb_maxarea = minarea<MAXDEBAREA?minarea:MAXDEBAREA; /* 3 or less */
  plistint = plistin;
  stx = objparent->xmin;
  sty = objparent->ymin;
  enx = objparent->xmax;
  eny = objparent->ymax;
  thresh = objlist->thresh;
  cn = 0;

  iscan = objrootsubmap + (sty-suby)*subw + (stx-subx);

  /* As we only analyse a fraction of the map, a step occurs between lines */
  step = subw - (++enx-stx);
  eny++;

  /*------Allocate memory to store object data */
  free(objlist->obj);

  if (!(obj=objlist->obj=malloc(nobjm*sizeof(objstruct))))
    {
      out = MEMORY_ALLOC_ERROR;
      plist = NULL;			/* To avoid gcc -Wall warnings */
      goto exit_lutz;
    }

  /*------Allocate memory for the pixel list */
  free(objlist->plist);
  if (!(objlist->plist
	= malloc((eny-sty)*(enx-stx)*plistsize)))
    {
      out = MEMORY_ALLOC_ERROR;
      plist = NULL;			/* To avoid gcc -Wall warnings */
      goto exit_lutz;
    }

  pixel = plist = objlist->plist;

  /*----------------------------------------*/
  for (xl=stx; xl<=enx; xl++)
    buffers->marker[xl] = 0;

  objlist->nobj = 0;
  co = pstop = 0;
  curpixinfo.pixnb = 1;

  for (yl=sty; yl<=eny; yl++, iscan += step)
    {
      ps = COMPLETE;
      cs = NONOBJECT;
      trunflag = (yl==0 || yl==buffers->ymax) ? SEP_OBJ_TRUNC : 0;
      if (yl==eny)
	iscan = buffers->discan;

      for (xl=stx; xl<=enx; xl++)
	{
	  newmarker = buffers->marker[xl];
	  buffers->marker[xl] = 0;
	  if ((inewsymbol = (xl!=enx)?*(iscan++):-1) < 0)
	    luflag = 0;
	  else
	    {
	      curpixinfo.flag = trunflag;
	      plistint = plistin+inewsymbol;
	      luflag = (PLISTPIX(plistint, cdvalue) > thresh?1:0);
	    }
	  if (luflag)
	    {
	      if (xl==0 || xl==buffers->xmax)
		curpixinfo.flag |= SEP_OBJ_TRUNC;
	      memcpy(pixel, plistint, (size_t)plistsize);
	      PLIST(pixel, nextpix) = -1;
	      curpixinfo.lastpix = curpixinfo.firstpix = cn;
	      cn += plistsize;
	      pixel += plistsize;

	      /*----------------- Start Segment -----------------------------*/
	      if (cs != OBJECT)
		{
		  cs = OBJECT;
		  if (ps == OBJECT)
		    {
		      if (buffers->start[co] == UNKNOWN)
			{
			  buffers->marker[xl] = 'S';
			  buffers->start[co] = xl;
			}
		      else  buffers->marker[xl] = 's';
		    }
		  else
		    {
		      buffers->psstack[pstop++] = ps;
		      buffers->marker[xl] = 'S';
		      buffers->start[++co] = xl;
		      ps = COMPLETE;
		      buffers->info[co] = initinfo;
		    }
		}
	    }

	  /*-------------------Process New Marker ---------------------------*/
	  if (newmarker)
	    {
	      if (newmarker == 'S')
		{
		  buffers->psstack[pstop++] = ps;
		  if (cs == NONOBJECT)
		    {
		      buffers->psstack[pstop++] = COMPLETE;
		      buffers->info[++co] = buffers->store[xl];
		      buffers->start[co] = UNKNOWN;
		    }
		  else
		    update(&buffers->info[co], &buffers->store[xl], plist);
		  ps = OBJECT;
		}

	      else if (newmarker == 's')
		{
		  if ((cs == OBJECT) && (ps == COMPLETE))
		    {
		      pstop--;
		      xl2 = buffers->start[co];
		      update(&buffers->info[co-1], &buffers->info[co], plist);
		      if (buffers->start[--co] == UNKNOWN)
			buffers->start[co] = xl2;
		      else
			buffers->marker[xl2] = 's';
		    }
		  ps = OBJECT;
		}
	      else if (newmarker == 'f')
		ps = INCOMPLETE;
	      else if (newmarker == 'F')
		{
		  ps = buffers->psstack[--pstop];
		  if ((cs == NONOBJECT) && (ps == COMPLETE))
		    {
		      if (buffers->start[co] == UNKNOWN)
			{
			  if ((int)buffers->info[co].pixnb >= deb_maxarea)
			    {
			      if (objlist->nobj>=nobjm)
				if (!(obj = objlist->obj = (objstruct *)
				      realloc(obj, (nobjm+=nobjm/2)*
					      sizeof(objstruct))))
				  {
				    out = MEMORY_ALLOC_ERROR;
				    goto exit_lutz;
				  }
			      lutzsort(&buffers->info[co], objlist);
			    }
			}
		      else
			{
			  buffers->marker[buffers->end[co]] = 'F';
			  buffers->store[buffers->start[co]] = buffers->info[co];
			}
		      co--;
		      ps = buffers->psstack[--pstop];
		    }
		}
	    }
	  /* end process new marker -----------------------------------------*/

	  if (luflag)
	    update (&buffers->info[co],&curpixinfo, plist);
	  else
	    {
	      /* ----------------- End Segment ------------------------------*/
	      if (cs == OBJECT)
		{
		  cs = NONOBJECT;
		  if (ps != COMPLETE)
		    {
		      buffers->marker[xl] = 'f';
		      buffers->end[co] = xl;
		    }
		  else
		    {
		      ps = buffers->psstack[--pstop];
		      buffers->marker[xl] = 'F';
		      buffers->store[buffers->start[co]] = buffers->info[co];
		      co--;
		    }
		}
	    }
	}
    }

 exit_lutz:

  if (objlist->nobj && out == RETURN_OK)
    {
      if (!(objlist->obj=
	    realloc(obj, objlist->nobj*sizeof(objstruct))))
	out = MEMORY_ALLOC_ERROR;
    }
  else
    {
      free(obj);
      objlist->obj = NULL;
    }

  if (cn && out == RETURN_OK)
    {
      if (!(objlist->plist=realloc(plist,cn)))
	out = MEMORY_ALLOC_ERROR;
    }
  else
    {
      free(objlist->plist);
      objlist->plist = NULL;
    }

  return out;
}

/********************************* lutzsort ***********************************/
/*
Add an object to the object list based on info (pixel info)
*/
void  lutzsort(infostruct *info, objliststruct *objlist)
{
  objstruct *obj = objlist->obj+objlist->nobj;

  memset(obj, 0, (size_t)sizeof(objstruct));
  obj->firstpix = info->firstpix;
  obj->lastpix = info->lastpix;
  obj->flag = info->flag;
  objlist->npix += info->pixnb;

  preanalyse(objlist->nobj, objlist);

  objlist->nobj++;
}

/********************************* update ************************************/
/*
update object's properties each time one of its pixels is scanned by lutz()
*/
void  update(infostruct *infoptr1, infostruct *infoptr2, pliststruct *pixel)
{
  infoptr1->pixnb += infoptr2->pixnb;
  infoptr1->flag |= infoptr2->flag;
  if (infoptr1->firstpix == -1)
    {
      infoptr1->firstpix = infoptr2->firstpix;
      infoptr1->lastpix = infoptr2->lastpix;
    }
  else if (infoptr2->lastpix != -1)
    {
      PLIST(pixel+infoptr1->lastpix, nextpix) = infoptr2->firstpix;
      infoptr1->lastpix = infoptr2->lastpix;
    }
}
