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
#include "extract.h"

#ifndef	RAND_MAX
#define	RAND_MAX 2147483647
#endif
#define	NSONMAX	1024  /* max. number per level */
#define	NBRANCH	16    /* starting number per branch */

int addobj(int, objliststruct *, objliststruct *);
int belong(int, objliststruct *, int, objliststruct *);
int gatherup(objliststruct *, objliststruct *);

static objliststruct *objlist=NULL;
static short	     *son=NULL, *ok=NULL;

/******************************** parcelout **********************************
PURPOSE Divide a list of isophotal detections in several parts (deblending).
INPUT   input objlist,
        output objlist,
OUTPUT  RETURN_OK if success, RETURN_ERROR otherwise (memory overflow).
NOTES   Even if the object is not deblended, the output objlist threshold is
        recomputed if a variable threshold is used.
 ***/
int parcelout(objliststruct *objlistin, objliststruct *objlistout,
	      int deblend_nthresh, double deblend_mincont, int minarea)
{
  objstruct		*obj;
  static objliststruct	debobjlist, debobjlist2;
  double		dthresh, dthresh0, value0;
  int			h,i,j,k,l,m,
                        xn,
			nbm = NBRANCH,
			status;

  status = RETURN_OK;
  xn = deblend_nthresh;

  /* ---- initialize lists of objects */
  debobjlist.obj = debobjlist2.obj =  NULL;
  debobjlist.plist = debobjlist2.plist = NULL;
  debobjlist.nobj = debobjlist2.nobj = 0;
  debobjlist.npix = debobjlist2.npix = 0;
  objlistout->thresh = debobjlist2.thresh = objlistin->thresh;
  memset(objlist, 0, (size_t)xn*sizeof(objliststruct));

  for (l=0; l<objlistin->nobj && status==RETURN_OK; l++)
      {
	dthresh0 = objlistin->obj[l].dthresh;
	
	objlistout->dthresh = debobjlist2.dthresh = dthresh0;
	if ((status = addobj(l, objlistin, &objlist[0])) != RETURN_OK)
	  goto exit;
	if ((status = addobj(l, objlistin, &debobjlist2)) != RETURN_OK)
	  goto exit;
	value0 = objlist[0].obj[0].fdflux*deblend_mincont;
	ok[0] = (short)1;
	for (k=1; k<xn; k++)
	  {
	    /*------ Calculate threshold */
	    dthresh = objlistin->obj[l].fdpeak;
	    debobjlist.dthresh = dthresh > 0.0? 
	      dthresh0*pow(dthresh/dthresh0,(double)k/xn) : dthresh0;
	    
	    /*--------- Build tree (bottom->up) */
	    if (objlist[k-1].nobj>=NSONMAX)
	      {
		status = RETURN_ERROR;
		goto exit;
	      }

	    for (i=0; i<objlist[k-1].nobj; i++)
	      {
		if ((status=lutz(objlistin, l, &objlist[k-1].obj[i],
				 &debobjlist, minarea)) != RETURN_OK)
		  goto exit;

		for (j=h=0; j<debobjlist.nobj; j++)
		  if (belong(j, &debobjlist, i, &objlist[k-1]))
		    {
		      debobjlist.obj[j].dthresh = debobjlist.dthresh;
		      if ((status = addobj(j, &debobjlist, &objlist[k]))
			  != RETURN_OK)
			goto exit;
		      m = objlist[k].nobj - 1;
		      if (m>=NSONMAX)
			{
			  status = RETURN_ERROR;
			  goto exit;
			}
		      if (h>=nbm-1)
			if (!(son = (short *)
			      realloc(son,xn*NSONMAX*(nbm+=16)*sizeof(short))))
			  {
			    status = RETURN_ERROR;
			    goto exit;
			  }
		      son[k-1+xn*(i+NSONMAX*(h++))] = (short)m;
		      ok[k+xn*m] = (short)1;
		    }
		son[k-1+xn*(i+NSONMAX*h)] = (short)-1;
	      }
	  }

	/*------- cut the right branches (top->down) */
	for (k = xn-2; k>=0; k--)
	  {
	    obj = objlist[k+1].obj;
	    for (i=0; i<objlist[k].nobj; i++)
	      {
		for (m=h=0; (j=(int)son[k+xn*(i+NSONMAX*h)])!=-1; h++)
		  {
		    if (obj[j].fdflux - obj[j].dthresh*obj[j].fdnpix > value0)
		      m++;
		    ok[k+xn*i] &= ok[k+1+xn*j];
		  }
		if (m>1)	
		  {
		    for (h=0; (j=(int)son[k+xn*(i+NSONMAX*h)])!=-1; h++)
		      if (ok[k+1+xn*j] &&
			  obj[j].fdflux-obj[j].dthresh*obj[j].fdnpix > value0)
			{
			  objlist[k+1].obj[j].flag |= OBJ_MERGED
			    | ((OBJ_ISO_PB|OBJ_APERT_PB|OBJ_OVERFLOW)
			       &debobjlist2.obj[0].flag);
			  if ((status = addobj(j, &objlist[k+1], &debobjlist2))
			      == RETURN_ERROR)
			    goto exit;
			}
		    ok[k+xn*i] = (short)0;
		  }
	      }
	  }

	if (ok[0])
	  status = addobj(0, &debobjlist2, objlistout);
	else
	  status = gatherup(&debobjlist2, objlistout);

      exit:
	free(debobjlist2.obj);
	free(debobjlist2.plist);
	
	for (k=0; k<xn; k++)
	  {
	    free(objlist[k].obj);
	    free(objlist[k].plist);
	  }
      }

  free(debobjlist.obj);
  free(debobjlist.plist);
  
  return status;
}


/******************************* allocparcelout ******************************/
/*
Allocate the memory allocated by global pointers in refine.c
*/
int allocparcelout(int deblend_nthresh)
{
  int status=RETURN_OK;
  QMALLOC(son, short,  deblend_nthresh*NSONMAX*NBRANCH, status);
  QMALLOC(ok, short,  deblend_nthresh*NSONMAX, status);
  QMALLOC(objlist, objliststruct, deblend_nthresh, status);

  return status;
 exit:
  freeparcelout();
  return status;
}

/******************************* freeparcelout *******************************/
/*
Free the memory allocated by global pointers in refine.c
*/
void freeparcelout(void)
{
  free(son);
  son = NULL;
  free(ok);
  ok = NULL;
  free(objlist);
  objlist = NULL;
  return;
}

/********************************* gatherup **********************************/
/*
Collect faint remaining pixels and allocate them to their most probable
progenitor.
*/
int gatherup(objliststruct *objlistin, objliststruct *objlistout)
{
  char        *bmp;
  float       *amp, *p, dx,dy, drand, dist, distmin;
  objstruct   *objin = objlistin->obj, *objout, *objt;

  pliststruct *pixelin = objlistin->plist, *pixelout, *pixt,*pixt2;

  int         i,k,l, *n, iclst, npix, bmwidth,
              nobj = objlistin->nobj, xs,ys, x,y, status;

  bmp = NULL;
  amp = p = NULL;
  n = NULL;
  status = RETURN_OK;

  objlistout->dthresh = objlistin->dthresh;
  objlistout->thresh = objlistin->thresh;

  QMALLOC(amp, float, nobj, status);
  QMALLOC(p, float, nobj, status);
  QMALLOC(n, int, nobj, status);

  for (i=1; i<nobj; i++)
    preanalyse(i, objlistin, ANALYSE_FULL);

  p[0] = 0.0;
  bmwidth = objin->xmax - (xs=objin->xmin) + 1;
  npix = bmwidth * (objin->ymax - (ys=objin->ymin) + 1);
  if (!(bmp = (char *)calloc(1, npix*sizeof(char))))
    {
      bmp = NULL;
      status = RETURN_ERROR;
      goto exit;
    }
  
  for (objt = objin+(i=1); i<nobj; i++, objt++)
    {
      /*-- Now we have passed the deblending section, reset thresholds */
      objt->dthresh = objlistin->dthresh;
      objt->thresh = objlistin->thresh;

      /* ------------	flag pixels which are already allocated */
      for (pixt=pixelin+objin[i].firstpix; pixt>=pixelin;
	   pixt=pixelin+PLIST(pixt,nextpix))
	bmp[(PLIST(pixt,x)-xs) + (PLIST(pixt,y)-ys)*bmwidth] = '\1';
      
      status = addobj(i, objlistin, objlistout);
      if (status == RETURN_ERROR)
	goto exit;
      n[i] = objlistout->nobj - 1;

      dist = objt->fdnpix/(2*PI*objt->abcor*objt->a*objt->b);
      amp[i] = dist<70.0? objt->dthresh*expf(dist) : 4.0*objt->fdpeak;

      /* ------------ limitate expansion ! */
      if (amp[i]>4.0*objt->fdpeak)
	amp[i] = 4.0*objt->fdpeak;
    }

  objout = objlistout->obj;		/* DO NOT MOVE !!! */

  if (!(pixelout=(pliststruct *)realloc(objlistout->plist,
					(objlistout->npix + npix)*plistsize)))
    {
      status = RETURN_ERROR;
      goto exit;
    }
  
  objlistout->plist = pixelout;
  k = objlistout->npix;
  iclst = 0;				/* To avoid gcc -Wall warnings */
  for (pixt=pixelin+objin->firstpix; pixt>=pixelin;
       pixt=pixelin+PLIST(pixt,nextpix))
    {
      x = PLIST(pixt,x);
      y = PLIST(pixt,y);
      if (!bmp[(x-xs) + (y-ys)*bmwidth])
	{
	  pixt2 = pixelout + (l=(k++*plistsize));
	  memcpy(pixt2, pixt, (size_t)plistsize);
	  PLIST(pixt2, nextpix) = -1;
	  distmin = 1e+31;
	  for (objt = objin+(i=1); i<nobj; i++, objt++)
	    {
	      dx = x - objt->mx;
	      dy = y - objt->my;
	      dist=0.5*(objt->cxx*dx*dx+objt->cyy*dy*dy+objt->cxy*dx*dy)/objt->abcor;
	      p[i] = p[i-1] + (dist<70.0?amp[i]*expf(-dist) : 0.0);
	      if (dist<distmin)
		{
		  distmin = dist;
		  iclst = i;
		}
	    }			
	  if (p[nobj-1] > 1.0e-31)
	    {
	      drand = p[nobj-1]*rand()/RAND_MAX;
	      for (i=1; i<nobj && p[i]<drand; i++);
	      if (i==nobj)
		i=iclst;
	    }
	  else
	    i = iclst;
	  objout[n[i]].lastpix=PLIST(pixelout+objout[n[i]].lastpix,nextpix)=l;
	}
    }

  objlistout->npix = k;
  if (!(objlistout->plist = (pliststruct *)realloc(pixelout,
						   objlistout->npix*plistsize)))
    status = GATHERUP_MEMORY_ERROR;

 exit:
  free(bmp);
  free(amp);
  free(p);
  free(n);

  return status;
}

/**************** belong (originally in manobjlist.c) ************************/
/*
say if an object is "included" in another.
*/

int belong(int corenb, objliststruct *coreobjlist,
	   int shellnb, objliststruct *shellobjlist)
{
  objstruct	*cobj = &(coreobjlist->obj[corenb]),
                *sobj = &(shellobjlist->obj[shellnb]);
  pliststruct	*cpl = coreobjlist->plist, *spl = shellobjlist->plist, *pixt;

  int		xc=PLIST(cpl+cobj->firstpix,x), yc=PLIST(cpl+cobj->firstpix,y);

  for (pixt = spl+sobj->firstpix; pixt>=spl; pixt = spl+PLIST(pixt,nextpix))
    if ((PLIST(pixt,x) == xc) && (PLIST(pixt,y) == yc))
      return 1;

  return 0;
}


/********** addobj (originally in manobjlist.c) ******************************/
/*
 *  Add object number `objnb` from list `objl1` to list `objl2`.
 */
int addobj(int objnb, objliststruct *objl1, objliststruct *objl2)
{
  objstruct	*objl2obj;
  pliststruct	*plist1 = objl1->plist, *plist2 = objl2->plist;
  int		fp, i, j, npx, objnb2;
  
  j = (fp = objl2->npix)*plistsize; /* fp = 2nd lists's plist size in pix */
                                    /* j = 2nd list's plist size in bytes */
  objnb2 = objl2->nobj;             /* # of objects currently in 2nd list*/

  /* Allocate space in `objl2` for the new object */
  if (objl2->nobj)
    objl2obj = (objstruct *)realloc(objl2->obj, (++objl2->nobj) *
				    sizeof(objstruct));
  else
    objl2obj = (objstruct *)malloc((++objl2->nobj)*sizeof(objstruct));
  if (!objl2obj)
    goto exit_addobj;
  objl2->obj = objl2obj;

  /* Allocate space for the new object's pixels in 2nd list's plist */
  npx = objl1->obj[objnb].fdnpix;
  if (fp)
    plist2 = (pliststruct *)realloc(plist2, (objl2->npix+=npx) * plistsize);
  else
    plist2 = (pliststruct *)malloc((objl2->npix=npx)*plistsize);
  if (!plist2)
    goto exit_addobj;
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
 exit_addobj:
  objl2->nobj--;
  objl2->npix = fp;
  return RETURN_ERROR;
}



