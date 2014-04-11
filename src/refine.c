/*
*				refine.c
*
* Deblend sources based on their pixel lists.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <math.h>
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

static objliststruct *objlist;
static short	     *son, *ok;

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
			out;

  out = RETURN_OK;

  xn = deblend_nthresh;

  /* ---- initialize lists of objects */
  debobjlist.obj = debobjlist2.obj =  NULL;
  debobjlist.plist = debobjlist2.plist = NULL;
  debobjlist.nobj = debobjlist2.nobj = 0;
  debobjlist.npix = debobjlist2.npix = 0;
  objlistout->thresh = debobjlist2.thresh = objlistin->thresh;
  memset(objlist, 0, (size_t)xn*sizeof(objliststruct));

  for (l=0; l<objlistin->nobj && out==RETURN_OK; l++)
      {
	dthresh0 = objlistin->obj[l].dthresh;
	
	objlistout->dthresh = debobjlist2.dthresh = dthresh0;
	if ((out = addobj(l, objlistin, &objlist[0])) == RETURN_ERROR)
	  goto exit_parcelout;
	if ((out = addobj(l, objlistin, &debobjlist2)) == RETURN_ERROR)
	  goto exit_parcelout;
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
		out = RETURN_ERROR;
		goto exit_parcelout;
	      }

	    for (i=0; i<objlist[k-1].nobj; i++)
	      {
		if ((out=lutz(objlistin, l, &objlist[k-1].obj[i],
			      &debobjlist, minarea)) == RETURN_ERROR)
		  goto exit_parcelout;

		for (j=h=0; j<debobjlist.nobj; j++)
		  if (belong(j, &debobjlist, i, &objlist[k-1]))
		    {
		      debobjlist.obj[j].dthresh = debobjlist.dthresh;
		      m = addobj(j, &debobjlist, &objlist[k]);
		      if (m==RETURN_ERROR || m>=NSONMAX)
			{
			  out = RETURN_ERROR;
			  goto exit_parcelout;
			}
		      if (h>=nbm-1)
			if (!(son = (short *)
			      realloc(son,xn*NSONMAX*(nbm+=16)*sizeof(short))))
			  {
			    out = RETURN_ERROR;
			    goto exit_parcelout;
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
			  if ((out = addobj(j, &objlist[k+1], &debobjlist2))
			      == RETURN_ERROR)
			    goto exit_parcelout;
			}
		    ok[k+xn*i] = (short)0;
		  }
	      }
	  }

	if (ok[0])
	  out = addobj(0, &debobjlist2, objlistout);
	else
	  out = gatherup(&debobjlist2, objlistout);

      exit_parcelout:
	
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
  
  return out;
}


/******************************* allocparcelout ******************************/
/*
Allocate the memory allocated by global pointers in refine.c
*/
void	allocparcelout(int deblend_nthresh)
  {
  QMALLOC(son, short,  deblend_nthresh*NSONMAX*NBRANCH);
  QMALLOC(ok, short,  deblend_nthresh*NSONMAX);
  QMALLOC(objlist, objliststruct,  deblend_nthresh);

  return;
  }

/******************************* freeparcelout *******************************/
/*
Free the memory allocated by global pointers in refine.c
*/
void	freeparcelout(void)
  {
  QFREE(son);
  QFREE(ok);
  QFREE(objlist);
  return;
  }

/********************************* gatherup **********************************/
/*
Collect faint remaining pixels and allocate them to their most probable
progenitor.
*/
int gatherup(objliststruct *objlistin, objliststruct *objlistout)
{
  char		*bmp;
  float	        *amp, *p, dx,dy, drand, dist, distmin;
  objstruct	*objin = objlistin->obj, *objout, *objt;

  pliststruct	*pixelin = objlistin->plist, *pixelout, *pixt,*pixt2;

  int		i,k,l, *n, iclst, npix, bmwidth,
		nobj = objlistin->nobj, xs,ys, x,y, out;

  out = RETURN_OK;

  objlistout->dthresh = objlistin->dthresh;
  objlistout->thresh = objlistin->thresh;

  QMALLOC(amp, float, nobj);
  QMALLOC(p, float, nobj);
  QMALLOC(n, int, nobj);

  for (i=1; i<nobj; i++)
    preanalyse(i, objlistin, ANALYSE_FULL);

  p[0] = 0.0;
  bmwidth = objin->xmax - (xs=objin->xmin) + 1;
  npix = bmwidth * (objin->ymax - (ys=objin->ymin) + 1);
  if (!(bmp = (char *)calloc(1, npix*sizeof(char))))
    {
      bmp = 0;
      out = RETURN_ERROR;
      goto exit_gatherup;
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
      
      if ((n[i] = addobj(i, objlistin, objlistout)) == RETURN_ERROR)
	{
	  out = RETURN_ERROR;
	  goto exit_gatherup;
	}

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
      out = RETURN_ERROR;
      goto exit_gatherup;
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
    out = GATHERUP_MEMORY_ERROR;

 exit_gatherup:

  free(bmp);
  free(amp);
  free(p);
  free(n);

  return out;
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
Add an object to an objlist.
*/
int	addobj(int objnb, objliststruct *objl1, objliststruct *objl2)
{
  objstruct	*objl2obj;
  pliststruct	*plist1 = objl1->plist, *plist2 = objl2->plist;
  int		fp, i, j, npx, objnb2;
  
  j = (fp = objl2->npix)*plistsize;
  objnb2 = objl2->nobj;

  /* Update the object list */
  if (objl2->nobj)
    {
      if (!(objl2obj = (objstruct *)realloc(objl2->obj, (++objl2->nobj) *
					    sizeof(objstruct))))
	goto exit_addobj;
    }
  else
    if (!(objl2obj = (objstruct *)malloc((++objl2->nobj)*sizeof(objstruct))))
      goto exit_addobj;

  /* Update the pixel list */
  npx = objl1->obj[objnb].fdnpix;
  if (fp)
    {
      if (!(plist2 = (pliststruct *)realloc(plist2,
					    (objl2->npix+=npx) * plistsize)))
	goto exit_addobj;
    }
  else
    if (!(plist2=(pliststruct *)malloc((objl2->npix=npx)*plistsize)))
      goto exit_addobj;

  objl2->obj = objl2obj;
  objl2->plist = plist2;
  
  plist2 += j;
  for(i=objl1->obj[objnb].firstpix; i!=-1; i=PLIST(plist1+i,nextpix))
    {
      memcpy(plist2, plist1+i, (size_t)plistsize);
      PLIST(plist2,nextpix) = (j+=plistsize);
      plist2 += plistsize;
    }
  
  PLIST(plist2-=plistsize, nextpix) = -1;
  
  objl2->obj[objnb2] = objl1->obj[objnb];
  objl2->obj[objnb2].firstpix = fp*plistsize;
  objl2->obj[objnb2].lastpix = j-plistsize;
  return objnb2;
  
 exit_addobj:

  objl2->nobj--;
  objl2->npix = fp;
  return RETURN_ERROR;
}



