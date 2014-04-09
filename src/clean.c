/*
*				clean.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include	<math.h>
#include	<stdlib.h>
#include	<string.h>
#include	"sep.h"

/*------------------------------- variables ---------------------------------*/

static LONG		*cleanvictim;

/******************************* initclean **********************************
initialize cleanvictim and cleanobjlist.
 ***/
void	initclean(void)
  {
  if (CLEAN_FLAG)
    QMALLOC(cleanvictim, LONG, CLEAN_STACKSIZE);
  QMALLOC(cleanobjlist, objliststruct, 1);
  cleanobjlist->obj = NULL;
  cleanobjlist->plist = NULL;
  cleanobjlist->nobj = cleanobjlist->npix = 0;

  return;
  }


/******************************** endclean **********************************
PROTO   void endclean(void)
PURPOSE End things related to CLEANing.
INPUT   -.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP & Leiden & ESO)
VERSION 03/07/97
 ***/
void	endclean(void)
  {
  if (prefs.clean_flag)
    free(cleanvictim);
  free(cleanobjlist);
  return;
  }


/********************************** clean ***********************************
PURPOSE Remove object from frame -buffer and put it in the "CLEANlist".
INPUT   Object number,
        Object list (source).
OUTPUT  0 if the object was CLEANed, 1 otherwise.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden & ESO)
VERSION 08/02/2001
 ***/
int	clean(picstruct *field, picstruct *dfield, int objnb,
		objliststruct *objlistin)
  {
   objstruct		*objin, *obj;
   int			i,j,k;
   double		amp,ampin,alpha,alphain, unitarea,unitareain,beta,val;
   float	       	dx,dy,rlim;

  objin = objlistin->obj+objnb;
  beta = prefs.clean_param;
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
        val = 1+alphain*(objin->cxx*dx*dx+objin->cyy*dy*dy+objin->cxy*dx*dy);
        if (val>1.0
	    && ((float)(val<1e10?ampin*pow(val,-beta) : 0.0) > obj->mthresh))
/*------- the newcomer puts that object in its menu! */
          cleanvictim[j++] = i;
        }
      else
        {
        unitarea = PI*obj->a*obj->b;
        amp = obj->fdflux/(2*unitarea*obj->abcor);
        alpha = (pow(amp/obj->dthresh, 1.0/beta) - 1)*unitarea/obj->fdnpix;
        val = 1+alpha*(obj->cxx*dx*dx+obj->cyy*dy*dy+obj->cxy*dx*dy);
        if (val>1.0
	    && ((float)(val<1e10?amp*pow(val,-beta) : 0.0) > objin->mthresh))
          {
/*------- the newcomer is eaten!! */
          mergeobject(objin, obj);
          if (prefs.blank_flag)
            {
/*---------- Paste back ``CLEANed'' object pixels before forgetting them */
            if (objin->blank)
              {
              pasteimage(field, objin->blank, objin->subw, objin->subh,
			objin->subx, objin->suby);
              free(objin->blank);
              }
            if (objin->dblank)
              {
              pasteimage(dfield, objin->dblank, objin->subw, objin->subh,
			objin->subx, objin->suby);
              free(objin->dblank);
              }
            }

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
    if (prefs.blank_flag)
      {
/*---- Paste back ``CLEANed'' object pixels before forgetting them */
      if (obj->blank)
        {
        pasteimage(field, obj->blank, obj->subw, obj->subh,
		obj->subx, obj->suby);
        free(obj->blank);
        }
      if (obj->dblank)
        {
        pasteimage(dfield, obj->dblank, obj->subw, obj->subh,
		obj->subx, obj->suby);
        free(obj->dblank);
        }
      }
    subcleanobj(k);
    }

  return 1;
  }


/******************************* addcleanobj ********************************/
/*
Add an object to the "cleanobjlist".
*/
void	addcleanobj(objstruct *objin)

  {
   int		margin, y;
   float	hh1,hh2;

/*Update the object list */
  if (cleanobjlist->nobj)
    {
    if (!(cleanobjlist->obj = (objstruct *)realloc(cleanobjlist->obj,
		    (++cleanobjlist->nobj) * sizeof(objstruct))))
      error(EXIT_FAILURE, "Not enough memory for ", "CLEANing");
    }
  else
    {
    if (!(cleanobjlist->obj = (objstruct *)malloc((++cleanobjlist->nobj)
		    * sizeof(objstruct))))
      error(EXIT_FAILURE, "Not enough memory for ", "CLEANing");
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
  if ((y=(int)(objin->my+0.49999)+prefs.cleanmargin)>objin->ycmax)
    objin->ycmax = y;
  objin->ycmin = objin->ymin-margin;
  if ((y=(int)(objin->my+0.49999)-prefs.cleanmargin)<objin->ycmin)
    objin->ycmin = y;

  cleanobjlist->obj[cleanobjlist->nobj-1] = *objin;

  return;
  }


/******************************** mergeobject *******************************/
/*
Merge twos objects from "objlist".
*/
void	mergeobject(objstruct *objslave,objstruct *objmaster)

  {
   checkstruct	*check;

  if ((check = prefs.check[CHECK_SEGMENTATION]))
    {
     ULONG	*pix;
     ULONG	colorslave = objslave->number,
		colormaster = objmaster->number;
     int	dx,dx0,dy,dpix;

    pix = (ULONG *)check->pix+check->width*objslave->ymin + objslave->xmin;
    dx0 = objslave->xmax-objslave->xmin+1;
    dpix = check->width-dx0;
    for (dy=objslave->ymax-objslave->ymin+1; dy--; pix += dpix)
      for (dx=dx0; dx--; pix++)
        if (*pix==colorslave)
          *pix = colormaster;
    }

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
  mergeflags(objmaster, objslave);

  return;
  }


/******************************* subcleanobj ********************************/
/* 
remove an object from a "cleanobjlist".
*/
void subcleanobj(int objnb)
{
  /* Update the object list */
  if (objnb>=cleanobjlist->nobj)
    error(EXIT_FAILURE, "*Internal Error*: no CLEAN object to remove ",
	  "in subcleanobj()");
  
  if (--cleanobjlist->nobj)
    {
      if (cleanobjlist->nobj != objnb)
	cleanobjlist->obj[objnb] = cleanobjlist->obj[cleanobjlist->nobj];
      if (!(cleanobjlist->obj = (objstruct *)realloc(cleanobjlist->obj,
						     cleanobjlist->nobj * sizeof(objstruct))))
	error(EXIT_FAILURE, "Not enough memory for ", "CLEANing");
    }
  else
    free(cleanobjlist->obj);
  
  return;
}

