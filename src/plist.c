/*
*				plist.c
*
* Manage pixel lists.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include	<stdio.h>
#include	<stdlib.h>
#include "sep.h"


/******************************** createsubmap *******************************
PURPOSE Create pixel-index submap for deblending.
OUTPUT  RETURN_OK if success, RETURN_FATAL_ERROR otherwise (memory overflow).
*/
int	createsubmap(objliststruct *objlist, int no)
  {
    objstruct	*obj;
    pliststruct	*pixel, *pixt;
    int		i, n, xmin,ymin, w, *pix, *pt;

    obj = objlist->obj+no;
    pixel = objlist->plist;

    obj->subx = xmin = obj->xmin;
    obj->suby = ymin = obj->ymin;
    obj->subw = w = obj->xmax - xmin + 1;
    obj->subh = obj->ymax - ymin + 1;
    n = w*obj->subh;
    if (!(obj->submap = pix = (int *)malloc(n*sizeof(int))))
      return RETURN_FATAL_ERROR;
    pt = pix;
    for (i=n; i--;)
      *(pt++) = -1;
    
    for (i=obj->firstpix; i!=-1; i=PLIST(pixt,nextpix))
      {
	pixt = pixel+i;
	*(pix+(PLIST(pixt,x)-xmin) + (PLIST(pixt,y)-ymin)*w) = i;
      }
    
    return RETURN_OK;
  }


/****************************** init_plist ************************************
PURPOSE	initialize a pixel-list and its components.
 ***/
void	init_plist(PIXTYPE *conv, PIXTYPE *cdwfield)

  {
   pbliststruct	*pbdum = NULL;
   int		i;

  plistsize = sizeof(pbliststruct);
  plistoff_value = (char *)&pbdum->value - (char *)pbdum;

  if (conv)
    {
    plistexist_cdvalue = 1;
    plistoff_cdvalue = plistsize;
    plistsize += sizeof(PIXTYPE);
    }
  else
    {
    plistexist_cdvalue = 0;
    plistoff_cdvalue = plistoff_dvalue;
    }

  if (cdwfield)
    {
    plistexist_var = 1;
    plistoff_var = plistsize;
    plistsize += sizeof(PIXTYPE);
    }
  else
    plistexist_var = 0;

  if (cdwfield)
    {
    plistexist_dthresh = 1;
    plistoff_dthresh = plistsize;
    plistsize += sizeof(PIXTYPE);
    }
  else
    plistexist_dthresh = 0;

  return;

  /* can we remove these? */
  plistexist_dvalue = 0;
  plistoff_dvalue = plistoff_value;
  plistexist_flag = 0;
  plistexist_wflag = 0;

  return;
  }

