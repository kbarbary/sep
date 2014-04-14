#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sep.h"
#include "extract.h"

#define DETECT_MAXAREA 0        /* replaces prefs.ext_maxarea */
#define MEMORY_PIXSTACK 300000  /* number of pixels in stack */
                                /* (replaces prefs.mem_pixstack) */
#define	WTHRESH_CONVFAC	1e-4    /* Factor to apply to weights when */
			        /* thresholding filtered weight-maps */

void convolve(PIXTYPE *, int, int, int, float *, int, int, PIXTYPE *);
int  sortit(infostruct *, objliststruct *, PIXTYPE *, PIXTYPE *, int,
	    int, double, objliststruct *, LONG *, int, double);
int  createsubmap(objliststruct *, int);
void plistinit(PIXTYPE *, PIXTYPE *);

/****************************** extract **************************************/
objliststruct *extract(PIXTYPE *cfield, PIXTYPE *cdwfield, int w, int h,
		       PIXTYPE dthresh, PIXTYPE athresh, PIXTYPE cdwthresh,
		       int threshabsolute, int minarea,
		       float *conv, int convw, int convh,
		       int deblend_nthresh, double deblend_mincont,
		       int clean_flag, double clean_param,
		       int *status)
{
  static infostruct	curpixinfo, *info, *store, initinfo, freeinfo, *victim;
  objliststruct       	objlist, *cleanobjlist;
  pliststruct		*pixel, *pixt; 
  char			*marker, newmarker;
  int			co, i,j, flag, luflag,pstop, xl,xl2,yl, cn,
			nposize, stacksize, maxpixnb,
                        varthreshflag, convn;
  short	       	        trunflag;
  PIXTYPE		thresh, relthresh, cdnewsymbol;
  PIXTYPE               *scan,*cdscan,*cdwscan,*wscan,*dumscan;
  float                 sum, *convnorm;
  pixstatus		cs, ps, *psstack;
  int			*start, *end;
  LONG                  *cleanvictim;

  *status = RETURN_OK;

  /*----- Beginning of the main loop: Initialisations  */
  pixel = NULL;
  convnorm = NULL;
  scan = wscan = cdscan = cdwscan = dumscan = NULL;
  victim = NULL;
  cleanvictim = NULL;
  info = NULL;
  store = NULL;
  marker = NULL;
  psstack = NULL;
  start = end = NULL;
  cleanobjlist = NULL;
  convn = 0;
  sum = 0.0;

  /* cdwfield is the detection weight-field if available (can be null)*/
  if (cdwthresh>BIG*WTHRESH_CONVFAC)
    cdwthresh = BIG*WTHRESH_CONVFAC;

  /* If WEIGHTing and no absolute thresholding, activate threshold scaling */
  varthreshflag = (cdwfield && threshabsolute==0);
  relthresh = varthreshflag ? thresh : 0.0;/* To avoid gcc warnings*/
  objlist.dthresh = dthresh;
  objlist.thresh = athresh;

  /*Allocate memory for buffers */
  stacksize = w+1;
  QMALLOC(info, infostruct, stacksize, *status);
  QCALLOC(store, infostruct, stacksize, *status);
  QMALLOC(marker, char, stacksize, *status);
  QMALLOC(dumscan, PIXTYPE, stacksize, *status);
  QMALLOC(psstack, pixstatus, stacksize, *status);
  QCALLOC(start, int, stacksize, *status);
  QMALLOC(end, int, stacksize, *status);
  if ((*status = lutzalloc(w, h)) != RETURN_OK)
    goto exit;
  if ((*status = allocparcelout(deblend_nthresh)) != RETURN_OK)
    goto exit;

  /* Some initializations */
  thresh = objlist.dthresh;
  initinfo.pixnb = 0;
  initinfo.flag = 0;
  initinfo.firstpix = initinfo.lastpix = -1;

  for (xl=0; xl<stacksize; xl++)
    {
    marker[xl]  = 0 ;
    dumscan[xl] = -BIG ;
    }

  co = pstop = 0;
  objlist.nobj = 1;
  curpixinfo.pixnb = 1;

  /* Init cleaning procedure */
  if (clean_flag)
    QMALLOC(cleanvictim, LONG, CLEAN_STACKSIZE, *status);
  QMALLOC(cleanobjlist, objliststruct, 1, *status);
  cleanobjlist->obj = NULL;
  cleanobjlist->plist = NULL;
  cleanobjlist->nobj = cleanobjlist->npix = 0;


  /* Allocate memory for the pixel list */
  plistinit(conv, cdwfield);
  if (!(pixel = objlist.plist = malloc(nposize=MEMORY_PIXSTACK*plistsize)))
    {
      *status = MEMORY_PIXSTACK_ERROR;
      goto exit;
    }

  /*----- at the beginning, "free" object fills the whole pixel list */
  freeinfo.firstpix = 0;
  freeinfo.lastpix = nposize-plistsize;
  pixt = pixel;
  for (i=plistsize; i<nposize; i += plistsize, pixt += plistsize)
    PLIST(pixt, nextpix) = i;
  PLIST(pixt, nextpix) = -1;

  if (conv)
    {
      /* allocate memory for convolved buffers */
      QMALLOC(cdscan, PIXTYPE, stacksize, *status);
      if (cdwfield)
	QCALLOC(cdwscan, PIXTYPE, stacksize, *status);

      /* normalize the filter */
      convn = convw * convh;
      QMALLOC(convnorm, PIXTYPE, convn, *status);
      for (i=0; i<convn; i++)
	sum += fabs(conv[i]);
      for (i=0; i<convn; i++)
	convnorm[i] = conv[i] / sum;
    }

  /*----- MAIN LOOP ------ */
  for (yl=0; yl<=h;)
    {

      ps = COMPLETE;
      cs = NONOBJECT;
    
      /* Need an empty line for Lutz' algorithm to end gracely */
      if (yl==h)
	{
	  if (conv)
	    {
	      free(cdscan);
	      cdscan = NULL;
	      if (cdwfield)
		{
		  free(cdwscan);
		  cdwscan = NULL;
		}
	    }
	  cdwscan = cdscan = dumscan;
	}

      else
	{
	  scan = cfield + yl*w;
	  if (cdwfield)
	    wscan = cdwfield + yl*w;

	  /* filter the lines */
	  if (conv)
	    {
	      convolve(cfield, w, h, yl, conv, convw, convh, cdscan);
	      if (cdwfield)
		convolve(cdwfield, w, h, yl, conv, convw, convh, cdwscan);
	    }
	  else
	    {
	      cdscan = scan;
	      cdwscan = wscan;
	    }	  
	}
      
      trunflag = (yl==0 || yl==h-1)? OBJ_TRUNC:0;
      
      for (xl=0; xl<=w; xl++)
	{

	  if (xl == w)
	    cdnewsymbol = -BIG;
	  else
	    cdnewsymbol = cdscan[xl];

	  newmarker = marker[xl];
	  marker[xl] = 0;

	  curpixinfo.flag = trunflag;
	  if (varthreshflag)
	    thresh = relthresh*sqrt((xl==w || yl==h)? 0.0:cdwscan[xl]);
	  luflag = cdnewsymbol > thresh?1:0;  /* is pixel above threshold? */

	  if (luflag)
	    {
	      /* flag the current object if we're near the image bounds */
	      if (xl==0 || xl==w-1)
		curpixinfo.flag |= OBJ_TRUNC;
	      
	      /* point pixt to first free pixel in pixel list */
	      /* and increment the "first free pixel" */
	      pixt = pixel + (cn=freeinfo.firstpix);
	      freeinfo.firstpix = PLIST(pixt, nextpix);
	      
	      /* Running out of pixels, the largest object becomes a "victim" */
	      if (freeinfo.firstpix==freeinfo.lastpix)
		{
		  sprintf(errdetail, "Pixel stack overflow at position %d,%d",
			  xl+1, yl+1);
		  *status = PIXSTACK_OVERFLOW_ERROR;
		  goto exit;
		  
		  /* NOTE: The above error was originally just a warning.
		     with the change to an error, the following code in this
		     if block is never executed. TODO: should this just
		     be a warning (or nothing?)
		  */

		  /* loop over pixels in row to find largest object */
		  maxpixnb = 0;
		  for (i=0; i<=w; i++)
		    if (store[i].pixnb>maxpixnb)
		      if (marker[i]=='S' || (newmarker=='S' && i==xl))
			{
			  flag = 0;
			  if (i<xl)
			    for (j=0; j<=co; j++)
			      flag |= (start[j]==i);
			  if (!flag)
			    maxpixnb = (victim = &store[i])->pixnb;
			}
		  for (j=1; j<=co; j++)
		    if (info[j].pixnb>maxpixnb)
		      maxpixnb = (victim = &info[j])->pixnb;
		  
		  if (!maxpixnb)
		    {
		      *status = FATAL_ERROR;
		      goto exit;
		    }
		  if (maxpixnb <= 1)
		    {
		      *status = PIXSTACK_OVERFLOW_ERROR;
		      goto exit;
		    }
		  freeinfo.firstpix = PLIST(pixel+victim->firstpix, nextpix);
		  PLIST(pixel+victim->lastpix, nextpix) = freeinfo.lastpix;
		  PLIST(pixel+(victim->lastpix=victim->firstpix), nextpix) = -1;
		  victim->pixnb = 1;
		  victim->flag |= OBJ_OVERFLOW;
		}
	      /*------------------------------------------------------------*/

	      curpixinfo.lastpix = curpixinfo.firstpix = cn;
	      PLIST(pixt, nextpix) = -1;
	      PLIST(pixt, x) = xl;
	      PLIST(pixt, y) = yl;
	      PLIST(pixt, value) = scan[xl];
	      if (PLISTEXIST(cdvalue))
		PLISTPIX(pixt, cdvalue) = cdnewsymbol;

	      /* Detect pixels with a low weight ----------------------------*/
	      if (PLISTEXIST(var))
		PLISTPIX(pixt, var) = wscan[xl];
	      
	      if (cs != OBJECT)
/*------------------------------- Start Segment -----------------------------*/
		{
		  cs = OBJECT;
		  if (ps == OBJECT)
		    {
		      if (start[co] == UNKNOWN)
			{
			  marker[xl] = 'S';
			  start[co] = xl;
			}
		      else
			marker[xl] = 's';
		    }
		  else
		    {
		      psstack[pstop++] = ps;
		      marker[xl] = 'S';
		      start[++co] = xl;
		      ps = COMPLETE;
		      info[co] = initinfo;
		    }
		}

	    } /* closes if pixel above threshold */

	  /* process new marker ---------------------------------------------*/
	  /* newmarker is marker[ ] at this pixel position. We'll only       */
	  /* enter this if marker[ ] was set on a previous loop iteration.   */
	  if (newmarker)
	    {
	      if (newmarker == 'S')
		{
		  psstack[pstop++] = ps;
		  if (cs == NONOBJECT)
		    {
		      psstack[pstop++] = COMPLETE;
		      info[++co] = store[xl];
		      start[co] = UNKNOWN;
		    }
		  else
		    update(&info[co], &store[xl], pixel);
		  ps = OBJECT;
		}

	      else if (newmarker == 's')
		{
		  if ((cs == OBJECT) && (ps == COMPLETE))
		    {
		      pstop--;
		      xl2 = start[co];
		      update (&info[co-1],&info[co], pixel);
		      if (start[--co] == UNKNOWN)
			start[co] = xl2;
		      else
			marker[xl2] = 's';
		    }
		  ps = OBJECT;
		}

	      else if (newmarker == 'f')
		ps = INCOMPLETE;

	      else if (newmarker == 'F')
		{
		  ps = psstack[--pstop];
		  if ((cs == NONOBJECT) && (ps == COMPLETE))
		    {
		      if (start[co] == UNKNOWN)
			{
			  if ((int)info[co].pixnb >= minarea)
			    {
			      *status = sortit(&info[co], &objlist,
					       cdwscan, wscan, minarea,
					       clean_flag, clean_param,
					       cleanobjlist, cleanvictim,
					       deblend_nthresh,
					       deblend_mincont);
			      if (*status != RETURN_OK)
				goto exit;
			    }

			  /* free the chain-list */
			  PLIST(pixel+info[co].lastpix, nextpix) =
			    freeinfo.firstpix;
			  freeinfo.firstpix = info[co].firstpix;
			}
		      else
			{
			  marker[end[co]] = 'F';
			  store[start[co]] = info[co];
			}
		      co--;
		      ps = psstack[--pstop];
		    }
		}
	    }
	  /* end of if (newmarker) ------------------------------------------*/

	  if (luflag)
	    update(&info[co], &curpixinfo, pixel);

	  /*---------------------- End Segment ------------------------------*/
	  else if (cs == OBJECT)
	    {
	      cs = NONOBJECT;
	      if (ps != COMPLETE)
		{
		  marker[xl] = 'f';
		  end[co] = xl;
		}
	      else
		{
		  ps = psstack[--pstop];
		  marker[xl] = 'F';
		  store[start[co]] = info[co];
		  co--;
		}
	    }

	} /*------------ End of the loop over the x's -----------------------*/

      /*-- Prepare markers for the next line */
      yl++;

    } /*--------------------- End of the loop over the y's ------------------*/

  
  /* Now that all "detected" pixels have been removed, analyse detections */
  /* removed this!
  ontotal = 0;
  for (j=cleanobjlist->nobj; j--;)
    {
      ontotal = thecat.ntotal;
      endobject(field, dfield, wfield, cdwfield, 0, cleanobjlist);
      subcleanobj(0);
    }
  */

 exit:
  if (clean_flag)
    free(cleanvictim);
  freeparcelout();
  free(pixel);
  lutzfree();
  free(info);
  free(store);
  free(marker);
  free(dumscan);
  free(psstack);
  free(start);
  free(end);
  if (conv)
    free(convnorm);

  if (*status != RETURN_OK)
    {
      free(cdscan);   /* only need to free these in case of early exit */
      free(cdwscan);
      return NULL;
    }

  return cleanobjlist;
}


/********************************* sortit ************************************/
/*
build the object structure.
*/
int sortit(infostruct *info, objliststruct *objlist,
	   PIXTYPE *cdwscan, PIXTYPE *wscan, int minarea,
	   int clean_flag, double clean_param, objliststruct *cleanobjlist,
	   LONG *cleanvictim, int deblend_nthresh, double deblend_mincont)
{
  objliststruct	        objlistout, *objlist2;
  static objstruct	obj;
  objstruct		*cobj;
  int 			i,j,n;
  int status=RETURN_OK;

  objlistout.obj = NULL;
  objlistout.plist = NULL;
  objlistout.nobj = objlistout.npix = 0;

  /*----- Allocate memory to store object data */
  objlist->obj = &obj;
  objlist->nobj = 1;

  memset(&obj, 0, (size_t)sizeof(objstruct));
  objlist->npix = info->pixnb;
  obj.firstpix = info->firstpix;
  obj.lastpix = info->lastpix;
  obj.flag = info->flag;
  obj.dthresh = objlist->dthresh;
  obj.thresh = objlist->thresh;

  preanalyse(0, objlist, ANALYSE_FAST);

  /*----- Check if the current strip contains the lower isophote
    (it always should since the "current strip" is the entire image!) */
  if ((int)obj.ymin < 0)
    obj.flag |= OBJ_ISO_PB;

  if (!(obj.flag & OBJ_OVERFLOW) && (createsubmap(objlist, 0) == RETURN_OK))
    {
      if (parcelout(objlist, &objlistout, deblend_nthresh, deblend_mincont,
		    minarea) == RETURN_OK)
	objlist2 = &objlistout;
      else
	{
	  objlist2 = objlist;
	  for (i=0; i<objlist2->nobj; i++)
	    objlist2->obj[i].flag |= OBJ_DOVERFLOW;
	  sprintf(errdetail, "Deblending overflow for detection at %.0f,%.0f", obj.mx+1, obj.my+1);
	  status = DEBLEND_OVERFLOW_ERROR;
	  goto exit_sortit;
	}
      free(obj.submap);
    }
  else
    objlist2 = objlist;
  
  for (i=0; i<objlist2->nobj; i++)
    {
      preanalyse(i, objlist2, ANALYSE_FULL|ANALYSE_ROBUST);
      if (DETECT_MAXAREA && objlist2->obj[i].fdnpix > DETECT_MAXAREA)
	continue;

      /* removed analyse (defined in analyse.c in sextractor) */
      /* analyse(field, dfield, i, objlist2); */
      
      cobj = objlist2->obj + i;

      if ((n=cleanobjlist->nobj) >= CLEAN_STACKSIZE)
	{
	  objstruct	*cleanobj;
	  int		ymin, victim=0;

	  ymin = 2000000000;	/* No image is expected to be that tall ! */
	  cleanobj = cleanobjlist->obj;
	  for (j=0; j<n; j++, cleanobj++)
	    if (cleanobj->ycmax < ymin)
	      {
		victim = j;
		ymin = cleanobj->ycmax;
	      }
	  
	  /* removed (defined in analyse.c in sextractor) */
	  /* endobject(field, dfield, wfield, dwfield,victim,cleanobjlist); */

	  /* TODO don't think I should be removing this here!! */
	  subcleanobj(victim, cleanobjlist);
	}

      /* Only add the object if it is not swallowed by cleaning */
      if (!clean_flag || clean(i, objlist2, cleanobjlist, cleanvictim,
			       clean_param, &status))
	  status = addcleanobj(cobj, cleanobjlist);
      if (status != RETURN_OK)
	    goto exit_sortit;
    }

 exit_sortit:
  free(objlistout.plist);
  free(objlistout.obj);
  return status;
}


/******************************** convolve ***********************************/
/* (originally in filter.c in sextractor)
Convolve a scan line with an array.
*/
void	convolve(PIXTYPE *im,                    /* full image (was field) */
		 int w, int h,                   /* image size */
		 int y,                          /* line in image */
		 float *conv,                    /* convolution mask */
		 int convw, int convh,           /* mask size */
		 PIXTYPE *mscan)                 /* convolved line */
{
  int		convw2,m0,me,m,mx,dmx, y0, dy;
  float	        *mask;
  PIXTYPE	*mscane, *s,*s0, *d,*de, mval;

  convw2 = convw/2;
  mscane = mscan+w; /* limit of scanline */
  y0 = y - (convh/2); /* starting y for convolution */

  /* check if start extends beyond image */
  if (y0 < 0)
    {
      m0 = convw*(-y0);
      y0 = 0;
    }
  else
    m0 = 0;

  if ((dy = h - y0) < convh)
    me = convw*dy;
  else
    me = convw*convh;

  memset(mscan, 0, w*sizeof(PIXTYPE));
  s0 = NULL;				/* To avoid gcc -Wall warnings */
  mask = conv+m0;
  for (m = m0, mx = 0; m<me; m++, mx++) /* loop over pixels in conv mask */
                                        /* mx is x position in mask */
    {
      if (mx==convw) 
	mx = 0;
      if (!mx)
	s0 = im + w*((y0++)%h);  /* every time mx goes to 0, increment */
				 /* start line in the image */

      if ((dmx = mx-convw2)>=0)  /* dmx is x-offset in mask */
	{
	  s = s0 + dmx;
	  d = mscan;
	  de = mscane - dmx;
	}
      else
	{
	  s = s0;
	  d = mscan - dmx;
	  de = mscane;
	}

      mval = *(mask++);
      while (d<de)
	*(d++) += mval**(s++);
    }

  return;
}


/******************************** createsubmap *******************************
PURPOSE Create pixel-index submap for deblending.
OUTPUT  RETURN_OK if success, RETURN_ERROR otherwise (memory overflow).
*/
int createsubmap(objliststruct *objlist, int no)
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
    return RETURN_ERROR;
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

/****************************** plistinit ************************************
 * (originally init_plist() in sextractor)
PURPOSE	initialize a pixel-list and its components.
 ***/
void plistinit(PIXTYPE *conv, PIXTYPE *cdwfield)
{
  pbliststruct	*pbdum = NULL;

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
