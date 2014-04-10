#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	"sep.h"

/****************************** scanimage ************************************

 ***/

/* sextractor defaults
   -------------------
   threshabsolute 0 (No; relative threholding)
   dthresh 1.5
   athresh 1.5
   minarea 5
   deblend_nthresh 32
   deblend_mincont 0.005
   clean_flag      1 (Yes)
   clean_param     1.0
   conv = [1 2 1; 2 4 2; 1 2 1]
*/
int scanimage(PIXTYPE *cfield, PIXTYPE *cdwfield, int w, int h,
	      PIXTYPE dthresh, PIXTYPE athresh, PIXTYPE cdwthresh,
	      int threshabsolute, int minarea,
	      float *conv, int convw, int convh,
	      int deblend_nthresh, double deblend_mincont,
	      int clean_flag, double clean_param)
{
  static infostruct	curpixinfo, *info, *store, initinfo, freeinfo, *victim;
  objliststruct       	objlist, *cleanobjlist;
  objstruct		*cleanobj;
  pliststruct		*pixel, *pixt; 
  char			*marker, newmarker, *blankpad, *bpt,*bpt0;
  int			co, i,j, flag, luflag,pstop, xl,xl2,yl, cn,
			nposize, stacksize, w, h, blankh, maxpixnb,
                        varthreshflag, convn;
  short	       	        trunflag;
  PIXTYPE		thresh, relthresh, cdnewsymbol, cdwthresh,wthresh,
			*scan,*dscan,*cdscan,*dwscan,*dwscanp,*dwscann,
			*cdwscan,*cdwscanp,*cdwscann,*wscand,
                        *scant, *wscan,*wscann,*wscanp;
  float                 sum, *convnorm;
  FLAGTYPE		*pfscan[MAXFLAG];
  status		cs, ps, *psstack;
  int			*start, *end, ymax;
  LONG                  *cleanvictim;
  
int yblank, stripy, y, ymin, stripylim, stripsclim;  /* from picstruct */
  int wyblank, wstripy, wy, wymin, wstripylim, wstripsclim;

  /*----- Beginning of the main loop: Initialisations  */
  convn = 0;
  sum = 0.0;
  convnorm = NULL;
  cdscan = cdwscan = NULL;              /* Avoid gcc -Wall warnings */
  victim = NULL;			/* Avoid gcc -Wall warnings */
  cleanvictim = NULL;
  blankh = 0;				/* Avoid gcc -Wall warnings */

  /* cdwfield is the detection weight-field if available (can be null)*/
  if (cdwthresh>BIG*WTHRESH_CONVFAC)
    cdwthresh = BIG*WTHRESH_CONVFAC;
  wthresh = 0.0;

  /* If WEIGHTing and no absolute thresholding, activate threshold scaling */
  varthreshflag = (cdwfield && threshabsolute==0);
  relthresh = varthreshflag ? thresh : 0.0;/* To avoid gcc warnings*/
  objlist.dthresh = dthresh;
  objlist.thresh = athresh;
  yblank = 1;
  y = stripy = 0;
  ymin = stripylim = 0;
  stripysclim = 0;

  /* init regardless of whether weightmap is NULL; avoid gcc -Wall warnings */
  wy = wstripy = 0;
  wymin = wstripylim = 0;
  wstripysclim = 0;

  /*Allocate memory for buffers */
  stacksize = w+1;
  QMALLOC(info, infostruct, stacksize);
  QCALLOC(store, infostruct, stacksize);
  QMALLOC(marker, char, stacksize);
  QMALLOC(dumscan, PIXTYPE, stacksize);
  QMALLOC(psstack, status, stacksize);
  QCALLOC(start, int, stacksize);
  QMALLOC(end, int, stacksize);
  blankpad = bpt = NULL;
  lutzalloc(w,h);
  allocparcelout(deblend_nthresh);

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
    QMALLOC(cleanvictim, LONG, CLEAN_STACKSIZE);
  QMALLOC(cleanobjlist, objliststruct, 1);
  cleanobjlist->obj = NULL;
  cleanobjlist->plist = NULL;
  cleanobjlist->nobj = cleanobjlist->npix = 0;


  /* Allocate memory for the pixel list */
  init_plist(conv, cdwfield);
  if (!(pixel = objlist.plist = malloc(nposize=MEMORY_PIXSTACK*plistsize)))
    return MEMORY_PIXSTACK_ERROR;

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
      QMALLOC(cdscan, PIXTYPE, stacksize);
      if (cdwfield)
	  QCALLOC(cdwscan, PIXTYPE, stacksize);

      /* normalize the filter */
      convn = convw * convh;
      QMALLOC(convnorm, PIXTYPE, convn);
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
	      if (cdwfield)
		free(cdwscan);
	    }
	  cdwscan = cdscan = dumscan;
	}

      else
	{
	  scan = cfield + stripy*w;
	  if (cdwfield)
	    wscan = cdwfield + stripy*w;

	  /* no separate detection image */
	  dscan = scan;
	  dwscan = wscan;

	  /* filter the lines */
	  if (conv)
	    {
	      convolve(cfield, w, h, stripy, conv, convw, convh, cdscan);
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
		  sprintf(gstr, "%d,%d", xl+1, yl+1);
		  warning("Pixel stack overflow at position ", gstr);
		  
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
		      warning("*Fatal Error*: something is badly bugged in ",
			      "scanimage()!");
		      return FATAL_ERROR;
		    }
		  if (maxpixnb <= 1)
		    {
		      warning("Pixel stack overflow in ", "scanimage()");
		      return PIXSTACK_OVERFLOW_ERROR;
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
			    sortit(&info[co], &objlist,
				   cdwscan, wscan, minarea,
				   clean_flag, clean_param, cleanobjlist,
				   cleanvictim);

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
      stripy = y = yl;
      if (cdwfield)
	wstripy = wy = yl;

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

  if (clean_flag)
    free(cleanvictim);
  free(cleanobjlist); /* TODO don't free this! return it */

  /* TODO free mem if early exit (goto here)!!!! */

  /*Free memory */
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
  return;

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

  return;
}

/********************************* sortit ************************************/
/*
build the object structure.
*/
void sortit(infostruct *info, objliststruct *objlist,
	    PIXTYPE *cdwscan, PIXTYPE *wscan, int minarea,
	    int clean_flag, double clean_param, objliststruct *cleanobjlist, LONG *cleanvictim)
{
  objliststruct	        objlistout, *objlist2;
  static objstruct	obj;
  objstruct		*cobj;
  pliststruct		*pixel;
  int 			i,j,n;

  pixel = objlist->plist;
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
	  sprintf(gstr, "%.0f,%.0f", obj.mx+1, obj.my+1);
	  warning("Deblending overflow for detection at ", gstr);
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
	  int		ymin, ymax, victim=0;

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
      if (!clean_flag ||
	  clean(i, objlist2, cleanobjlist, cleanvictim, clean_param))
	addcleanobj(cobj, cleanobjlist);
    }

  free(objlistout.plist);
  free(objlistout.obj);
  
  return;
}


/******************************** preanalyse *********************************
PROTO   void preanalyse(int no, objliststruct *objlist, int analyse_type)
PURPOSE Compute basic image parameters from the pixel-list for each detection.
INPUT   objlist number,
        objlist pointer,
        analysis switch flag.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP & Leiden & ESO)
VERSION 28/11/2003
 ***/
void  preanalyse(int no, objliststruct *objlist, int analyse_type)
{
  objstruct	*obj = &objlist->obj[no];
  pliststruct	*pixel = objlist->plist, *pixt;
  PIXTYPE	peak, cpeak, val, cval, minthresh, thresht;
  double	thresh,thresh2, t1t2,darea,
                mx,my, mx2,my2,mxy, rv, tv,
		xm,ym, xm2,ym2,xym,
		temp,temp2, theta,pmx2,pmy2;
  int		x, y, xmin,xmax, ymin,ymax,area2, fdnpix, dnpix;
  
  /*-----  initialize stacks and bounds */
  thresh = obj->dthresh;
  minthresh = 0.0;
  fdnpix = dnpix = 0;
  rv = 0.0;
  peak = cpeak = -BIG;
  ymin = xmin = 2*MAXPICSIZE;    /* to be really sure!! */
  ymax = xmax = 0;

  /*-----  integrate results */
  for (pixt=pixel+obj->firstpix; pixt>=pixel; pixt=pixel+PLIST(pixt,nextpix))
    {
      x = PLIST(pixt, x);
      y = PLIST(pixt, y);
      val=PLISTPIX(pixt, dvalue);
      if (cpeak < (cval=PLISTPIX(pixt, cdvalue)))
	cpeak = cval;
      if (peak < val)
	peak = val;
      rv += cval;
      if (xmin > x)
	xmin = x;
      if (xmax < x)
	xmax = x;
      if (ymin > y)
	ymin = y;
      if (ymax < y)
	ymax = y;
      fdnpix++;
    }    
  
  /* copy some data to "obj" structure */
  obj->fdnpix = (LONG)fdnpix;
  obj->fdflux = (float)rv;
  obj->fdpeak = cpeak;
  obj->dpeak = peak;
  obj->xmin = xmin;
  obj->xmax = xmax;
  obj->ymin = ymin;
  obj->ymax = ymax;

  if (analyse_type & ANALYSE_FULL)
    {
      mx = my = tv = 0.0;
      mx2 = my2 = mxy = 0.0;
      thresh2 = (thresh + peak)/2.0;
      area2 = 0;
      for (pixt=pixel+obj->firstpix; pixt>=pixel;
	   pixt=pixel+PLIST(pixt,nextpix))
	{
	  x = PLIST(pixt,x)-xmin;  /* avoid roundoff errors on big images */
	  y = PLIST(pixt,y)-ymin;  /* avoid roundoff errors on big images */
	  cval = PLISTPIX(pixt, cdvalue);
	  tv += (val = PLISTPIX(pixt, dvalue));
	  if (val>thresh)
	    dnpix++;
	  if (val > thresh2)
	    area2++;
	  mx += cval * x;
	  my += cval * y;
	  mx2 += cval * x*x;
	  my2 += cval * y*y;
	  mxy += cval * x*y;
	}

      /* compute object's properties */
      xm = mx / rv;    /* mean x */
      ym = my / rv;    /* mean y */

      /* In case of blending, use previous barycenters */
      if ((analyse_type & ANALYSE_ROBUST) && (obj->flag & OBJ_MERGED))
	{
	  double	xn,yn;

	  xn = obj->mx-xmin;
	  yn = obj->my-ymin;
	  xm2 = mx2 / rv + xn*xn - 2*xm*xn;
	  ym2 = my2 / rv + yn*yn - 2*ym*yn;
	  xym = mxy / rv + xn*yn - xm*yn - xn*ym;
	  xm = xn;
	  ym = yn;
	}
      else
	{
	  xm2 = mx2 / rv - xm * xm;	 /* variance of x */
	  ym2 = my2 / rv - ym * ym;	 /* variance of y */
	  xym = mxy / rv - xm * ym;	 /* covariance */
	}

      /* Handle fully correlated x/y (which cause a singularity...) */
      if ((temp2=xm2*ym2-xym*xym)<0.00694)
	{
	  xm2 += 0.0833333;
	  ym2 += 0.0833333;
	  temp2 = xm2*ym2-xym*xym;
	  obj->singuflag = 1;
	}
      else
	obj->singuflag = 0;

      if ((fabs(temp=xm2-ym2)) > 0.0)
	theta = atan2(2.0 * xym,temp) / 2.0;
      else
	theta = PI/4.0;

    temp = sqrt(0.25*temp*temp+xym*xym);
    pmy2 = pmx2 = 0.5*(xm2+ym2);
    pmx2+=temp;
    pmy2-=temp;

    obj->dnpix = (obj->flag & OBJ_OVERFLOW)? obj->fdnpix:(LONG)dnpix;
    obj->dflux = tv;
    obj->mx = xm+xmin;	/* add back xmin */
    obj->my = ym+ymin;	/* add back ymin */
    obj->mx2 = xm2;
    obj->my2 = ym2;
    obj->mxy = xym;
    obj->a = (float)sqrt(pmx2);
    obj->b = (float)sqrt(pmy2);
    obj->theta = theta*180.0/PI;

    obj->cxx = (float)(ym2/temp2);
    obj->cyy = (float)(xm2/temp2);
    obj->cxy = (float)(-2*xym/temp2);

    darea = (double)area2 - dnpix;
    t1t2 = thresh/thresh2;
    if (t1t2>0.0 && !plistexist_dthresh)  /* was: prefs.dweight_flag */
      {
	obj->abcor = (darea<0.0?darea:-1.0)/(2*PI*log(t1t2<1.0?t1t2:0.99)
					     *obj->a*obj->b);
	if (obj->abcor>1.0)
	  obj->abcor = 1.0;
      }
    else
      obj->abcor = 1.0;
    }
  
  return;
}
