#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	"sep.h"

/****************************** scanimage ************************************
PROTO   void scanimage(picstruct *field, picstruct *dfield, picstruct *ffield,
        picstruct *wfield, picstruct *dwfield)
PURPOSE Scan of the large pixmap(s). Main loop and heart of the program.
INPUT   Measurement field pointer,
        Detection field pointer,
        Flag field pointer,
        Measurement weight-map field pointer,
        Detection weight-map field pointer,
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 21/12/2011
 ***/
int	scanimage(PIXTYPE *cfield, PIXTYPE *cdwfield, int w, int h,
		  PIXTYPE dthresh, PIXTYPE athresh, PIXTYPE cdwthresh,
		  int threshabsolute, 
		  float *conv, int convw, int convh)
{
  static infostruct	curpixinfo, *info, *store, initinfo, freeinfo, *victim;
  checkstruct		*check;
  objliststruct       	objlist;
  objstruct		*cleanobj;
  pliststruct		*pixel, *pixt; 
  char			*marker, newmarker, *blankpad, *bpt,*bpt0;
  int			co, i,j, flag, luflag,pstop, xl,xl2,yl, cn,
			nposize, stacksize, w, h, blankh, maxpixnb,
                        varthreshflag, ontotal, convn;
  short	       	        trunflag;
  PIXTYPE		thresh, relthresh, cdnewsymbol, cdwthresh,wthresh,
			*scan,*dscan,*cdscan,*dwscan,*dwscanp,*dwscann,
			*cdwscan,*cdwscanp,*cdwscann,*wscand,
                        *scant, *wscan,*wscann,*wscanp;
  float                 sum, *convnorm;
  FLAGTYPE		*pfscan[MAXFLAG];
  status		cs, ps, *psstack;
  int			*start, *end, ymax;

  sexcatstruct         thecat; /* from globals.h */
  int yblank, stripy, y, ymin, stripylim, stripsclim;  /* from picstruct */
  int wyblank, wstripy, wy, wymin, wstripylim, wstripsclim;

  convn = 0;
  sum = 0.0;
  convnorm = NULL;
  cdscan = cdwscan = NULL;              /* Avoid gcc -Wall warnings */
  victim = NULL;			/* Avoid gcc -Wall warnings */
  blankh = 0;				/* Avoid gcc -Wall warnings */

  /*----- Beginning of the main loop: Initialisations  */
  thecat.ntotal = thecat.ndetect = 0;

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
  allocparcelout();

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
  initclean();

  /* Allocate memory for the pixel list */
  init_plist();
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
	  scan = cdfield + stripy*w;
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



/* remember to free convnorm if conv */
