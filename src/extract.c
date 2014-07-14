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

/* Note: was scan.c in SExtractor. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sep.h"
#include "sepcore.h"
#include "extract.h"

#define DETECT_MAXAREA 0        /* replaces prefs.ext_maxarea */
#define MEMORY_PIXSTACK 300000  /* number of pixels in stack */
                                /* (replaces prefs.mem_pixstack) */
#define	WTHRESH_CONVFAC	1e-4    /* Factor to apply to weights when */
			        /* thresholding filtered weight-maps */

/* globals */
int plistexist_cdvalue, plistexist_dthresh, plistexist_var;
int plistoff_value, plistoff_cdvalue, plistoff_dthresh, plistoff_var;
int plistsize;

typedef void (*convolver)(void *image, int w, int h, int y,
			  float *conv, int convw, int convh, PIXTYPE *buf);
int get_convolver(int dtype, convolver *f);

int  sortit(infostruct *, objliststruct *, int,
	    objliststruct *, int, double);
void plistinit(void *, void *);
void clean(objliststruct *objlist, double clean_param, int *survives);
int convertobj(int l, objliststruct *objlist, sepobj *objout, int w);


/****************************** extract **************************************/
int sep_extract(void *image, void *noise,
		int dtype, int ndtype, short noise_flag, int w, int h,
	        float thresh, int minarea, float *conv, int convw, int convh,
		int deblend_nthresh, double deblend_cont,
		int clean_flag, double clean_param,
		sepobj **objects, int *nobj)
{
  static infostruct curpixinfo, initinfo, freeinfo;
  objliststruct     objlist;
  char              newmarker;
  int               co, i, j, flag, luflag, pstop, xl, xl2, yl, cn;
  int               nposize, stacksize, maxpixnb, convn, status;
  int               elsize_im, elsize_noise;
  short             trunflag;
  PIXTYPE           relthresh, cdnewsymbol;
  float             sum;
  pixstatus         cs, ps;

  static infostruct *info, *store, *victim;
  objliststruct     *finalobjlist;
  pliststruct	    *pixel, *pixt;
  char              *marker;
  PIXTYPE           *scan, *cdscan, *cdwscan, *wscan, *dumscan;
  float             *convnorm;
  int               *start, *end, *survives;
  pixstatus         *psstack;
  BYTE              *imageline, *noiseline;
  convolver         convolve_im, convolve_noise;
  array_converter   convert_im, convert_noise;
  char              errtext[80];

  status = RETURN_OK;
  pixel = NULL;
  convnorm = NULL;
  scan = wscan = cdscan = cdwscan = dumscan = NULL;
  victim = NULL;
  info = NULL;
  store = NULL;
  marker = NULL;
  psstack = NULL;
  start = end = NULL;
  finalobjlist = NULL; /* final return value */
  convn = 0;
  sum = 0.0;
  convolve_im = NULL;
  convolve_noise = NULL;
  convert_im = NULL;
  convert_noise = NULL;
  imageline = (BYTE *)image;
  noiseline = (BYTE *)noise;

  /* If we have a noise array, set relative threshold */
  relthresh = noise? thresh : 0.0; /* To avoid gcc warnings*/

  objlist.dthresh = thresh;
  objlist.thresh = thresh;

  /*Allocate memory for buffers */
  stacksize = w+1;
  QMALLOC(info, infostruct, stacksize, status);
  QCALLOC(store, infostruct, stacksize, status);
  QMALLOC(marker, char, stacksize, status);
  QMALLOC(dumscan, PIXTYPE, stacksize, status);
  QMALLOC(psstack, pixstatus, stacksize, status);
  QCALLOC(start, int, stacksize, status);
  QMALLOC(end, int, stacksize, status);
  if ((status = lutzalloc(w, h)) != RETURN_OK)
    goto exit;
  if ((status = allocdeblend(deblend_nthresh)) != RETURN_OK)
    goto exit;

  /* allocate scan buffer(s) and get array converter function(s) */
  QMALLOC(scan, PIXTYPE, stacksize, status);
  status = get_array_converter(dtype, &convert_im, &elsize_im);
  if (status != RETURN_OK)
    goto exit;
  if (noise)
    {
      QMALLOC(wscan, PIXTYPE, stacksize, status);
      status = get_array_converter(ndtype, &convert_noise, &elsize_noise);
      if (status != RETURN_OK)
	goto exit;
    }

  /* More initializations */
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

  /* Init finalobjlist (the return catalog) */
  QMALLOC(finalobjlist, objliststruct, 1, status);
  finalobjlist->obj = NULL;
  finalobjlist->plist = NULL;
  finalobjlist->nobj = finalobjlist->npix = 0;


  /* Allocate memory for the pixel list */
  plistinit(conv, noise);
  if (!(pixel = objlist.plist = malloc(nposize=MEMORY_PIXSTACK*plistsize)))
    {
      status = MEMORY_ALLOC_ERROR;
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
      QMALLOC(cdscan, PIXTYPE, stacksize, status);
      if (noise)
	QCALLOC(cdwscan, PIXTYPE, stacksize, status);

      /* normalize the filter */
      convn = convw * convh;
      QMALLOC(convnorm, PIXTYPE, convn, status);
      for (i=0; i<convn; i++)
	sum += fabs(conv[i]);
      for (i=0; i<convn; i++)
	convnorm[i] = conv[i] / sum;

      /* get the right convolve function for the image & noise data types */
      status = get_convolver(dtype, &convolve_im);
      if (status != RETURN_OK)
	goto exit;
      if (noise)
	{
	  status = get_convolver(ndtype, &convolve_noise);
	  if (status != RETURN_OK)
	    goto exit;
	}

    }

  /*----- MAIN LOOP ------ */
  for (yl=0; yl<=h; yl++)
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
	      if (noise)
		{
		  free(cdwscan);
		  cdwscan = NULL;
		}
	    }
	  cdwscan = cdscan = dumscan;
	}

      else
	{
	  /* read the current (yl) line of input arrays into PIXTYPE buffers*/
	  convert_im(imageline, w, scan);
	  if (noise)
	    convert_noise(noiseline, w, wscan);

	  /* filter the lines */
	  if (conv)
	    {
	      convolve_im(image, w, h, yl, convnorm, convw, convh, cdscan);
	      if (noise)
		convolve_noise(noise, w, h, yl, convnorm, convw, convh,
			       cdwscan);
	    }
	  else
	    {
	      cdscan = scan;
	      cdwscan = wscan;
	    }	  
	}
      
      trunflag = (yl==0 || yl==h-1)? SEP_OBJ_TRUNC:0;
      
      for (xl=0; xl<=w; xl++)
	{
	  if (xl == w)
	    cdnewsymbol = -BIG;
	  else
	    cdnewsymbol = cdscan[xl];

	  newmarker = marker[xl];  /* marker at this pixel */
	  marker[xl] = 0;

	  curpixinfo.flag = trunflag;
	  if (noise)
	    thresh = relthresh * ((xl==w || yl==h)? 0.0: cdwscan[xl]);
	  luflag = cdnewsymbol > thresh? 1: 0;  /* is pixel above thresh? */

	  if (luflag)
	    {
	      /* flag the current object if we're near the image bounds */
	      if (xl==0 || xl==w-1)
		curpixinfo.flag |= SEP_OBJ_TRUNC;
	      
	      /* point pixt to first free pixel in pixel list */
	      /* and increment the "first free pixel" */
	      pixt = pixel + (cn=freeinfo.firstpix);
	      freeinfo.firstpix = PLIST(pixt, nextpix);
	      curpixinfo.lastpix = curpixinfo.firstpix = cn;

	      /* set values for the new pixel */ 
	      PLIST(pixt, nextpix) = -1;
	      PLIST(pixt, x) = xl;
	      PLIST(pixt, y) = yl;
	      PLIST(pixt, value) = scan[xl];
	      if (PLISTEXIST(cdvalue))
		PLISTPIX(pixt, cdvalue) = cdnewsymbol;
	      if (PLISTEXIST(var))
		PLISTPIX(pixt, var) = wscan[xl];

	      /* Check if we are running out of free pixels in objlist.plist */
	      /* (previously, the largest object became a "victim") */
	      if (freeinfo.firstpix==freeinfo.lastpix)
		{
		  status = SEP_INTERNAL_ERROR;
		  sprintf(errtext, "Pixel stack overflow at position %d,%d.",
			  xl+1, yl+1);
		  put_errdetail(errtext);
		  goto exit;
		  
		  /* NOTE: The above error was originally just a warning.
		     with the change to an error, the following code in this
		     if block is never executed.
		     TODO: should this just be a warning (or nothing?)
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
		  
		  if ((!maxpixnb) || (maxpixnb <= 1))
		    {
		      status = SEP_INTERNAL_ERROR;
		      goto exit;
		    }
		  freeinfo.firstpix = PLIST(pixel+victim->firstpix, nextpix);
		  PLIST(pixel+victim->lastpix, nextpix) = freeinfo.lastpix;
		  PLIST(pixel+(victim->lastpix=victim->firstpix), nextpix) = -1;
		  victim->pixnb = 1;
		  victim->flag |= SEP_OBJ_OVERFLOW;
		}
	      /*------------------------------------------------------------*/

	      /* if the current status on this line is not already OBJECT... */
	      /* start segment */
	      if (cs != OBJECT)
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
	  /* newmarker is marker[ ] at this pixel position before we got to
	     it. We'll only enter this if marker[ ] was set on a previous
	     loop iteration.   */
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
			      status = sortit(&info[co], &objlist, minarea,
					      finalobjlist,
					      deblend_nthresh,deblend_cont);
			      if (status != RETURN_OK)
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

	  /* update the info or end segment */
	  if (luflag)
	    {
	      update(&info[co], &curpixinfo, pixel);
	    }
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

      /* increment array pointers */ 
      imageline += w*elsize_im;
      if (noise)
	noiseline += w*elsize_noise;

    } /*---------------- End of the loop over the y's -----------------------*/

  /* convert `finalobjlist` to an array of `sepobj` structs */
  if (clean_flag)
    {
      /* Calculate mthresh for all objects in the list (needed for cleaning) */
      for (i=0; i<finalobjlist->nobj; i++)
	{
	  status = analysemthresh(i, finalobjlist, minarea, thresh);
	  if (status != RETURN_OK)
	    goto exit;
	}

      QMALLOC(survives, int, finalobjlist->nobj, status);
      clean(finalobjlist, clean_param, survives);

      /* count surviving objects and allocate space accordingly*/
      *nobj = 0;
      for (i=0; i<finalobjlist->nobj; i++)
	*nobj += survives[i];
      QMALLOC(*objects, sepobj, *nobj, status);

      /* fill */
      j=0;
      for (i=0; i<finalobjlist->nobj; i++)
	if (survives[i])
	    convertobj(i, finalobjlist, (*objects) + j++, w);
    }
  else
    {
      *nobj = finalobjlist->nobj;
      QMALLOC(*objects, sepobj, *nobj, status);
      for (i=0; i<finalobjlist->nobj; i++)
	convertobj(i, finalobjlist, (*objects) + i, w);
    }

 exit:
  free(finalobjlist->obj);
  free(finalobjlist->plist);
  free(finalobjlist);
  freedeblend();
  free(pixel);
  lutzfree();
  free(info);
  free(store);
  free(marker);
  free(dumscan);
  free(psstack);
  free(start);
  free(end);
  free(scan);
  free(wscan);
  if (conv)
    free(convnorm);

  if (status != RETURN_OK)
    {
      free(cdscan);   /* only need to free these in case of early exit */
      free(cdwscan);
      *objects = NULL;
      *nobj = 0;
    }

  return status;
}


/********************************* sortit ************************************/
/*
build the object structure.
*/
int sortit(infostruct *info, objliststruct *objlist, int minarea,
	   objliststruct *finalobjlist,
	   int deblend_nthresh, double deblend_mincont)
{
  objliststruct	        objlistout, *objlist2;
  static objstruct	obj;
  int 			i, status;
  PIXTYPE thresh;

  status=RETURN_OK;  
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
  obj.dthresh = thresh = objlist->dthresh;
  obj.thresh = objlist->thresh;

  preanalyse(0, objlist);

  /*----- Check if the current strip contains the lower isophote
    (it always should since the "current strip" is the entire image!) */
  if ((int)obj.ymin < 0)
    obj.flag |= SEP_OBJ_ISO_PB;

  if (!(obj.flag & SEP_OBJ_OVERFLOW))
    {
      status = deblend(objlist, 0, &objlistout, deblend_nthresh,
		       deblend_mincont, minarea);
      if (status == RETURN_OK)
	{
	  objlist2 = &objlistout;
	}
      else
	{
	  /* formerly, this wasn't a fatal error, so a flag was set for
	     the object and we continued. I'm leaving the flag-setting here
	     in case we want to change this to a non-fatal error in the
	     future, but currently the flag setting is irrelevant. */
	  objlist2 = objlist;
	  for (i=0; i<objlist2->nobj; i++)
	    objlist2->obj[i].flag |= SEP_OBJ_DOVERFLOW;
	  goto exit;
	}
    }
  else
    objlist2 = objlist;
  
  /* Analyze the deblended objects and add to the final list */
  for (i=0; i<objlist2->nobj; i++)
    {
      analyse(i, objlist2, 1);

      /* this does nothing if DETECT_MAXAREA is 0 (and it currently is) */
      if (DETECT_MAXAREA && objlist2->obj[i].fdnpix > DETECT_MAXAREA)
	continue;

      /* add the object to the final list */
      status = addobjdeep(i, objlist2, finalobjlist);
      if (status != RETURN_OK)
	goto exit;
    }

 exit:
  free(objlistout.plist);
  free(objlistout.obj);
  return status;
}


/********** addobjdeep (originally in manobjlist.c) **************************/
/*
Add object number `objnb` from list `objl1` to list `objl2`.
Unlike `addobjshallow` this also copies plist pixels to the second list.
*/

int addobjdeep(int objnb, objliststruct *objl1, objliststruct *objl2)
{
  objstruct	*objl2obj;
  pliststruct	*plist1 = objl1->plist, *plist2 = objl2->plist;
  int		fp, i, j, npx, objnb2;
  
  fp = objl2->npix;      /* 2nd list's plist size in pixels */
  j = fp*plistsize;      /* 2nd list's plist size in bytes */
  objnb2 = objl2->nobj;  /* # of objects currently in 2nd list*/

  /* Allocate space in `objl2` for the new object */
  if (objnb2)
    objl2obj = (objstruct *)realloc(objl2->obj,
				    (++objl2->nobj)*sizeof(objstruct));
  else
    objl2obj = (objstruct *)malloc((++objl2->nobj)*sizeof(objstruct));

  if (!objl2obj)
    goto earlyexit;
  objl2->obj = objl2obj;

  /* Allocate space for the new object's pixels in 2nd list's plist */
  npx = objl1->obj[objnb].fdnpix;
  if (fp)
    plist2 = (pliststruct *)realloc(plist2, (objl2->npix+=npx)*plistsize);
  else
    plist2 = (pliststruct *)malloc((objl2->npix=npx)*plistsize);

  if (!plist2)
    goto earlyexit;
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
 earlyexit:
  objl2->nobj--;
  objl2->npix = fp;
  return MEMORY_ALLOC_ERROR;
}



/******************************** convolve ***********************************/
/* Convolve functions.
 * 
 * Convolve an image with a convolution kernel (at a single line).
 *
 * There is a separate function for each type of image.
 * In sep_extract(), we select the correct function for the given array type
 * at runtime.
 *
 * (originally in filter.c in sextractor)
 *
 * image : full input array
 * w, h : width and height of image
 * y : line to convolve in image
 * conv : convolution kernel
 * convw, convh : width and height of conv
 * buf : output convolved line (`w` elements long)
 */
void convolve_flt(void *image, int w, int h, int y,
		  float *conv, int convw, int convh, PIXTYPE *buf)
{
  int	  convw2, cx, dcx, y0, dy;
  float	  *convend;
  PIXTYPE *bufend, *buft, *buftend, cval;

  float *im, *line, *imt;
  im = (float *)image;

  line = NULL;        /* To avoid gcc -Wall warnings */
  convw2 = convw/2;
  bufend = buf+w;     /* pointer to end of output buffer */
  y0 = y - (convh/2); /* starting y line for convolution */

  /* reduce height of convolution kernel if it extends beyond image */
  dy = h - y0;        /* distance from starting line to top of image */
  if (convh > dy)
    convh = dy;
  convend = conv + (convw*convh);

  /* Set start position in convolution kernel (and start line in image) */
  if (y0 < 0)
    {
      conv += convw*(-y0);
      y0 = 0;
    }

  /* initialize output buffer to zero */
  memset(buf, 0, w*sizeof(PIXTYPE));

  /* loop over pixels in the mask */
  for (cx=0; conv<convend; conv++, cx++)
    {
      cval = *conv;

      /* get the x position in the mask */
      if (cx==convw)
	cx = 0;

      /* when cx goes to zero, increment the start line in the image */
      if (!cx)
	line = im + w*((y0++)%h);

      /* get start and end positions in the source and target line */
      if ((dcx = cx-convw2)>=0)
	{
	  imt = line + dcx;
	  buft = buf;
	  buftend = bufend - dcx;
	}
      else
	{
	  imt = line;
	  buft = buf - dcx;
	  buftend = bufend;
	}

      while (buft < buftend)
	*(buft++) += cval * *(imt++);
    }

  return;
}

void convolve_dbl(void *image, int w, int h, int y,
		  float *conv, int convw, int convh, PIXTYPE *buf)
{
  int	  convw2, cx, dcx, y0, dy;
  float	  *convend;
  PIXTYPE *bufend, *buft, *buftend, cval;

  double *im, *line, *imt;
  im = (double *)image;

  line = NULL;        /* To avoid gcc -Wall warnings */
  convw2 = convw/2;
  bufend = buf+w;     /* pointer to end of output buffer */
  y0 = y - (convh/2); /* starting y line for convolution */

  /* reduce height of convolution kernel if it extends beyond image */
  dy = h - y0;        /* distance from starting line to top of image */
  if (convh > dy)
    convh = dy;
  convend = conv + (convw*convh);

  /* Set start position in convolution kernel (and start line in image) */
  if (y0 < 0)
    {
      conv += convw*(-y0);
      y0 = 0;
    }

  /* initialize output buffer to zero */
  memset(buf, 0, w*sizeof(PIXTYPE));

  /* loop over pixels in the mask */
  for (cx=0; conv<convend; conv++, cx++)
    {
      cval = *conv;

      /* get the x position in the mask */
      if (cx==convw)
	cx = 0;

      /* when cx goes to zero, increment the start line in the image */
      if (!cx)
	line = im + w*((y0++)%h);

      /* get start and end positions in the source and target line */
      if ((dcx = cx-convw2)>=0)
	{
	  imt = line + dcx;
	  buft = buf;
	  buftend = bufend - dcx;
	}
      else
	{
	  imt = line;
	  buft = buf - dcx;
	  buftend = bufend;
	}

      while (buft < buftend)
	*(buft++) += cval * *(imt++);
    }

  return;
}

void convolve_int(void *image, int w, int h, int y,
		  float *conv, int convw, int convh, PIXTYPE *buf)
{
  int	  convw2, cx, dcx, y0, dy;
  float	  *convend;
  PIXTYPE *bufend, *buft, *buftend, cval;

  int *im, *line, *imt;
  im = (int *)image;

  line = NULL;        /* To avoid gcc -Wall warnings */
  convw2 = convw/2;
  bufend = buf+w;     /* pointer to end of output buffer */
  y0 = y - (convh/2); /* starting y line for convolution */

  /* reduce height of convolution kernel if it extends beyond image */
  dy = h - y0;        /* distance from starting line to top of image */
  if (convh > dy)
    convh = dy;
  convend = conv + (convw*convh);

  /* Set start position in convolution kernel (and start line in image) */
  if (y0 < 0)
    {
      conv += convw*(-y0);
      y0 = 0;
    }

  /* initialize output buffer to zero */
  memset(buf, 0, w*sizeof(PIXTYPE));

  /* loop over pixels in the mask */
  for (cx=0; conv<convend; conv++, cx++)
    {
      cval = *conv;

      /* get the x position in the mask */
      if (cx==convw)
	cx = 0;

      /* when cx goes to zero, increment the start line in the image */
      if (!cx)
	line = im + w*((y0++)%h);

      /* get start and end positions in the source and target line */
      if ((dcx = cx-convw2)>=0)
	{
	  imt = line + dcx;
	  buft = buf;
	  buftend = bufend - dcx;
	}
      else
	{
	  imt = line;
	  buft = buf - dcx;
	  buftend = bufend;
	}

      while (buft < buftend)
	*(buft++) += cval * *(imt++);
    }

  return;
}

/* return the correct converter depending on the datatype code */
int get_convolver(int dtype, convolver *f)
{
  int status = RETURN_OK;
  char errtext[80];

  if (dtype == SEP_TFLOAT)
    {
      *f = convolve_flt;
    }
  else if (dtype == SEP_TDOUBLE)
    {
      *f = convolve_dbl;
    }
  else if (dtype == SEP_TINT)
    {
      *f = convolve_int;
    }
  else
    {
      *f = NULL;
      status = ILLEGAL_DTYPE;
      sprintf(errtext, "unsupported data type in get_convolver(): %d",
	      dtype);
      put_errdetail(errtext);
    }
  return status;
}

/****************************** plistinit ************************************
 * (originally init_plist() in sextractor)
PURPOSE	initialize a pixel-list and its components.
 ***/
void plistinit(void *conv, void *var)
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
      plistoff_cdvalue = plistoff_value;
    }

  if (var)
    {
      plistexist_var = 1;
      plistoff_var = plistsize;
      plistsize += sizeof(PIXTYPE);
    }
  else
    plistexist_var = 0;

  if (var)
    {
      plistexist_dthresh = 1;
      plistoff_dthresh = plistsize;
      plistsize += sizeof(PIXTYPE);
    }
  else
    plistexist_dthresh = 0;

  return;

}


/************************** clean an objliststruct ***************************/
/*
Fill a list with whether each object in the list survived the cleaning 
(assumes that mthresh has already been calculated for all objects in the list)
*/

void clean(objliststruct *objlist, double clean_param, int *survives)
{
  objstruct     *obj1, *obj2;
  int	        i,j;
  double        amp,ampin,alpha,alphain, unitarea,unitareain,beta,val;
  float	       	dx,dy,rlim;

  beta = clean_param;

  /* initialize to all surviving */
  for (i=0; i<objlist->nobj; i++)
    survives[i] = 1;

  obj1 = objlist->obj;
  for (i=0; i<objlist->nobj; i++, obj1++)
    {
      if (!survives[i])
	continue;

      /* parameters for test object */
      unitareain = PI*obj1->a*obj1->b;
      ampin = obj1->fdflux/(2*unitareain*obj1->abcor);
      alphain = (pow(ampin/obj1->dthresh, 1.0/beta)-1)*unitareain/obj1->fdnpix;

      /* loop over remaining objects in list*/
      obj2 = obj1 + 1;
      for (j=i+1; j<objlist->nobj; j++, obj2++)
	{
	  if (!survives[j])
	    continue;

	  dx = obj1->mx - obj2->mx;
	  dy = obj1->my - obj2->my;
	  rlim = obj1->a + obj2->a;
	  rlim *= rlim;
	  if (dx*dx + dy*dy > rlim*CLEAN_ZONE*CLEAN_ZONE)
	    continue;

	  /* if obj1 is bigger, see if it eats obj2 */
	  if (obj2->fdflux < obj1->fdflux)
	    {
	      val = 1 + alphain*(obj1->cxx*dx*dx + obj1->cyy*dy*dy +
				 obj1->cxy*dx*dy);
	      if (val>1.0 && ((float)(val<1e10?ampin*pow(val,-beta):0.0) >
			      obj2->mthresh))
		  survives[j] = 0; /* the test object eats this one */
	    }

	  /* if obj2 is bigger, see if it eats obj1 */
	  else
	    {
	      unitarea = PI*obj2->a*obj2->b;
	      amp = obj2->fdflux/(2*unitarea*obj2->abcor);
	      alpha = (pow(amp/obj2->dthresh, 1.0/beta) - 1) *
		unitarea/obj2->fdnpix;
	      val = 1 + alpha*(obj2->cxx*dx*dx + obj2->cyy*dy*dy +
			       obj2->cxy*dx*dy);
	      if (val>1.0 && ((float)(val<1e10?amp*pow(val,-beta):0.0) >
			      obj1->mthresh))
		survives[i] = 0;  /* this object eats the test object */
	    }

	} /* inner loop over objlist (obj2) */
    } /* outer loop of objlist (obj1) */
}


/*****************************************************************************/
/*
Convert to an output object.
*/

int convertobj(int l, objliststruct *objlist, sepobj *objout, int w)
{
  int j, status = RETURN_OK;
  objstruct *obj = objlist->obj + l;
  pliststruct *pixt, *pixel;

  pixel = objlist->plist;

  objout->thresh = obj->dthresh;  /* these change names */
  objout->npix = obj->fdnpix;
  objout->tnpix = obj->dnpix;

  objout->xmin = obj->xmin;
  objout->xmax = obj->xmax;
  objout->ymin = obj->ymin;
  objout->ymax = obj->ymax;
  objout->x = obj->mx;
  objout->y = obj->my;
  objout->x2 = obj->mx2;
  objout->y2 = obj->my2;
  objout->xy = obj->mxy;

  objout->a = obj->a;
  objout->b = obj->b;
  objout->theta = obj->theta;

  objout->cxx = obj->cxx;
  objout->cyy = obj->cyy;
  objout->cxy = obj->cxy;

  objout->cflux = obj->fdflux; /* these change names */
  objout->flux = obj->dflux;
  objout->cpeak = obj->fdpeak;
  objout->peak = obj->dpeak;

  objout->xpeak = obj->xpeak;
  objout->ypeak = obj->ypeak;
  objout->xcpeak = obj->xcpeak;
  objout->ycpeak = obj->ycpeak;

  objout->flag = obj->flag;
  if (obj->singuflag)
    objout->flag |= SEP_OBJ_SINGU;
  
  /* Allocate object's pixel list */
  QMALLOC(objout->pix, int, objout->npix, status);

  /* fill it */
  for (pixt=pixel+obj->firstpix, j=0; pixt>=pixel;
       pixt=pixel+PLIST(pixt,nextpix), j++)
      objout->pix[j] = PLIST(pixt,x) + w*PLIST(pixt,y);

 exit:
  return status;
}

void sep_freeobjarray(sepobj *objects, int nobj)
/* free memory associated with an array of sepobj, including pixel lists */
{
  while (nobj > 0)
    free(objects[--nobj].pix);
  free(objects);
}
