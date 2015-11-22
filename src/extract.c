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

#define DETECT_MAXAREA 0             /* replaces prefs.ext_maxarea */
#define	WTHRESH_CONVFAC	1e-4         /* Factor to apply to weights when */
			             /* thresholding filtered weight-maps */

/* globals */
int plistexist_cdvalue, plistexist_thresh, plistexist_var;
int plistoff_value, plistoff_cdvalue, plistoff_thresh, plistoff_var;
int plistsize;
size_t extract_pixstack = 300000;

/* get and set pixstack */
void sep_set_extract_pixstack(size_t val)
{
  extract_pixstack = val;
}

size_t sep_get_extract_pixstack()
{
  return extract_pixstack;
}

int  sortit(infostruct *, objliststruct *, int,
	    objliststruct *, int, double);
void plistinit(void *, void *);
void clean(objliststruct *objlist, double clean_param, int *survives);
int convertobj(int l, objliststruct *objlist, sepobj *objout, int w);

int arraybuffer_init(arraybuffer *buf, void *arr, int dtype, int w, int h,
                     int bufw, int bufh);
void arraybuffer_readline(arraybuffer *buf);
void arraybuffer_free(arraybuffer *buf);

/********************* array buffer functions ********************************/

/* initialize buffer */
/* bufw must be less than or equal to w */
int arraybuffer_init(arraybuffer *buf, void *arr, int dtype, int w, int h,
                     int bufw, int bufh)
{
  int status, yl;
  status = RETURN_OK;

  /* data info */
  buf->dptr = arr;
  buf->dw = w;
  buf->dh = h;

  /* buffer array info */
  buf->bptr = NULL;
  QMALLOC(buf->bptr, PIXTYPE, bufw*bufh, status);
  buf->bw = bufw;
  buf->bh = bufh;

  /* pointers to within buffer */
  buf->midline = buf->bptr + bufw*(bufh/2);  /* ptr to middle buffer line */
  buf->lastline = buf->bptr + bufw*(bufh-1);  /* ptr to last buffer line */

  status = get_array_converter(dtype, &(buf->readline), &(buf->elsize));
  if (status != RETURN_OK)
    goto exit;

  /* initialize yoff */
  buf->yoff = -bufh;

  /* read in lines until the first data line is one line above midline */
  for (yl=0; yl < bufh - bufh/2 - 1; yl++)
    arraybuffer_readline(buf);

  return status;

 exit:
  free(buf->bptr);
  buf->bptr = NULL;
  return status;
}

/* read a line into the buffer at the top, shifting all lines down one */
void arraybuffer_readline(arraybuffer *buf)
{
  PIXTYPE *line;
  int y;

  /* shift all lines down one */
  for (line = buf->bptr; line < buf->lastline; line += buf->bw)
    memcpy(line, line + buf->bw, sizeof(PIXTYPE) * buf->bw);

  /* which image line now corresponds to the last line in buffer? */
  buf->yoff++;
  y = buf->yoff + buf->bh - 1;

  if (y < buf->dh)
    buf->readline(buf->dptr + buf->elsize * buf->dw * y, buf->dw,
                  buf->lastline);

  return;
}

void arraybuffer_free(arraybuffer *buf)
{
  free(buf->bptr);
  buf->bptr = NULL;
}

/* apply_mask_line: Apply the mask to the image and noise buffers.
 *
 * If convolution is off, masked values should simply be not
 * detected. For this, would be sufficient to either set data to zero or
 * set noise (if present) to infinity.
 *
 * If convolution is on, strictly speaking, a masked (unknown) pixel
 * should "poison" the convolved value whenever it is present in the
 * convolution kernel (e.g., NaN behavior). However, practically we'd
 * rather use a "best guess" for the value. Without doing
 * interpolation from neighbors, 0 is the best guess (assuming image
 * is background subtracted).
 *
 * For the purpose of the full matched filter, we should set noise = infinity.
 *
 * So, this routine sets masked pixels to zero in the image buffer and
 * infinity in the noise buffer (if present). It affects the first 
 */
void apply_mask_line(arraybuffer *mbuf, arraybuffer *imbuf, arraybuffer *nbuf)
{
  int i;
  
  for (i=0; i<mbuf->bw; i++)
    {
      if (mbuf->lastline[i] > 0.0)
        {
          imbuf->lastline[i] = 0.0;
          if (nbuf)
            nbuf->lastline[i] = BIG;
        }
    }
}

/****************************** extract **************************************/
int sep_extract(void *image, void *noise, void *mask,
                int dtype, int ndtype, int mdtype, int w, int h,
	        float thresh, int minarea, float *conv, int convw, int convh,
		int filter_type, int deblend_nthresh, double deblend_cont,
		int clean_flag, double clean_param,
		sepobj **objects, int *nobj)
{
  arraybuffer       imbuf, nbuf, mbuf;
  infostruct        curpixinfo, initinfo, freeinfo;
  objliststruct     objlist;
  char              newmarker;
  size_t            mem_pixstack;
  int               nposize, oldnposize;
  int               co, i, j, luflag, pstop, xl, xl2, yl, cn;
  int               stacksize, convn, status;
  int               bufh;
  short             trunflag;
  PIXTYPE           relthresh, cdnewsymbol;
  float             sum;
  pixstatus         cs, ps;

  infostruct        *info, *store;
  objliststruct     *finalobjlist;
  pliststruct	    *pixel, *pixt;
  char              *marker;
  PIXTYPE           *scan, *cdscan, *wscan, *dummyscan;
  PIXTYPE           *sigscan, *workscan;
  float             *convnorm;
  int               *start, *end, *survives;
  pixstatus         *psstack;
  char              errtext[512];

  status = RETURN_OK;
  pixel = NULL;
  convnorm = NULL;
  scan = wscan = cdscan = dummyscan = NULL;
  sigscan = workscan = NULL;
  info = NULL;
  store = NULL;
  marker = NULL;
  psstack = NULL;
  start = end = NULL;
  finalobjlist = NULL; /* final return value */
  convn = 0;
  sum = 0.0;

  mem_pixstack = sep_get_extract_pixstack();

  /* seed the random number generator consistently on each call to get
     consistent results. rand() is used in deblending. */
  srand(1);

  /* If we have a noise array, thresh is actually relthresh; we'll use
   * relthresh below to set the threshold for each pixel. */
  relthresh = thresh;

  objlist.thresh = thresh;

  /*Allocate memory for buffers */
  stacksize = w+1;
  QMALLOC(info, infostruct, stacksize, status);
  QCALLOC(store, infostruct, stacksize, status);
  QMALLOC(marker, char, stacksize, status);
  QMALLOC(dummyscan, PIXTYPE, stacksize, status);
  QMALLOC(psstack, pixstatus, stacksize, status);
  QCALLOC(start, int, stacksize, status);
  QMALLOC(end, int, stacksize, status);
  if ((status = lutzalloc(w, h)) != RETURN_OK)
    goto exit;
  if ((status = allocdeblend(deblend_nthresh)) != RETURN_OK)
    goto exit;

  /* Initialize buffers for input array(s).
   * The buffer size depends on whether or not convolution is active.
   * If not convolving, the buffer size is just a single line. If convolving,
   * the buffer height equals the height of the convolution kernel.
   */
  bufh = conv ? convh : 1;
  status = arraybuffer_init(&imbuf, image, dtype, w, h, stacksize, bufh);
  if (status != RETURN_OK)
    goto exit;
  if (noise)
    {
      status = arraybuffer_init(&nbuf, noise, ndtype, w, h, stacksize, bufh);
      if (status != RETURN_OK)
        goto exit;
    }
  if (mask)
    {
      status = arraybuffer_init(&mbuf, mask, mdtype, w, h, stacksize, bufh);
      if (status != RETURN_OK)
        goto exit;
    }

  /* `scan` (or `wscan`) is always a pointer to the current line being
   * processed. It might be the only line in the buffer, or it might be the 
   * middle line. */
  scan = imbuf.midline;
  if (noise)
    wscan = nbuf.midline;

  /* More initializations */
  initinfo.pixnb = 0;
  initinfo.flag = 0;
  initinfo.firstpix = initinfo.lastpix = -1;

  for (xl=0; xl<stacksize; xl++)
    {
    marker[xl]  = 0 ;
    dummyscan[xl] = -BIG ;
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
  if (!(pixel = objlist.plist = malloc(nposize=mem_pixstack*plistsize)))
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

  /* can only use a matched filter when convolving and when there is a noise
   * array */
  if (!(conv && noise))
    filter_type = SEP_FILTER_CONV;

  if (conv)
    {
      /* allocate memory for convolved buffers */
      QMALLOC(cdscan, PIXTYPE, stacksize, status);
      if (filter_type == SEP_FILTER_MATCHED)
        {
          QMALLOC(sigscan, PIXTYPE, stacksize, status);
          QMALLOC(workscan, PIXTYPE, stacksize, status);
        }

      /* normalize the filter */
      convn = convw * convh;
      QMALLOC(convnorm, PIXTYPE, convn, status);
      for (i=0; i<convn; i++)
	sum += fabs(conv[i]);
      for (i=0; i<convn; i++)
	convnorm[i] = conv[i] / sum;
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
	    }
	  cdscan = dummyscan;
	}

      else
	{
          arraybuffer_readline(&imbuf);
          if (noise)
            arraybuffer_readline(&nbuf);
          if (mask)
            {
              arraybuffer_readline(&mbuf);
              apply_mask_line(&mbuf, &imbuf, (noise? &nbuf: NULL));
            }

	  /* filter the lines */
	  if (conv)
	    {
              status = convolve(&imbuf, yl, convnorm, convw, convh, cdscan);
              if (status != RETURN_OK)
                goto exit;

	      if (filter_type == SEP_FILTER_MATCHED)
                {
                  status = matched_filter(&imbuf, &nbuf, yl, convnorm, convw,
                                          convh, workscan, sigscan);

                  if (status != RETURN_OK)
                    goto exit;
                }
            }
	  else
	    {
	      cdscan = scan;
	    }	  
	}
      
      trunflag = (yl==0 || yl==h-1)? SEP_OBJ_TRUNC: 0;
      
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
            thresh = relthresh * ((xl==w || yl==h)? 0.0: wscan[xl]);

          /* luflag: is pixel above thresh (Y/N)? */
          if (filter_type == SEP_FILTER_MATCHED)
            luflag = ((xl != w) && (sigscan[xl] > relthresh))? 1: 0;
          else
            luflag = cdnewsymbol > thresh? 1: 0;

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
		PLISTPIX(pixt, var) = wscan[xl];  /* not currently used */
	      if (PLISTEXIST(thresh))
		PLISTPIX(pixt, thresh) = thresh;

	      /* Check if we are running out of free pixels in objlist.plist */
	      if (freeinfo.firstpix==freeinfo.lastpix)
		{
		  status = PIXSTACK_FULL;
		  sprintf(errtext,
			  "The limit of %d active object pixels over the "
			  "detection threshold was reached. Check that "
			  "the image is background subtracted and the "
			  "detection threshold is not too low. If you "
                          "need to increase the limit, use "
                          "sep.set_extract_pixstack.",
			  (int)mem_pixstack);
		  put_errdetail(errtext);
		  goto exit;

		  /* The code in the rest of this block increases the
		   * stack size as needed. Currently, it is never
		   * executed. This is because it isn't clear that
		   * this is a good idea: most times when the stack
		   * overflows it is due to user error: too-low
		   * threshold or image not background subtracted. */

		  /* increase the stack size */
		  oldnposize = nposize;
 		  mem_pixstack = (int)(mem_pixstack * 2);
		  nposize = mem_pixstack * plistsize;
		  pixel = (pliststruct *)realloc(pixel, nposize);
		  objlist.plist = pixel;
		  if (!pixel)
		    {
		      status = MEMORY_ALLOC_ERROR;
		      goto exit;
		    }

		  /* set next free pixel to the start of the new block 
		   * and link up all the pixels in the new block */
		  PLIST(pixel+freeinfo.firstpix, nextpix) = oldnposize;
		  pixt = pixel + oldnposize;
		  for (i=oldnposize + plistsize; i<nposize;
		       i += plistsize, pixt += plistsize)
		    PLIST(pixt, nextpix) = i;
		  PLIST(pixt, nextpix) = -1;

		  /* last free pixel is now at the end of the new block */
		  freeinfo.lastpix = nposize - plistsize;
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
			      /* update threshold before object is processed */
			      objlist.thresh = thresh;

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
  free(dummyscan);
  free(psstack);
  free(start);
  free(end);
  arraybuffer_free(&imbuf);
  if (noise)
    arraybuffer_free(&nbuf);
  if (mask)
    arraybuffer_free(&mbuf);
  if (conv)
    free(convnorm);
  if (filter_type == SEP_FILTER_MATCHED)
    {
      free(sigscan);
      free(workscan);
    }

  if (status != RETURN_OK)
    {
      free(cdscan);   /* only need to free these in case of early exit */
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
  obj.thresh = objlist->thresh;

  preanalyse(0, objlist);

  status = deblend(objlist, 0, &objlistout, deblend_nthresh, deblend_mincont,
		   minarea);
  if (status)
    {
      /* formerly, this wasn't a fatal error, so a flag was set for
       * the object and we continued. I'm leaving the flag-setting
       * here in case we want to change this to a non-fatal error in
       * the future, but currently the flag setting is irrelevant. */
      objlist2 = objlist;
      for (i=0; i<objlist2->nobj; i++)
	objlist2->obj[i].flag |= SEP_OBJ_DOVERFLOW;
      goto exit;
    }
  else
    objlist2 = &objlistout;

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

      plistexist_thresh = 1;
      plistoff_thresh = plistsize;
      plistsize += sizeof(PIXTYPE);
    }
  else
    {
      plistexist_var = 0;
      plistexist_thresh = 0;
    }

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
      alphain = (pow(ampin/obj1->thresh, 1.0/beta)-1)*unitareain/obj1->fdnpix;

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
	      alpha = (pow(amp/obj2->thresh, 1.0/beta) - 1) *
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

  objout->thresh = obj->thresh;
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
