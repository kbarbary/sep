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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sep.h"
#include "sepcore.h"

#define	BACK_MINGOODFRAC   0.5   /* min frac with good weights*/
#define	QUANTIF_NSIGMA     5     /* histogram limits */
#define	QUANTIF_NMAXLEVELS 4096  /* max nb of quantif. levels */
#define	QUANTIF_AMIN       4     /* min nb of "mode pixels" */

/* Background info in a single mesh*/
typedef struct
{
  float	 mode, mean, sigma;	/* Background mode, mean and sigma */
  LONG	 *histo;	       	/* Pointer to a histogram */
  int	 nlevels;		/* Nb of histogram bins */
  float	 qzero, qscale;		/* Position of histogram */
  float	 lcut, hcut;		/* Histogram cuts */
  int	 npix;			/* Number of pixels involved */
} backstruct;

void backhisto(backstruct *, PIXTYPE *, PIXTYPE *,
	       int, int, int, int, PIXTYPE);
void backstat(backstruct *, PIXTYPE *, PIXTYPE *,
	      int, int, int, int, PIXTYPE);
int filterback(sepbackmap *bkmap, int fw, int fh, float fthresh);
float backguess(backstruct *, float *, float *);
int makebackspline(sepbackmap *, float *, float *);

int sep_makeback(void *im, void *mask, int dtype, int mdtype, int w, int h,
		 int bw, int bh, float mthresh, int fw, int fh,
		 float fthresh, sepbackmap **bkm)
{
  BYTE *imt, *maskt;
  int npix;                   /* size of image */
  int nx, ny, nb;             /* number of background boxes in x, y, total */
  int bufsize;                /* size of a "row" of boxes in pixels (w*bh) */
  int elsize;                 /* size (in bytes) of an image array element */
  int melsize;                /* size (in bytes) of a mask array element */
  PIXTYPE *buf, *buft, *mbuf, *mbuft;
  PIXTYPE maskthresh;
  array_converter convert, mconvert;
  backstruct *backmesh, *bm;  /* info about each background "box" */
  sepbackmap *bkmap;          /* output */
  int j,k,m, status;

  status = RETURN_OK;
  npix = w*h;
  bufsize = w*bh;
  maskthresh = mask? mthresh: 0.0;

  backmesh = bm = NULL;
  bkmap = NULL;
  buf = mbuf = buft = mbuft = NULL;
  convert = mconvert = NULL;

  /* determine number of background boxes */
  if ((nx = (w-1)/bw + 1) < 1)
    nx = 1;
  if ((ny = (h-1)/bh + 1) < 1)
    ny = 1;
  nb = nx*ny;

  /* Allocate temp memory & initialize */
  QMALLOC(backmesh, backstruct, nx, status);
  bm = backmesh;
  for (m=nx; m--; bm++)
    bm->histo=NULL;

  /* Allocate the returned struct */
  QMALLOC(bkmap, sepbackmap, 1, status);
  bkmap->w = w;
  bkmap->h = h;
  bkmap->nx = nx;
  bkmap->ny = ny;
  bkmap->n = nb;
  bkmap->bw = bw;
  bkmap->bh = bh;
  bkmap->back = NULL;
  bkmap->sigma = NULL;
  bkmap->dback = NULL;
  bkmap->dsigma = NULL;
  QMALLOC(bkmap->back, float, nb, status);
  QMALLOC(bkmap->sigma, float, nb, status);
  QMALLOC(bkmap->dback, float, nb, status);
  QMALLOC(bkmap->dsigma, float, nb, status);

  /* cast input array pointers. These are used to step through the arrays. */
  imt = (BYTE *)im;
  maskt = (BYTE *)mask;

  /* get the correct array converter and element size, based on dtype code */
  status = get_array_converter(dtype, &convert, &elsize);
  if (status != RETURN_OK)
    goto exit;
  if (mask)
    {
      status = get_array_converter(mdtype, &mconvert, &melsize);
      if (status != RETURN_OK)
	goto exit;
    }

  /* If the input array type is not PIXTYPE, allocate a buffer to hold
     converted values */
  if (dtype != PIXDTYPE)
    {
      QMALLOC(buf, PIXTYPE, bufsize, status);
      buft = buf;
      if (status != RETURN_OK)
	goto exit;
    }
  if (mask && (mdtype != PIXDTYPE))
    {
      QMALLOC(mbuf, PIXTYPE, bufsize, status);
      mbuft = mbuf;
      if (status != RETURN_OK)
	goto exit;
    }

  /* loop over rows of background boxes.
   * (here, we could loop over individual boxes rather than entire
   * rows, but this is convenient for converting the image and mask
   * arrays.  This is also how it is originally done in SExtractor,
   * because the pixel buffers are only read in from disk in
   * increments of a row of background boxes at a time.)
   */
  for (j=0; j<ny; j++)
    {
      /* if the last row, modify the width appropriately*/
      if (j == ny-1 && npix%bufsize)
        bufsize = npix%bufsize;

      /* convert this row to PIXTYPE and store in buffer(s)*/
      if (dtype != PIXDTYPE)
	convert(imt, bufsize, buft);
      else
	buft = (PIXTYPE *)imt;

      if (mask)
	{
	  if (mdtype != PIXDTYPE)
	    mconvert(maskt, bufsize, mbuft);
	  else
	    mbuft = (PIXTYPE *)maskt;
	}

      /* Get clipped mean, sigma for all boxes in the row */
      backstat(backmesh, buft, mbuft, bufsize, nx, w, bw, maskthresh);

      /* Allocate histograms in each box in this row. */
      bm = backmesh;
      for (m=nx; m--; bm++)
	if (bm->mean <= -BIG)
	  bm->histo=NULL;
	else
	  QCALLOC(bm->histo, LONG, bm->nlevels, status);
      backhisto(backmesh, buft, mbuft, bufsize, nx, w, bw, maskthresh);

      /* Compute background statistics from the histograms */
      bm = backmesh;
      for (m=0; m<nx; m++, bm++)
	{
	  k = m+nx*j;
	  backguess(bm, bkmap->back+k, bkmap->sigma+k);
	  free(bm->histo);
	  bm->histo = NULL;
	}

      /* increment array pointers to next row of background boxes */
      imt += elsize * bufsize;
      if (mask)
	maskt += melsize * bufsize;
    }

  /* free memory */
  free(buf);
  buf = NULL;
  free(mbuf);
  mbuf = NULL;
  free(backmesh);
  backmesh = NULL;

  /* Median-filter and check suitability of the background map */
  if ((status = filterback(bkmap, fw, fh, fthresh)) != RETURN_OK)
    goto exit;

  /* Compute 2nd derivatives along the y-direction */
  if ((status = makebackspline(bkmap, bkmap->back, bkmap->dback)) !=
      RETURN_OK)
    goto exit;
  if ((status = makebackspline(bkmap, bkmap->sigma, bkmap->dsigma)) !=
      RETURN_OK)
    goto exit;

  *bkm = bkmap;
  return status;

  /* If we encountered a problem, clean up any allocated memory */
 exit:
  free(buf);
  free(mbuf);
  if (backmesh)
    {
      bm = backmesh;
      for (m=0; m<nx; m++, bm++)
	free(bm->histo);
    }
  free(backmesh);
  sep_freeback(bkmap);
  *bkm = NULL;
  return status;
}

/******************************** backstat **********************************/
/*
Compute robust statistical estimators in a row of meshes.
*/
void backstat(backstruct *backmesh,
	      PIXTYPE *buf, PIXTYPE *wbuf, int bufsize,
	      int n, int w, int bw, PIXTYPE maskthresh)
{
  backstruct	*bm;
  double	pix, wpix, sig, mean, sigma, step;
  PIXTYPE	*buft,*wbuft;
  PIXTYPE       lcut,hcut;
  int		m,h,x,y, npix,wnpix, offset, lastbite;
  
  h = bufsize/w;  /* height of background boxes in this row */
  bm = backmesh;
  offset = w - bw;
  step = sqrt(2/PI)*QUANTIF_NSIGMA/QUANTIF_AMIN;

  for (m = n; m--; bm++,buf+=bw)
    {
      if (!m && (lastbite=w%bw))
	{
	  bw = lastbite;
	  offset = w-bw;
	}

      mean = sigma = 0.0;
      buft=buf;
      npix = 0;

      /* We separate the weighted case at this level to avoid penalty in CPU */
      if (wbuf)
	{
	  wbuft = wbuf;
	  for (y=h; y--; buft+=offset,wbuft+=offset)
	    for (x=bw; x--;)
	      {
		pix = *(buft++);
		if ((wpix = *(wbuft++)) <= maskthresh && pix > -BIG)
		  {
		    mean += pix;
		    sigma += pix*pix;
		    npix++;
		  }
	      }
	}
      else
	for (y=h; y--; buft+=offset)
	  for (x=bw; x--;)
	    if ((pix = *(buft++)) > -BIG)
	      {
		mean += pix;
		sigma += pix*pix;
		npix++;
	      }

      /*-- If not enough valid pixels, discard this mesh */
      if ((float)npix < (float)(bw*h*BACK_MINGOODFRAC))
	{
	  bm->mean = bm->sigma = -BIG;
	  if (wbuf)
	    wbuf += bw;
	  continue;
	}

      mean /= (double)npix;
      sigma = (sig = sigma/npix - mean*mean)>0.0? sqrt(sig):0.0;
      lcut = bm->lcut = (PIXTYPE)(mean - 2.0*sigma);
      hcut = bm->hcut = (PIXTYPE)(mean + 2.0*sigma);
      mean = sigma = 0.0;
      npix = wnpix = 0;
      buft = buf;
      
      /* do statistics for this mesh again, with cuts */
      if (wbuf)
	{
	  wbuft=wbuf;
	  for (y=h; y--; buft+=offset, wbuft+=offset)
	    for (x=bw; x--;)
	      {
		pix = *(buft++);
		if ((wpix = *(wbuft++))<=maskthresh && pix<=hcut && pix>=lcut)
		  {
		    mean += pix;
		    sigma += pix*pix;
		    npix++;
		  }
	      }
	}
      else
	for (y=h; y--; buft+=offset)
	  for (x=bw; x--;)
	    {
	      pix = *(buft++);
	      if (pix<=hcut && pix>=lcut)
		{
		  mean += pix;
		  sigma += pix*pix;
		  npix++;
		}
	    }

      bm->npix = npix;
      mean /= (double)npix;
      sig = sigma/npix - mean*mean;
      sigma = sig>0.0 ? sqrt(sig):0.0;
      bm->mean = mean;
      bm->sigma = sigma;
      if ((bm->nlevels = (int)(step*npix+1)) > QUANTIF_NMAXLEVELS)
	bm->nlevels = QUANTIF_NMAXLEVELS;
      bm->qscale = sigma>0.0? 2*QUANTIF_NSIGMA*sigma/bm->nlevels : 1.0;
      bm->qzero = mean - QUANTIF_NSIGMA*sigma;

      if (wbuf)
	wbuf += bw;
    }

  return;

}

/******************************** backhisto *********************************/
/*
Fill histograms in a row of meshes.
*/
void backhisto(backstruct *backmesh,
	       PIXTYPE *buf, PIXTYPE *wbuf, int bufsize,
	       int n, int w, int bw, PIXTYPE maskthresh)
{
  backstruct	*bm;
  PIXTYPE	*buft,*wbuft;
  float	        qscale, cste, wpix;
  LONG		*histo;
  int		h,m,x,y, nlevels, lastbite, offset, bin;

  h = bufsize/w;
  bm = backmesh;
  offset = w - bw;
  for (m=0; m++<n; bm++ , buf+=bw)
    {
      if (m==n && (lastbite=w%bw))
	{
	  bw = lastbite;
	  offset = w-bw;
	}

      /*-- Skip bad meshes */
      if (bm->mean <= -BIG)
	{
	  if (wbuf)
	    wbuf += bw;
	  continue;
	}

      nlevels = bm->nlevels;
      histo = bm->histo;
      qscale = bm->qscale;
      cste = 0.499999 - bm->qzero/qscale;
      buft = buf;

      if (wbuf)
	{
	  wbuft = wbuf;
	  for (y=h; y--; buft+=offset, wbuft+=offset)
	    for (x=bw; x--;)
	      {
		bin = (int)(*(buft++)/qscale + cste);
		if ((wpix = *(wbuft++))<=maskthresh && bin<nlevels && bin>=0)
		  (*(histo+bin))++;
	      }
	  wbuf += bw;
	}
      else
	for (y=h; y--; buft += offset)
	  for (x=bw; x--;)
	    {
	      bin = (int)(*(buft++)/qscale + cste);
	      
	      if (bin>=0 && bin<nlevels)
		(*(histo+bin))++;
	    }
    }
  return;
}

/******************************* backguess **********************************/
/*
Estimate the background from a histogram;
*/
float	backguess(backstruct *bkg, float *mean, float *sigma)

#define	EPS	(1e-4)	/* a small number */

{
  LONG		*histo, *hilow, *hihigh, *histot;
  unsigned long lowsum, highsum, sum;
  double	ftemp, mea, sig, sig1, med, dpix;
  int		i, n, lcut,hcut, nlevelsm1, pix;

  /* Leave here if the mesh is already classified as `bad' */
  if (bkg->mean<=-BIG)
    {
      *mean = *sigma = -BIG;
      return -BIG;
    }

  histo = bkg->histo;
  hcut = nlevelsm1 = bkg->nlevels-1;
  lcut = 0;

  sig = 10.0*nlevelsm1;
  sig1 = 1.0;
  mea = med = bkg->mean;

  /* iterate until sigma converges or drops below 0.1 (up to 100 iterations) */
  for (n=100; n-- && (sig>=0.1) && (fabs(sig/sig1-1.0)>EPS);)
    {
      sig1 = sig;
      sum = mea = sig = 0.0;
      lowsum = highsum = 0;
      histot = hilow = histo+lcut;
      hihigh = histo+hcut;

      for (i=lcut; i<=hcut; i++)
	{
	  if (lowsum<highsum)
	    lowsum += *(hilow++);
	  else
	    highsum +=  *(hihigh--);
	  sum += (pix = *(histot++));
	  mea += (dpix = (double)pix*i);
	  sig += dpix*i;
	}

      med = hihigh>=histo?((hihigh-histo) + 0.5 +
			   ((double)highsum-lowsum)/(2.0*(*hilow>*hihigh?
							   *hilow:*hihigh)))
	: 0.0;
      if (sum)
	{
	  mea /= (double)sum;
	  sig = sig/sum - mea*mea;
	}

      sig = sig>0.0?sqrt(sig):0.0;
      lcut = (ftemp=med-3.0*sig)>0.0 ?(int)(ftemp>0.0?ftemp+0.5:ftemp-0.5):0;
      hcut = (ftemp=med+3.0*sig)<nlevelsm1 ?(int)(ftemp>0.0?ftemp+0.5:ftemp-0.5)
	: nlevelsm1;

    }
  *mean = fabs(sig)>0.0? (fabs(bkg->sigma/(sig*bkg->qscale)-1) < 0.0 ?
			  bkg->qzero+mea*bkg->qscale
			  :(fabs((mea-med)/sig)< 0.3 ?
			    bkg->qzero+(2.5*med-1.5*mea)*bkg->qscale
			    :bkg->qzero+med*bkg->qscale))
    :bkg->qzero+mea*bkg->qscale;
  
  *sigma = sig*bkg->qscale;

  return *mean;
}

/****************************************************************************/

int filterback(sepbackmap *bkmap, int fw, int fh, float fthresh)
/* Median filterthe background map to remove the contribution
 * from bright sources. */
{
  float	*back, *sigma, *back2, *sigma2, *bmask, *smask, *sigmat;
  float d2, d2min, med, val, sval;
  int i, j, px, py, np, nx, ny, npx, npx2, npy, npy2, dpx, dpy, x, y, nmin;
  int status;
  
  status = RETURN_OK;
  bmask = smask = back2 = sigma2 = NULL;

  nx = bkmap->nx;
  ny = bkmap->ny;
  np = bkmap->n;
  npx = fw/2;
  npy = fh/2;
  npy *= nx;

  QMALLOC(bmask, float, (2*npx+1)*(2*npy+1), status);
  QMALLOC(smask, float, (2*npx+1)*(2*npy+1), status);
  QMALLOC(back2, float, np, status);
  QMALLOC(sigma2, float, np, status);

  back = bkmap->back;
  sigma = bkmap->sigma;
  val = sval = 0.0;  /* to avoid gcc -Wall warnings */

/* Look for `bad' meshes and interpolate them if necessary */
  for (i=0,py=0; py<ny; py++)
    for (px=0; px<nx; px++,i++)
      if ((back2[i]=back[i])<=-BIG)
        {
/*------ Seek the closest valid mesh */
        d2min = BIG;
        nmin = 0.0;
        for (j=0,y=0; y<ny; y++)
          for (x=0; x<nx; x++,j++)
            if (back[j]>-BIG)
              {
              d2 = (float)(x-px)*(x-px)+(y-py)*(y-py);
              if (d2<d2min)
                {
                val = back[j];
                sval = sigma[j];
                nmin = 1;
                d2min = d2;
                }
              else if (d2==d2min)
                {
                val += back[j];
                sval += sigma[j];
                nmin++;
                }
              }
        back2[i] = nmin? val/nmin: 0.0;
        sigma[i] = nmin? sval/nmin: 1.0;
        }
  memcpy(back, back2, (size_t)np*sizeof(float));

/* Do the actual filtering */
  for (py=0; py<np; py+=nx)
    {
    npy2 = np - py - nx;
    if (npy2>npy)
      npy2 = npy;
    if (npy2>py)
      npy2 = py;
    for (px=0; px<nx; px++)
      {
      npx2 = nx - px - 1;
      if (npx2>npx)
        npx2 = npx;
      if (npx2>px)
        npx2 = px;
      i=0;
      for (dpy = -npy2; dpy<=npy2; dpy+=nx)
        {
        y = py+dpy;
        for (dpx = -npx2; dpx <= npx2; dpx++)
          {
          x = px+dpx;
          bmask[i] = back[x+y];
          smask[i++] = sigma[x+y];
          }
        }
      if (fabs((med=fqmedian(bmask, i))-back[px+py])>=fthresh)
        {
        back2[px+py] = med;
        sigma2[px+py] = fqmedian(smask, i);
        }
      else
        {
        back2[px+py] = back[px+py];
        sigma2[px+py] = sigma[px+py];
        }
      }
    }

  free(bmask);
  free(smask);
  bmask = smask = NULL;
  memcpy(back, back2, np*sizeof(float));
  bkmap->globalback = fqmedian(back2, np);
  free(back2);
  back2 = NULL;
  memcpy(sigma, sigma2, np*sizeof(float));
  bkmap->globalrms = fqmedian(sigma2, np);

  if (bkmap->globalrms <= 0.0)
    {
    sigmat = sigma2+np;
    for (i=np; i-- && *(--sigmat)>0.0;);
    if (i>=0 && i<(np-1))
      bkmap->globalrms = fqmedian(sigmat+1, np-1-i);
    else
      bkmap->globalrms = 1.0;
    }

  free(sigma2);
  sigma2 = NULL;

  return status;

 exit:
  free(bmask);
  free(smask);
  free(back2);
  free(sigma2);
  return status;
}


/******************************* makebackspline ******************************/
/*
 * Pre-compute 2nd derivatives along the y direction at background nodes.
 */
int makebackspline(sepbackmap *bkmap, float *map, float *dmap)
{
  int   x, y, nbx, nby, nbym1, status;
  float *dmapt, *mapt, *u, temp;
  u = NULL;
  status = RETURN_OK;

  nbx = bkmap->nx;
  nby = bkmap->ny;
  nbym1 = nby - 1;
  for (x=0; x<nbx; x++)
    {
      mapt = map+x;
      dmapt = dmap+x;
      if (nby>1)
	{
	  QMALLOC(u, float, nbym1, status); /* temporary array */
	  *dmapt = *u = 0.0;	/* "natural" lower boundary condition */
	  mapt += nbx;
	  for (y=1; y<nbym1; y++, mapt+=nbx)
	    {
	      temp = -1/(*dmapt+4);
	      *(dmapt += nbx) = temp;
	      temp *= *(u++) - 6*(*(mapt+nbx)+*(mapt-nbx)-2**mapt);
	      *u = temp;
	    }
	  *(dmapt+=nbx) = 0.0;	/* "natural" upper boundary condition */
	  for (y=nby-2; y--;)
	    {
	      temp = *dmapt;
	      dmapt -= nbx;
	      *dmapt = (*dmapt*temp+*(u--))/6.0;
	    }
	  free(u);
	  u = NULL;
	}
      else
	*dmapt = 0.0;
    }

  return status;

 exit: 
  free(u);
  return status;
}

/*****************************************************************************/

float sep_globalback(sepbackmap *bkmap)
{
  return bkmap->globalback;
}

float sep_globalrms(sepbackmap *bkmap)
{
  return bkmap->globalrms;
}


/*****************************************************************************/

float sep_backpix_linear(sepbackmap *bkmap, int x, int y)
/*
 * return background at position x,y.
 * (linear interpolation between background map vertices).
 */
{
  int    nx, ny, xl, yl, pos;
  double dx, dy, cdx;
  float	 *b;
  float  b0, b1, b2, b3;

  b = bkmap->back;
  nx = bkmap->nx;
  ny = bkmap->ny;

  dx = (double)x/bkmap->bw - 0.5;
  dy = (double)y/bkmap->bh - 0.5;
  dx -= (xl = (int)dx);
  dy -= (yl = (int)dy);

  if (xl<0)
    {
    xl = 0;
    dx -= 1.0;
    }
  else if (xl>=nx-1)
    {
    xl = nx<2 ? 0 : nx-2;
    dx += 1.0;
    }

  if (yl<0)
    {
    yl = 0;
    dy -= 1.0;
    }
  else if (yl>=ny-1)
    {
    yl = ny<2 ? 0 : ny-2;
    dy += 1.0;
    }

  pos = yl*nx + xl;
  cdx = 1 - dx;

  b0 = *(b+=pos);		/* consider when nbackx or nbacky = 1 */
  b1 = nx<2? b0:*(++b);
  b2 = ny<2? *b:*(b+=nx);
  b3 = nx<2? *b:*(--b);

  return (float)((1-dy)*(cdx*b0 + dx*b1) + dy*(dx*b2 + cdx*b3));
}


/*****************************************************************************/

int sep_backline_flt(sepbackmap *bkmap, int y, float *line)
/* Interpolate background at line y (bicubic spline interpolation between
 * background map vertices) and save to line */
{
  int i,j,x,yl, nbx,nbxm1,nby, nx,width, ystep, changepoint, status;
  float	dx,dx0,dy,dy3, cdx,cdy,cdy3, temp, xstep;
  float *nodebuf, *dnodebuf, *u;
  float *node, *nodep, *dnode, *blo, *bhi, *dblo, *dbhi;

  status = RETURN_OK;
  nodebuf = node = NULL;
  dnodebuf = dnode = NULL;
  u = NULL;

  width = bkmap->w;
  nbx = bkmap->nx;
  nbxm1 = nbx - 1;
  nby = bkmap->ny;
  if (nby > 1)
    {
      dy = (float)y/bkmap->bh - 0.5;
      dy -= (yl = (int)dy);
      if (yl<0)
	{
	  yl = 0;
	  dy -= 1.0;
	}
      else if (yl>=nby-1)
	{
	  yl = nby<2 ? 0 : nby-2;
	  dy += 1.0;
	}
      /*-- Interpolation along y for each node */
      cdy = 1 - dy;
      dy3 = (dy*dy*dy-dy);
      cdy3 = (cdy*cdy*cdy-cdy);
      ystep = nbx*yl;
      blo = bkmap->back + ystep;
      bhi = blo + nbx;
      dblo = bkmap->dback + ystep;
      dbhi = dblo + nbx;
      QMALLOC(nodebuf, float, nbx, status);  /* Interpolated background */
      nodep = node = nodebuf;
      for (x=nbx; x--;)
	*(nodep++) = cdy**(blo++) + dy**(bhi++) + cdy3**(dblo++) +
	  dy3**(dbhi++);

      /*-- Computation of 2nd derivatives along x */
      QMALLOC(dnodebuf, float, nbx, status);  /* 2nd derivative along x */
      dnode = dnodebuf;
      if (nbx>1)
	{
	  QMALLOC(u, float, nbxm1, status);  /* temporary array */
	  *dnode = *u = 0.0;	/* "natural" lower boundary condition */
	  nodep = node+1;
	  for (x=nbxm1; --x; nodep++)
	    {
	      temp = -1/(*(dnode++)+4);
	      *dnode = temp;
	      temp *= *(u++) - 6*(*(nodep+1)+*(nodep-1)-2**nodep);
	      *u = temp;
	    }
	  *(++dnode) = 0.0;	/* "natural" upper boundary condition */
	  for (x=nbx-2; x--;)
	    {
	      temp = *(dnode--);
	      *dnode = (*dnode*temp+*(u--))/6.0;
	    }
	  free(u);
	  u = NULL;
	  dnode--;
	}
    }
  else
    {
      /*-- No interpolation and no new 2nd derivatives needed along y */
      node = bkmap->back;
      dnode = bkmap->dback;
    }

  /*-- Interpolation along x */
  if (nbx>1)
    {
      nx = bkmap->bw;
      xstep = 1.0/nx;
      changepoint = nx/2;
      dx  = (xstep - 1)/2;	/* dx of the first pixel in the row */
      dx0 = ((nx+1)%2)*xstep/2;	/* dx of the 1st pixel right to a bkgnd node */
      blo = node;
      bhi = node + 1;
      dblo = dnode;
      dbhi = dnode + 1;
      for (x=i=0,j=width; j--; i++, dx += xstep)
	{
	  if (i==changepoint && x>0 && x<nbxm1)
	    {
	      blo++;
	      bhi++;
	      dblo++;
	      dbhi++;
	      dx = dx0;
	    }
	  cdx = 1 - dx;

	  *(line++) = (float)(cdx*(*blo+(cdx*cdx-1)**dblo) 
			      + dx*(*bhi+(dx*dx-1)**dbhi));

	  if (i==nx)
	    {
	      x++;
	      i = 0;
	    }
	}
    }
  else
    for (j=width; j--;)
      {
	*(line++) = (float)*node;
      }

 exit:
  free(nodebuf);
  free(dnodebuf);
  free(u);
  return status;
}


/*****************************************************************************/

int sep_backrmsline_flt(sepbackmap *bkmap, int y, float *line)
/* Bicubic-spline interpolation of the background noise along the current
 * scanline (y). NOTE: Most of the code is a copy of subbackline(), for
 * optimization reasons.
 */
{
  int i,j,x,yl, nbx,nbxm1,nby, nx, width, ystep, changepoint, status;
  float	dx,dx0,dy,dy3, cdx,cdy,cdy3, temp, xstep;
  float *nodebuf, *dnodebuf, *u;
  float *node, *nodep, *dnode, *blo, *bhi, *dblo, *dbhi;
  status = RETURN_OK;
  nodebuf = node = NULL;
  dnodebuf = dnode = NULL;
  u = NULL;

  nbx = bkmap->nx;
  nbxm1 = nbx - 1;
  nby = bkmap->ny;
  if (nby > 1)
    {
      dy = (float)y/bkmap->bh - 0.5;
      dy -= (yl = (int)dy);
      if (yl<0)
	{
	  yl = 0;
	  dy -= 1.0;
	}
      else if (yl>=nby-1)
	{
	  yl = nby<2 ? 0 : nby-2;
	  dy += 1.0;
	}
      /*-- Interpolation along y for each node */
      cdy = 1 - dy;
      dy3 = (dy*dy*dy-dy);
      cdy3 = (cdy*cdy*cdy-cdy);
      ystep = nbx*yl;
      blo = bkmap->sigma + ystep;
      bhi = blo + nbx;
      dblo = bkmap->dsigma + ystep;
      dbhi = dblo + nbx;
      QMALLOC(nodebuf, float, nbx, status); /* Interpolated background */
      nodep = node = nodebuf;
      for (x=nbx; x--;)
	*(nodep++) = cdy**(blo++)+dy**(bhi++)+cdy3**(dblo++)+dy3**(dbhi++);

      /*-- Computation of 2nd derivatives along x */
      QMALLOC(dnodebuf, float, nbx, status); /* 2nd derivative along x */
      dnode = nodebuf;
      if (nbx>1)
	{
	  QMALLOC(u, float, nbxm1, status);	/* temporary array */
	  *dnode = *u = 0.0;	/* "natural" lower boundary condition */
	  nodep = node+1;
	  for (x=nbxm1; --x; nodep++)
	    {
	      temp = -1/(*(dnode++)+4);
	      *dnode = temp;
	      temp *= *(u++) - 6*(*(nodep+1)+*(nodep-1)-2**nodep);
	      *u = temp;
	    }
	  *(++dnode) = 0.0;	/* "natural" upper boundary condition */
	  for (x=nbx-2; x--;)
	    {
	      temp = *(dnode--);
	      *dnode = (*dnode*temp+*(u--))/6.0;
	    }
	  free(u);
	  u = NULL;
	  dnode--;
	}
    }
  else
    {
      /*-- No interpolation and no new 2nd derivatives needed along y */
      node = bkmap->sigma;
      dnode = bkmap->dsigma;
    }
  
  /*-- Interpolation along x */
  width = bkmap->w;
  if (nbx>1)
    {
      nx = bkmap->bw;
      xstep = 1.0/nx;
      changepoint = nx/2;
      dx  = (xstep - 1)/2;	/* dx of the first pixel in the row */
      dx0 = ((nx+1)%2)*xstep/2;	/* dx of the 1st pixel right to a bkgnd node */
      blo = node;
      bhi = node + 1;
      dblo = dnode;
      dbhi = dnode + 1;
      for (x=i=0,j=width; j--; i++, dx += xstep)
	{
	  if (i==changepoint && x>0 && x<nbxm1)
	    {
	      blo++;
	      bhi++;
	      dblo++;
	      dbhi++;
	      dx = dx0;
	    }
	  cdx = 1 - dx;
	  *(line++) = (float)(cdx*(*blo+(cdx*cdx-1)**dblo)
			      + dx*(*bhi+(dx*dx-1)**dbhi));
	  if (i==nx)
	    {
	      x++;
	      i = 0;
	    }
	}
    }
  else
    for (j=width; j--;)
      *(line++) = (float)*node;

 exit:
  free(node);
  free(dnode);
  free(u);
  return status;
}

/*****************************************************************************/
/* Multiple dtype functions and convenience functions.
 * These mostly wrap the two "line" functions above. */

int sep_backline(sepbackmap *bkmap, int y, void *line, int dtype)
{
  array_writer write_array;
  int size, status;
  float *tmpline; 

  if (dtype == SEP_TFLOAT)
    return sep_backline_flt(bkmap, y, (float *)line);
   
  tmpline = NULL;

  status = get_array_writer(dtype, &write_array, &size);
  if (status != RETURN_OK)
    goto exit;

  QMALLOC(tmpline, float, bkmap->w, status);
  status = sep_backline_flt(bkmap, y, tmpline);
  if (status != RETURN_OK)
    goto exit;

  /* write to desired output type */
  write_array(tmpline, bkmap->w, line);

 exit:
  free(tmpline);
  return status;
}

int sep_backrmsline(sepbackmap *bkmap, int y, void *line, int dtype)
{
  array_writer write_array;
  int size, status;
  float *tmpline; 

  if (dtype == SEP_TFLOAT)
    return sep_backrmsline_flt(bkmap, y, (float *)line);
   
  tmpline = NULL;

  status = get_array_writer(dtype, &write_array, &size);
  if (status != RETURN_OK)
    goto exit;

  QMALLOC(tmpline, float, bkmap->w, status);
  status = sep_backrmsline_flt(bkmap, y, tmpline);
  if (status != RETURN_OK)
    goto exit;

  /* write to desired output type */
  write_array(tmpline, bkmap->w, line);

 exit:
  free(tmpline);
  return status;
}

int sep_backarray(sepbackmap *bkmap, void *arr, int dtype)
{
  int y, width, size, status;
  array_writer write_array;
  float *tmpline;
  BYTE *line;

  tmpline = NULL;
  status = RETURN_OK;
  width = bkmap->w;

  if (dtype == SEP_TFLOAT)
    {
      tmpline = (float *)arr;
      for (y=0; y<bkmap->h; y++, tmpline+=width)
	if ((status = sep_backline_flt(bkmap, y, tmpline)) != RETURN_OK)
	  return status;
      return status;
    }
  
  if ((status = get_array_writer(dtype, &write_array, &size)) != RETURN_OK)
    goto exit;

  QMALLOC(tmpline, float, width, status);

  line = (BYTE *)arr;
  for (y=0; y<bkmap->h; y++, line += size*width)
    {
      if ((status = sep_backline_flt(bkmap, y, tmpline)) != RETURN_OK)
	goto exit;
      write_array(tmpline, width, line);
    }

 exit:
  free(tmpline);
  return status;
}

int sep_backrmsarray(sepbackmap *bkmap, void *arr, int dtype)
{
  int y, width, size, status;
  array_writer write_array;
  float *tmpline;
  BYTE *line;

  tmpline = NULL;
  status = RETURN_OK;
  width = bkmap->w;

  if (dtype == SEP_TFLOAT)
    {
      tmpline = (float *)arr;
      for (y=0; y<bkmap->h; y++, tmpline+=width)
	if ((status = sep_backrmsline_flt(bkmap, y, tmpline)) != RETURN_OK)
	  return status;
      return status;
    }
  
  if ((status = get_array_writer(dtype, &write_array, &size)) != RETURN_OK)
    goto exit;

  QMALLOC(tmpline, float, width, status);

  line = (BYTE *)arr;
  for (y=0; y<bkmap->h; y++, line += size*width)
    {
      if ((status = sep_backrmsline_flt(bkmap, y, tmpline)) != RETURN_OK)
	goto exit;
      write_array(tmpline, width, line);
    }

 exit:
  free(tmpline);
  return status;
}

int sep_subbackline(sepbackmap *bkmap, int y, void *line, int dtype)
{
  array_writer subtract_array;
  int status, size;
  PIXTYPE *tmpline;

  tmpline = NULL;
  status = RETURN_OK;

  QMALLOC(tmpline, PIXTYPE, bkmap->w, status);

  status = sep_backline_flt(bkmap, y, tmpline);
  if (status != RETURN_OK)
    goto exit;

  status = get_array_subtractor(dtype, &subtract_array, &size);
  if (status != RETURN_OK)
    goto exit;

  subtract_array(tmpline, bkmap->w, line);

 exit:
  free(tmpline);
  return status;
}

int sep_subbackarray(sepbackmap *bkmap, void *arr, int dtype)
{
  array_writer subtract_array;
  int y, status, size, width;
  PIXTYPE *tmpline;
  BYTE *arrt;

  tmpline = NULL;
  status = RETURN_OK;
  width = bkmap->w;
  arrt = (BYTE *)arr;

  QMALLOC(tmpline, PIXTYPE, width, status);

  status = get_array_subtractor(dtype, &subtract_array, &size);
  if (status != RETURN_OK)
    goto exit;

  for (y=0; y<bkmap->h; y++, arrt+=(width*size))
    { 
      if ((status = sep_backline_flt(bkmap, y, tmpline)) != RETURN_OK)
	goto exit;
      subtract_array(tmpline, width, arrt);
    }

 exit:
  free(tmpline);
  return status;
}

/*****************************************************************************/

void sep_freeback(sepbackmap *bkmap)
{
  if (bkmap)
    {
      free(bkmap->back);
      free(bkmap->dback);
      free(bkmap->sigma);
      free(bkmap->dsigma);
    }
  free(bkmap);
  
  return;
}
