
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sep.h"

/* weight >= wthresh implies that pixel will be used. */
/* w, h is image size in pixels */
/* bw, bh is size of a single background tile in pixels */
backsplstruct *makeback(PIXTYPE *im, PIXTYPE *weight,
                        int w, int h, int bw, int bh, PIXTYPE wthresh)
{
  int npix;                   /* size of image */
  int nx, ny, nb;             /* number of background boxes in x, y, total */
  int bufsize;                /* size of a "row" of boxes in pixels (w*bh) */
  backstruct *backmesh, *bm;  /* info about each background "box" */
  backsplstruct *backspline;  /* output */
  int j,k,m;

  npix = w*h;
  bufsize = w*bh;

  /* determine number of background boxes */
  if ((nx = (w-1)/bw + 1) < 1)
    nx = 1;
  if ((ny = (h-1)/bh + 1) < 1)
    ny = 1;
  nb = nx*ny;

  /* Allocate memory */
  QMALLOC(backmesh, backstruct, nx);		/* background information */

  /* Allocate the returned struct */
  QMALLOC(backspline, backsplstruct, 1);
  backspline->imnx = w;
  backspline->imny = h;
  backspline->nx = nx;
  backspline->ny = ny;
  backspline->n = nx*ny;
  backspline->bw = bw;
  backspline->bh = bh;
  QMALLOC(backspline->back, float, nb);
  QMALLOC(backspline->sigma, float, nb);

  /* loop over rows of background boxes.
     (here, we could loop over all boxes, but this is how its done in
     SExtractor, because the pixel buffers are only read in from disk
     in increments of a row of background boxes at a time.) */
  for (j=0; j<ny; j++, im+=bufsize, weight+=bufsize)
    {
      /* if the last row, modify the width appropriately*/
      if (j == ny-1 && npix%bufsize)
        bufsize = npix%bufsize;

      /* Get clipped mean, sigma for all boxes in the row */
      backstat(backmesh, im, weight, bufsize, nx, w, bw, weight?wthresh:0.0);

      /* Allocate histograms in each box in this row. */
      bm = backmesh;
      for (m=nb; m--; bm++)
	if (bm->mean <= -BIG)
	  bm->histo=NULL;
	else
	  QCALLOC(bm->histo, LONG, bm->nlevels);
      backhisto(backmesh, im, weight, bufsize, nx, w, bw, weight?wthresh:0.0);

      /*-- Compute background statistics from the histograms */
      bm = backmesh;
      for (m=0; m<nx; m++, bm++)
	{
	  k = m+nx*j;
	  backguess(bm, backspline->back+k, backspline->sigma+k);
	  free(bm->histo);
	}

    }

  /* free memory */
  free(backmesh);

  /* Median-filter and check suitability of the background map */
  filterback(backspline, 3, 0.0);

  /* Compute 2nd derivatives along the y-direction */
  backspline->dback = makebackspline(backspline, backspline->back);
  backspline->dsigma = makebackspline(backspline, backspline->sigma);

  return backspline;
}


/******************************** backstat **********************************/
/*
Compute robust statistical estimators in a row of meshes.
*/
void backstat(backstruct *backmesh,
	      PIXTYPE *buf, PIXTYPE *wbuf, int bufsize,
	      int n, int w, int bw, PIXTYPE wthresh)
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
		if ((wpix = *(wbuft++)) < wthresh && pix > -BIG)
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
		if ((wpix = *(wbuft++))<wthresh && pix<=hcut && pix>=lcut)
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
	       int n, int w, int bw, PIXTYPE wthresh)
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
	if (wbuf)
	  wbuf += bw;
	continue;

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
		if ((wpix = *(wbuft++))<wthresh && bin<nlevels && bin>=0)
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

      med = hihigh>=histo?
	((hihigh-histo)+0.5+((double)highsum-lowsum)/(2.0*(*hilow>*hihigh?
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

/******************************* filterback *********************************/
/*
Median filtering of the background map to remove the contribution from bright
sources.
*/
void	filterback(backsplstruct *backspline, int filtersize,
		   float filterthresh)
{
   float	*back,*sigma, *back2,*sigma2, *bmask,*smask, *sigmat,
		d2,d2min, fthresh, med, val,sval;
   int		i,j,px,py, np, nx,ny, npx,npx2, npy,npy2, dpx,dpy, x,y, nmin;

  fthresh = filterthresh;
  nx = backspline->nx;
  ny = backspline->ny;
  np = backspline->n;
  npx = filtersize/2;
  npy = filtersize/2;
  npy *= nx;

  QMALLOC(bmask, float, (2*npx+1)*(2*npy+1));
  QMALLOC(smask, float, (2*npx+1)*(2*npy+1));
  QMALLOC(back2, float, np);
  QMALLOC(sigma2, float, np);

  back = backspline->back;
  sigma = backspline->sigma;
  val = sval = 0.0;			/* to avoid gcc -Wall warnings */

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
  memcpy(back, back2, np*sizeof(float));
  backspline->backmean = fqmedian(back2, np);
  free(back2);
  memcpy(sigma, sigma2, np*sizeof(float));
  backspline->backsig = fqmedian(sigma2, np);

  if (backspline->backsig<=0.0)
    {
    sigmat = sigma2+np;
    for (i=np; i-- && *(--sigmat)>0.0;);
    if (i>=0 && i<(np-1))
      backspline->backsig = fqmedian(sigmat+1, np-1-i);
    else
      backspline->backsig = 1.0;
    }

  free(sigma2);


  return;
  }


/******************************* makebackspline ******************************/
/*
Pre-compute 2nd derivatives along the y direction at background nodes.
*/
float *makebackspline(backsplstruct *backspline, float *map)

  {
   int		x,y, nbx,nby,nbym1;
   float	*dmap,*dmapt,*mapt, *u, temp;

  nbx = backspline->nx;
  nby = backspline->ny;
  nbym1 = nby - 1;
  QMALLOC(dmap, float, backspline->n);
  for (x=0; x<nbx; x++)
    {
    mapt = map+x;
    dmapt = dmap+x;
    if (nby>1)
      {
      QMALLOC(u, float, nbym1);	/* temporary array */
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
      }
    else
      *dmapt = 0.0;
    }

  return dmap;
  }
