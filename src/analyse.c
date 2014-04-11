#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sep.h"
#include "extract.h"

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
  PIXTYPE	peak, cpeak, val, cval;
  double	thresh,thresh2, t1t2,darea,
                mx,my, mx2,my2,mxy, rv, tv,
		xm,ym, xm2,ym2,xym,
		temp,temp2, theta,pmx2,pmy2;
  int		x, y, xmin,xmax, ymin,ymax,area2, fdnpix, dnpix;
  
  /*-----  initialize stacks and bounds */
  thresh = obj->dthresh;
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
	  double xn, yn;

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
