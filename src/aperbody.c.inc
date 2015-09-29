/*
	Adding (void *) pointers is a GNU C extension, not part of standard C. 
	When compiling on Windows with MS Visual C compiler need to cast the
	(void *) to something the size of one byte.
*/
#if defined(_MSC_VER)
	#define MSVC_VOID_CAST (char *)
#else
  	#define MSVC_VOID_CAST
#endif

int APER_NAME(void *data, void *error, void *mask,
	      int dtype, int edtype, int mdtype, int w, int h,
	      double maskthresh, double gain, short inflag,
	      double x, double y, APER_ARGS, int subpix,
	      double *sum, double *sumerr, double *area, short *flag)
{
  PIXTYPE pix, varpix;
  double dx, dy, dx1, dy2, offset, scale, scale2, tmp;
  double tv, sigtv, totarea, maskarea, overlap, rpix2;
  int ix, iy, xmin, xmax, ymin, ymax, sx, sy, status, size, esize, msize;
  long pos;
  short errisarray, errisstd;
  BYTE *datat, *errort, *maskt;
  converter convert, econvert, mconvert;
  APER_DECL;

  /* input checks */
  APER_CHECKS;
  if (subpix < 0)
    return ILLEGAL_SUBPIX;

  /* initializations */
  size = esize = msize = 0;
  tv = sigtv = 0.0;
  overlap = totarea = maskarea = 0.0;
  datat = maskt = NULL;
  errort = error;
  *flag = 0;
  varpix = 0.0;
  scale = 1.0/subpix;
  scale2 = scale*scale;
  offset = 0.5*(scale-1.0);

  APER_INIT;

  /* get data converter(s) for input array(s) */
  if ((status = get_converter(dtype, &convert, &size)))
    return status;
  if (error && (status = get_converter(edtype, &econvert, &esize)))
    return status;
  if (mask && (status = get_converter(mdtype, &mconvert, &msize)))
    return status;

  /* get options */
  errisarray = inflag & SEP_ERROR_IS_ARRAY;
  if (!error)
    errisarray = 0; /* in case user set flag but error is NULL */
  errisstd = !(inflag & SEP_ERROR_IS_VAR);

  /* If error exists and is scalar, set the pixel variance now */
  if (error && !errisarray)
    {
      varpix = econvert(errort);
      if (errisstd)
	varpix *= varpix;
    }

  /* get extent of box */
  APER_BOXEXTENT;

  /* loop over rows in the box */
  for (iy=ymin; iy<ymax; iy++)
    {
      /* set pointers to the start of this row */
      pos = (iy%h) * w + xmin;
      datat = MSVC_VOID_CAST data + pos*size;
      if (errisarray)
	errort = MSVC_VOID_CAST error + pos*esize;
      if (mask)
	maskt = MSVC_VOID_CAST mask + pos*msize;

      /* loop over pixels in this row */
      for (ix=xmin; ix<xmax; ix++)
	{
	  dx = ix - x;
	  dy = iy - y;
	  rpix2 = APER_RPIX2;
	  if (APER_COMPARE1)
	    {
	      if (APER_COMPARE2)  /* might be partially in aperture */
		{
		  if (subpix == 0)
		    overlap = APER_EXACT;
		  else
		    {
		      dx += offset;
		      dy += offset;
		      overlap = 0.0;
		      for (sy=subpix; sy--; dy+=scale)
			{
			  dx1 = dx;
			  dy2 = dy*dy;
			  for (sx=subpix; sx--; dx1+=scale)
			    {
			      rpix2 = APER_RPIX2_SUBPIX;
			      if (APER_COMPARE3)
				overlap += scale2;
			    }
			}
		    }
		}
	      else
		/* definitely fully in aperture */
		overlap = 1.0;
	      
	      pix = convert(datat);

	      if (errisarray)
		{
		  varpix = econvert(errort);
		  if (errisstd)
		    varpix *= varpix;
		}

	      if (mask && (mconvert(maskt) > maskthresh))
		{ 
		  *flag |= SEP_APER_HASMASKED;
		  maskarea += overlap;
		}
	      else
		{
		  tv += pix*overlap;
		  sigtv += varpix*overlap;
		}

	      totarea += overlap;

	    } /* closes "if pixel might be within aperture" */
	  
	  /* increment pointers by one element */
	  datat += size;
	  if (errisarray)
	    errort += esize;
	  maskt += msize;
	}
    }

  /* correct for masked values */
  if (mask)
    {
      if (inflag & SEP_MASK_IGNORE)
	totarea -= maskarea;
      else
	{
	  tv *= (tmp = totarea/(totarea-maskarea));
	  sigtv *= tmp;
	}
    }

  /* add poisson noise, only if gain > 0 */
  if (gain > 0.0 && tv>0.0)
    sigtv += tv/gain;

  *sum = tv;
  *sumerr = sqrt(sigtv);
  *area = totarea;

  return status;
}
