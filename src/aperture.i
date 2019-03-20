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

int APER_NAME(sep_image *im,
	      double x, double y, APER_ARGS, int id, int subpix, short inflag,
	      double *sum, double *sumerr, double *area, short *flag)
{
  PIXTYPE pix, varpix;
  double dx, dy, dx1, dy2, offset, scale, scale2, tmp;
  double tv, sigtv, totarea, maskarea, overlap, rpix2;
  int ix, iy, xmin, xmax, ymin, ymax, sx, sy, status, size, esize, msize, ssize;
  int ismasked;
  long pos;
  short errisarray, errisstd;
  BYTE *datat, *errort, *maskt, *segt;
  converter convert, econvert, mconvert, sconvert;
  APER_DECL;

  /* input checks */
  APER_CHECKS;
  if (subpix < 0)
    return ILLEGAL_SUBPIX;

  /* initializations */
  size = esize = msize = ssize = 0;
  tv = sigtv = 0.0;
  overlap = totarea = maskarea = 0.0;
  datat = maskt = segt = NULL;
  errort = im->noise;
  *flag = 0;
  varpix = 0.0;
  scale = 1.0/subpix;
  scale2 = scale*scale;
  offset = 0.5*(scale-1.0);
  errisarray = 0;
  errisstd = 0;

  APER_INIT;

  /* get data converter(s) for input array(s) */
  if ((status = get_converter(im->dtype, &convert, &size)))
    return status;
  if (im->mask && (status = get_converter(im->mdtype, &mconvert, &msize)))
    return status;

  if (im->segmap && (status = get_converter(im->sdtype, &sconvert, &ssize)))
    return status;
      
  /* get image noise */
  if (im->noise_type != SEP_NOISE_NONE)
    {
      errisstd = (im->noise_type == SEP_NOISE_STDDEV);
      if (im->noise)
        {
          errisarray = 1;
          if ((status = get_converter(im->ndtype, &econvert, &esize)))
             return status;
        }
      else
        {
          varpix = (errisstd)?  im->noiseval * im->noiseval: im->noiseval;
        }
    }

  /* get extent of box */
  APER_BOXEXTENT;
  
  /* loop over rows in the box */
  for (iy=ymin; iy<ymax; iy++)
    {
      /* set pointers to the start of this row */
      pos = (iy%im->h) * im->w + xmin;
      datat = MSVC_VOID_CAST im->data + pos*size;
      if (errisarray)
	errort = MSVC_VOID_CAST im->noise + pos*esize;
      if (im->mask)
	maskt = MSVC_VOID_CAST im->mask + pos*msize;
      if (im->segmap)
  	segt = MSVC_VOID_CAST im->segmap + pos*ssize;
  	
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
              
              ismasked = 0;
	      if (im->mask && (mconvert(maskt) > im->maskthresh))
	        {
	          ismasked = 1;
	        }
	      
	      /* Segmentation image:  
	           
	           If `id` is negative, require segmented pixels within the 
	           aperture.
	           
	           If `id` is positive, mask pixels with nonzero segment ids
	           not equal to `id`.
	           
	      */ 
	      if (im->segmap)
  	        {
  	          if (id > 0) 
  	            {
  	              if ((sconvert(segt) > 0.) & (sconvert(segt) != id))
  	                {
  	                  ismasked = 1;
  	                }
  	            } else {
	              if (sconvert(segt) != -1*id)
	                {
	                  ismasked = 1;
	                }  	            
  	            }
  	        }
  	      
	      if (ismasked > 0) 
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
	  segt += ssize;
	}
    }

  /* correct for masked values */
  if (im->mask)
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
  if (im->gain > 0.0 && tv>0.0)
    sigtv += tv / im->gain;

  *sum = tv;
  *sumerr = sqrt(sigtv);
  *area = totarea;

  return status;
}
