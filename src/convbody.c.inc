
/* Convolve one line of an image with a given kernel.
 *
 * image : full input array
 * w, h : width and height of image
 * y : line to convolve in image
 * conv : convolution kernel
 * convw, convh : width and height of conv
 * buf : output convolved line (`w` elements long)
 */
void CONVOLVE_FN(void *image, int w, int h, int y,
	         float *conv, int convw, int convh, PIXTYPE *buf)
{
  int	  convw2, cx, dcx, y0, dy;
  float	  *convend;
  PIXTYPE *bufend, *buft, *buftend, cval;

  CONVOLVE_TYPE *im, *line, *imt;
  im = (CONVOLVE_TYPE *)image;

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
