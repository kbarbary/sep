
/* Apply a matched filter to one line of an image with a given kernel.
 *
 * image : full input array
 * noise : full input noise array
 * w, h : width and height of image
 * y : line to apply the matched filter to in an image
 * conv : convolution kernel
 * convw, convh : width and height of conv
 * num_buf : output numerator line (`w` elements long)
 * denom_buf : output denominator line (`w` elements long)
 */
void MATCHED_FILTER_FN(void *image, void *noise, int w, int h, int y,
                       float *conv, int convw, int convh, PIXTYPE *num_buf,
                       PIXTYPE *denom_buf)
{
  int	  convw2, cx, dcx, y0, dy;
  float	  *convend;
  PIXTYPE *num_bufend, *num_buft, *num_buftend, *denom_bufend, *denom_buft,
          *denom_buftend, cval;

  MATCHED_FILTER_IMAGE_TYPE image_val;
  MATCHED_FILTER_NOISE_TYPE noise_val;

  MATCHED_FILTER_IMAGE_TYPE *image_buf, *image_line, *image_target;
  MATCHED_FILTER_NOISE_TYPE *noise_buf, *noise_line, *noise_target;

  image_buf = (MATCHED_FILTER_IMAGE_TYPE *)image;
  noise_buf = (MATCHED_FILTER_NOISE_TYPE *)noise;

  image_line = NULL;              /* To avoid gcc -Wall warnings */
  noise_line = NULL;              /* To avoid gcc -Wall warnings */
  image_val = 0;                  /* To avoid gcc -Wall warnings */
  noise_val = 0;                  /* To avoid gcc -Wall warnings */
  convw2 = convw/2;
  num_bufend = num_buf+w;         /* pointer to end of output buffer */
  denom_bufend = denom_buf+w;     /* pointer to end of output buffer */
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

  /* initialize output buffers to zero */
  memset(num_buf, 0, w*sizeof(PIXTYPE));
  memset(denom_buf, 0, w*sizeof(PIXTYPE));

  /* loop over pixels in the mask */
  for (cx=0; conv<convend; conv++, cx++)
    {
      cval = *conv;

      /* get the x position in the mask */
      if (cx==convw)
	cx = 0;

      /* when cx goes to zero, increment the start line in the image */
      if (!cx)
        {
          image_line = image_buf + w*(y0%h);
          noise_line = noise_buf + w*(y0%h);
          y0++;
        }

      /* get start and end positions in the source and target line */
      if ((dcx = cx-convw2)>=0)
	{
	  image_target = image_line + dcx;
	  noise_target = noise_line + dcx;
	  num_buft = num_buf;
	  num_buftend = num_bufend - dcx;
	  denom_buft = denom_buf;
	  denom_buftend = denom_bufend - dcx;
	}
      else
	{
	  image_target = image_line;
	  noise_target = noise_line;
	  num_buft = num_buf - dcx;
	  num_buftend = num_bufend;
	  denom_buft = denom_buf - dcx;
	  denom_buftend = denom_bufend;
	}

      while (num_buft < num_buftend)
        {
          image_val = *(image_target++);
          noise_val = *(noise_target++);
          if (noise_val != 0)
            {
              *num_buft += cval * image_val / noise_val / noise_val;
              *denom_buft += cval * cval / noise_val / noise_val;
            }
          
          num_buft++;
          denom_buft++;
        }
    }

  /* take the square root of the denominator */
  denom_buft = denom_buf;
  while (denom_buft < denom_bufend)
    {
      *denom_buft = sqrt(*denom_buft);
      denom_buft++;
    }

  return;
}
