#include <stdio.h>
#include <stdlib.h>
#include <fitsio.h>
#include "sep.h"

int main(int argc, char **argv)
{
  char *fname;
  float *im;
  fitsfile *file;
  int status;
  /*backspline *bkspl;*/
  int i, npix;
  long fpixel[2] = {1, 1};
  long naxes[2];

  /* Get filename */
  if (argc != 2)
    {
      puts("Usage: runtest FILENAME");
      exit(1);
    }
  fname = argv[1];

  /* Read in image */
  puts(fname);
  fits_open_file(&file, fname, READONLY, &status); 
  puts("opened file");
  fits_get_img_size(file, 2, naxes, &status);
  npix = naxes[0]*naxes[1];
  printf("npix = %d\n", npix);
  im = (float*)malloc(npix*sizeof(float));
  fits_read_pix(file, FLOAT_IMG, fpixel, npix, 0, im, 0, &status);

  /*for (i=0; i<n; i++)
    im[i] = (float)rand()/RAND_MAX; */
  for (i=0; i<20; i++)
      printf("%f  ", im[i]);

  /*puts("Before makeback...");
  bkspl = makeback(im, NULL, 100, 100, 20, 20, 0.0, 3, 3, 0.0);
  puts("After makeback...");
  freeback(bkspl); */

  fits_close_file(file, &status);

  return 0;
}
