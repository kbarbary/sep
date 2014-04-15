#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/time.h>
#include <fitsio.h>
#include "sep.h"

uint64_t gettime_ns()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (uint64_t)tv.tv_sec * 1000000000ULL + tv.tv_usec * 1000ULL;
}

int main(int argc, char **argv)
{
  char fname[] = "DECam_00179121_24.fits[0]";
  float *im;
  fitsfile *file;
  int status=0;
  int npix;
  uint64_t t0, t1;
  long fpixel[2] = {1, 1};
  long naxes[2];
  int *anynull=0;

  /* Get filename
  if (argc != 2)
    {
      puts("Usage: runtest FILENAME");
      exit(1);
    }
  fname = argv[1];
  */

  /* Read in image */
  puts(fname);
  fits_open_file(&file, fname, READONLY, &status); 
  puts("opened file");
  fits_get_img_size(file, 2, naxes, &status);
  npix = naxes[0]*naxes[1];
  printf("image size = %ld, %ld\n", naxes[0], naxes[1]);
  im = (float*)malloc(npix*sizeof(float));
  ffgpxv(file, TFLOAT, fpixel, npix, 0, im, anynull, &status);
  fits_close_file(file, &status);
  if (status)
    goto exit;

  printf("Running makeback... ");
  backmap *bkmap;
  t0 = gettime_ns();
  bkmap = makeback(im, NULL, naxes[0], naxes[1], 64, 64, 0.0, 3, 3, 0.0,
		   &status);
  if (status)
    goto exit;
  t1 = gettime_ns();
  printf("done in %.1f ms.\n", (double)(t1 - t0)/1000000.);

  /* write background image out */
  float *bkim;
  bkim = (float*)malloc(npix*sizeof(float));

  printf("evaluting background map...");
  t0 = gettime_ns();
  status = backim(bkmap, bkim);
  if (status)
    goto exit;
  t1 = gettime_ns();
  printf("done in %.1f ms.\n", (double)(t1 - t0)/1000000.);


  /* skip this for now */
  if (0)
    {
      printf("writing to file: sepback.fits...");
      fitsfile *f;
      ffinit(&f, "!sepback.fits", &status); /* open new image */
      ffcrim(f, FLOAT_IMG, 2, naxes, &status); /* create image extension */
      ffppx(f, TFLOAT, fpixel, npix, bkim, &status); 
      fits_close_file(f, &status);
    }
  free(bkim);
  bkim = NULL;

  /* subtract background map from image data */
  int x, y;
  float *imsub = (float*)malloc(npix*sizeof(float));
  float *buf = (float*)malloc(naxes[0]*sizeof(float));
  float *imsubt = imsub;
  float *imt = im;
  float *buft = buf;
  for (y=0; y<naxes[1]; y++)
    {
      backline(bkmap, y, buf);
      for (x=0, buft=buf; x<naxes[0]; x++, imsubt++, buft++, imt++)
	{
	  *imsubt = *imt - *buft;
	}
    }

  float *bkgrms = (float*)malloc(npix*sizeof(float));
  backrmsim(bkmap, bkgrms);
  for (x=0; x<npix; x++)
    bkgrms[x] = bkgrms[x]*bkgrms[x];

  freeback(bkmap);

  /* find objects */
  float conv[] = {1,2,1, 2,4,2, 1,2,1};
  objliststruct *catalog = NULL;

  printf("extracting...");
  t0 = gettime_ns();
  status = extract(imsub, bkgrms, naxes[0], naxes[1], 1.5, 1.5, 0.0,
		   0, 5, conv, 3, 3, 32, 0.005, 1, 1.0, &catalog);
  if (status)
    goto exit;
  t1 = gettime_ns();
  printf("done in %.1f ms.\n", (double)(t1 - t0)/1000000.);
  printf("%d objects\n", catalog->nobj);


 exit:
  printf("status: %d\n", status);
  puts(errdetail);
  return status;
}
