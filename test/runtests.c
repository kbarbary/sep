#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/time.h>
#include "sep.h"

uint64_t gettime_ns()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (uint64_t)tv.tv_sec * 1000000000ULL + tv.tv_usec * 1000ULL;
}


int main(int argc, char **argv)
{
  int status, npix, nx, ny, i;
  float *im;
  uint64_t t0, t1;
  backmap *bkmap = NULL;

  status = 0;

  nx = 2048;
  ny = 4096;

  npix = nx * ny;

  /* Create a "simulated" image. */
  im = (float*)malloc(npix*sizeof(float));

  /* fill with noise and "objects" */



  printf("Running makeback... ");
  t0 = gettime_ns();
  status = makeback(im, NULL, nx, ny, 64, 64, 0.0, 3, 3, 0.0, &bkmap);
  if (status)
    goto exit;
  t1 = gettime_ns();
  printf("done in %.1f ms.\n", (double)(t1 - t0)/1000000.);


  /* subtract background */
  status = subbackarray(bkmap, im);

  /* background variance */
  float *bkgvar = (float*)malloc(npix*sizeof(float));
  status = backvararray(bkmap, bkgvar);

  /* find objects */
  float conv[] = {1,2,1, 2,4,2, 1,2,1};
  objliststruct *catalog = NULL;

  printf("extracting...");
  t0 = gettime_ns();
  float thresh = 1.5 * bkmap->backsig;
  status = extract(im, NULL, nx, ny, thresh, 5,
		   conv, 3, 3, 32, 0.005, 1, 1.0, &catalog);
  if (status)
    goto exit;
  t1 = gettime_ns();
  printf("done in %.1f ms.\n", (double)(t1 - t0)/1000000.);
  printf("%d objects\n", catalog->nobj);


 exit:
  freeback(bkmap);
  free(im);
  if (status)
    {
      printf("FAILED with status: %d\n", status);
      puts(errdetail);
    }
  else
    {
      printf("Tests passed\n");
    }
  return status;
}
