#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <sys/time.h>
#include "sep.h"

uint64_t gettime_ns()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (uint64_t)tv.tv_sec * 1000000000ULL + tv.tv_usec * 1000ULL;
}

double gaussrand()
{
  static double V1, V2, S;
  static int phase = 0;
  double X;
  
  if (phase == 0)
    {
      do
	{
	  double U1 = (double)rand() / RAND_MAX;
	  double U2 = (double)rand() / RAND_MAX;
	  
	  V1 = 2 * U1 - 1;
	  V2 = 2 * U2 - 1;
	  S = V1 * V1 + V2 * V2;
	} 
      while (S >= 1 || S == 0);
      
      X = V1 * sqrt(-2 * log(S) / S);
    }
  else
    {
      X = V2 * sqrt(-2 * log(S) / S);
    }

  phase = 1 - phase;
  
  return X;
}

float *makenoiseim(nx, ny, mean, sigma)
{
  int i, npix;
  float *im, *imt;

  im = (float*)malloc((npix=nx*ny)*sizeof(float));

  /* fill with noise and "objects" */
  imt = im;
  for (i=0; i<npix; i++, imt++)
    *imt = mean + sigma*(float)gaussrand();

  return im;
}

float *ones(nx, ny)
{
  int i, npix;
  float *im, *imt;

  im = (float *)malloc((npix = nx*ny)*sizeof(float));
  imt = im;
  for (i=0; i<npix; i++, imt++)
    *imt = 1.0;

  return im;
}

void addbox(float *im, int w, int h, float xc, float yc, float r, float val)
/* n = sersic index */
{
  int xmin, xmax, ymin, ymax;
  int x, y;

  int rmax = (int)r;

  xmin = (int)xc - rmax;
  xmin = (xmin < 0)? 0: xmin;
  xmax = (int)xc + rmax;
  xmax = (xmax > w)? w: xmax;
  ymin = (int)yc - rmax;
  ymin = (ymin < 0)? 0: ymin;
  ymax = (int)yc + rmax;
  ymax = (ymax > h)? h: ymax;

  for (y=ymin; y<ymax; y++)
    for (x=xmin; x<xmax; x++)
      im[x+w*y] += val;

  return;
}

int main(int argc, char **argv)
{
  int status, nx, ny;
  float trueback, truesigma;
  double flux, fluxerr;
  short flag;
  float *im;
  uint64_t t0, t1;
  sepbackmap *bkmap = NULL;
  float conv[] = {1,2,1, 2,4,2, 1,2,1};
  int nobj = 0;
  sepobj *objects = NULL;

  status = 0;

  nx = 2048;
  ny = 4096;
  trueback = 10.;
  truesigma = 1.;

  printf("sep version: %s\n", sep_version_string);

  /* create a blank image */
  im = makenoiseim(nx, ny, trueback, truesigma);
  addbox(im, nx, ny, 1000., 1000., 5., 20.);

  printf("sep_makeback... ");
  t0 = gettime_ns();
  status = sep_makeback(im, NULL, nx, ny, 64, 64, 0.0, 3, 3, 0.0, &bkmap);
  if (status)
    goto exit;
  t1 = gettime_ns();
  printf("done in %.1f ms.\n", (double)(t1 - t0)/1000000.);

  printf("sep_subbackarray... ");
  t0 = gettime_ns();
  status = sep_subbackarray(bkmap, im);
  if (status)
    goto exit;
  t1 = gettime_ns();
  printf("done in %.1f ms.\n", (double)(t1 - t0)/1000000.);

  printf("sep_extract...");
  t0 = gettime_ns();
  status = sep_extract(im, NULL, nx, ny, 1.5*bkmap->backsig, 5,
		       conv, 3, 3, 32, 0.005, 1, 1.0, &nobj, &objects);
  if (status)
    goto exit;
  t1 = gettime_ns();
  printf("done in %.1f ms.\n", (double)(t1 - t0)/1000000.);

  printf("%d objects\n", nobj);

  free(im);

  /* aperture photometry */
  im = ones(nx, ny);
  
  flux = fluxerr = 0.0;
  flag = 0;
  sep_apercirc(im, NULL, nx, ny, 0.0, 0.0,
	       10.0, 10.0, 5.0, 5, &flux, &fluxerr, &flag);
  printf("flux = %.3f, fluxerr = %.3f\n", flux, fluxerr);

 exit:
  sep_freeback(bkmap);
  free(im);
  if (status)
    {
      printf("FAILED with status: %d\n", status);
      char errtext[512];
      sep_get_errdetail(errtext);
      puts(errtext);
    }
  else
    {
      printf("All tests pass.\n");
    }
  return status;
}
