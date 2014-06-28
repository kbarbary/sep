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

float *uniformf(float a, float b, int n)
/* an array of n random numbers from the uniform interval (a, b) */
{
  int i;
  float *result;
    
  result = (float*)malloc(n*sizeof(float));
  for (i=0; i<n; i++)
      result[i] = rand() / ((double)RAND_MAX);

  return result;
}

double gaussrand()
{
  static double V1, V2, S;
  static int phase = 0;
  double x, u1, u2;
  
  if (phase == 0)
    {
      do
	{
	  u1 = (double)rand() / RAND_MAX;
	  u2 = (double)rand() / RAND_MAX;
	  
	  V1 = 2.*u1 - 1.;
	  V2 = 2.*u2 - 1.;
	  S = V1*V1 + V2*V2;
	} 
      while (S >= 1. || S == 0.);
      
      x = V1 * sqrt(-2. * log(S) / S);
    }
  else
    {
      x = V2 * sqrt(-2. * log(S) / S);
    }

  phase = 1 - phase;
  
  return x;
}

float *makenoiseim(int nx, int ny, float mean, float sigma)
{
  int i, npix;
  float *im, *imt;

  im = (float*)malloc((npix=nx*ny)*sizeof(float));

  /* fill with noise */
  imt = im;
  for (i=0; i<npix; i++, imt++)
    *imt = mean + sigma*(float)gaussrand();

  return im;
}

float *ones(int nx, int ny)
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

}

void printbox(float *im, int w, int h, int xmin, int xmax, int ymin, int ymax)
/* print image values to the screen in a grid

   Don't make box size too big!!
*/
{
  int i, j, npix;

  npix = w*h;

  for (j = ymin; j < ymax && j < h; j++)
    {
      for (i = xmin; i < xmax && i < w; i++)
	printf("%6.2f ", im[w*j+i]);
      printf("\n");
    }
}



int main(int argc, char **argv)
{
  int i, status, nx, ny;
  float trueback, truesigma, xc, yc, r;
  double flux, fluxerr;
  short flag;
  float *im;
  uint64_t t0, t1;
  sepbackmap *bkmap = NULL;
  float conv[] = {1,2,1, 2,4,2, 1,2,1};
  int nobj = 0;
  int nsources = 0;
  sepobj *objects = NULL;

  status = 0;

  nx = 2048;
  ny = 4096;
  trueback = 10.;
  truesigma = 1.;

  printf("sep version: %s\n", sep_version_string);

  /* create a blank image with noise */
  im = makenoiseim(nx, ny, trueback, truesigma);

  /* add fake sources */
  for (yc = 100.; yc < ny - 100.; yc += 100.)
    for (xc = 100.; xc < nx - 100.; xc += 100.)
      {
	nsources += 1;
	addbox(im, nx, ny, xc, yc, 5., 10.*truesigma);
      }

  /***************************************************************************/
  /* Background map */

  printf("sep_makeback()            ");
  t0 = gettime_ns();
  status = sep_makeback(im, NULL, nx, ny, 64, 64, 0.0, 3, 3, 0.0, &bkmap);
  t1 = gettime_ns();
  if (status)
    goto exit;
  printf("%6.1f ms\n", (double)(t1 - t0)/1000000.);

  printf("sep_subbackarray()        ");
  t0 = gettime_ns();
  status = sep_subbackarray(bkmap, im);
  t1 = gettime_ns();
  if (status)
    goto exit;
  printf("%6.1f ms\n", (double)(t1 - t0)/1000000.);

  /***************************************************************************/
  /* extract */

  printf("sep_extract()             ");
  t0 = gettime_ns();
  status = sep_extract(im, NULL, nx, ny, 1.5*bkmap->backsig, 5,
		       conv, 3, 3, 32, 0.005, 1, 1.0, &nobj, &objects);
  t1 = gettime_ns();
  if (status)
    goto exit;
  printf("%6.1f ms", (double)(t1 - t0)/1000000.);
  printf("  [%d input, %d output]\n", nsources, nobj);


  free(im);

  /***************************************************************************/
  /* aperture photometry */
  int naper, j;
  float *xcs, *ycs;

  im = ones(nx, ny);
  naper = 1000;
  flux = fluxerr = 0.0;
  flag = 0;

  float rs[] = {3., 5., 10., 20.};
  for (j=0; j<4; j++)
    {
      r = rs[j];
      xcs = uniformf(2.*r, nx - 2.*r, naper);
      ycs = uniformf(2.*r, ny - 2.*r, naper);

      printf("sep_apercirc() [r=%4.1f]   ", r);
      t0 = gettime_ns();
      for (i=0; i<naper; i++)
	sep_apercirc(im, NULL, nx, ny, 0.0, 0.0,
		     xcs[i], ycs[i], r, 5, &flux, &fluxerr, &flag);
      t1 = gettime_ns();
      printf("%6.3f us/aperture\n", (double)(t1 - t0) / 1000. / naper);
      free(xcs);
      free(ycs);
    }


  /***************************************************************************/
  /* clean-up */

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
