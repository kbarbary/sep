#include <stdio.h>
#include <stdlib.h>
#include "sep.h"

int main(int argc, char **argv)
{
  float *im;
  backspline *bkspl;
  int i;
  int n = 10000;
  im = (float*)calloc(n, sizeof(float));
  for (i=0; i<n; i++)
    im[i] = (float)rand()/RAND_MAX;
  /*  for (i=0; i<20; i++)
      printf("%f  ", im[i]); */

  puts("Before makeback...");
  bkspl = makeback(im, NULL, 100, 100, 20, 20, 0.0, 3, 3, 0.0);
  puts("After makeback...");
  freeback(bkspl);

  return 0;
}
