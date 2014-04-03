/*
*				extract.c
*
* Extract connected pixels from an image raster.
*
*/

#include	<math.h>
#include	<stdlib.h>
#include	<string.h>
#include	"sep.h"

/*------------------------- Static buffers for lutz() -----------------------*/

static infostruct	*info, *store;
static char		*marker;
static status		*psstack;
static int		*start, *end, *discan, xmin,ymin,xmax,ymax;


/******************************* lutzalloc ***********************************/
/*
Allocate once for all memory space for buffers used by lutz().
*/
void	lutzalloc(int width, int height)
  {
   int	*discant,
	stacksize, i;

  stacksize = width+1;
  xmin = ymin = 0;
  xmax = width-1;
  ymax = height-1;
  QMALLOC(info, infostruct, stacksize);
  QMALLOC(store, infostruct, stacksize);
  QMALLOC(marker, char, stacksize);
  QMALLOC(psstack, status, stacksize);
  QMALLOC(start, int, stacksize);
  QMALLOC(end, int, stacksize);
  QMALLOC(discan, int, stacksize);
  discant = discan;
  for (i=stacksize; i--;)
    *(discant++) = -1;

  return;
  }
