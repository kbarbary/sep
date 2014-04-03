/*
*				refine.c
*
* Deblend sources based on their pixel lists.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include	<math.h>
#include	<stdlib.h>
#include	<string.h>

#ifndef	RAND_MAX
#define	RAND_MAX	2147483647
#endif
#define	NSONMAX			1024	/* max. number per level */
#define	NBRANCH			16	/* starting number per branch */

/* was in prefs... TODO make this an argument to the scan funciton */
#define DEBLEND_NTHRESH 32  /* sextractor default is 32 */


static objliststruct	*objlist;
static short		*son, *ok;



/******************************* allocparcelout ******************************/
/*
Allocate the memory allocated by global pointers in refine.c
*/
void	allocparcelout(void)
  {
  QMALLOC(son, short,  DEBLEND_NTHRESH*NSONMAX*NBRANCH);
  QMALLOC(ok, short,  DEBLEND_NTHRESH*NSONMAX);
  QMALLOC(objlist, objliststruct,  DEBLEND_NTHRESH);

  return;
  }
