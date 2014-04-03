/*
*				clean.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include	<math.h>
#include	<stdlib.h>
#include	<string.h>
#include	"sep.h"

/*------------------------------- variables ---------------------------------*/

static LONG		*cleanvictim;

/******************************* initclean **********************************
initialize cleanvictim and cleanobjlist.
 ***/
void	initclean(void)
  {
  if (CLEAN_FLAG)
    QMALLOC(cleanvictim, LONG, CLEAN_STACKSIZE);
  QMALLOC(cleanobjlist, objliststruct, 1);
  cleanobjlist->obj = NULL;
  cleanobjlist->plist = NULL;
  cleanobjlist->nobj = cleanobjlist->npix = 0;

  return;
  }

