/*
*				plist.c
*
* Manage pixel lists.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <stdio.h>
#include <stdlib.h>
#include "sep.h"
#include "defs.h"
#include "plist.h"

/****************************** init_plist ************************************
PURPOSE	initialize a pixel-list and its components.
 ***/
void init_plist(PIXTYPE *conv, PIXTYPE *cdwfield)
{
  pbliststruct	*pbdum = NULL;
  int		i;

  plistsize = sizeof(pbliststruct);
  plistoff_value = (char *)&pbdum->value - (char *)pbdum;

  if (conv)
    {
      plistexist_cdvalue = 1;
      plistoff_cdvalue = plistsize;
      plistsize += sizeof(PIXTYPE);
    }
  else
    {
      plistexist_cdvalue = 0;
      plistoff_cdvalue = plistoff_dvalue;
    }

  if (cdwfield)
    {
      plistexist_var = 1;
      plistoff_var = plistsize;
      plistsize += sizeof(PIXTYPE);
    }
  else
    plistexist_var = 0;

  if (cdwfield)
    {
      plistexist_dthresh = 1;
      plistoff_dthresh = plistsize;
      plistsize += sizeof(PIXTYPE);
    }
  else
    plistexist_dthresh = 0;

  return;

  /* can we remove these? */
  plistexist_dvalue = 0;
  plistoff_dvalue = plistoff_value;
  plistexist_flag = 0;
  plistexist_wflag = 0;

  return;
}
