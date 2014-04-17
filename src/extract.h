/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
* This file is part of SEP
*
* Copyright 2014 Kyle Barbary
*
* SEP is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* SEP is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with SEP.  If not, see <http://www.gnu.org/licenses/>.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
* This file part of: SExtractor
*
* Copyright:         (C) 1993-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
*
* License:           GNU General Public License
*
* SExtractor is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* SExtractor is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* You should have received a copy of the GNU General Public License
* along with SExtractor. If not, see <http://www.gnu.org/licenses/>.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#define	UNKNOWN	        -1  /* flag for LUTZ */
#define	ANALYSE_FAST	 0  /* flags for preanalyse */
#define	ANALYSE_FULL	 1
#define	ANALYSE_ROBUST	 2
#define	CLEAN_ZONE      10.0  /* zone (in sigma) to consider for processing */
#define CLEAN_STACKSIZE 3000  /* replaces prefs.clean_stacksize  */
                              /* (MEMORY_OBJSTACK in sextractor inputs) */
#define CLEAN_MARGIN    0  /* replaces prefs.cleanmargin which was set based */
                           /* on stuff like apertures and vignet size */
#define	MARGIN_SCALE   2.0 /* Margin / object height */ 
#define	MARGIN_OFFSET  4.0 /* Margin offset (pixels) */ 
#define	MAXDEBAREA     3   /* max. area for deblending (must be >= 1)*/
#define	MAXFLAG	       4   /* max. # of FLAG-images (TODO: remove)*/
#define	MAXPICSIZE     1048576 /* max. image size in any dimension */

/* plist-related macros */
#define	PLIST(ptr, elem)	(((pbliststruct *)(ptr))->elem)
#define	PLISTEXIST(elem)	(plistexist_##elem)
#define	PLISTPIX(ptr, elem)	(*((PIXTYPE *)((ptr)+plistoff_##elem)))
#define	PLISTFLAG(ptr, elem)	(*((FLAGTYPE *)((ptr)+plistoff_##elem)))

/* Extraction status */
typedef	enum {COMPLETE, INCOMPLETE, NONOBJECT, OBJECT} pixstatus;

/* Temporary object parameters during extraction */
typedef struct structinfo
{
  LONG	pixnb;	    /* Number of pixels included */
  LONG	firstpix;   /* Pointer to first pixel of pixlist */
  LONG	lastpix;    /* Pointer to last pixel of pixlist */
  short	flag;	    /* Extraction flag */
} infostruct;

typedef struct
{
  int     nextpix;
  int     x, y;
  PIXTYPE value;
} pbliststruct;

/* globals */
int plistexist_value, plistexist_dvalue, plistexist_cdvalue,
  plistexist_flag, plistexist_wflag, plistexist_dthresh, plistexist_var,
  plistoff_value, plistoff_dvalue, plistoff_cdvalue,
  plistoff_flag[MAXFLAG], plistoff_wflag, plistoff_dthresh, plistoff_var,
  plistsize;

int objmthresh(int objnb, objliststruct *objlist, int minarea,
	       PIXTYPE dthresh);
void preanalyse(int, objliststruct *, int);
int  addcleanobj(objstruct *, objliststruct *);
int  subcleanobj(int, objliststruct *);
int  clean(int, objliststruct *, objliststruct *, double, int *);
int  lutzalloc(int, int);
void lutzfree(void);
int  lutz(objliststruct *, int, objstruct *, objliststruct *, int);
void update(infostruct *, infostruct *, pliststruct *);
int  allocparcelout(int);
void freeparcelout(void);
int  parcelout(objliststruct *, objliststruct *, int, double, int);
