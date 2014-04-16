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

#define	SEP_VERSION  "0.1.dev"
#define	SEP_DATE     "2014-03-14"

/*------------------------- global typedefs ---------------------------------*/

typedef float PIXTYPE;   /* type of image arrays */

/*-------------------- error codes & messages -------------------------------*/
#define	RETURN_OK		0
#define	RETURN_ERROR		(-1) /* general unspecified error */
#define MEMORY_PIXSTACK_ERROR   1
#define PIXSTACK_OVERFLOW_ERROR 2
#define FATAL_ERROR             3
#define MEMORY_CLEAN_ERROR      4
#define NO_CLEAN_OBJ_ERROR      5 /* Internal Error: no CLEAN object to
                                     remove in subcleanobj()*/
#define LUTZ_REALLOC_ERROR      6 /* problem with mem. realloc. in lutz() */
#define GATHERUP_MEMORY_ERROR   7 /* Not enough memory to update pixel list in
				     gatherup()" */
#define MEMORY_ALLOC_ERROR      8 /* Could not allocate memory for.. */ 
#define DEBLEND_OVERFLOW_ERROR  9

/*------------------------------ background ---------------------------------*/
typedef struct
{
  int imnx, imny;          /* original image width, height */
  int bw, bh;              /* single tile width, height */
  int nx, ny;              /* number of tiles in x, y */
  int n;                   /* nx*ny */
  float backmean, backsig; /* global mean, sigma */
  float *back;             /* background map */
  float *dback;            /* background map */
  float *sigma;
  float *dsigma;
} backmap;

/* w, h is image size in pixels */
/* bw, bh is size of a single background tile in pixels */
/* var > varthresh will be ignored. */
int makeback(PIXTYPE *im, PIXTYPE *mask, int w, int h,
	     int bw, int bh, PIXTYPE maskthresh, int fbx, int fby,
	     float fthresh, backmap **bkm);
PIXTYPE	backpixlinear(backmap *, int, int);
int backline(backmap *, int, PIXTYPE *);
int backrmsline(backmap *, int, PIXTYPE *);
int backvarline(backmap *, int, PIXTYPE *);
int backarray(backmap *, PIXTYPE *);
int backrmsarray(backmap *, PIXTYPE *);
int backvararray(backmap *, PIXTYPE *);
int subbackline(backmap *, int, PIXTYPE *);
int subbackarray(backmap *, PIXTYPE *);
void freeback(backmap *);


/*-------------------------- source extraction ------------------------------*/

#define	OBJ_CROWDED   0x0001
#define	OBJ_MERGED    0x0002
#define	OBJ_SATUR     0x0004
#define	OBJ_TRUNC     0x0008
#define OBJ_APERT_PB  0x0010
#define	OBJ_ISO_PB    0x0020
#define	OBJ_DOVERFLOW 0x0040
#define	OBJ_OVERFLOW  0x0080
#define	NISO          8	     /* number of isophotes */

typedef	unsigned char	BYTE;  /* a byte */
typedef	char pliststruct;      /* Dummy type for plist */

typedef struct
{
  /* ---- basic parameters */
  int	   number;			/* ID */
  int	   fdnpix;		       	/* nb of extracted pix */
  int	   dnpix;	       		/* nb of pix above thresh  */
  int	   npix;       			/* "" in measured frame */
  int	   nzdwpix;			/* nb of zero-dweights around */
  int	   nzwpix;		       	/* nb of zero-weights inside */
  float	   fdflux;	       		/* integrated ext. flux */
  float	   dflux;      			/* integrated det. flux */
  float	   flux;       			/* integrated mes. flux */
  float	   fluxerr;			/* integrated variance */
  PIXTYPE  fdpeak;	       		/* peak intensity (ADU) */
  PIXTYPE  dpeak;      			/* peak intensity (ADU) */
  PIXTYPE  peak;       			/* peak intensity (ADU) */
  /* ---- astrometric data */
  int	   peakx,peaky;                       /* pos of brightest pix */
  double   mx, my;	       		      /* barycenter */
  double   poserr_mx2,poserr_my2,poserr_mxy;  /* Error ellips moments */
  /* ---- morphological data */			
  int	   xmin,xmax,ymin,ymax,ycmin,ycmax;  /* x,y limits */
  PIXTYPE  *blank, *dblank;                  /* BLANKing sub-images  */
  int	   *submap;                          /* Pixel-index sub-map */
  int	   subx,suby, subw,subh;	     /* sub-image pos. and size */
  short	   flag;			     /* extraction flags */
  /* BYTE	   wflag; */	       	     /* weighted extraction flags */
  /* FLAGTYPE imaflag[MAXFLAG]; */     	     /* flags from FLAG-images */
  BYTE	   singuflag;			     /* flags for singularities */
  /* int      imanflag[MAXFLAG]; */    	     /* number of MOST flags */
  double   mx2,my2,mxy;			     /* variances and covariance */
  float	   a, b, theta, abcor;		     /* moments and angle */
  float	   cxx,cyy,cxy;			     /* ellipse parameters */
  int	   firstpix;			     /* ptr to first pixel */
  int	   lastpix;			     /* ptr to last pixel */
  float	   bkg, dbkg, sigbkg, dsigbkg;	     /* Background stats (ADU) */
  float	   thresh;		             /* measur. threshold (ADU) */
  float	   dthresh;		       	     /* detect. threshold (ADU) */
  float	   mthresh;		             /* max. threshold (ADU) */
  int	   iso[NISO];			     /* isophotal areas */
  float	   fwhm;			     /* IMAGE FWHM */
} objstruct;

typedef struct
{
  int           nobj;	  /* number of objects in list */
  objstruct     *obj;	  /* pointer to the object array */
  int           npix;	  /* number of pixels in pixel-list */
  pliststruct   *plist;	  /* pointer to the pixel-list */
  PIXTYPE       dthresh;  /* detection threshold */
  PIXTYPE       thresh;	  /* analysis threshold */
} objliststruct;

int extract(PIXTYPE *cfield, PIXTYPE *cdwfield, int w, int h,
	    PIXTYPE dthresh, PIXTYPE athresh, PIXTYPE cdwthresh,
	    int threshabsolute, int minarea,
	    float *conv, int convw, int convh,
	    int deblend_nthresh, double deblend_mincont,
	    int clean_flag, double clean_param,
	    objliststruct **catalog);

/* sextractor defaults show in []                     */
/*----------------------------------------------------*/
/* image array                                        */
/* variance array (can be NULL)                       */
/* size of arrays (w, h)                              */
/* detection threshold                          [1.5] */
/* analysis threshold                           [1.5] */
/* cdwthresh (???)                                    */
/* threshabsolute (0=relative)                    [0] */
/* minarea                                        [5] */
/* conv array (can be NULL)     [{1 2 1 2 4 2 1 2 1}] */
/* conv size (w, h)                            [3, 3] */
/* deblend_nthresh                               [32] */
/* deblend_mincont                             [0.005] */
/* clean_flag (1 = YES)                           [1] */
/* clean_param                                  [1.0] */

/*-------------------- global internal definitions --------------------------*/
/*
 * These are not actually needed by external code, but included here for ease
 * of #include statements: not worth having a separate file.
 */

#define	BIG 1e+30  /* a huge number (< biggest value a float can store) */
#define	PI  3.1415926535898 /* never met before? */
#define	DEG (PI/180.0)	    /* 1 deg in radians */

typedef	int	      LONG;
typedef	unsigned int  ULONG;

char errdetail[512];

#define	QCALLOC(ptr, typ, nel, status)				     	\
  {if (!(ptr = (typ *)calloc((size_t)(nel),sizeof(typ))))		\
      {									\
	sprintf(errdetail, #ptr " (" #nel "=%lu elements) "		\
		"at line %d in module " __FILE__ " !",			\
		(size_t)(nel)*sizeof(typ), __LINE__);			\
	status = MEMORY_ALLOC_ERROR;					\
	goto exit;							\
      };								\
  }

#define	QMALLOC(ptr, typ, nel, status)					\
  {if (!(ptr = (typ *)malloc((size_t)(nel)*sizeof(typ))))		\
      {									\
	sprintf(errdetail, #ptr " (" #nel "=%lu elements) "		\
		"at line %d in module " __FILE__ " !",			\
		(size_t)(nel)*sizeof(typ), __LINE__);			\
	status = MEMORY_ALLOC_ERROR;					\
	goto exit;							\
      };								\
  }
