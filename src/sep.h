

#define		VERSION		"0.1"
#define		DATE		"2014-03-14"

/*--------------------------- Internal constants ----------------------------*/
#define	BIG			1e+30		/* a huge number */
#define	LESSBIG			1e+25		/* a somewhat smaller number */
#define	DATA_BUFSIZE		262144		/* data buffer size */
#define	MARGIN_SCALE		2.0		/* Margin / object height */ 
#define	MARGIN_OFFSET		4.0		/* Margin offset (pixels) */ 
#define	MAXCHAR			512		/* max. number of characters */
#define	MAXCHARL		16384		/* max.nb of chars in strlist*/
#define	MAXDEBAREA		3		/* max. area for deblending */
#define	MAXFLAG			4		/* max. # of FLAG-images */
#define	MAXIMAGE		2		/* max. # of input images */
#define	MAXNAPER		32		/* max. number of apertures */
#define	MAXNASSOC		32		/* max. number of assoc. */
#define	MAXPICSIZE		1048576		/* max. image size */
#define	NISO			8		/* number of isophotes */
#define	OUTPUT			stderr		/* where all msgs are sent */
#define PSF_NPSFMAX		9		/* Max number of fitted PSFs */

#ifndef PI
#define	PI			3.1415926535898	/* never met before? */
#endif
#define	DEG			(PI/180.0)	/* 1 deg in radians */

/* NOTES:
 *
 *One must have:	BIG < the biggest element a float can store
 *			DATA_BUFSIZE >= 2880 with DATA_BUFSIZE%8 = 0
 *			MAXCHAR >= 16
 *			1 <= MAXCHECK <= MAXLIST (see prefs.h)
 *			1 <= MAXDEBAREA (see prefs.c & extract.c)
 *			1 <= MAXFLAG <= MAXLIST (see prefs.h)
 *			1 <= MAXIMAGE <= MAXLIST (see prefs.h)
 *			1 <= MAXNAPER <= MAXLIST (see prefs.h)
 *			1 <= MAXNASSOC <= MAXLIST (see prefs.h)
 *			MAXPICSIZE > size of any image!!
 *			NISO = 8 (otherwise need to change prefs.h)
 *			1 <= PSF_NPSFMAX
*/

/* definitions from define.h ------------------------------------------------*/
#define		OBJ_CROWDED	0x0001
#define		OBJ_MERGED	0x0002
#define		OBJ_SATUR	0x0004
#define		OBJ_TRUNC	0x0008
#define		OBJ_APERT_PB	0x0010
#define		OBJ_ISO_PB	0x0020
#define		OBJ_DOVERFLOW	0x0040
#define		OBJ_OVERFLOW	0x0080

/*--------------------- in case of missing constants ------------------------*/
#ifndef	EXIT_SUCCESS
#define	EXIT_SUCCESS		0
#endif
#ifndef	EXIT_FAILURE
#define	EXIT_FAILURE		-1
#endif

/*------------------------------- types ------------------------------------*/
typedef float PIXTYPE;
typedef	int		LONG;
typedef	unsigned int	ULONG;

/* globals -----------------------------------------------------------------*/
char			gstr[MAXCHAR];

/*------------------------------- Other Macros -----------------------------*/

#define	DEXP(x)		exp(2.30258509299*(x))		/* 10^x */
#define	DEXPF(x)	expf(2.30258509299f*(x))	/* 10^x */

#define	QCALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)calloc((size_t)(nel),sizeof(typ)))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lu elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QMALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lu elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QMALLOC16(ptr, typ, nel) \
		{if (posix_memalign((void **)&ptr, 16, (size_t)(nel)*sizeof(typ))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lu elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QREALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)realloc(ptr, (size_t)(nel)*sizeof(typ))))\
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lu elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define QMEMCPY(ptrin, ptrout, typ, nel) \
		{if (ptrin) \
                  {if (!(ptrout = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		     { \
		     sprintf(gstr, #ptrout " (" #nel "=%lu elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		     error(EXIT_FAILURE,"Could not allocate memory for ",gstr);\
                     }; \
                   memcpy(ptrout, ptrin, (size_t)(nel)*sizeof(typ)); \
                   }; \
                 }

#define	RINT(x)	(int)(floor(x+0.5))

#define FLAG(x)		(*((char *)&flag##x))

#define VECFLAG(x)	(*((char *)flag##x))


/*---------------------------------------------------------------------------*/
/* back.h                                                                    */
/*---------------------------------------------------------------------------*/

/*----------------------------- Internal constants --------------------------*/
#define	BACK_BUFSIZE		1048576		/* bkgnd buffer */
#define	BACK_MINGOODFRAC	0.5		/* min frac with good weights*/
#define	QUANTIF_NSIGMA		5		/* histogram limits */
#define	QUANTIF_NMAXLEVELS	4096		/* max nb of quantif. levels */
#define	QUANTIF_AMIN		4		/* min nb of "mode pixels" */

#define	BACK_WSCALE		1		/* Activate weight scaling */
#define	BACK_NOWSCALE		0		/* No weight scaling */

/*------------------------------- structures --------------------------------*/
typedef float PIXTYPE;

/* Background info */
typedef struct structback
{
  float		mode, mean, sigma;	/* Background mode, mean and sigma */
  LONG		*histo;			/* Pointer to a histogram */
  int		nlevels;		/* Nb of histogram bins */
  float		qzero, qscale;		/* Position of histogram */
  float		lcut, hcut;		/* Histogram cuts */
  int		npix;			/* Number of pixels involved */
} backstruct;

typedef struct structbackspline
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
} backspline;

/* Public API ---------------------------------------------------------------*/
backspline *makeback(PIXTYPE *, PIXTYPE *, int, int, int, int, PIXTYPE,
		     int, int, float);
PIXTYPE	backpixlinear(backspline *, int, int);
void backline(backspline *, int, PIXTYPE *);
void backim(backspline *, PIXTYPE *);
void backrmsline(backspline *, int, PIXTYPE *);
void backrmsim(backspline *, PIXTYPE *);
void freeback(backspline *);

/* priviate functions -------------------------------------------------------*/
void backhisto(backstruct *, PIXTYPE *, PIXTYPE *,
	       int, int, int, int, PIXTYPE);
void backstat(backstruct *, PIXTYPE *, PIXTYPE *,
	      int, int, int, int, PIXTYPE);
void filterback(backspline *, int, int, float);
float backguess(backstruct *, float *, float *);
float *makebackspline(backspline *, float *);

extern void error(int, char *, char *);

/*---------------------------------------------------------------------------*/
/* misc                                                                      */
/*---------------------------------------------------------------------------*/

float fqmedian(float *, int);

/*---------------------------------------------------------------------------*/
/* types.h                                                                   */
/*---------------------------------------------------------------------------*/


/*-------------------------------- catalog  ---------------------------------*/

typedef struct
  {
  int		ndetect;				/* nb of detections */
  int		ntotal;					/* Total object nb */
  int		nparam;					/* Nb of parameters */
/*----- Misc. strings defining the extraction */
  char		prefs_name[MAXCHAR];			/* Prefs filename*/
  char		image_name[MAXCHAR];			/* image filename*/
  char		psf_name[MAXCHAR];			/* PSF filename*/
  char		nnw_name[MAXCHAR];			/* NNW name */
  char		filter_name[MAXCHAR];			/* Filter name */
  char		soft_name[MAXCHAR];			/* Sextractor version*/
/*----- time */
  char		ext_date[16],ext_time[16];		/* date and time */
  double	ext_elapsed;				/* processing time */
/*----- MEF */
  int		currext;				/* current extension */
  int		next;					/* Nb of extensions */
  }		sexcatstruct;

/*---------------------------------------------------------------------------*/
/* weight.h                                                                  */
/*---------------------------------------------------------------------------*/


#define	WTHRESH_CONVFAC		1e-4	/* Factor to apply to weights when */
					/* thresholding filtered weight-maps */



/*---------------------------------------------------------------------------*/
/* extract.h                                                                 */
/*---------------------------------------------------------------------------*/


/*------------------------------ definitions --------------------------------*/

#define	NOBJ			256		/* starting number of obj. */
#define	UNKNOWN			-1		/* flag for LUTZ */

/*--------------------------------- typedefs --------------------------------*/
typedef	char		pliststruct;		/* Dummy type for plist */
typedef	enum		{COMPLETE, INCOMPLETE, NONOBJECT, OBJECT}
				status;	/* Extraction status */

/*--------------------------------- variables -------------------------------*/

PIXTYPE		*dumscan;

/*------------------------------- structures --------------------------------*/
/* Temporary object parameters during extraction */
typedef struct structinfo
  {
  LONG		pixnb;			/* Number of pixels included */
  LONG		firstpix;		/* Pointer to first pixel of pixlist */
  LONG		lastpix;		/* Pointer to last pixel of pixlist */
  short		flag;			/* Extraction flag */
  }       infostruct;

typedef struct
  {
  int		nobj;			/* number of objects in list */
  objstruct	*obj;			/* pointer to the object array */
  int		npix;			/* number of pixels in pixel-list */
  pliststruct	*plist;			/* pointer to the pixel-list */
  PIXTYPE	dthresh;		/* detection threshold */
  PIXTYPE	thresh;			/* analysis threshold */
  }	objliststruct;

/*--------------------------------- objects ---------------------------------*/
/* I: "PIXEL" parameters */

typedef struct
  {
/* ---- basic parameters */
  int		number;				/* ID */
  int		fdnpix;				/* nb of extracted pix */
  int		dnpix;				/* nb of pix above thresh  */
  int		npix;				/* "" in measured frame */
  int		nzdwpix;			/* nb of zero-dweights around */
  int		nzwpix;				/* nb of zero-weights inside */
  float		fdflux;				/* integrated ext. flux */
  float		dflux;				/* integrated det. flux */
  float		flux;				/* integrated mes. flux */
  float		fluxerr;			/* integrated variance */
  PIXTYPE	fdpeak;				/* peak intensity (ADU) */
  PIXTYPE	dpeak;				/* peak intensity (ADU) */
  PIXTYPE	peak;				/* peak intensity (ADU) */
/* ---- astrometric data */
  int		peakx,peaky;			/* pos of brightest pix */
  double       	mx, my;				/* barycenter */
  double	poserr_mx2, poserr_my2,
		poserr_mxy;			/* Error ellips moments */
/* ---- morphological data */			
  int		xmin,xmax,ymin,ymax,ycmin,ycmax;/* x,y limits */
  PIXTYPE	*blank, *dblank; 	       	/* BLANKing sub-images  */
  int		*submap;			/* Pixel-index sub-map */
  int		subx,suby, subw,subh;		/* sub-image pos. and size */
  short		flag;				/* extraction flags */
  BYTE		wflag;				/* weighted extraction flags */
  FLAGTYPE	imaflag[MAXFLAG];		/* flags from FLAG-images */
  BYTE		singuflag;			/* flags for singularities */
  int		imanflag[MAXFLAG];     		/* number of MOST flags */
  double	mx2,my2,mxy;			/* variances and covariance */
  float		a, b, theta, abcor;		/* moments and angle */
  float		cxx,cyy,cxy;			/* ellipse parameters */
  int		firstpix;			/* ptr to first pixel */
  int		lastpix;			/* ptr to last pixel */
  float		bkg, dbkg, sigbkg, dsigbkg;	/* Background stats (ADU) */
  float		thresh;				/* measur. threshold (ADU) */
  float		dthresh;		       	/* detect. threshold (ADU) */
  float		mthresh;		       	/* max. threshold (ADU) */
  int		iso[NISO];			/* isophotal areas */
  float		fwhm;				/* IMAGE FWHM */
  
  }	objstruct;

/*------------------------------- functions ---------------------------------*/
void		lutzalloc(int, int);

/*		lutzfree(void),
		lutzsort(infostruct *, objliststruct *),
		sortit(picstruct *, picstruct *, picstruct *, picstruct *,
			infostruct *, objliststruct *, PIXTYPE *, PIXTYPE *),
		update(infostruct *, infostruct *, pliststruct *);

int		lutz(objliststruct *, int, objstruct *, objliststruct *); 

*/

/*---------------------------------------------------------------------------*/
/* globals.h                                                                 */
/*---------------------------------------------------------------------------*/

void allocparcelout(void);

/*---------------------------------------------------------------------------*/
/* clean.h                                                                   */
/*---------------------------------------------------------------------------*/

/*------------------------------ definitions --------------------------------*/

#define		CLEAN_ZONE		10.0	/* zone (in sigma) to */
						/* consider for processing */

#define CLEAN_FLAG 1      /* replaces prefs.clean_flag (move to scan input?) */
#define CLEAN_STACKSIZE 3000  /* replaces prefs.clean_stacksize  */
                              /* (MEMORY_OBJSTACK in sextractor inputs) */

/*------------------------------- variables ---------------------------------*/

objliststruct	*cleanobjlist;		/* laconic, isn't it? */

/*------------------------------- functions ---------------------------------*/

extern void	addcleanobj(objstruct *),
		endclean(void),
		initclean(void),
		subcleanobj(int);

extern int	clean(picstruct *field, picstruct *dfield,
		      int, objliststruct *);

/*---------------------------------------------------------------------------*/
/* plist.h                                                                   */
/*---------------------------------------------------------------------------*/


/*------------------------------- definitions -------------------------------*/

#define	PLIST(ptr, elem)	(((pbliststruct *)(ptr))->elem)

#define	PLISTEXIST(elem)	(plistexist_##elem)

#define	PLISTPIX(ptr, elem)	(*((PIXTYPE *)((ptr)+plistoff_##elem)))

#define	PLISTFLAG(ptr, elem)	(*((FLAGTYPE *)((ptr)+plistoff_##elem)))

/*------------------------------- structures --------------------------------*/

typedef struct
  {
  int		nextpix;
  int		x, y;
  PIXTYPE       value;
  }	pbliststruct;

/*-------------------------------- globals ----------------------------------*/

int	plistexist_value, plistexist_dvalue, plistexist_cdvalue,
	plistexist_flag, plistexist_wflag, plistexist_dthresh, plistexist_var,
	plistoff_value, plistoff_dvalue, plistoff_cdvalue,
	plistoff_flag[MAXFLAG], plistoff_wflag, plistoff_dthresh, plistoff_var,
	plistsize;

/*------------------------------- functions ---------------------------------*/

void	init_plist(PIXTYPE *filter, PIXTYPE *cdwfield);

int	createblank(objliststruct *objlist, int n),
	createsubmap(objliststruct *objlist, int n);

#define MEMORY_PIXSTACK  300000         /* number of pixels in stack */
                                        /* replaces prefs.mem_pixstack */
#define MEMORY_BUFSIZE 1024             /* number of lines in buffer */
/*------------------------------- error codes -------------------------------*/

#define MEMORY_PIXSTACK_ERROR 1

/*---------------------------------------------------------------------------*/
/* filter.h                                                                  */
/*---------------------------------------------------------------------------*/

void	convolve(PIXTYPE *im,                      /* full image (was field) */
		 int w, int h,                     /* image size */
		 int y,                            /* line in image */
		 float *conv,                      /* convolution mask */
		 int convw, int convh, int convn,  /* size of conv */
		 PIXTYPE *mscan);                  /* convolved line */
