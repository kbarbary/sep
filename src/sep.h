

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

/*------------------------------- functions ---------------------------------*/
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
