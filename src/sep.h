/*---------------------------------------------------------------------------*/
/* sep public interface                                                      */
/*---------------------------------------------------------------------------*/

#define	SEP_VERSION  "0.1.0"
#define	SEP_DATE     "2014-03-14"

/*------------------------------- error codes -------------------------------*/
#define	RETURN_OK		0
#define	RETURN_ERROR		(-1) /* general unspecified error */
#define	RETURN_FATAL_ERROR	(-1)
#define MEMORY_PIXSTACK_ERROR   1
#define PIXSTACK_OVERFLOW_ERROR 2
#define FATAL_ERROR             3
#define MEMORY_CLEAN_ERROR      4
#define NO_CLEAN_OBJ_ERROR      5 /* Internal Error: no CLEAN object to
                                     remove in subcleanobj()*/
#define LUTZ_REALLOC_ERROR      6 /* problem with mem. realloc. in lutz() */
#define GATHERUP_MEMORY_ERROR   7 /* Not enough memory to update pixel list in
				     gatherup()" */
#define	EXIT_SUCCESS		0
#define	EXIT_FAILURE		-1

typedef float PIXTYPE;

/*------------------------------ background ---------------------------------*/
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

backspline *makeback(PIXTYPE *, PIXTYPE *, int, int, int, int, PIXTYPE,
		     int, int, float);
PIXTYPE	backpixlinear(backspline *, int, int);
void backline(backspline *, int, PIXTYPE *);
void backim(backspline *, PIXTYPE *);
void backrmsline(backspline *, int, PIXTYPE *);
void backrmsim(backspline *, PIXTYPE *);
void freeback(backspline *);


/*-------------------------- source extraction ------------------------------*/

#define		OBJ_CROWDED	0x0001
#define		OBJ_MERGED	0x0002
#define		OBJ_SATUR	0x0004
#define		OBJ_TRUNC	0x0008
#define		OBJ_APERT_PB	0x0010
#define		OBJ_ISO_PB	0x0020
#define		OBJ_DOVERFLOW	0x0040
#define		OBJ_OVERFLOW	0x0080

typedef struct
{
  int           nobj;	  /* number of objects in list */
  objstruct     *obj;	  /* pointer to the object array */
  int           npix;	  /* number of pixels in pixel-list */
  pliststruct   *plist;	  /* pointer to the pixel-list */
  PIXTYPE       dthresh;  /* detection threshold */
  PIXTYPE       thresh;	  /* analysis threshold */
} objliststruct;

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
  BYTE	   wflag;			     /* weighted extraction flags */
  FLAGTYPE imaflag[MAXFLAG];		     /* flags from FLAG-images */
  BYTE	   singuflag;			     /* flags for singularities */
  int      imanflag[MAXFLAG];     	     /* number of MOST flags */
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

/* below, sextractor defaults are shown in [] */
int extract(PIXTYPE *, /* image array                                        */
	    PIXTYPE *, /* variance array (can be NULL)                       */
	    int, int,  /* size of arrays (w, h)                              */
	    PIXTYPE,   /* detection threshold                          [1.5] */
	    PIXTYPE,   /* analysis threshold                           [1.5] */
	    PIXTYPE,   /* cdwthresh (???)                                    */
	    int,       /* threshabsolute (0=relative)                    [0] */
	    int,       /* minarea                                        [5] */
	    float *,   /* conv array (can be NULL)     [{1 2 1 2 4 2 1 2 1}] */
	    int, int,  /* conv size (w, h)                            [3, 3] */
	    int,       /* deblend_nthresh                               [32] */
	    double,    /*deblend_mincont                             [0.005] */
	    int,       /* clean_flag (1 = YES)                           [1] */
	    double)    /* clean_param                                  [1.0] */
