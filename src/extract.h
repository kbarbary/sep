
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
typedef	enum {COMPLETE, INCOMPLETE, NONOBJECT, OBJECT} status;

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

void preanalyse(int, objliststruct *, int);
extern int addcleanobj(objstruct *, objliststruct *);
extern int subcleanobj(int, objliststruct *);
extern int clean(int, objliststruct *, objliststruct *, LONG *, double);
void lutzalloc(int, int);
void lutzfree(void);
int  lutz(objliststruct *, int, objstruct *, objliststruct *, int);
void update(infostruct *, infostruct *, pliststruct *);
void allocparcelout(void);
void freeparcelout(void);
int  parcelout(objliststruct *, objliststruct *, int, int);
