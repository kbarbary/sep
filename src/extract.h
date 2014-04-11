#define	UNKNOWN			-1		/* flag for LUTZ */

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
