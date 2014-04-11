/*---------------------------------------------------------------------------*/
/* global internal definitions                                               */
/*---------------------------------------------------------------------------*/

#define	BIG	       1e+30		/* a huge number */
#define	LESSBIG	       1e+25		/* a somewhat smaller number */
#define	DATA_BUFSIZE   262144		/* data buffer size */
#define	MARGIN_SCALE   2.0		/* Margin / object height */ 
#define	MARGIN_OFFSET  4.0		/* Margin offset (pixels) */ 
#define	MAXCHAR	       512		/* max. number of characters */
#define	MAXCHARL       16384		/* max.nb of chars in strlist*/
#define	MAXDEBAREA     3		/* max. area for deblending */
#define	MAXFLAG	       4		/* max. # of FLAG-images */
#define	MAXPICSIZE     1048576		/* max. image size */
#define	NISO	       8		/* number of isophotes */
#define	OUTPUT	       stderr		/* where all msgs are sent */
#define PSF_NPSFMAX    9		/* Max number of fitted PSFs */
#ifndef PI
#define	PI	       3.1415926535898	/* never met before? */
#endif
#define	DEG	       (PI/180.0)	/* 1 deg in radians */

/* NOTE: One must have
 *
 * BIG < the biggest element a float can store
 * DATA_BUFSIZE >= 2880 with DATA_BUFSIZE%8 = 0
 * MAXCHAR >= 16
 * 1 <= MAXCHECK <= MAXLIST (see prefs.h)
 * 1 <= MAXDEBAREA (see prefs.c & extract.c)
 * 1 <= MAXFLAG <= MAXLIST (see prefs.h)
 *			1 <= MAXIMAGE <= MAXLIST (see prefs.h)
 *			1 <= MAXNAPER <= MAXLIST (see prefs.h)
 *			1 <= MAXNASSOC <= MAXLIST (see prefs.h)
 *			MAXPICSIZE > size of any image!!
 *			NISO = 8 (otherwise need to change prefs.h)
 *			1 <= PSF_NPSFMAX
*/

typedef	int		LONG;
typedef	unsigned int	ULONG;

char			gstr[MAXCHAR];

#define	QCALLOC(ptr, typ, nel)						\
  {if (!(ptr = (typ *)calloc((size_t)(nel),sizeof(typ))))		\
      {									\
	sprintf(gstr, #ptr " (" #nel "=%lu elements) "			\
		"at line %d in module " __FILE__ " !",			\
		(size_t)(nel)*sizeof(typ), __LINE__);			\
	error(EXIT_FAILURE, "Could not allocate memory for ", gstr);	\
      };								\
  }

#define	QMALLOC(ptr, typ, nel)						\
  {if (!(ptr = (typ *)malloc((size_t)(nel)*sizeof(typ))))		\
      {									\
	sprintf(gstr, #ptr " (" #nel "=%lu elements) "			\
		"at line %d in module " __FILE__ " !",			\
		(size_t)(nel)*sizeof(typ), __LINE__);			\
	error(EXIT_FAILURE, "Could not allocate memory for ", gstr);	\
      };								\
  }

#define	QMALLOC16(ptr, typ, nel)					\
  {if (posix_memalign((void **)&ptr, 16, (size_t)(nel)*sizeof(typ)))	\
      {									\
	sprintf(gstr, #ptr " (" #nel "=%lu elements) "			\
		"at line %d in module " __FILE__ " !",			\
		(size_t)(nel)*sizeof(typ), __LINE__);			\
	error(EXIT_FAILURE, "Could not allocate memory for ", gstr);	\
      };								\
  }

#define	QREALLOC(ptr, typ, nel)						\
  {if (!(ptr = (typ *)realloc(ptr, (size_t)(nel)*sizeof(typ))))		\
      {									\
	sprintf(gstr, #ptr " (" #nel "=%lu elements) "			\
		"at line %d in module " __FILE__ " !",			\
		(size_t)(nel)*sizeof(typ), __LINE__);			\
	error(EXIT_FAILURE, "Could not allocate memory for ", gstr);	\
      };								\
  }

#define QMEMCPY(ptrin, ptrout, typ, nel)	                        \
  {if (ptrin)								\
      {if (!(ptrout = (typ *)malloc((size_t)(nel)*sizeof(typ))))	\
	  {								\
	    sprintf(gstr, #ptrout " (" #nel "=%lu elements) "		\
		    "at line %d in module " __FILE__ " !",		\
		    (size_t)(nel)*sizeof(typ), __LINE__);		\
	    error(EXIT_FAILURE,"Could not allocate memory for ",gstr);	\
	  };								\
	memcpy(ptrout, ptrin, (size_t)(nel)*sizeof(typ));		\
      };								\
  }

/* can remove this after we get rid of the refs in this file */
/* TODO: make a public errormsg function that returns the message in gstr */
extern void error(int, char *, char *);
