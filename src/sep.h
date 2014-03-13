/*
*				define.h
*
* Global definitions
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	SExtractor is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*	SExtractor is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with SExtractor. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		12/04/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/* Check if we are using a configure script here */
#ifndef HAVE_CONFIG_H
#define		VERSION		"2.x"
#define		DATE		"2009-03-31"
#define		THREADS_NMAX	1024		/* max. number of threads */
#endif

/*------------------------ what, who, when and where ------------------------*/

#define		BANNER		"SExtractor"
#define		MYVERSION	VERSION
#define		EXECUTABLE	"sex"
#define         COPYRIGHT       "2012 IAP/CNRS/UPMC"
#define		DISCLAIMER	BANNER " comes with ABSOLUTELY NO WARRANTY\n" \
		"You may redistribute copies of " BANNER "\n" \
		"under the terms of the GNU General Public License."
#define		AUTHORS		"Emmanuel BERTIN <bertin@iap.fr>"
#define		WEBSITE		"http://astromatic.net/software/sextractor"
#define		INSTITUTE	"IAP  http://www.iap.fr"

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

/*---- Set defines according to machine's specificities and customizing -----*/

#if _LARGEFILE_SOURCE
#define	FSEEKO	fseeko
#define	FTELLO	ftello
#else
#define	FSEEKO	fseek
#define	FTELLO	ftell
#endif
/*--------------------- in case of missing constants ------------------------*/

#ifndef		SEEK_SET
#define		SEEK_SET	0
#endif
#ifndef		SEEK_CUR
#define		SEEK_CUR	1
#endif

#ifndef	EXIT_SUCCESS
#define	EXIT_SUCCESS		0
#endif
#ifndef	EXIT_FAILURE
#define	EXIT_FAILURE		-1
#endif

/*---------------------------- return messages ------------------------------*/

#define		RETURN_OK		0
#define		RETURN_ERROR		(-1)
#define		RETURN_FATAL_ERROR	(-2)

/*------------------- a few definitions to read FITS parameters ------------*/

#define	FBSIZE	2880L	/* size (in bytes) of one FITS block */

#define	FITSTOF(k, def) \
			(st[0]=0,((point = fitsnfind(buf, k, n))? \
				 fixexponent(point), \
				atof(strncat(st, &point[10], 70)) \
				:(def)))

#define	FITSTOI(k, def) \
			(st[0]=0,(point = fitsnfind(buf, k, n))? \
				 atoi(strncat(st, &point[10], 70)) \
				:(def))

#define	FITSTOS(k, str, def) \
                { if (fitsread(buf,k,str,H_STRING,T_STRING)!= RETURN_OK) \
                    strcpy(str, (def)); \
                }

/*------------------------------- Other Macros -----------------------------*/

#define	DEXP(x)		exp(2.30258509299*(x))		/* 10^x */
#define	DEXPF(x)	expf(2.30258509299f*(x))	/* 10^x */

#define QFREAD(ptr, size, afile, fname) \
		if (fread(ptr, (size_t)(size), (size_t)1, afile)!=1) \
		  error(EXIT_FAILURE, "*Error* while reading ", fname)

#define QFWRITE(ptr, size, afile, fname) \
		if (fwrite(ptr, (size_t)(size), (size_t)1, afile)!=1) \
		  error(EXIT_FAILURE, "*Error* while writing ", fname)

#define	QFSEEK(afile, offset, pos, fname) \
		if (FSEEKO(afile, (offset), pos)) \
		  error(EXIT_FAILURE,"*Error*: file positioning failed in ", \
			fname)

#define	QFTELL(afile, pos, fname) \
		if ((pos=FTELLO(afile))==-1) \
		  error(EXIT_FAILURE,"*Error*: file position unknown in ", \
			fname)

#define	QFREE(ptr) \
		{free(ptr); \
		ptr = NULL;}

#define	QCALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)calloc((size_t)(nel),sizeof(typ)))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QMALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QMALLOC16(ptr, typ, nel) \
		{if (posix_memalign((void **)&ptr, 16, (size_t)(nel)*sizeof(typ))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QREALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)realloc(ptr, (size_t)(nel)*sizeof(typ))))\
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define QMEMCPY(ptrin, ptrout, typ, nel) \
		{if (ptrin) \
                  {if (!(ptrout = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		     { \
		     sprintf(gstr, #ptrout " (" #nel "=%lld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		     error(EXIT_FAILURE,"Could not allocate memory for ",gstr);\
                     }; \
                   memcpy(ptrout, ptrin, (size_t)(nel)*sizeof(typ)); \
                   }; \
                 }

#define	RINT(x)	(int)(floor(x+0.5))

#define	PIX(pic, x, y)	pic->strip[(((int)y)%pic->stripheight) \
				*pic->width +(int)x]

#define	NPRINTF		if (prefs.verbose_type == NORM \
				|| prefs.verbose_type==WARN) fprintf

#define	NFPRINTF(w,x)	{if (prefs.verbose_type==NORM \
				|| prefs.verbose_type==WARN) \
				fprintf(w, "\33[1M> %s\n\33[1A",x); \
			else if (prefs.verbose_type == FULL) \
				fprintf(w, "%s.\n", x);}

#define	QPRINTF		if (prefs.verbose_type != QUIET)	fprintf

#define	FPRINTF		if (prefs.verbose_type == FULL)	fprintf

#define	QWARNING       	if (prefs.verbose_type==WARN \
				|| prefs.verbose_type==FULL)	warning

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

typedef struct backspl
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
} backsplstruct;

/*------------------------------- functions ---------------------------------*/
void backhisto(backstruct *, backstruct *, PIXTYPE *, PIXTYPE *,
	       size_t, int, int, int, PIXTYPE);
void backstat(backstruct *, backstruct *, PIXTYPE *, PIXTYPE *,
	      size_t, int, int, int, PIXTYPE);
/* void backrmsline(picstruct *, int, PIXTYPE *);
void copyback(picstruct *infield, picstruct *outfield);
void endback(picstruct *); */
void filterback(backsplstruct *);
/* void subbackline(picstruct *, int, PIXTYPE *); */

backsplstruct *makeback(float *, float *, int, int, int, int, float);

float backguess(backstruct *, float *, float *),
/* float localback(picstruct *, objstruct *), */
float *makebackspline(picstruct *, float *);

extern PIXTYPE	back(picstruct *, int, int);


//prototypes
float makeback(float *, int, int);

/*---------------------------------------------------------------------------*/
/* misc                                                                      */
/*---------------------------------------------------------------------------*/

float fqmedian(float *, int);
