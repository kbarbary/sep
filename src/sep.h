/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
* This file is part of SEP
*
* Copyright 1993-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
* Copyright 2014 SEP developers
*
* SEP is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* SEP is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with SEP.  If not, see <http://www.gnu.org/licenses/>.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

typedef float PIXTYPE;   /* type of image arrays */

/*------------------------- constant definitions ----------------------------*/

#define TFLOAT       42
#define TDOUBLE      82

#define SEP_OBJ_CROWDED   0x0001
#define SEP_OBJ_MERGED    0x0002
#define SEP_OBJ_SATUR     0x0004
#define SEP_OBJ_TRUNC     0x0008
#define SEP_OBJ_APERT_PB  0x0010
#define SEP_OBJ_ISO_PB    0x0020
#define SEP_OBJ_DOVERFLOW 0x0040
#define SEP_OBJ_OVERFLOW  0x0080
#define SEP_OBJ_SINGU     0x0100

/*--------------------- global background estimation ------------------------*/

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
} sepbackmap;

/* w, h is image size in pixels */
/* bw, bh is size of a single background tile in pixels */
/* var > varthresh will be ignored. */
int sep_makeback(PIXTYPE *im, PIXTYPE *mask, int w, int h,
		 int bw, int bh, PIXTYPE maskthresh, int fbx, int fby,
		 float fthresh, sepbackmap **bkm);
PIXTYPE	sep_backpixlinear(sepbackmap *, int, int);
int sep_backline(sepbackmap *, int, PIXTYPE *);
int sep_backrmsline(sepbackmap *, int, PIXTYPE *);
int sep_backvarline(sepbackmap *, int, PIXTYPE *);
int sep_backarray(sepbackmap *, PIXTYPE *);
int sep_backrmsarray(sepbackmap *, PIXTYPE *);
int sep_backvararray(sepbackmap *, PIXTYPE *);
int sep_subbackline(sepbackmap *, int, PIXTYPE *);
int sep_subbackarray(sepbackmap *, PIXTYPE *);
void sep_freeback(sepbackmap *);

/*-------------------------- source extraction ------------------------------*/

typedef struct
{
  float	   thresh;               /* threshold (ADU)                          */
  int	   npix;                 /* # pixels extracted (size of pix array)   */
  int      tnpix;                /* # pixels above thresh (unconvolved)      */
  int	   xmin,xmax,ymin,ymax;  /* x,y limits                               */
  double   mx, my;               /* barycenter                               */
  double   mx2,my2,mxy;		 /* variances and covariance                 */
  float	   a, b, theta, abcor;   /* moments and angle                        */
  float	   cxx,cyy,cxy;	         /* ellipse parameters                       */
  float	   cflux;                /* total flux of pixels (convolved im)      */
  float	   flux;      		 /* total flux of pixels (unconvolved)       */
  PIXTYPE  cpeak;                /* peak intensity (ADU) (convolved)         */
  PIXTYPE  peak;                 /* peak intensity (ADU) (unconvolved)       */
  short	   flag;                 /* extraction flags                         */
  int      *pix;                 /* pixel array (length is npix)             */
} sepobj;

int sep_extract(PIXTYPE *im,            /* image array                      */
		PIXTYPE *var,           /* variance array (can be NULL)     */
		int w, int h,           /* width and height of arrays       */
		PIXTYPE thresh,         /* detection threshold  [1.5*sigma] */
		int minarea,            /* minimum area in pixels       [5] */
		float *conv,            /* convolution array (can be NULL)  */
                                        /*            [{1 2 1 2 4 2 1 2 1}] */
		int convw, int convh,   /* w, h of convolution array  [3,3] */
		int deblend_nthresh,    /* deblending thresholds       [32] */
		double deblend_mincont, /* min. deblending contrast [0.005] */
		int clean_flag,         /* perform cleaning?            [1] */
		double clean_param,     /* clean parameter            [1.0] */
		int *nobj,              /* OUTPUT: number of objects        */
		sepobj **objects);      /* OUTPUT: object array             */


void sep_freeobjarray(sepobj *objects, int nobj);
/* free memory associated with an sepobj array, including pixel lists */

/*-------------------------- aperture photometry ----------------------------*/

int sep_apercirc(void *im, void *var, int dtype, int w, int h,
		  PIXTYPE gain, PIXTYPE varthresh,
		  double cx, double cy, double r, int subpix,
		  double *flux, double *fluxerr, short *flag);

/*----------------------- info & error messaging ----------------------------*/

extern char *sep_version_string;
void sep_get_errmsg(int status, char *errtext);
void sep_get_errdetail(char *errtext);
