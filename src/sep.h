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

/* datatype codes */
#define SEP_TBYTE        11  /* 8-bit unsigned byte */
#define SEP_TINT         31  /* native int type */
#define SEP_TFLOAT       42
#define SEP_TDOUBLE      82

/* object & aperture flags */
#define SEP_OBJ_MERGED       0x0001  /* object is result of deblending */
#define SEP_OBJ_TRUNC        0x0002  /* object truncated at image boundary */
#define SEP_OBJ_DOVERFLOW    0x0004  /* not currently used, but could be */
#define SEP_OBJ_SINGU        0x0008  /* x,y fully correlated */
#define SEP_APER_TRUNC       0x0010
#define SEP_APER_HASMASKED   0x0020
#define SEP_APER_ALLMASKED   0x0040
#define SEP_APER_NONPOSITIVE 0x0080

/* noise_type values in sep_image */
#define SEP_NOISE_NONE   0
#define SEP_NOISE_STDDEV 1
#define SEP_NOISE_VAR    2

/* input flags for aperture photometry */
#define SEP_MASK_IGNORE      0x0004

/* threshold interpretation for sep_extract */
#define SEP_THRESH_RELATIVE 0  /* in units of standard deviations (sigma) */
#define SEP_THRESH_ABSOLUTE 1  /* absolute data values */

/* filter types for sep_extract */
#define SEP_FILTER_CONV    0
#define SEP_FILTER_MATCHED 1

/*---------------------- image type -----------------------------------------*/

typedef struct {
  void *data;      /* data array                */
  void *noise;     /* noise array (can be NULL) */
  void *mask;      /* mask array (can be NULL)  */
  int dtype;       /* element type of image     */
  int ndtype;      /* element type of noise     */
  int mdtype;      /* element type of mask      */
  int w;           /* array width               */
  int h;           /* array height              */
  double noiseval; /* scalar noise value; used only if noise == NULL */
  short noise_type; /* interpretation of noise value */
  double maskthresh; /* only (mask<=maskthresh) counted     */
  double gain;   /* (poisson counts / data unit) */
} sep_image;


/*--------------------- global background estimation ------------------------*/

typedef struct {
  int w, h;          /* original image width, height */
  int bw, bh;        /* single tile width, height */
  int nx, ny;        /* number of tiles in x, y */
  int n;             /* nx*ny */
  float globalback;  /* global mean */
  float globalrms;   /* global sigma */
  float *back;       /* node data for interpolation */
  float *dback;
  float *sigma;    
  float *dsigma;
} sep_backmap;

/* sep_makeback()
 * 
 * Create representation of spatially varying image background and variance.
 *
 * Note that the returned pointer must eventually be freed by calling 
 * `sep_freeback()`.
 *
 * If a mask is supplied, only pixels with mask value <= mthresh are counted.
 * In addition to the mask, pixels <= -1e30 and NaN are ignored.
 * 
 * Source Extractor defaults:
 * 
 * - bw, bh = (64, 64)
 * - fw, fh = (3, 3)
 * - fthresh = 0.0
 */
int sep_makeback(sep_image *image,
		 int bw, int bh,      /* size of a single background tile    */
		 float mthresh,       
		 int fw, int fh,      /* filter size in tiles                */
		 float fthresh,       /* filter threshold                    */
		 sep_backmap **bkmap); /* OUTPUT                              */

/* Get the estimate of the global background "mean" or standard deviation */
float sep_globalback(sep_backmap *bkmap);
float sep_globalrms(sep_backmap *bkmap);

/* Return background at (x, y).
 * Unlike other routines, this uses simple linear interpolation. */
float sep_backpix_linear(sep_backmap *bkmap, int x, int y);

/* Evaluate the background, RMS, or variance at line y.
 * Uses bicubic spline interpolation between background map verticies.
 * The second function subtracts the background from the input array.
 * Line must be an array with same width as original image. */
int sep_backline(sep_backmap *bkmap, int y, void *line, int dtype);
int sep_subbackline(sep_backmap *bkmap, int y, void *line, int dtype);
int sep_backrmsline(sep_backmap *bkmap, int y, void *line, int dtype);

/* Evaluate the background, RMS, or variance for entire image.
 * Uses bicubic spline interpolation between background map verticies.
 * The second function subtracts the background from the input array.
 * Arr must be an array of the same size as original image. */
int sep_backarray(sep_backmap *bkmap, void *arr, int dtype);
int sep_subbackarray(sep_backmap *bkmap, void *arr, int dtype);
int sep_backrmsarray(sep_backmap *bkmap, void *arr, int dtype);

/* Free memory associated with bkmap */
void sep_freeback(sep_backmap *bkmap);

/*-------------------------- source extraction ------------------------------*/

typedef struct
{
  float	   thresh;               /* threshold (ADU)                          */
  int	   npix;                 /* # pixels extracted (size of pix array)   */
  int      tnpix;                /* # pixels above thresh (unconvolved)      */
  int	   xmin,xmax,ymin,ymax;  /* x,y limits                               */
  double   x, y;                 /* barycenter (first moments)               */
  double   x2,y2,xy;		 /* second moments                           */
  float	   a, b, theta;          /* ellipse parameters                       */
  float	   cxx,cyy,cxy;	         /* ellipse parameters (alternative)         */
  float	   cflux;                /* total flux of pixels (convolved im)      */
  float	   flux;      		 /* total flux of pixels (unconvolved)       */
  float    cpeak;                /* peak intensity (ADU) (convolved)         */
  float    peak;                 /* peak intensity (ADU) (unconvolved)       */
  int      xcpeak, ycpeak;       /* x, y coords of peak (convolved) pixel    */
  int      xpeak, ypeak;         /* x, y coords of peak (unconvolved) pixel  */
  short	   flag;                 /* extraction flags                         */
  int      *pix;                 /* pixel array (length is npix)             */
} sepobj;

/* Extract sources from an image.
 *
 * Source Extractor defaults are shown in [ ] above.
 *
 * Notes
 * -----
 * `dtype` and `ndtype` indicate the data type (float, int, double) of the 
 * image and noise arrays, respectively.
 *
 * If `noise` is NULL, thresh is interpreted as an absolute threshold.
 * If `noise` is not null, thresh is interpreted as a relative threshold
 * (the absolute threshold will be thresh*noise[i,j]).
 * 
 */
int sep_extract(sep_image *image,
		float thresh,         /* detection threshold     [1.5] */
                int thresh_type,    /* threshold units [SEP_THRESH_RELATIVE] */
		int minarea,          /* minimum area in pixels          [5] */
		float *conv,          /* convolution array (can be NULL)     */
                                      /*               [{1 2 1 2 4 2 1 2 1}] */
		int convw, int convh, /* w, h of convolution array     [3,3] */
                int filter_type,      /* convolution (0) or matched (1)  [0] */
		int deblend_nthresh,  /* deblending thresholds          [32] */
		double deblend_cont,  /* min. deblending contrast    [0.005] */
		int clean_flag,       /* perform cleaning?               [1] */
		double clean_param,   /* clean parameter               [1.0] */
		sepobj **objects,     /* OUTPUT: object array                */
		int *nobj);           /* OUTPUT: number of objects           */


/* set and get the size of the pixel stack used in extract() */
void sep_set_extract_pixstack(size_t val);
size_t sep_get_extract_pixstack(void);

/* free memory associated with an sepobj array, including pixel lists */
void sep_freeobjarray(sepobj *objects, int nobj);

/*-------------------------- aperture photometry ----------------------------*/


/* Sum array values within a circular aperture.
 * 
 * Notes
 * -----
 * error : Can be a scalar (default), an array, or NULL
 *         If an array, set the flag SEP_ERROR_IS_ARRAY in `inflags`.
 *         Can represent 1-sigma std. deviation (default) or variance.
 *         If variance, set the flag SEP_ERROR_IS_VARIANCE in `inflags`.
 *
 * gain : If 0.0, poisson noise on sum is ignored when calculating error.
 *        Otherwise, (sum / gain) is added to the variance on sum.
 *
 * area : Total pixel area included in sum. Includes masked pixels that were
 *        corrected. The area can differ from the exact area of a circle due
 *        to inexact subpixel sampling and intersection with array boundaries.
 */
int sep_sum_circle(sep_image *image,
		   double maskthresh, /* pixel masked if mask > maskthresh */
		   short inflags,     /* input flags (see below) */
		   double x,          /* center of aperture in x */
		   double y,          /* center of aperture in y */
		   double r,          /* radius of aperture */
		   int subpix,        /* subpixel sampling */
		   double *sum,       /* OUTPUT: sum */
		   double *sumerr,    /* OUTPUT: error on sum */
		   double *area,      /* OUTPUT: area included in sum */
		   short *flag);      /* OUTPUT: flags */


int sep_sum_circann(sep_image *image, double maskthresh, short inflags,
		    double x, double y, double rin, double rout, int subpix,
		    double *sum, double *sumerr, double *area, short *flag);

int sep_sum_ellipse(sep_image *image, double maskthresh, short inflags,
		    double x, double y, double a, double b, double theta,
		    double r, int subpix,
		    double *sum, double *sumerr, double *area, short *flag);

int sep_sum_ellipann(sep_image *image, double maskthresh, short inflags,
		     double x, double y, double a, double b, double theta,
		     double rin, double rout, int subpix,
		     double *sum, double *sumerr, double *area, short *flag);

/* sep_sum_circann_multi()
 *
 * Sum an array of circular annuli more efficiently (but with no exact mode).
 *
 * Notable parameters:
 * 
 * rmax:     Input radii are  [rmax/n, 2*rmax/n, 3*rmax/n, ..., rmax].
 * n:        Length of input and output arrays.
 * sum:      Preallocated array of length n holding sums in annuli. sum[0]
 *           corrresponds to r=[0, rmax/n], sum[n-1] to outermost annulus.
 * sumvar:   Preallocated array of length n holding variance on sums.
 * area:     Preallocated array of length n holding area summed in each annulus.
 * maskarea: Preallocated array of length n holding masked area in each
             annulus (if mask not NULL).
 * flag:     Output flag (non-array).
 */
int sep_sum_circann_multi(sep_image *im, double maskthresh, short inflag,
			  double x, double y, double rmax, int n, int subpix,
			  double *sum, double *sumvar, double *area,
			  double *maskarea, short *flag);

/* sep_flux_radius()
 *
 * Calculate the radii enclosing the requested fraction of flux relative
 * to radius rmax. 
 *
 * (see previous functions for most arguments)
 * rmax : maximum radius to analyze
 * fluxtot : scale requested flux fractions to this. (If NULL, flux within
             `rmax` is used.)
 * fluxfrac : array of requested fractions.
 * n : length of fluxfrac
 * r : (output) result array of length n.
 * flag : (output) scalar flag
 */
int sep_flux_radius(sep_image *im, double maskthresh, short inflag,
		    double x, double y, double rmax, int subpix,
		    double *fluxtot, double *fluxfrac, int n,
		    double *r, short *flag);

/* sep_kron_radius()
 *
 * Calculate Kron radius within an ellipse given by 
 *
 *     cxx*(x'-x)^2 + cyy*(y'-y)^2 + cxy*(x'-x)*(y'-y) < r^2
 *
 * The Kron radius is sum(r_i * v_i) / sum(v_i) where v_i is the value of pixel
 * i and r_i is the "radius" of pixel i, as given by the left hand side of
 * the above equation.
 *
 * Flags that might be set:
 * SEP_APER_HASMASKED - at least one of the pixels in the ellipse is masked.
 * SEP_APER_ALLMASKED - All pixels in the ellipse are masked. kronrad = 0.
 * SEP_APER_NONPOSITIVE - There was a nonpositive numerator or deminator.
 *                        kronrad = 0.
 */
int sep_kron_radius(sep_image *im, double maskthresh, double x, double y,
		    double cxx, double cyy, double cxy, double r,
		    double *kronrad, short *flag);


/* Calculate "windowed" position parameters.
 *
 * This is an iterative procedure.
 *
 * x, y       : initial center
 * sig        : sigma of Gaussian to use for weighting. The integration
 *              radius is 4 * sig.
 * subpix     : Subpixels to use in aperture-pixel overlap.
 *              SExtractor uses 11. 0 is supported for exact overlap.
 * xout, yout : output center.
 * niter      : number of iterations used.
 */
int sep_windowed(sep_image *im, double maskthresh, short inflag,
                 double x, double y, double sig, int subpix,
                 double *xout, double *yout, int *niter, short *flag);



void sep_set_ellipse(unsigned char *arr, int w, int h,
		     double x, double y, double cxx, double cyy, double cxy,
		     double r, unsigned char val);
/* Set array elements within an ellipitcal aperture to a given value.
 *
 * Ellipse: cxx*(x'-x)^2 + cyy*(y'-y)^2 + cxy*(x'-x)*(y'-y) = r^2  
 */

int sep_ellipse_axes(double cxx, double cyy, double cxy,
		     double *a, double *b, double *theta);
void sep_ellipse_coeffs(double a, double b, double theta,
			double *cxx, double *cyy, double *cxy);

/*----------------------- info & error messaging ----------------------------*/

extern char *sep_version_string;
/* library version string (e.g., "0.2.0") */

void sep_get_errmsg(int status, char *errtext);
/* Return a short descriptive error message that corresponds to the input
 * error status value.  The message may be up to 60 characters long, plus
 * the terminating null character. */

void sep_get_errdetail(char *errtext);
/* Return a longer error message with more specifics about the problem.
   The message may be up to 512 characters */
