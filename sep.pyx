# This wrapper licensed under an MIT license.

"""
Source Extraction and Photometry

This module is a wrapper of the SEP C library.
"""
import numpy as np
cimport numpy as np
from libc cimport limits
from libc.math cimport sqrt
cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.version cimport PY_MAJOR_VERSION
from warnings import warn

np.import_array()  # To access the numpy C-API.

__version__ = "0.6.0"

# -----------------------------------------------------------------------------
# Definitions from the SEP C library

# macro definitions from sep.h
DEF SEP_TBYTE = 11
DEF SEP_TINT = 31
DEF SEP_TFLOAT = 42
DEF SEP_TDOUBLE = 82

# input flag values (C macros)
DEF SEP_NOISE_NONE = 0
DEF SEP_NOISE_STDDEV = 1
DEF SEP_NOISE_VAR = 2

# filter types for sep_extract
DEF SEP_FILTER_CONV = 0
DEF SEP_FILTER_MATCHED = 1

# Threshold types
DEF SEP_THRESH_REL = 0
DEF SEP_THRESH_ABS = 1 

# input flags for aperture photometry
DEF SEP_MASK_IGNORE = 0x0004

# Output flag values accessible from python
OBJ_MERGED = np.short(0x0001)
OBJ_TRUNC = np.short(0x0002)
OBJ_DOVERFLOW = np.short(0x0004)
OBJ_SINGU = np.short(0x0008)
APER_TRUNC = np.short(0x0010)
APER_HASMASKED = np.short(0x0020)
APER_ALLMASKED = np.short(0x0040)
APER_NONPOSITIVE = np.short(0x0080)

# macro defintion from sepcore.h
# This is not part of the SEP API, but we pull it out because want to
# explicitly detect memory errors so that we can raise MemoryError().
DEF MEMORY_ALLOC_ERROR = 1

# header definitions
cdef extern from "sep.h":

    ctypedef struct sep_image:
        void *data
        void *noise
        void *mask
        int dtype
        int ndtype
        int mdtype
        int w
        int h
        double noiseval
        short noise_type
        double gain
        double maskthresh

    ctypedef struct sep_bkg:
        int w
        int h
        float globalback
        float globalrms
    
    ctypedef struct sep_catalog:
        int    nobj
        float  *thresh
        int    *npix
        int    *tnpix
        int    *xmin
        int    *xmax
        int    *ymin
        int    *ymax
        double *x
        double *y
        double *x2
        double *y2
        double *xy
        float  *a
        float  *b
        float  *theta
        float  *cxx
        float  *cyy
        float  *cxy
        float  *cflux
        float  *flux
        float  *cpeak
        float  *peak
        int    *xcpeak
        int    *ycpeak
        int    *xpeak
        int    *ypeak
        short  *flag
        int    **pix
        int    *objectspix

    int sep_background(sep_image *im,
                       int bw, int bh,
                       int fw, int fh,
                       double fthresh,
                       sep_bkg **bkg)

    float sep_bkg_global(sep_bkg *bkg)
    float sep_bkg_globalrms(sep_bkg *bkg)
    int sep_bkg_array(sep_bkg *bkg, void *arr, int dtype)
    int sep_bkg_rmsarray(sep_bkg *bkg, void *arr, int dtype)
    int sep_bkg_subarray(sep_bkg *bkg, void *arr, int dtype)
    void sep_bkg_free(sep_bkg *bkg)

    int sep_extract(sep_image *image,
                    float thresh,
                    int thresh_type,
                    int minarea,
                    float *conv,
                    int convw, int convh,
                    int filter_type,
                    int deblend_nthresh,
                    double deblend_cont,
                    int clean_flag,
                    double clean_param,
                    sep_catalog **catalog)

    void sep_catalog_free(sep_catalog *catalog)

    int sep_sum_circle(sep_image *image,
                       double x, double y, double r, int subpix, short inflags,
                       double *sum, double *sumerr, double *area, short *flag)

    int sep_sum_circann(sep_image *image,
                        double x, double y, double rin, double rout,
                        int subpix, short inflags,
                        double *sum, double *sumerr, double *area, short *flag)

    int sep_sum_ellipse(sep_image *image,
                        double x, double y, double a, double b, double theta,
                        double r, int subpix, short inflags,
                        double *sum, double *sumerr, double *area,
                        short *flag)

    int sep_sum_ellipann(sep_image *image,
                         double x, double y, double a, double b,
                         double theta, double rin, double rout, int subpix,
                         short inflags,
                         double *sum, double *sumerr, double *area,
                         short *flag)

    int sep_flux_radius(sep_image *image,
                        double x, double y, double rmax, int subpix,
                        short inflag,
                        double *fluxtot, double *fluxfrac, int n,
                        double *r, short *flag)

    int sep_kron_radius(sep_image *image,
                        double x, double y, double cxx, double cyy,
                        double cxy, double r,
                        double *kronrad, short *flag)

    int sep_windowed(sep_image *image,
                     double x, double y, double sig, int subpix, short inflag,
                     double *xout, double *yout, int *niter, short *flag)

    int sep_ellipse_axes(double cxx, double cyy, double cxy,
                         double *a, double *b, double *theta)

    void sep_ellipse_coeffs(double a, double b, double theta,
                            double *cxx, double *cyy, double *cxy)

    void sep_set_ellipse(unsigned char *arr, int w, int h,
                         double x, double y,
                         double cxx, double cyy, double cxy, double r,
                         unsigned char val)

    void sep_set_extract_pixstack(size_t val)
    size_t sep_get_extract_pixstack()

    void sep_get_errmsg(int status, char *errtext)
    void sep_get_errdetail(char *errtext)

# -----------------------------------------------------------------------------
# Utility functions

cdef int _get_sep_dtype(dtype) except -1:
    """Convert a numpy dtype to the corresponding SEP dtype integer code."""
    if not dtype.isnative:
        raise ValueError(
            "Input array with dtype '{0}' has non-native byte order. "
            "Only native byte order arrays are supported. "
            "To change the byte order of the array 'data', do "
            "'data = data.byteswap().newbyteorder()'".format(dtype))
    t = dtype.type
    if t is np.single:
        return SEP_TFLOAT
    elif t is np.bool_ or t is np.ubyte:
        return SEP_TBYTE
    elif dtype == np.double:
        return SEP_TDOUBLE
    elif dtype == np.intc:
        return SEP_TINT
    raise ValueError('input array dtype not supported: {0}'.format(dtype))


cdef int _check_array_get_dims(np.ndarray arr, int *w, int *h) except -1:
    """Check some things about an array and return dimensions"""

    # Raise an informative message if array is not C-contiguous
    if not arr.flags["C_CONTIGUOUS"]:
        raise ValueError("array is not C-contiguous")

    # Check that there are exactly 2 dimensions
    if arr.ndim != 2:
        raise ValueError("array must be 2-d")

    # ensure that arr dimensions are not too large for C ints.
    if arr.shape[0] <= <Py_ssize_t> limits.INT_MAX:
        h[0] = arr.shape[0]
    else:
        raise ValueError("array height  ({0:d}) greater than INT_MAX ({1:d})"
                         .format(arr.shape[0], limits.INT_MAX))
    if arr.shape[1] <= <Py_ssize_t> limits.INT_MAX:
       w[0] = arr.shape[1]
    else:
        raise ValueError("array width ({0:d}) greater than INT_MAX ({1:d})"
                         .format(arr.shape[1], limits.INT_MAX))
    return 0

cdef int _assert_ok(int status) except -1:
    """Get the SEP error message corresponding to status code"""
    cdef char *errmsg
    cdef char *errdetail

    if status == 0:
        return 0

    # First check if we have an out-of-memory error, so we don't try to
    # allocate more memory to hold the error message.
    if status == MEMORY_ALLOC_ERROR:
        raise MemoryError

    # Otherwise, get error message.
    errmsg = <char *>PyMem_Malloc(61 * sizeof(char))
    sep_get_errmsg(status, errmsg)
    pyerrmsg = <bytes> errmsg
    PyMem_Free(errmsg)

    # Get error detail.
    errdetail = <char *>PyMem_Malloc(512 * sizeof(char))
    sep_get_errdetail(errdetail)
    pyerrdetail = <bytes> errdetail
    PyMem_Free(errdetail)

    # If error detail is present, append it to the message.
    if pyerrdetail != b"":
        pyerrmsg = pyerrmsg + b": " + pyerrdetail

    # Convert string to unicode if on python 3
    if PY_MAJOR_VERSION == 3:
        msg = pyerrmsg.decode()
    else:
        msg = pyerrmsg

    raise Exception(msg)


cdef int _parse_arrays(np.ndarray data, err, var, mask,
                       sep_image *im) except -1:
    """Helper function for functions accepting data, error & mask arrays.
    Fills in an sep_image struct."""

    cdef int ew, eh, mw, mh
    cdef np.uint8_t[:,:] buf, ebuf, mbuf

    # Clear im fields we might not touch (everything besides data, dtype, w, h)
    im.noise = NULL
    im.mask = NULL
    im.ndtype = 0
    im.mdtype = 0
    im.noiseval = 0.0
    im.noise_type = SEP_NOISE_NONE
    im.gain = 0.0
    im.maskthresh = 0.0

    # Get main image info
    _check_array_get_dims(data, &(im.w), &(im.h))
    im.dtype = _get_sep_dtype(data.dtype)
    buf = data.view(dtype=np.uint8)
    im.data = <void*>&buf[0, 0]

    # Check if noise is error or variance.
    if err is not None:
        if var is not None:
            raise ValueError("Cannot specify both err and var")
        im.noise_type = SEP_NOISE_STDDEV
    elif var is not None:
        err = var
        im.noise_type = SEP_NOISE_VAR

    # parse noise
    if err is None:
        im.noise = NULL
        im.noise_type = SEP_NOISE_NONE
        im.noiseval = 0.0
    elif isinstance(err, np.ndarray):
        if err.ndim == 0:
            im.noise = NULL
            im.noiseval = err
        elif err.ndim == 2:
            _check_array_get_dims(err, &ew, &eh)
            if ew != im.w or eh != im.h:
                raise ValueError("size of error/variance array must match"
                                 " data")
            im.ndtype = _get_sep_dtype(err.dtype)
            ebuf = err.view(dtype=np.uint8)
            im.noise = <void*>&ebuf[0, 0]
        else:
            raise ValueError("error/variance array must be 0-d or 2-d")
    else:
        im.noise = NULL
        im.noiseval = err

    # Optional input: mask
    if mask is None:
        im.mask = NULL
    else:
        _check_array_get_dims(mask, &mw, &mh)
        if mw != im.w or mh != im.h:
            raise ValueError("size of mask array must match data")
        im.mdtype = _get_sep_dtype(mask.dtype)
        mbuf = mask.view(dtype=np.uint8)
        im.mask = <void*>&mbuf[0, 0]

# -----------------------------------------------------------------------------
# Background Estimation

cdef class Background:
    """
    Background(data, mask=None, maskthresh=0.0, bw=64, bh=64,
               fw=3, fh=3, fthresh=0.0)

    Representation of spatially variable image background and noise.

    Parameters
    ----------
    data : 2-d `~numpy.ndarray`
        Data array.
    mask : 2-d `~numpy.ndarray`, optional
        Mask array, optional
    maskthresh : float, optional
        Mask threshold. This is the inclusive upper limit on the mask value
        in order for the corresponding pixel to be unmasked. For boolean
        arrays, False and True are interpreted as 0 and 1, respectively.
        Thus, given a threshold of zero, True corresponds to masked and
        False corresponds to unmasked.
    bw, bh : int, optional
        Size of background boxes in pixels. Default is 64.
    fw, fh : int, optional
        Filter width and height in boxes. Default is 3.
    fthresh : float, optional
        Filter threshold. Default is 0.0.
    """

    cdef sep_bkg *ptr      # pointer to C struct
    cdef np.dtype orig_dtype  # dtype code of original image

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __cinit__(self, np.ndarray data not None, np.ndarray mask=None,
                  float maskthresh=0.0, int bw=64, int bh=64,
                  int fw=3, int fh=3, float fthresh=0.0):

        cdef int status
        cdef sep_image im

        _parse_arrays(data, None, None, mask, &im)
        im.maskthresh = maskthresh
        status = sep_background(&im, bw, bh, fw, fh, fthresh, &self.ptr)
        _assert_ok(status)

        self.orig_dtype = data.dtype

    # Note: all initialization work is done in __cinit__. This is just here
    # for the docstring.
    def __init__(self, np.ndarray data not None, np.ndarray mask=None,
                 float maskthresh=0.0, int bw=64, int bh=64,
                 int fw=3, int fh=3, float fthresh=0.0):
        """Background(data, mask=None, maskthresh=0.0, bw=64, bh=64,
                      fw=3, fh=3, fthresh=0.0)"""
        pass

    property globalback:
        """Global background level."""
        def __get__(self):
            return sep_bkg_global(self.ptr)

    property globalrms:
        """Global background RMS."""
        def __get__(self):
            return sep_bkg_globalrms(self.ptr)

    def back(self, dtype=None):
        """back(dtype=None)

        Create an array of the background.

        Parameters
        ----------
        dtype : `~numpy.dtype`, optional
             Data type of output array. Default is the dtype of the original
             data.

        Returns
        -------
        back : `~numpy.ndarray`
            Array with same dimensions as original data.
        """
        cdef int sep_dtype
        cdef np.uint8_t[:, :] buf

        if dtype is None:
            dtype = self.orig_dtype
        else:
            dtype = np.dtype(dtype)
        sep_dtype = _get_sep_dtype(dtype)

        result = np.empty((self.ptr.h, self.ptr.w), dtype=dtype)
        buf = result.view(dtype=np.uint8)        
        status = sep_bkg_array(self.ptr, &buf[0, 0], sep_dtype)
        _assert_ok(status)

        return result

    def rms(self, dtype=None):
        """rms(dtype=None)

        Create an array of the background rms.

        Parameters
        ----------
        dtype : `~numpy.dtype`, optional
             Data type of output array. Default is the dtype of the original
             data.

        Returns
        -------
        rms : `~numpy.ndarray`
            Array with same dimensions as original data.
        """
        cdef int sep_dtype
        cdef np.uint8_t[:, :] buf

        if dtype is None:
            dtype = self.orig_dtype
        else:
            dtype = np.dtype(dtype)
        sep_dtype = _get_sep_dtype(dtype)

        result = np.empty((self.ptr.h, self.ptr.w), dtype=dtype)
        buf = result.view(dtype=np.uint8)
        status = sep_bkg_rmsarray(self.ptr, &buf[0, 0], sep_dtype)
        _assert_ok(status)

        return result


    def subfrom(self, np.ndarray data not None):
        """subfrom(data)

        Subtract the background from an existing array.

        Like ``data = data - bkg``, but avoids making a copy of the data.

        Parameters
        ----------
        data : `~numpy.ndarray`
            Input array, which will be updated in-place. Shape must match
            that of the original image used to measure the background. 
        """

        cdef int w, h, status, sep_dtype
        cdef np.uint8_t[:, :] buf

        assert self.ptr is not NULL

        _check_array_get_dims(data, &w, &h)
        sep_dtype = _get_sep_dtype(data.dtype)
        buf = data.view(dtype=np.uint8)

        # ensure dimensions match original image
        if (w != self.ptr.w or h != self.ptr.h):
            raise ValueError("Data dimensions do not match background "
                             "dimensions")

        status = sep_bkg_subarray(self.ptr, &buf[0, 0], sep_dtype)
        _assert_ok(status)

    def __array__(self, dtype=None):
        return self.back(dtype=dtype)

    def __rsub__(self, np.ndarray data not None):
        data = np.copy(data)
        self.subfrom(data)
        return data

    def __dealloc__(self):
        if self.ptr is not NULL:
            sep_bkg_free(self.ptr)

# -----------------------------------------------------------------------------
# Source Extraction

# This needs to match the result from extract
cdef packed struct Object:
    np.float64_t thresh
    np.int_t npix
    np.int_t tnpix
    np.int_t xmin
    np.int_t xmax
    np.int_t ymin
    np.int_t ymax
    np.float64_t x
    np.float64_t y
    np.float64_t x2
    np.float64_t y2
    np.float64_t xy
    np.float64_t a
    np.float64_t b
    np.float64_t theta
    np.float64_t cxx
    np.float64_t cyy
    np.float64_t cxy
    np.float64_t cflux
    np.float64_t flux
    np.float64_t cpeak
    np.float64_t peak
    np.int_t xcpeak
    np.int_t ycpeak
    np.int_t xpeak
    np.int_t ypeak
    np.int_t flag

default_kernel = np.array([[1.0, 2.0, 1.0],
                           [2.0, 4.0, 2.0],
                           [1.0, 2.0, 1.0]], dtype=np.float32)

def extract(np.ndarray data not None, float thresh, np.ndarray err=None,
            np.ndarray mask=None, int minarea=5,
            np.ndarray filter_kernel=default_kernel, filter_type='matched',
            int deblend_nthresh=32, double deblend_cont=0.005,
            bint clean=True, double clean_param=1.0,
            segmentation_map=False, np.ndarray conv=default_kernel):
    """extract(data, thresh, err=None, mask=None, minarea=5,
               filter_kernel=default_kernel, filter_type='matched',
               deblend_nthresh=32, deblend_cont=0.005, clean=True,
               clean_param=1.0, segmentation_map=False)

    Extract sources from an image.

    Parameters
    ----------
    data : `~numpy.ndarray`
        Data array (2-d).
    thresh : float
        Threshold pixel value for detection. If an ``err`` array is not given,
        this is interpreted as an absolute threshold. If ``err`` is
        given, this is interpreted as a relative threshold: the absolute
        threshold at pixel (j, i) will be ``thresh * err[j, i]``.
    err : `~numpy.ndarray`, optional
        Noise array for specifying a pixel-by-pixel detection threshold.
    mask : `~numpy.ndarray`, optional
        Mask array. ``True`` values, or numeric values greater than 0,
        are considered masked. Masking a pixel is equivalent to setting data
        to zero and noise (if present) to infinity.
    minarea : int, optional
        Minimum number of pixels required for an object. Default is 5.
    filter_kernel : `~numpy.ndarray` or None, optional
        Filter kernel used for on-the-fly filtering (used to
        enhance detection). Default is a 3x3 array:
        [[1,2,1], [2,4,2], [1,2,1]]. Set to ``None`` to skip
        convolution.
    filter_type : {'matched', 'conv'}, optional
        Filter treatment. This affects filtering behavior when a noise
        array is supplied. ``'matched'`` (default) accounts for
        pixel-to-pixel noise in the filter kernel. ``'conv'`` is
        simple convolution of the data array, ignoring pixel-to-pixel
        noise across the kernel.  ``'matched'`` should yield better
        detection of faint sources in areas of rapidly varying noise
        (such as found in coadded images made from semi-overlapping
        exposures).  The two options are equivalent when noise is
        constant.
    deblend_nthresh : int, optional
        Number of thresholds used for object deblending. Default is 32.
    deblend_cont : float, optional
        Minimum contrast ratio used for object deblending. Default is 0.005.
        To entirely disable deblending, set to 1.0.
    clean : bool, optional
        Perform cleaning? Default is True.
    clean_param : float, optional
        Cleaning parameter (see SExtractor manual). Default is 1.0.
    segmentation_map : bool, optional
        If True, also return a "segmentation map" giving the member
        pixels of each object. Default is False.

    Returns
    -------
    objects : `~numpy.ndarray`
        Extracted object parameters (structured array). Available fields are:

        * ``thresh`` (float) Threshold at object location.
        * ``npix`` (int) Number of pixels belonging to the object.
        * ``tnpix`` (int) Number of pixels above threshold (unconvolved data).
        * ``xmin``, ``xmax`` (int) Minimum, maximum x coordinates of pixels.
        * ``ymin``, ``ymax`` (int) Minimum, maximum y coordinates of pixels.
        * ``x``, ``y`` (float) object barycenter (first moments).
        * ``x2``, ``y2``, ``xy`` (float) Second moments.
        * ``a``, ``b``, ``theta`` (float) Ellipse parameters, scalled as
            described by Section 8.4.2 in "The Source Extractor Guide" or 
            Section 10.1.5-6 of v2.13 of SExtractor's User Manual.
        * ``cxx``, ``cyy``, ``cxy`` (float) Alternative ellipse parameters.
        * ``cflux`` (float) Sum of member pixels in convolved data.
        * ``flux`` (float) Sum of member pixels in unconvolved data.
        * ``cpeak`` (float) Peak value in convolved data.
        * ``peak`` (float) Peak value in unconvolved data.
        * ``xcpeak``, ``ycpeak`` (int) Coordinate of convolved peak pixel.
        * ``xpeak``, ``ypeak`` (int) Coordinate of convolved peak pixel.
        * ``flag`` (int) Extraction flags.

    segmap : `~numpy.ndarray`, optional
        Array of integers with same shape as data. Pixels not belonging to
        any object have value 0. All pixels belonging to the ``i``-th object
        (e.g., ``objects[i]``) have value ``i+1``. Only returned if
        ``segmentation_map=True``.
    """

    cdef int kernelw, kernelh, status, i, j
    cdef int filter_typecode, thresh_type
    cdef sep_catalog *catalog = NULL
    cdef np.ndarray[Object] result
    cdef float[:, :] kernelflt
    cdef float *kernelptr
    cdef np.int32_t[:, :] segmap_buf
    cdef np.int32_t *segmap_ptr
    cdef int *objpix
    cdef sep_image im

    # parse arrays
    _parse_arrays(data, err, None, mask, &im)

    # 'conv' has been renamed to filter_kernel. If the user has set it
    # explicitly, issue a warning. Don't use DeprecationWarning: no one will
    # ever see it.
    if conv is not default_kernel:
        warn("The 'conv' keyword argument is deprecated. Use the "
             "'filter_kernel' keyword argument instead.")
        if filter_kernel is default_kernel:
            filter_kernel = conv

    # Parse filter input
    if filter_kernel is None:
        kernelptr = NULL
        kernelw = 0
        kernelh = 0
    else:
        kernelflt = filter_kernel.astype(np.float32)
        kernelptr = &kernelflt[0, 0]
        kernelw = kernelflt.shape[1]
        kernelh = kernelflt.shape[0]

    if filter_type == 'matched':
        filter_typecode = SEP_FILTER_MATCHED
    elif filter_type == 'conv':
        filter_typecode = SEP_FILTER_CONV
    else:
        raise ValueError("unknown filter_type: {!r}".format(filter_type))

    # If there is an error array, the threshold is relative, otherwise
    # it is absolute
    if im.noise == NULL:
        thresh_type = SEP_THRESH_ABS
    else:
        thresh_type = SEP_THRESH_REL

    status = sep_extract(&im,
                         thresh, thresh_type, minarea,
                         kernelptr, kernelw, kernelh, filter_typecode,
                         deblend_nthresh, deblend_cont, clean, clean_param,
                         &catalog)
    _assert_ok(status)

    # Allocate result record array and fill it
    result = np.empty(catalog.nobj,
                      dtype=np.dtype([('thresh', np.float64),
                                      ('npix', np.int),
                                      ('tnpix', np.int),
                                      ('xmin', np.int),
                                      ('xmax', np.int),
                                      ('ymin', np.int),
                                      ('ymax', np.int),
                                      ('x', np.float64),
                                      ('y', np.float64),
                                      ('x2', np.float64),
                                      ('y2', np.float64),
                                      ('xy', np.float64),
                                      ('a', np.float64),
                                      ('b', np.float64),
                                      ('theta', np.float64),
                                      ('cxx', np.float64),
                                      ('cyy', np.float64),
                                      ('cxy', np.float64),
                                      ('cflux', np.float64),
                                      ('flux', np.float64),
                                      ('cpeak', np.float64),
                                      ('peak', np.float64),
                                      ('xcpeak', np.int),
                                      ('ycpeak', np.int),
                                      ('xpeak', np.int),
                                      ('ypeak', np.int),
                                      ('flag', np.int)]))

    for i in range(catalog.nobj):
        result['thresh'][i] = catalog.thresh[i]
        result['npix'][i] = catalog.npix[i]
        result['tnpix'][i] = catalog.tnpix[i]
        result['xmin'][i] = catalog.xmin[i]
        result['xmax'][i] = catalog.xmax[i]
        result['ymin'][i] = catalog.ymin[i]
        result['ymax'][i] = catalog.ymax[i]
        result['x'][i] = catalog.x[i]
        result['y'][i] = catalog.y[i]
        result['x2'][i] = catalog.x2[i]
        result['y2'][i] = catalog.y2[i]
        result['xy'][i] = catalog.xy[i]
        result['a'][i] = catalog.a[i]
        result['b'][i] = catalog.b[i]
        result['theta'][i] = catalog.theta[i]
        result['cxx'][i] = catalog.cxx[i]
        result['cyy'][i] = catalog.cyy[i]
        result['cxy'][i] = catalog.cxy[i]
        result['cflux'][i] = catalog.cflux[i]
        result['flux'][i] = catalog.flux[i]
        result['cpeak'][i] = catalog.cpeak[i]
        result['peak'][i] = catalog.peak[i]
        result['xcpeak'][i] = catalog.xcpeak[i]
        result['ycpeak'][i] = catalog.ycpeak[i]
        result['xpeak'][i] = catalog.xpeak[i]
        result['ypeak'][i] = catalog.ypeak[i]
        result['flag'][i] = catalog.flag[i]

    # construct a segmentation map, if it was requested.
    if segmentation_map:
        # Note: We have to write out `(data.shape[0], data.shape[1])` because
        # because Cython turns `data.shape` later into an int pointer when
        # the function argument is typed as np.ndarray.
        segmap = np.zeros((data.shape[0], data.shape[1]), dtype=np.int32)
        segmap_buf = segmap
        segmap_ptr = &segmap_buf[0, 0]
        for i in range(catalog.nobj):
            objpix = catalog.pix[i]
            for j in range(catalog.npix[i]):
                segmap_ptr[objpix[j]] = i + 1

    # Free the C catalog
    sep_catalog_free(catalog)

    if segmentation_map:
        return result, segmap
    else:
        return result

# -----------------------------------------------------------------------------
# Aperture Photometry

@cython.boundscheck(False)
@cython.wraparound(False)
def sum_circle(np.ndarray data not None, x, y, r,
               var=None, err=None, gain=None, np.ndarray mask=None,
               double maskthresh=0.0, bkgann=None, int subpix=5):
    """sum_circle(data, x, y, r, err=None, var=None, mask=None, maskthresh=0.0,
                  bkgann=None, gain=None, subpix=5)

    Sum data in circular aperture(s).

    Parameters
    ----------
    data : `~numpy.ndarray`
        2-d array to be summed.

    x, y, r : array_like
        Center coordinates and radius (radii) of aperture(s). 
        ``x`` corresponds to the second ("fast") axis of the input array
        and ``y`` corresponds to the first ("slow") axis.
        ``x, y = (0.0, 0.0)`` corresponds to the center of the first
        element of the array. These inputs obey numpy broadcasting rules.

    err, var : float or `~numpy.ndarray`
        Error *or* variance (specify at most one).

    mask : `~numpy.ndarray`, optional
        Mask array. If supplied, a given pixel is masked if its value
        is greater than ``maskthresh``.

    maskthresh : float, optional
        Threshold for a pixel to be masked. Default is ``0.0``.

    bkgann : tuple, optional
        Length 2 tuple giving the inner and outer radius of a
        "background annulus". If supplied, the background is estimated
        by averaging unmasked pixels in this annulus. If supplied, the inner
        and outer radii obey numpy broadcasting rules along with ``x``,
        ``y`` and ``r``.

    gain : float, optional
        Conversion factor between data array units and poisson counts,
        used in calculating poisson noise in aperture sum. If ``None``
        (default), do not add poisson noise.

    subpix : int, optional
        Subpixel sampling factor. If 0, exact overlap is calculated.
        Default is 5.

    Returns
    -------
    sum : `~numpy.ndarray`
        The sum of the data array within the aperture.

    sumerr : `~numpy.ndarray`
        Error on the sum.

    flags : `~numpy.ndarray`
        Integer giving flags. (0 if no flags set.)
    """

    cdef double flux1, fluxerr1, area1,
    cdef double bkgflux, bkgfluxerr, bkgarea
    cdef short flag1, bkgflag
    cdef size_t i
    cdef int status
    cdef np.broadcast it
    cdef sep_image im

    _parse_arrays(data, err, var, mask, &im)
    im.maskthresh = maskthresh
    if gain is not None:
        im.gain = gain

    # Require that inputs are float64 arrays. This has to be done because we
    # are using a broadcasting iterator below, where we need to know the type
    # in advance. There are other ways to do this, e.g., using NpyIter_Multi
    # in the numpy C-API. However, the best way to use this from cython
    # is not clear to me at this time.
    #
    # docs.scipy.org/doc/numpy/reference/c-api.iterator.html#NpyIter_MultiNew
    dt = np.dtype(np.double)
    x = np.require(x, dtype=dt)
    y = np.require(y, dtype=dt)
    r = np.require(r, dtype=dt)

    if bkgann is None:

        # allocate ouput arrays
        shape = np.broadcast(x, y, r).shape
        sum = np.empty(shape, dt)
        sumerr = np.empty(shape, dt)
        flag = np.empty(shape, np.short)

        it = np.broadcast(x, y, r, sum, sumerr, flag)
        while np.PyArray_MultiIter_NOTDONE(it):
            status = sep_sum_circle(
                &im,
                (<double*>np.PyArray_MultiIter_DATA(it, 0))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 1))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 2))[0],
                subpix, 0,
                <double*>np.PyArray_MultiIter_DATA(it, 3),
                <double*>np.PyArray_MultiIter_DATA(it, 4),
                &area1,
                <short*>np.PyArray_MultiIter_DATA(it, 5))
            _assert_ok(status)

            # Advance the iterator
            np.PyArray_MultiIter_NEXT(it)

        return sum, sumerr, flag

    else:
        rin, rout = bkgann

        # Require float arrays (see note above)
        rin = np.require(rin, dtype=dt)
        rout = np.require(rout, dtype=dt)

        # allocate ouput arrays
        shape = np.broadcast(x, y, r, rin, rout).shape
        sum = np.empty(shape, dt)
        sumerr = np.empty(shape, dt)
        flag = np.empty(shape, np.short)

        it = np.broadcast(x, y, r, rin, rout, sum, sumerr, flag)
        while np.PyArray_MultiIter_NOTDONE(it):
            status = sep_sum_circle(
                &im,
                (<double*>np.PyArray_MultiIter_DATA(it, 0))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 1))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 2))[0],
                subpix, 0, &flux1, &fluxerr1, &area1, &flag1)
            _assert_ok(status)
                
            # background subtraction
            # Note that background output flags are not used.
            status = sep_sum_circann(
                &im,
                (<double*>np.PyArray_MultiIter_DATA(it, 0))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 1))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 3))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 4))[0],
                1, SEP_MASK_IGNORE, &bkgflux, &bkgfluxerr, &bkgarea, &bkgflag)
            _assert_ok(status)

            flux1 -= bkgflux / bkgarea * area1
            bkgfluxerr = bkgfluxerr / bkgarea * area1
            fluxerr1 = sqrt(fluxerr1*fluxerr1 + bkgfluxerr*bkgfluxerr)
            (<double*>np.PyArray_MultiIter_DATA(it, 5))[0] = flux1
            (<double*>np.PyArray_MultiIter_DATA(it, 6))[0] = fluxerr1
            (<short*>np.PyArray_MultiIter_DATA(it, 7))[0] = flag1

            np.PyArray_MultiIter_NEXT(it)

        return sum, sumerr, flag

@cython.boundscheck(False)
@cython.wraparound(False)
def sum_circann(np.ndarray data not None, x, y, rin, rout,
                var=None, err=None, gain=None, np.ndarray mask=None,
                double maskthresh=0.0, int subpix=5):
    """sum_circann(data, x, y, rin, rout, err=None, var=None, mask=None,
                   maskthresh=0.0, gain=None, subpix=5)

    Sum data in circular annular aperture(s).

    Parameters
    ----------
    data : `~numpy.ndarray`
        2-d array to be summed.

    x, y, rin, rout : array_like
        Center coordinates and inner and outer radii of aperture(s). 
        ``x`` corresponds to the second ("fast") axis of the input array
        and ``y`` corresponds to the first ("slow") axis.
        ``x, y = (0.0, 0.0)`` corresponds to the center of the first
        element of the array. These inputs obey numpy broadcasting rules.
        It is required that ``rout >= rin >= 0.0``.

    err, var : float or ndarray
        Error *or* variance (specify at most one).

    mask : `~numpy.ndarray`, optional
        Mask array. If supplied, a given pixel is masked if its value
        is greater than ``maskthresh``.

    maskthresh : float, optional
        Threshold for a pixel to be masked. Default is ``0.0``.

    gain : float, optional
        Conversion factor between data array units and poisson counts,
        used in calculating poisson noise in aperture sum. If ``None``
        (default), do not add poisson noise.

    subpix : int, optional
        Subpixel sampling factor. Default is 5.

    Returns
    -------
    sum : `~numpy.ndarray`
        The sum of the data array within the aperture.

    sumerr : `~numpy.ndarray`
        Error on the sum.

    flags : `~numpy.ndarray`
        Integer giving flags. (0 if no flags set.)
    """

    cdef double area1
    cdef size_t i
    cdef int status
    cdef np.broadcast it
    cdef sep_image im

    _parse_arrays(data, err, var, mask, &im)
    im.maskthresh = maskthresh
    if gain is not None:
        im.gain = gain

    # convert inputs to double arrays
    dt = np.dtype(np.double)
    x = np.require(x, dtype=dt)
    y = np.require(y, dtype=dt)
    rin = np.require(rin, dtype=dt)
    rout = np.require(rout, dtype=dt)

    # allocate ouput arrays
    shape = np.broadcast(x, y, rin, rout).shape
    sum = np.empty(shape, dt)
    sumerr = np.empty(shape, dt)
    flag = np.empty(shape, np.short)

    it = np.broadcast(x, y, rin, rout, sum, sumerr, flag)
    while np.PyArray_MultiIter_NOTDONE(it):
        status = sep_sum_circann(
            &im,
            (<double*>np.PyArray_MultiIter_DATA(it, 0))[0],
            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0],
            (<double*>np.PyArray_MultiIter_DATA(it, 2))[0],
            (<double*>np.PyArray_MultiIter_DATA(it, 3))[0],
            subpix, 0,
            <double*>np.PyArray_MultiIter_DATA(it, 4),
            <double*>np.PyArray_MultiIter_DATA(it, 5),
            &area1,
            <short*>np.PyArray_MultiIter_DATA(it, 6))

        _assert_ok(status)

        np.PyArray_MultiIter_NEXT(it)

    return sum, sumerr, flag


def sum_ellipse(np.ndarray data not None, x, y, a, b, theta, r=1.0,
                var=None, err=None, gain=None, np.ndarray mask=None,
                double maskthresh=0.0, bkgann=None, int subpix=5):
    """sum_ellipse(data, x, y, a, b, theta, r, err=None, var=None, mask=None,
                   maskthresh=0.0, bkgann=None, gain=None, subpix=5)

    Sum data in elliptical aperture(s).

    Parameters
    ----------
    data : `~numpy.ndarray`
        2-d array to be summed.

    x, y : array_like
        Center coordinates and radius (radii) of aperture(s). 
        ``x`` corresponds to the second ("fast") axis of the input array
        and ``y`` corresponds to the first ("slow") axis.
        ``x, y = (0.0, 0.0)`` corresponds to the center of the first
        element of the array. These inputs obey numpy broadcasting rules.

    a, b, theta : array_like
        Ellipse parameters. These inputs, along with ``x``, ``y``, and ``r``,
        obey numpy broadcasting rules. ``a`` is the semi-major axis,
        ``b`` is the semi-minor axis and ``theta`` is angle in radians between
        the positive x axis and the major axis. It must be in the range
        ``[-pi/2, pi/2]``. It is also required that ``a >= b >= 0.0``. 

    r : array_like, optional
        Scaling factor for the semi-minor and semi-major axes. The
        actual ellipse used will have semi-major axis ``a * r`` and
        semi-minor axis ``b * r``. Setting this parameter to a value
        other than 1.0 is exactly equivalent to scaling both ``a`` and
        ``b`` by the same value. Default is 1.0.

    err, var : float or `~numpy.ndarray`
        Error *or* variance (specify at most one).

    mask : `~numpy.ndarray`, optional
        Mask array. If supplied, a given pixel is masked if its value
        is greater than ``maskthresh``.

    maskthresh : float, optional
        Threshold for a pixel to be masked. Default is ``0.0``.

    bkgann : tuple, optional
        Length 2 tuple giving the inner and outer radius of a
        "background annulus". If supplied, the background is estimated
        by averaging unmasked pixels in this annulus. If supplied, the inner
        and outer radii obey numpy broadcasting rules, along with ``x``,
        ``y``, and ellipse parameters.

    gain : float, optional
        Conversion factor between data array units and poisson counts,
        used in calculating poisson noise in aperture sum. If ``None``
        (default), do not add poisson noise.

    subpix : int, optional
        Subpixel sampling factor. Default is 5.

    Returns
    -------
    sum : `~numpy.ndarray`
        The sum of the data array within the aperture.

    sumerr : `~numpy.ndarray`
        Error on the sum.

    flags : `~numpy.ndarray`
        Integer giving flags. (0 if no flags set.)

    """

    cdef double flux1, fluxerr1, x1, y1, r1, area1, rin1, rout1
    cdef double bkgflux, bkgfluxerr, bkgarea
    cdef short flag1, bkgflag
    cdef size_t i
    cdef int status
    cdef np.broadcast it
    cdef sep_image im

    _parse_arrays(data, err, var, mask, &im)
    im.maskthresh = maskthresh
    if gain is not None:
        im.gain = gain

    # Require that inputs are float64 arrays. See note in circular aperture.
    dt = np.dtype(np.double)
    x = np.require(x, dtype=dt)
    y = np.require(y, dtype=dt)
    a = np.require(a, dtype=dt)
    b = np.require(b, dtype=dt)
    theta = np.require(theta, dtype=dt)
    r = np.require(r, dtype=dt)

    if bkgann is None:

        # allocate ouput arrays
        shape = np.broadcast(x, y, a, b, theta, r).shape
        sum = np.empty(shape, dt)
        sumerr = np.empty(shape, dt)
        flag = np.empty(shape, np.short)

        it = np.broadcast(x, y, a, b, theta, r, sum, sumerr, flag)
        while np.PyArray_MultiIter_NOTDONE(it):
            status = sep_sum_ellipse(
                &im,
                (<double*>np.PyArray_MultiIter_DATA(it, 0))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 1))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 2))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 3))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 4))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 5))[0],
                subpix, 0,
                <double*>np.PyArray_MultiIter_DATA(it, 6),
                <double*>np.PyArray_MultiIter_DATA(it, 7),
                &area1,
                <short*>np.PyArray_MultiIter_DATA(it, 8))
            _assert_ok(status)

            np.PyArray_MultiIter_NEXT(it)

        return sum, sumerr, flag

    else:
        rin, rout = bkgann

        # Require float arrays (see note above)
        rin = np.require(rin, dtype=dt)
        rout = np.require(rout, dtype=dt)

        # allocate ouput arrays
        shape = np.broadcast(x, y, a, b, theta, r, rin, rout).shape
        sum = np.empty(shape, dt)
        sumerr = np.empty(shape, dt)
        flag = np.empty(shape, np.short)

        it = np.broadcast(x, y, a, b, theta, r, rin, rout, sum, sumerr, flag)
        while np.PyArray_MultiIter_NOTDONE(it):
            status = sep_sum_ellipse(
                &im,
                (<double*>np.PyArray_MultiIter_DATA(it, 0))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 1))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 2))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 3))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 4))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 5))[0],
                subpix, 0, &flux1, &fluxerr1, &area1, &flag1)
            _assert_ok(status)

            status = sep_sum_ellipann(
                &im,
                (<double*>np.PyArray_MultiIter_DATA(it, 0))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 1))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 2))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 3))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 4))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 6))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 7))[0],
                subpix, 0, &bkgflux, &bkgfluxerr, &bkgarea, &bkgflag)
            _assert_ok(status)

            flux1 -= bkgflux / bkgarea * area1
            bkgfluxerr = bkgfluxerr / bkgarea * area1
            fluxerr1 = sqrt(fluxerr1*fluxerr1 + bkgfluxerr*bkgfluxerr)

            (<double*>np.PyArray_MultiIter_DATA(it, 8))[0] = flux1
            (<double*>np.PyArray_MultiIter_DATA(it, 9))[0] = fluxerr1
            (<short*>np.PyArray_MultiIter_DATA(it, 10))[0] = flag1

            #PyArray_MultiIter_NEXT is used to advance the iterator
            np.PyArray_MultiIter_NEXT(it)

        return sum, sumerr, flag


@cython.boundscheck(False)
@cython.wraparound(False)
def sum_ellipann(np.ndarray data not None, x, y, a, b, theta, rin, rout,
                 var=None, err=None, gain=None, np.ndarray mask=None,
                 double maskthresh=0.0, int subpix=5):
    """sum_ellipann(data, x, y, a, b, theta, rin, rout, err=None, var=None,
                    mask=None, maskthresh=0.0, gain=None, subpix=5)

    Sum data in elliptical annular aperture(s).

    Parameters
    ----------
    data : `~numpy.ndarray`
        2-d array to be summed.

    x, y : array_like
        Center coordinates and radius (radii) of aperture(s). 
        ``x`` corresponds to the second ("fast") axis of the input array
        and ``y`` corresponds to the first ("slow") axis.
        ``x, y = (0.0, 0.0)`` corresponds to the center of the first
        element of the array. These inputs obey numpy broadcasting rules.

    a, b, theta, rin, rout : array_like
        Elliptical annulus parameters. These inputs, along with ``x`` and ``y``,
        obey numpy broadcasting rules. ``a`` is the semi-major axis,
        ``b`` is the semi-minor axis and ``theta`` is angle in radians between
        the positive x axis and the major axis. It must be in the range
        ``[-pi/2, pi/2]``. It is also required that ``a >= b >= 0.0`` and
        ``rout >= rin >= 0.0``

    err, var : float or `~numpy.ndarray`
        Error *or* variance (specify at most one).

    mask : `~numpy.ndarray`, optional
        Mask array. If supplied, a given pixel is masked if its value
        is greater than ``maskthresh``.

    maskthresh : float, optional
        Threshold for a pixel to be masked. Default is ``0.0``.

    gain : float, optional
        Conversion factor between data array units and poisson counts,
        used in calculating poisson noise in aperture sum. If ``None``
        (default), do not add poisson noise.

    subpix : int, optional
        Subpixel sampling factor. Default is 5.

    Returns
    -------
    sum : `~numpy.ndarray`
        The sum of the data array within the aperture(s).

    sumerr : `~numpy.ndarray`
        Error on the sum.

    flags : `~numpy.ndarray`
        Integer giving flags. (0 if no flags set.)
    """

    cdef double flux1, fluxerr1, x1, y1, r1, area1, rin1, rout1
    cdef double bkgflux, bkgfluxerr, bkgarea
    cdef short flag1, bkgflag
    cdef size_t i
    cdef int status
    cdef np.broadcast it
    cdef sep_image im

    _parse_arrays(data, err, var, mask, &im)
    im.maskthresh = maskthresh
    if gain is not None:
        im.gain = gain

    # Require that inputs are float64 arrays. See note in circular aperture.
    dt = np.dtype(np.double)
    x = np.require(x, dtype=dt)
    y = np.require(y, dtype=dt)
    a = np.require(a, dtype=dt)
    b = np.require(b, dtype=dt)
    theta = np.require(theta, dtype=dt)
    rin = np.require(rin, dtype=dt)
    rout = np.require(rout, dtype=dt)

    # allocate ouput arrays
    shape = np.broadcast(x, y, a, b, theta, rin, rout).shape
    sum = np.empty(shape, dt)
    sumerr = np.empty(shape, dt)
    flag = np.empty(shape, np.short)

    it = np.broadcast(x, y, a, b, theta, rin, rout, sum, sumerr, flag)
    while np.PyArray_MultiIter_NOTDONE(it):
        status = sep_sum_ellipann(
            &im,
            (<double*>np.PyArray_MultiIter_DATA(it, 0))[0],
            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0],
            (<double*>np.PyArray_MultiIter_DATA(it, 2))[0],
            (<double*>np.PyArray_MultiIter_DATA(it, 3))[0],
            (<double*>np.PyArray_MultiIter_DATA(it, 4))[0],
            (<double*>np.PyArray_MultiIter_DATA(it, 5))[0],
            (<double*>np.PyArray_MultiIter_DATA(it, 6))[0],
            subpix, 0,
            <double*>np.PyArray_MultiIter_DATA(it, 7),
            <double*>np.PyArray_MultiIter_DATA(it, 8),
            &area1,
            <short*>np.PyArray_MultiIter_DATA(it, 9))
        _assert_ok(status)
        np.PyArray_MultiIter_NEXT(it)

    return sum, sumerr, flag

@cython.boundscheck(False)
@cython.wraparound(False)
def flux_radius(np.ndarray data not None, x, y, rmax, frac, normflux=None,
                np.ndarray mask=None, double maskthresh=0.0, int subpix=5):
    """flux_radius(data, x, y, rmax, frac, normflux=None, mask=None,
                   maskthresh=0.0, subpix=5)

    Return radius of a circle enclosing requested fraction of total flux.

    Parameters
    ----------
    data : `~numpy.ndarray`
        2-d array to be summed.

    x, y : array_like
        Center coordinates and radius (radii) of aperture(s). 
        ``x`` corresponds to the second ("fast") axis of the input array
        and ``y`` corresponds to the first ("slow") axis.
        ``x, y = (0.0, 0.0)`` corresponds to the center of the first
        element of the array. Shapes must match.

    rmax : array_like
        Maximum radius to analyze. Used as normalizing flux if ``normflux``
        is None. Shape must match x and y.

    frac : array_like
        Requested fraction of light (in range 0 to 1). Can be scalar or array.

    normflux : array_like, optional
        Normalizing flux for each position. If not given, the sum
        within ``rmax`` is used as the normalizing flux. If given,
        shape must match x, y and rmax.

    mask : `~numpy.ndarray`, optional
        Mask array. If supplied, a given pixel is masked if its value
        is greater than ``maskthresh``.

    maskthresh : float, optional
        Threshold for a pixel to be masked. Default is ``0.0``.

    subpix : int, optional
        Subpixel sampling factor. Default is 5.

    Returns
    -------
    radius : `~numpy.ndarray`
        The sum of the data array within the aperture(s). Shape is
        same as ``x``, except if ``frac`` is an array; then the
        dimension of ``frac`` will be appended. For example, if ``x``
        and ``frac`` are both 1-d arrays, the result will be a 2-d
        array with the trailing dimension corresponding to ``frac``.

    flags : `~numpy.ndarray`
        Integer giving flags. Same shape as ``x``. (0 if no flags set.)

    """

    cdef double flux1, fluxerr1, x1, y1, r1, area1, rin1, rout1
    cdef double bkgflux, bkgfluxerr, bkgarea
    cdef short flag1, bkgflag
    cdef size_t i
    cdef int status, fracn
    cdef short[:] flag
    cdef double[:, :] radius
    cdef double[:] fractmp
    cdef double[:] xtmp
    cdef double[:] ytmp
    cdef double[:] rtmp
    cdef double[:] normfluxbuf
    cdef double *normfluxptr
    cdef sep_image im

    _parse_arrays(data, None, None, mask, &im)
    im.maskthresh = maskthresh

    # Require that inputs are float64 arrays with same shape. See note in
    # circular aperture.
    # Also require that frac is a contiguous array.
    dt = np.dtype(np.double)
    x = np.require(x, dtype=dt)
    y = np.require(y, dtype=dt)
    rmax = np.require(rmax, dtype=dt)
    frac = np.require(frac, dtype=dt)
    inshape = x.shape
    infracshape = frac.shape
    if (y.shape != inshape or rmax.shape != inshape):
        raise ValueError("shape of x, y, and r must match")

    # Convert input arrays to 1-d for correct looping and indexing.
    xtmp = np.ravel(x)
    ytmp = np.ravel(y)
    rtmp = np.ravel(rmax)
    fractmp = np.ravel(np.ascontiguousarray(frac))
    fracn = len(fractmp)

    # optional `normflux` input.
    normfluxptr = NULL
    if normflux is not None:
        normflux = np.require(normflux, dtype=dt)
        if normflux.shape != inshape:
            raise ValueError("shape of normflux must match shape of "
                             "x, y and r")
        normfluxbuf = np.ravel(normflux)
        normfluxptr = &normfluxbuf[0]

    # Allocate ouput arrays. (We'll reshape these later to match the
    # input shapes.)
    flag = np.empty(len(xtmp), np.short)
    radius = np.empty((len(xtmp), len(fractmp)), dt)

    for i in range(len(xtmp)):
        if normfluxptr != NULL:
            normfluxptr = &normfluxbuf[i]
        status = sep_flux_radius(&im,
                                 xtmp[i], ytmp[i], rtmp[i], subpix, 0,
                                 normfluxptr, &fractmp[0], fracn,
                                 &radius[i, 0], &flag[i])
        _assert_ok(status)

    return (np.asarray(radius).reshape(inshape + infracshape),
            np.asarray(flag).reshape(inshape))


@cython.boundscheck(False)
@cython.wraparound(False)
def mask_ellipse(np.ndarray arr not None, x, y, a=None, b=None, theta=None,
                 r=1.0, **kwargs):
    """mask_ellipse(arr, x, y, a, b, theta, r=1.0)

    Mask ellipse(s) in an array.

    Set array elements to True (or 1) if they fall within the given
    ellipse.  The ``r`` keyword can be used to scale the ellipse.
    Equivalently, after converting ``a``, ``b``, ``theta`` to a
    coefficient ellipse representation (``cxx``, ``cyy``, ``cxy``),
    pixels that fulfill the condition

    .. math::

       cxx(x_i - x)^2 + cyy(y_i - y)^2 + cxx(x_i - x)(y_i - y) < r^2

    will be masked.

    Parameters
    ----------
    arr : `~numpy.ndarray`
        Input array to be masked. Array is updated in-place.
    x, y : array_like
        Center of ellipse(s).
    a, b, theta : array_like
        Parameters defining the extent of the ellipe(s).
    r : array_like, optional
        Scale factor of ellipse(s). Default is 1.
    """

    cdef int w, h
    cdef np.uint8_t[:,:] buf
    cdef double cxx, cyy, cxy

    dt = np.dtype(np.double)

    # only boolean arrays supported
    if not (arr.dtype.type is np.bool_ or arr.dtype.type is np.ubyte):
        raise ValueError("Array data type not supported: {0:s}"
                         .format(arr.dtype))
    _check_array_get_dims(arr, &w, &h)
    buf = arr.view(dtype=np.uint8)

    x = np.require(x, dtype=dt)
    y = np.require(y, dtype=dt)

    # New Behavior:
    if (a is not None and b is not None and theta is not None):
        r = np.require(r, dtype=dt)
        a = np.require(a, dtype=dt)
        b = np.require(b, dtype=dt)
        theta = np.require(theta, dtype=dt)

        it = np.broadcast(x, y, a, b, theta, r)
        while np.PyArray_MultiIter_NOTDONE(it):
            sep_ellipse_coeffs((<double*>np.PyArray_MultiIter_DATA(it, 2))[0],
                               (<double*>np.PyArray_MultiIter_DATA(it, 3))[0],
                               (<double*>np.PyArray_MultiIter_DATA(it, 4))[0],
                               &cxx, &cyy, &cxy)
            sep_set_ellipse(<unsigned char *>&buf[0, 0], w, h,
                            (<double*>np.PyArray_MultiIter_DATA(it, 0))[0],
                            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0],
                            cxx, cyy, cxy,
                            (<double*>np.PyArray_MultiIter_DATA(it, 5))[0],
                            1)
            np.PyArray_MultiIter_NEXT(it)

    # deprecated behavior
    elif ("cxx" in kwargs and "cyy" in kwargs and "cxy" in kwargs):
        if "scale" in kwargs:
            r = kwargs["scale"]
        r = np.require(r, dtype=dt)
        cxx_ = np.require(kwargs["cxx"], dtype=dt)
        cyy_ = np.require(kwargs["cyy"], dtype=dt)
        cxy_ = np.require(kwargs["cxy"], dtype=dt)

        it = np.broadcast(x, y, cxx_, cyy_, cxy_, r)
        while np.PyArray_MultiIter_NOTDONE(it):
            sep_set_ellipse(<unsigned char *>&buf[0, 0], w, h,
                            (<double*>np.PyArray_MultiIter_DATA(it, 0))[0],
                            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0],
                            (<double*>np.PyArray_MultiIter_DATA(it, 2))[0],
                            (<double*>np.PyArray_MultiIter_DATA(it, 3))[0],
                            (<double*>np.PyArray_MultiIter_DATA(it, 4))[0],
                            (<double*>np.PyArray_MultiIter_DATA(it, 5))[0],
                            1)
            np.PyArray_MultiIter_NEXT(it)
    else:
        raise ValueError("Must specify a, b and theta")


def kron_radius(np.ndarray data not None, x, y, a, b, theta, r,
                np.ndarray mask=None, double maskthresh=0.0):
    """kron_radius(data, x, y, a, b, theta, r, mask=None, maskthresh=0.0)

    Calculate Kron "radius" within an ellipse.

    The Kron radius is given by 

    .. math::

       \sum_i r_i I(r_i) / \sum_i I(r_i)

    where the sum is over all pixels in the aperture and the radius is given
    in units of ``a`` and ``b``: ``r_i`` is the distance to the pixel relative
    to the distance to the ellipse specified by ``a``, ``b``, ``theta``.
    Equivalently, after converting the ellipse parameters to their coefficient
    representation, ``r_i`` is given by

    .. math::

       r_i^2 = cxx(x_i - x)^2 + cyy(y_i - y)^2 + cxx(x_i - x)(y_i - y)

    Parameters
    ----------
    data : `~numpy.ndarray`
        Data array.

    x, y : array_like
        Ellipse center(s).

    a, b, theta : array_like
        Ellipse parameters.

    r : array_like
        "Radius" of ellipse over which to integrate. If the ellipse
        extent correponds to second moments of an object, this is the
        number of "isophotal radii" in Source Extractor parlance. A
        Fixed value of 6 is used in Source Extractor.

    mask : `numpy.ndarray`, optional
        An optional mask.

    maskthresh : float, optional
        Pixels with mask > maskthresh will be ignored.

    Returns
    -------
    kronrad : array_like
        The Kron radius.

    flag : array_like
        Integer value indicating conditions about the aperture or how
        many masked pixels it contains.

    """

    cdef double cxx, cyy, cxy
    cdef sep_image im

    _parse_arrays(data, None, None, mask, &im)
    im.maskthresh = maskthresh

    # See note in apercirc on requiring specific array type
    dt = np.dtype(np.double)
    x = np.require(x, dtype=dt)
    y = np.require(y, dtype=dt)
    a = np.require(a, dtype=dt)
    b = np.require(b, dtype=dt)
    theta = np.require(theta, dtype=dt)
    r = np.require(r, dtype=dt)

    # allocate output arrays
    shape = np.broadcast(x, y, a, b, theta, r).shape
    kr = np.empty(shape, np.float)
    flag = np.empty(shape, np.short)

    it = np.broadcast(x, y, a, b, theta, r, kr, flag)
    while np.PyArray_MultiIter_NOTDONE(it):
        sep_ellipse_coeffs((<double*>np.PyArray_MultiIter_DATA(it, 2))[0],
                           (<double*>np.PyArray_MultiIter_DATA(it, 3))[0],
                           (<double*>np.PyArray_MultiIter_DATA(it, 4))[0],
                           &cxx, &cyy, &cxy)
        sep_kron_radius(&im,
                        (<double*>np.PyArray_MultiIter_DATA(it, 0))[0],
                        (<double*>np.PyArray_MultiIter_DATA(it, 1))[0],
                        cxx, cyy, cxy,
                        (<double*>np.PyArray_MultiIter_DATA(it, 5))[0],
                        <double*>np.PyArray_MultiIter_DATA(it, 6),
                        <short*>np.PyArray_MultiIter_DATA(it, 7))
        np.PyArray_MultiIter_NEXT(it)

    return kr, flag 

def winpos(np.ndarray data not None, xinit, yinit, sig,
           np.ndarray mask=None, double maskthresh=0.0, int subpix=11):
    """winpos(data, xinit, yinit, sig, mask=None, maskthresh=0.0, subpix=11)

    Calculate more accurate object centroids using 'windowed' algorithm.

    Starting from the supplied initial center position, an iterative
    algorithm is used to determine a better object centroid. On each
    iteration, the centroid is calculated using all pixels within a
    circular aperture of ``4*sig`` from the current position,
    weighting pixel positions by their flux and the amplitude of a 2-d
    Gaussian with sigma ``sig``. Iteration stops when the change in
    position falls under some threshold or a maximum number of
    iterations is reached. This is equivalent to ``XWIN_IMAGE`` and
    ``YWIN_IMAGE`` parameters in Source Extractor (for the correct choice
    of sigma for each object).

    Parameters
    ----------
    data : `~numpy.ndarray`
        Data array.

    xinit, yinit : array_like
        Initial center(s).

    sig : array_like
        Gaussian sigma used for weighting pixels. Pixels within a circular
        aperture of radius 4*sig are included.

    mask : `numpy.ndarray`, optional
        An optional mask.

    maskthresh : float, optional
        Pixels with mask > maskthresh will be ignored.

    subpix : int
        Subpixel sampling used to determine pixel overlap with
        aperture.  11 is used in Source Extractor. For exact overlap
        calculation, use 0.

    Returns
    -------
    x, y : np.ndarray
        New x and y position(s).

    flag : np.ndarray
        Flags.

    """

    cdef int status
    cdef double cxx, cyy, cxy
    cdef int niter = 0  # not currently returned
    cdef sep_image im

    _parse_arrays(data, None, None, mask, &im)
    im.maskthresh = maskthresh

    # See note in apercirc on requiring specific array type
    dt = np.dtype(np.double)
    xinit = np.require(xinit, dtype=dt)
    yinit = np.require(yinit, dtype=dt)
    sig = np.require(sig, dtype=dt)

    # allocate output arrays
    shape = np.broadcast(xinit, yinit, sig).shape
    x = np.empty(shape, np.float)
    y = np.empty(shape, np.float)
    flag = np.empty(shape, np.short)

    it = np.broadcast(xinit, yinit, sig, x, y, flag)
    while np.PyArray_MultiIter_NOTDONE(it):
        status = sep_windowed(&im,
                              (<double*>np.PyArray_MultiIter_DATA(it, 0))[0],
                              (<double*>np.PyArray_MultiIter_DATA(it, 1))[0],
                              (<double*>np.PyArray_MultiIter_DATA(it, 2))[0],
                              subpix, 0,
                              <double*>np.PyArray_MultiIter_DATA(it, 3),
                              <double*>np.PyArray_MultiIter_DATA(it, 4),
                              &niter,
                              <short*>np.PyArray_MultiIter_DATA(it, 5))
        _assert_ok(status)
        np.PyArray_MultiIter_NEXT(it)

    return x, y, flag 



def ellipse_coeffs(a, b, theta):
    """ellipse_coeffs(a, b, theta)

    Convert from ellipse axes and angle to coefficient representation.

    Parameters
    ----------
    a, b, theta : array_like
        Ellipse(s) semi-major, semi-minor axes and position angle
        respectively.  Position angle is radians counter clockwise
        from positive x axis to major axis, and lies in range
        ``[-pi/2, pi/2]``

    Returns
    -------
    cxx, cyy, cxy : `~numpy.ndarray`
        Describes the ellipse(s) ``cxx * x^2 + cyy * y^2 + cxy * xy = 1``
    """

    dt = np.dtype(np.double)
    a = np.require(a, dtype=dt)
    b = np.require(b, dtype=dt)
    theta = np.require(theta, dtype=dt)

    shape = np.broadcast(a, b, theta).shape
    cxx = np.empty(shape, dt)
    cyy = np.empty(shape, dt)
    cxy = np.empty(shape, dt)

    it = np.broadcast(a, b, theta, cxx, cyy, cxy)
    while np.PyArray_MultiIter_NOTDONE(it):
        sep_ellipse_coeffs((<double*>np.PyArray_MultiIter_DATA(it, 0))[0],
                           (<double*>np.PyArray_MultiIter_DATA(it, 1))[0],
                           (<double*>np.PyArray_MultiIter_DATA(it, 2))[0],
                           <double*>np.PyArray_MultiIter_DATA(it, 3),
                           <double*>np.PyArray_MultiIter_DATA(it, 4),
                           <double*>np.PyArray_MultiIter_DATA(it, 5))
        np.PyArray_MultiIter_NEXT(it)

    return cxx, cyy, cxy

def ellipse_axes(cxx, cyy, cxy):
    """ellipse_axes(cxx, cyy, cxy)

    Convert from coefficient ellipse representation to ellipse axes and angle.

    Parameters
    ----------
    cxx, cyy, cxy : array_like
        Describes the ellipse(s) ``cxx * x**2 + cyy * y**2 + cxy * x * y = 1``

    Returns
    -------
    a, b, theta : `~numpy.ndarray`
        Ellipse(s) semi-major, semi-minor axes and position angle
        respectively.  Position angle is radians counter clockwise
        from positive x axis to major axis, and lies in range
        ``(-pi/2, pi/2)``

    Raises
    ------
    ValueError
        If input parameters do not describe an ellipse.

    """

    cdef int status

    dt = np.dtype(np.double)
    cxx = np.require(cxx, dtype=dt)
    cyy = np.require(cyy, dtype=dt)
    cxy = np.require(cxy, dtype=dt)

    shape = np.broadcast(cxx, cyy, cxy).shape
    a = np.empty(shape, dt)
    b = np.empty(shape, dt)
    theta = np.empty(shape, dt)

    status = 0
    it = np.broadcast(cxx, cyy, cxy, a, b, theta)
    while np.PyArray_MultiIter_NOTDONE(it):
        status = sep_ellipse_axes(
            (<double*>np.PyArray_MultiIter_DATA(it, 0))[0],
            (<double*>np.PyArray_MultiIter_DATA(it, 1))[0],
            (<double*>np.PyArray_MultiIter_DATA(it, 2))[0],
            <double*>np.PyArray_MultiIter_DATA(it, 3),
            <double*>np.PyArray_MultiIter_DATA(it, 4),
            <double*>np.PyArray_MultiIter_DATA(it, 5))
        if status:
            break

        np.PyArray_MultiIter_NEXT(it)

    if status:
        raise ValueError(
            "parameters do not describe ellipse: "
            "cxx={0:f}, cyy={1:f}, cxy={2:f}".format(
                (<double*>np.PyArray_MultiIter_DATA(it, 0))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 1))[0],
                (<double*>np.PyArray_MultiIter_DATA(it, 2))[0]))

    return a, b, theta

# -----------------------------------------------------------------------------
# Utility functions

def set_extract_pixstack(size_t size):
    """set_extract_pixstack(size)

    Set the size in pixels of the internal pixel buffer used in extract().

    The current value can be retrieved with get_extract_pixstack. The
    initial default is 300000.
    """
    sep_set_extract_pixstack(size)

def get_extract_pixstack():
    """get_extract_pixstack()

    Get the size in pixels of the internal pixel buffer used in extract().
    """
    return sep_get_extract_pixstack()

# -----------------------------------------------------------------------------
# deprecated stuff

def istruncated(np.ndarray flag not None):
    """True where 'aperture truncated' flag is set."""
    return (flag & APER_TRUNC) != 0

def hasmasked(np.ndarray flag not None):
    """True where 'aperture has masked pixel(s)' flag is set."""
    return (flag & APER_HASMASKED) != 0

apercirc = sum_circle
