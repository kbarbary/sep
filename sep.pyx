# This wrapper licensed under an MIT license.

"""
Source Extraction and Photometry

This module is a wrapper of the SEP C library.
"""
import numpy as np
cimport numpy as np
from libc cimport limits
cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Free

np.import_array()  # To access the numpy C-API.

__version__ = "0.2.dev"

# -----------------------------------------------------------------------------
# Definitions from the SEP C library

# macro definitions from sep.h
DEF SEP_TBYTE = 11
DEF SEP_TINT = 31
DEF SEP_TFLOAT = 42
DEF SEP_TDOUBLE = 82

# Flag values (both input and output)
DEF SEP_ERROR_IS_VAR = 0x0001
DEF SEP_ERROR_IS_ARRAY = 0x0002
DEF SEP_MASK_IGNORE = 0x0004
DEF SEP_APER_TRUNC = 0x0010
DEF SEP_APER_HASMASKED = 0x0200

# macro defintion from sepcore.h
# This is not part of the SEP API, but we pull it out because want to
# explicitly detect memory errors so that we can raise MemoryError().
DEF MEMORY_ALLOC_ERROR = 1

# header definitions
cdef extern from "sep.h":

    ctypedef struct sepbackmap:
        int w
        int h
        float globalback
        float globalrms
    
    ctypedef struct sepobj:
        float  thresh
        int    npix
        int    tnpix
        int    xmin,xmax,ymin,ymax
        double x, y
        double x2, y2, xy
        float  a, b, theta
        float  cxx, cyy, cxy
        float  cflux
        float  flux
        float  cpeak
        float  peak
        int    xpeak, ypeak
        int    xcpeak, ycpeak
        short  flag
        int    *pix

    int sep_makeback(void* im, void *mask,
                     int dtype,
                     int mdtype,
                     int w, int h,
                     int bw, int bh,
                     float mthresh,
                     int fw, int fh,
                     float fthresh,
                     sepbackmap **bkmap)

    int sep_backarray(sepbackmap *bkmap, void *arr, int dtype)
    int sep_backrmsarray(sepbackmap *bkmap, void *arr, int dtype)
    int sep_subbackarray(sepbackmap *bkmap, void *arr, int dtype)
    void sep_freeback(sepbackmap *bkmap)

    int sep_extract(void *image,
                    void *noise,
                    int dtype,
                    int ndtype,
                    short noise_flag,
                    int w, int h,
                    float thresh,
                    int minarea,
                    float *conv,
                    int convw, int convh,
                    int deblend_nthresh,
                    double deblend_cont,
                    int clean_flag,
                    double clean_param,
                    sepobj **objects,
                    int *nobj)

    void sep_freeobjarray(sepobj *objects, int nobj)

    int sep_apercirc(void *data, void *error, void *mask,
                     int dtype, int edtype, int mdtype, int w, int h,
                     double maskthresh, double gain, short inflags,
                     double x, double y, double r, int subpix,
                     double *sum, double *sumerr, double *area, short *flag)

    int sep_apercircann(void *data, void *error, void *mask,
                        int dtype, int edtype, int mdtype, int w, int h,
                        double maskthresh, double gain, short inflags,
                        double x, double y, double rin, double rout, int subpix,
                        double *sum, double *sumerr, double *area, short *flag)

    void sep_setellipse_ucc(unsigned char *arr, int w, int h,
                            float x, float y, float cxx, float cyy, float cxy,
                            float r, unsigned char val)

    void sep_get_errmsg(int status, char *errtext)
    void sep_get_errdetail(char *errtext)

# -----------------------------------------------------------------------------
# Utility functions

cdef int _get_sep_dtype(dtype) except -1:
    """Convert a numpy dtype to the corresponding SEP dtype integer code."""
    if not dtype.isnative:
        raise ValueError(
            "Input array with dtype {0} has non-native byte order. "
            "Only native byte order arrays are supported. "
            "Arrays can be converted to native byte order in-place "
            "with 'a = a.byteswap(True).newbyteorder()'".format(dtype))
    t = dtype.type
    if t is np.single:
        return SEP_TFLOAT
    elif t is np.bool_ or t is np.ubyte:
        return SEP_TBYTE
    elif t == np.double:
        return SEP_TDOUBLE
    elif t == np.intc:
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
    cdef bytes pyerrmsg
    cdef bytes pyerrdetail

    if status == 0:
        return 0
    if status == MEMORY_ALLOC_ERROR:
        raise MemoryError

    # otherwise, get error message
    errmsg = <char *>PyMem_Malloc(61 * sizeof(char))
    sep_get_errmsg(status, errmsg)
    pyerrmsg = errmsg
    PyMem_Free(errmsg)

    errdetail = <char *>PyMem_Malloc(512 * sizeof(char))
    sep_get_errdetail(errdetail)
    pyerrdetail = errdetail
    PyMem_Free(errdetail)

    raise Exception(pyerrmsg + " " + pyerrdetail)    

# -----------------------------------------------------------------------------
# Background Estimation

cdef class Background:
    """Background(data, mask=None, maskthresh=0.0, bw=64, bh=64, fw=3, fh=3,
                  fthresh=0.0)

    Estimate the global image background and noise.

    Parameters
    ----------
    data : 2-d numpy.ndarray
        Data array.
    mask : 2-d numpy.ndarray, optional
        Mask array, optional
    maskthresh : float, optional
        Mask threshold. This is the inclusive upper limit on the mask value
        in order for the corresponding pixel to be unmasked. For boolean
        arrays, False and True are interpreted as 0 and 1, respectively.
        Thus, given a threshold of zero, True corresponds to masked and
        False corresponds to unmasked.
    bw, bh : int, optional
        Size of background boxes in pixels. Default is 3.
    fw, fh : int, optional
        Filter width and height in boxes. Default is 3.
    fthresh : float, optional
        Filter threshold.
    """

    cdef sepbackmap *ptr      # pointer to C struct
    cdef np.dtype orig_dtype  # dtype code of original image

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __cinit__(self, np.ndarray data not None, np.ndarray mask=None,
                  float maskthresh=0.0, int bw=64, int bh=64,
                  int fw=3, int fh=3, float fthresh=0.0):

        cdef int w, h, w2, h2, status, sep_dtype, sep_dtype2 
        cdef np.uint8_t[:, :] buf
        cdef np.uint8_t[:, :] buf2
        cdef void *maskptr

        _check_array_get_dims(data, &w, &h)
        sep_dtype = _get_sep_dtype(data.dtype)
        self.orig_dtype = data.dtype
        buf = data.view(dtype=np.uint8)

        if mask is not None:
            _check_array_get_dims(mask, &w2, &h2)
            if (w != w2) or (h != h2):
                raise ValueError("dimensions of data and mask must match")
            sep_dtype2 = _get_sep_dtype(mask.dtype)
            buf2 = mask.view(dtype=np.uint8)
            maskptr = &buf2[0, 0]
        else:
            sep_dtype2 = 0  # ignored if mask is NULL
            maskptr = NULL

        status = sep_makeback(&buf[0, 0], maskptr, sep_dtype, sep_dtype2,
                              w, h, bw, bh, maskthresh, fw, fh, fthresh,
                              &self.ptr)
        _assert_ok(status)

    property globalback:
        """Average background level."""
        def __get__(self):
            return self.ptr.globalback

    property globalrms:
        """Average background RMS"""
        def __get__(self):
            return self.ptr.globalrms

    def back(self, dtype=None):
        """back(dtype=None)

        Create an array of the background.
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
        status = sep_backarray(self.ptr, &buf[0, 0], sep_dtype)
        _assert_ok(status)

        return result

    def rms(self, dtype=None):
        """rms(dtype=None)

        Create an array of the background rms.
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
        status = sep_backrmsarray(self.ptr, &buf[0, 0], sep_dtype)
        _assert_ok(status)

        return result


    def subfrom(self, np.ndarray data not None):
        """subfrom(data)

        Subtract the background from an existing array.

        Parameters
        ----------
        data : 2-d numpy.ndarray
            Shape must match that of the original image used to measure the
            background.
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

        status = sep_subbackarray(self.ptr, &buf[0, 0], sep_dtype)
        _assert_ok(status)

    def __dealloc__(self):
        if self.ptr is not NULL:
            sep_freeback(self.ptr)

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

default_conv = np.array([[1.0, 2.0, 1.0],
                         [2.0, 4.0, 2.0],
                         [1.0, 2.0, 1.0]], dtype=np.float32)

def extract(np.ndarray data not None, float thresh, int minarea=5,
            np.ndarray conv=default_conv, int deblend_nthresh=32,
            double deblend_cont=0.005, bint clean=True,
            double clean_param=1.0):
    """extract(data, thresh, minarea=5, conv=default_conv, deblend_nthresh=32,
               deblend_cont=0.005, clean=True, clean_param=1.0)

    Extract sources from an image.

    Parameters
    ----------
    data : np.ndarray
        Data array (2-d).
    thresh : float
        Threshold pixel value for detection.
    minarea : int, optional
        Minimum number of pixels required for an object. Default is 5.
    conv : np.ndarray or None
        Convolution kernel used for on-the-fly image convolution (used to
        enhance detection). Default is a 3x3 array. Set to ``None`` to skip
        convolution.
    deblend_nthresh : int, optional
        Number of thresholds used for object deblending. Default is 32.
    deblend_cont : float, optional
        Minimum contrast ratio used for object deblending. Default is 0.005.
    clean : bool, optional
        Perform cleaning? Default is True.
    clean_param : float, optional
        Cleaning parameter (see SExtractor manual). Default is 1.0.

    Returns
    -------
    objects : np.ndarray
        Extracted object parameters (structured array).
        See `objects.dtype.names` for available fields.
    """

    cdef int w, h, convw, convh, status, sep_dtype, nobj, i
    cdef np.uint8_t[:, :] buf
    cdef sepobj *objects
    cdef np.ndarray[Object] result
    cdef float[:, :] convflt
    #cdef float *default_conv = [1., 2., 1., 2., 4., 2., 1., 2., 1.]
    #cdef int default_convw = 3
    #cdef int default_convh = 3
    cdef float *convptr

    _check_array_get_dims(data, &w, &h)
    sep_dtype = _get_sep_dtype(data.dtype)
    buf = data.view(dtype=np.uint8)

    # Parse convolution input
    #if conv == 'default':
    #    convptr = default_conv
    #    convw = default_convw
    #    convh = default_convh
    if conv is None:
        convptr = NULL
        convw = 0
        convh = 0
    else:
        convflt = conv.astype(np.float32)
        convptr = &convflt[0, 0]
        convw = convflt.shape[1]
        convh = convflt.shape[0]

    status = sep_extract(&buf[0,0], NULL, sep_dtype, 0, 0, w, h,
                         thresh, minarea, convptr, convw, convh,
                         deblend_nthresh, deblend_cont, clean, clean_param,
                         &objects, &nobj)
    _assert_ok(status)

    # Allocate result record array and fill it
    result = np.empty(nobj, dtype=np.dtype([('thresh', np.float64),
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

    for i in range(nobj):
        result['thresh'][i] = objects[i].thresh
        result['npix'][i] = objects[i].npix
        result['tnpix'][i] = objects[i].tnpix
        result['xmin'][i] = objects[i].xmin
        result['xmax'][i] = objects[i].xmax
        result['ymin'][i] = objects[i].ymin
        result['ymax'][i] = objects[i].ymax
        result['x'][i] = objects[i].x
        result['y'][i] = objects[i].y
        result['x2'][i] = objects[i].x2
        result['y2'][i] = objects[i].y2
        result['xy'][i] = objects[i].xy
        result['a'][i] = objects[i].a
        result['b'][i] = objects[i].b
        result['theta'][i] = objects[i].theta
        result['cxx'][i] = objects[i].cxx
        result['cyy'][i] = objects[i].cyy
        result['cxy'][i] = objects[i].cxy
        result['cflux'][i] = objects[i].cflux
        result['flux'][i] = objects[i].flux
        result['cpeak'][i] = objects[i].cpeak
        result['peak'][i] = objects[i].peak
        result['xcpeak'][i] = objects[i].xcpeak
        result['ycpeak'][i] = objects[i].ycpeak
        result['xpeak'][i] = objects[i].xpeak
        result['ypeak'][i] = objects[i].ypeak
        result['flag'][i] = objects[i].flag

    # Free C array
    sep_freeobjarray(objects, nobj)

    return result

# -----------------------------------------------------------------------------
# Aperture Photometry

@cython.boundscheck(False)
@cython.wraparound(False)
def apercirc(np.ndarray data not None, x, y, r,
             var=None, err=None, gain=None, np.ndarray mask=None,
             double maskthresh=0.0, bkgann=None, int subpix=5):
    """apercirc(data, x, y, r, err=None, var=None, mask=None, maskthresh=0.0,
                bkgann=None, gain=None, subpix=5)

    Sum data in circular apertures.

    Parameters
    ----------
    data : np.ndarray
        2-d array to be summed.
    x, y, r : array_like
        Center coordinates and radius (radii) of aperture(s). 
        ``x`` corresponds to the second ("fast") axis of the input array
        and ``y`` corresponds to the first ("slow") axis.
        ``x, y = (0.0, 0.0)`` corresponds to the center of the first
        element of the array. These inputs obey numpy broadcasting rules.
    err, var : float or ndarray
        Error *or* variance (specify at most one).
    mask : ndarray, optional
        Mask array. If supplied, a given pixel is masked if its value
        is greater than ``maskthresh``.
    maskthresh : float, optional
        Threshold for a pixel to be masked. Default is ``0.0``.
    bkgann : tuple, optional
        Length 2 tuple giving the inner and outer radius of a
        "background annulus". If supplied, the background is estimated
        by averaging unmasked pixels in this annulus. If supplied, the inner
        and outer radii obey numpy broadcasting rules, along with ``x``,
        ``y`` and ``r``.
    gain : float, optional
        Conversion factor between data array units and poisson counts,
        used in calculating poisson noise in aperture sum. If ``None``
        (default), do not add poisson noise.
    subpix : int, optional
        Subpixel sampling factor. Default is 5.

    Returns
    -------
    sum : ndarray
        The sum of the data array within the aperture.
    sumerr : ndarray
        Error on the sum.
    flags : ndarray
        Integer giving flags. (0 if no flags set.)
    """

    cdef double flux1, fluxerr1, x1, y1, r1, area1, rin1, rout1
    cdef double bkgflux, bkgfluxerr, bkgarea, gain_
    cdef float scalarerr
    cdef short flag1, bkgflag
    cdef short inflag
    cdef size_t i
    cdef int w, h, ew, eh, mw, mh, dtype, edtype, mdtype, status
    cdef np.uint8_t[:,:] buf, ebuf, mbuf
    cdef void *ptr
    cdef void *eptr
    cdef void *mptr
    cdef np.broadcast it

    inflag = 0
    eptr = NULL
    mptr = NULL
    edtype = 0
    mdtype = 0

    # Get main image info
    _check_array_get_dims(data, &w, &h)
    dtype = _get_sep_dtype(data.dtype)
    buf = data.view(dtype=np.uint8)
    ptr = <void*>&buf[0, 0]

    # Optional inputs: error or variance
    if err is not None and var is not None:
        raise ValueError("Cannot specify both err and var")
    elif var is not None:
        err = var
        inflag |= SEP_ERROR_IS_VAR
    if err is not None:
        if isinstance(err, np.ndarray):
            if err.ndim == 0:
                scalarerr = err
                eptr = <void*>&scalarerr
                edtype = SEP_TFLOAT
            elif err.ndim == 2:
                _check_array_get_dims(err, &ew, &eh)
                if ew != w or eh != h:
                    raise ValueError("size of error/variance array must match"
                                     " data")
                edtype = _get_sep_dtype(err.dtype)
                ebuf = err.view(dtype=np.uint8)
                eptr = <void*>&ebuf[0, 0]
                inflag |= SEP_ERROR_IS_ARRAY
            else:
                raise ValueError("error/variance array must be 0-d or 2-d")
        else:
            scalarerr = err
            eptr = <void*>&scalarerr
            edtype = SEP_TFLOAT

    # Optional input: mask
    if mask is not None:
        _check_array_get_dims(mask, &mw, &mh)
        if mw != w or mh != h:
            raise ValueError("size of mask array must match data")
        mdtype = _get_sep_dtype(mask.dtype)
        mbuf = mask.view(dtype=np.uint8)
        mptr = <void*>&mbuf[0, 0]

    # Optional input: gain
    gain_ = 0.0
    if gain is not None:
        gain_ = gain

    # Require that inputs are float64 arrays. This has to be done because we
    # are using a broadcasting iterator below, where we need to know the type
    # in advance. There are other ways to do this, e.g., using NpyIter_Multi
    # in the numpy C-API. However, the best way to use this from cython
    # is not clear.
    #
    # docs.scipy.org/doc/numpy/reference/c-api.iterator.html#NpyIter_MultiNew
    # 
    x = np.require(x, dtype=np.float)
    y = np.require(y, dtype=np.float)
    r = np.require(r, dtype=np.float)

    if bkgann is None:

        # allocate ouput arrays
        shape = np.broadcast(x, y, r).shape
        sum = np.empty(shape, np.float)
        sumerr = np.empty(shape, np.float)
        flag = np.empty(shape, np.short)

        # broadcasting iterator over x, y, r
        it = np.broadcast(x, y, r, sum, sumerr, flag)

        while np.PyArray_MultiIter_NOTDONE(it):

                x1 = (<double*>np.PyArray_MultiIter_DATA(it, 0))[0]
                y1 = (<double*>np.PyArray_MultiIter_DATA(it, 1))[0]
                r1 = (<double*>np.PyArray_MultiIter_DATA(it, 2))[0]

                status = sep_apercirc(ptr, eptr, mptr,
                                      dtype, edtype, mdtype, w, h,
                                      maskthresh, gain_, inflag,
                                      x1, y1, r1, subpix,
                                      &flux1, &fluxerr1, &area1, &flag1)

                _assert_ok(status)

                (<double*>np.PyArray_MultiIter_DATA(it, 3))[0] = flux1
                (<double*>np.PyArray_MultiIter_DATA(it, 4))[0] = fluxerr1
                (<short*>np.PyArray_MultiIter_DATA(it, 5))[0] = flag1

                #PyArray_MultiIter_NEXT is used to advance the iterator
                np.PyArray_MultiIter_NEXT(it)

        return sum, sumerr, flag

    else:
        rin, rout = bkgann

        # Require float arrays (see note above)
        rin = np.require(rin, dtype=np.float)
        rout = np.require(rout, dtype=np.float)

        # allocate ouput arrays
        shape = np.broadcast(x, y, r, rin, rout).shape
        sum = np.empty(shape, np.float)
        sumerr = np.empty(shape, np.float)
        flag = np.empty(shape, np.short)

        it = np.broadcast(x, y, r, rin, rout, sum, sumerr, flag)
        while np.PyArray_MultiIter_NOTDONE(it):

                x1 = (<double*>np.PyArray_MultiIter_DATA(it, 0))[0]
                y1 = (<double*>np.PyArray_MultiIter_DATA(it, 1))[0]
                r1 = (<double*>np.PyArray_MultiIter_DATA(it, 2))[0]
                rin1 = (<double*>np.PyArray_MultiIter_DATA(it, 3))[0]
                rout1 = (<double*>np.PyArray_MultiIter_DATA(it, 4))[0]

                status = sep_apercirc(ptr, eptr, mptr,
                                      dtype, edtype, mdtype, w, h,
                                      maskthresh, gain_, inflag,
                                      x1, y1, r1, subpix,
                                      &flux1, &fluxerr1, &area1, &flag1)
                _assert_ok(status)

                # background subtraction
                # Note that background output flags are not used.
                status = sep_apercircann(ptr, eptr, mptr,
                                         dtype, edtype, mdtype, w, h,
                                         maskthresh, gain_,
                                         inflag | SEP_MASK_IGNORE,
                                         x1, y1, rin1, rout1, 1,
                                         &bkgflux, &bkgfluxerr, &bkgarea,
                                         &bkgflag)

                flux1 -= bkgflux / bkgarea * area1
                bkgfluxerr = bkgfluxerr / bkgarea * area1
                fluxerr1 = fluxerr1*fluxerr1 + bkgfluxerr*bkgfluxerr

                (<double*>np.PyArray_MultiIter_DATA(it, 5))[0] = flux1
                (<double*>np.PyArray_MultiIter_DATA(it, 6))[0] = fluxerr1
                (<short*>np.PyArray_MultiIter_DATA(it, 7))[0] = flag1

                #PyArray_MultiIter_NEXT is used to advance the iterator
                np.PyArray_MultiIter_NEXT(it)

        return sum, sumerr, flag

def istruncated(np.ndarray flag not None):
    """True where 'aperture truncated' flag is set."""
    return (flag & SEP_APER_TRUNC) != 0

def hasmasked(np.ndarray flag not None):
    """True where 'aperture has masked pixel(s)' flag is set."""
    return (flag & SEP_APER_HASMASKED) != 0

def mask_ellipse(np.ndarray arr not None, x, y, scale=1.0, **kwargs):
    """mask_ellipse(arr, x, y, scale=1.0, **kwargs)

    Mask ellipse(s) in an array.

    Set array elements to True (1) if they fall within the given ellipse,
    defined by ``cxx*(x'-x)^2 + cyy*(y'-y)^2 + cxy*(x'-x')*(y'-y) = scale^2``

    Parameters
    ----------
    arr : np.ndarray
        Input array to be masked. Array is altered.
    x, y : array_like
        Center of ellipse(s).
    scale : array_like, optional
        Scale factor of ellipse(s). Default is 1.

    Additional Parameters
    ---------------------
    cxx, cyy, cxy : array_like
        Parameters defining the extent of the ellipse(s).
    """

    cdef int w, h
    cdef np.uint8_t[:,:] buf

    # only boolean arrays supported
    if not (arr.dtype == np.bool_ or arr.dtype == np.uint8):
        raise ValueError("Array data type not supported: {0:s}"
                         .format(arr.dtype))
    _check_array_get_dims(arr, &w, &h)
    buf = arr.view(dtype=np.uint8)

    # See note in apercirc on requiring specific array type
    x = np.require(x, dtype=np.float)
    y = np.require(y, dtype=np.float)
    scale = np.require(scale, dtype=np.float)

    if ('cxx' in kwargs and 'cyy' in kwargs and 'cxy' in kwargs):
        cxx = np.require(kwargs['cxx'], dtype=np.float)
        cyy = np.require(kwargs['cyy'], dtype=np.float)
        cxy = np.require(kwargs['cxy'], dtype=np.float)

        it = np.broadcast(x, y, cxx, cyy, cxy, scale)
        while np.PyArray_MultiIter_NOTDONE(it):
            sep_setellipse_ucc(<unsigned char *>&buf[0, 0], w, h,
                               (<double*>np.PyArray_MultiIter_DATA(it, 0))[0],
                               (<double*>np.PyArray_MultiIter_DATA(it, 1))[0],
                               (<double*>np.PyArray_MultiIter_DATA(it, 2))[0],
                               (<double*>np.PyArray_MultiIter_DATA(it, 3))[0],
                               (<double*>np.PyArray_MultiIter_DATA(it, 4))[0],
                               (<double*>np.PyArray_MultiIter_DATA(it, 5))[0],
                               1)
            np.PyArray_MultiIter_NEXT(it)

    else:
        raise ValueError("Must specify cxx, cyy and cxy keywords")
