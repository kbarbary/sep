SEP
===

*Source Extraction and Photometry in Python*

SEP makes available some of the algorithms in Source Extractor as
stand-alone functions and classes. These operate directly on in-memory
numpy arrays (no FITS files, configuration files, etc). It's derived
directly from (and tested against) the Source Extractor code base.

**Some features:**

- spatially variable background and noise estimation
- source extraction, with on-the-fly convolution and source deblending
- circular and elliptical aperture photometry
- extremely fast: implemented in C with Python bindings via Cython

**Additional features not in Source Extractor:**

- circular annulus and elliptical annulus apertures
- Local background subtraction in shape consistent with aperture
- exact pixel overlap mode in all aperture photometry functions
- ellipse masking


Installation
------------

**Requirements:**

- Tested on Python 2.6, 2.7, 3.3, 3.4
- numpy

**Install release version with pip:** ``pip install sep``

**Install release version with conda (64-bit Linux only):**
``conda install -c https://conda.binstar.org/kbarbary sep``

**Development version / source code:** http://github.com/kbarbary/sep


Usage Guide
-----------

.. toctree::
   :maxdepth: 1

   background
   detection
   apertures

Reference/API
-------------

**Background estimation & source detection**

.. autosummary::
   :toctree: api
   
   sep.Background
   sep.extract

**Aperture photometry**

.. autosummary::
   :toctree: api
   
   sep.sum_circle
   sep.sum_circann
   sep.sum_ellipse
   sep.sum_ellipann

**Aperture utilities**

.. autosummary::
   :toctree: api

   sep.kron_radius
   sep.flux_radius
   sep.mask_ellipse
   sep.ellipse_axes
   sep.ellipse_coeffs

.. note::

   The coordinate convention in SEP is that (0, 0) corresponds to the
   center of the first element of the data array. This agrees with the
   0-based indexing in Python and C.  However, note that
   this differs from the FITS convention where the center of the first
   element is at coordinates (1, 1). As Source Extractor deals with
   FITS files, its outputs follow the FITS convention. Thus, the
   coordinates from SEP will be offset from Source Extractor
   coordinates by -1 in x and y.

**Flags**

========================  ===========================================
Flag                      Description
========================  ===========================================
``sep.OBJ_MERGED``        object is result of deblending
``sep.OBJ_TRUNC``         object is truncated at image boundary
``sep.OBJ_SINGU``         x, y fully correlated in object
``sep.APER_TRUNC``        aperture truncated at image boundary
``sep.APER_HASMASKED``    aperture contains one or more masked pixels
``sep.APER_ALLMASKED``    aperture contains only masked pixels
``sep.APER_NONPOSITIVE``  aperture sum is negative in ``kron_radius``
========================  ===========================================

To see if a given flag is set in ``flags``::

    is_merged = (flags & sep.OBJ_MERGED) != 0
