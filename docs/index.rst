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

**Install:** ``pip install sep``

**Source code:** http://github.com/kbarbary/sep

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
   sep.mask_ellipse
   sep.ellipse_axes
   sep.ellipse_coeffs

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
