SEP
===

*Source Extraction and Photometry in Python*

SEP makes available some of the algorithms in Source Extractor as
stand-alone functions and classes. These operate directly on in-memory
numpy arrays (no FITS files, configuration files, etc).

**Some features:**

- spatially variable background and noise estimation
- optional mask in background estimation
- source extraction, with on-the-fly convolution, and source deblending
- variable detection threshold in source extraction
- circular and elliptical aperture photometry
- Extremely fast

**Features not in Source Extractor:**

- circular annulus and elliptical annulus photometry
- Local background subtraction in shape consistent with aperture
- exact pixel overlap mode in circular and elliptical apertures
- ellipse masking

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


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

