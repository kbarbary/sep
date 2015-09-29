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

SEP supports both Python 3 and Legacy Python (Python 2) and requires
only numpy. To install::

    pip install --no-deps sep

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
   sep.winpos
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

License and Citation
--------------------

The license for SEP is the Lesser GNU Public License (LGPL), granted
with the permission from the original author of Source Extractor.

If you use SEP in a publication, please cite the DOI
`10.5281/zenodo.15669 <http://dx.doi.org/10.5281/zenodo.15669>`_. The
link provides a variety of citation styles and BibTeX export. For example::

    @misc{kyle_barbary_2015_15669,
      author       = {Kyle Barbary and
                      Kyle Boone and
                      Christoph Deil},
      title        = {sep: v0.3.0},
      month        = feb,
      year         = 2015,
      doi          = {10.5281/zenodo.15669},
      url          = {http://dx.doi.org/10.5281/zenodo.15669}
    }


You may also wish to cite the original SExtractor paper: `Bertin &
Arnouts 1996 <http://adsabs.harvard.edu/abs/1996A%26AS..117..393B>`_.
