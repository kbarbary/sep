Reference/API
=============


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

**Low-level utilities**

.. autosummary::
   :toctree: api

   sep.get_extract_pixstack
   sep.set_extract_pixstack

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

.. note::

   The coordinate convention in SEP is that (0, 0) corresponds to the
   center of the first element of the data array. This agrees with the
   0-based indexing in Python and C.  However, note that
   this differs from the FITS convention where the center of the first
   element is at coordinates (1, 1). As Source Extractor deals with
   FITS files, its outputs follow the FITS convention. Thus, the
   coordinates from SEP will be offset from Source Extractor
   coordinates by -1 in x and y.
