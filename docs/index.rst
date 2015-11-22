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
- fast: implemented in C with Python bindings via Cython

**Additional features not in Source Extractor:**

- Optimized matched filter for variable noise in source extraction
- circular annulus and elliptical annulus apertures
- Local background subtraction in shape consistent with aperture
- exact pixel overlap mode in all aperture photometry functions
- ellipse masking


Installation
------------

SEP supports both Python 3 and Python 2 and requires only numpy.
To install::

    pip install --no-deps sep

**Development version / source code:** http://github.com/kbarbary/sep


Usage Guide
-----------

.. toctree::
   :maxdepth: 1

   detection
   apertures

.. toctree::
   :hidden:

   reference

For complete API documentation, see :doc:`reference`.


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
