---
title: 'SEP: Source Extractor as a library'
tags:
  - astronomy
authors:
 - name: Kyle Barbary
   orcid: 0000-0002-2532-3696
   affiliation: University of California, Berkeley
date: 29 August 2016
bibliography: paper.bib
---

# Summary

Source Extractor [@bertin96; @sextractor] is a widely used
command-line program for segmentation and analysis of astronomical
images. It reads in FITS format files, performs a configurable series
of tasks, including background estimation, source detection,
deblending and a wide array of source measurements, and finally
outputs a FITS format catalog file.

While Source Extractor is highly useful, the fact that it can only be
used as an executable -- reading input files, producing output files
and controlled by a limited set of configuration options specified in
another file -- can limit its applicability or lead to awkward
workflows. There is often a desire to have programmatic access to
perform one or more of the above tasks on in-memory images as part of
a larger custom analysis.

SEP makes available the core algorithms of Source Extractor in a
library of stand-alone functions and classes. These operate directly
on in-memory arrays (no FITS files or configuration files).  The code
is derived from the Source Extractor code base (written in C) and aims
to produce results compatible with Source Extractor whenever possible.
SEP consists of a C library with no dependencies outside the standard
library, and a Python module that wraps the C library in a Pythonic
API. The Python wrapper operates on NumPy arrays with NumPy as its
only dependency. It is generated using Cython.

From Source Extractor, SEP includes background estimation, image
segmentation (including on-the-fly filtering and source deblending),
aperture photometry in circular and elliptical apertures, and source
measurements such as Kron radius, "windowed" position fitting, and
half-light radius. Additionally, several features not in Source
Extractor have been added:

- Optimized matched filter for variable noise in source extraction.
- Circular annulus and elliptical annulus aperture photometry functions.
- Local background subtraction in shape consistent with aperture in aperture
  photometry functions.
- Exact pixel overlap mode in all aperture photometry functions.
- Masking of elliptical regions on images.

Finally, note that SEP is essentially a fork of Source Extractor that
has already diverged significantly from the original code base. One
might ask why SEP is not part of Source Extractor itself: the
command-line interface in Source Extractor could be built on top of
such a library. The answer is that a vast array of changes were
necessary in order to expose the functionality as stand-alone C
functions. It would be a lot of work to rewrite Source Extractor
itself in this way, with little gain for the executable itself.


# References
