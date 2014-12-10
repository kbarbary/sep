Python module versions

v0.2.0 (unreleased)
===================

* **[API change]** Change `mask_ellipse()` parameters from ellipse
  coefficients to ellipse axes and position angle, to match aperture
  functions.

* **[API change]** Change `apercirc()` to `sum_circle()`.
  (`apercirc()` still available as undocumented alias temporarily.)

* **[API change]** `theta` field in `extract()` output is now in
  radians rather than degrees, for compatibility with new ellipse
  aperture functions.

* Add `sum_circann()`, `sum_ellipse()`, `sum_ellipann()`,
  `kron_radius()`, `ellipse_coeffs()`, `ellipse_axes()` functions.

* Exact mode aperture photometery in all functions, with `subpix=0`.

* Enable variable thresholding in `sep.extract`. [#11]

* Fix bug in background masking. This bug impacted masking in all
  functions that used masking. Also affected C library.

* More detail in error messages coming from within the C library.
  More helpful error message for non-native byteorder arrays.

* Add ability to change pixel stack size used in `extract()`, with
  `set_extract_pixstack()` function

v0.1.0 (11 August 2014)
=======================

This is the first official release.
