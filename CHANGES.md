Python module versions

v0.5.3 (unreleased)
===================

- Fix segfault in `Background.back()` when box size is same as image size.
- Fix bug in creating long error messages on Python 3.

v0.5.2 (4 January 2015)
=======================

Adds OS X and Windows support.

v0.5.1 (30 November 2015)
=========================

Bugfix release for problem in setup.py in packaged code.

v0.5.0 (22 November 2015)
=========================

* `sep.extract()` now uses a more correct matched filter algorithm in the
  presence of a noise array, rather than simple convolution. The `conv`
  keyword has been changed to `filter_kernel` to reflect this, and a
  `filter_type` keyword has been added to allow selecting the old behavior
  of simple convolution.

* `sep.extract()` now accepts a `mask` keyword argument.

* `sep.extract()` can now return a segmentation map.

* Special methods added to allow `data - bkg` and `np.array(bkg)` where
  `bkg` is a Background object.

v0.4.1 (10 November 2015)
=========================

Bugfix release, fixing error estimate in `sep.sum_circle` and
`sep.sum_ellipse` when `bkgann` keyword argument is given.

v0.4.0 (1 June 2015)
====================

* New `sep.winpos()` function.

v0.3.0 (23 February 2015)
=========================

* New `sep.flux_radius()` function.

v0.2.0 (13 December 2014)
=========================

* **[breaking change]** `theta` field in `extract()` output is now in
  radians rather than degrees, for compatibility with new ellipse
  aperture functions.

* **[deprecation]** Change `mask_ellipse()` parameters from ellipse
  coefficients to ellipse axes and position angle, to match aperture
  functions. (Old behavior still works as well.)

* **[deprecation]** Change `apercirc()` to `sum_circle()`, to match
  new aperture functions. (Old name, `apercirc()`, still works.)

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
