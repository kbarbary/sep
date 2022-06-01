v1.2.1 (1 June 2022)
====================

* Same as v1.2.0 but with new wheels for Python 3.10 and AArch64.

v1.2.0 (1 May 2021)
===================

* Changed `numpy.float` and `numpy.int` types for deprecations in numpy 1.20 (#96).

* Make it possible to safely invoke C library from multiple threads on
  independent inputs.

  Global config functions such as `set_sub_object_limit()`
  and `set_extract_pixstack()` still configure global params
  (once for all threads), while other functions will retain their data
  in thread-local storages, so they can be invoked from multiple threads as
  long as they work on independent structures.

  Library compilation will now require a C11 compatible compiler, which should
  be nowadays available on all supported platforms.

* Mark some pointer parameters with `const *`. This is a backward-compatible
  change, but makes it easier to extract constants that can be safely shared
  between multiple threads and/or invocations.

v1.1.1 (6 January 2021)
=======================

* Same as v1.1.0 but with wheels built and uploaded to PyPI. Please report if you
  have problems with wheels.


v1.1.0 (3 January 2021)
=======================

* Add segmentation masking to the photometry and kron/auto functions (#69).

* Add functions `sep.set_sub_object_limit(limit)` and `sep.get_sub_object_limit()`
  for modifying and retrieving the sub-object deblending limit. Previously this
  parameter was hard-coded to 1024. 1024 is now the default value.

* This and future versions are now Python 3 only. Python 2 is no longer
  supported.

* Modernize setup.py with pyproject.toml


v1.0.3 (12 June 2018)
=====================

* Fix double-free bug in sep_extract() arising when an error status occurs
  and convolution is on. (#56)

* Work around numpy dependency in setup. (#59)


v1.0.2 (19 September 2017)
==========================

* Fix makefile so that `make install` works on OS X for the C library.
  Python module and C code are unchanged.


v1.0.1 (10 July 2017)
=====================

* Fix bug when using masked filter and noise array where objects with member
  pixels at end of image (maximum y coordinate) were erroneously missed.


v1.0.0 (30 September 2016)
==========================

* Remove features deprecated in previous versions.

* Fix bug in Background.rms() giving nonsensical results.

v0.6.0 (25 August 2016)
=======================

* New, more coherent C API. This change should be transparent to users
  of the Python module.

* Add variance uncertainty parameters `errx2`, `erry2` and `errxy` to
  output of `sep.extract()`.

* Add a minimum sigma to `sep.winpos()` to match Source Extractor
  behavior.

* Fix use of boolean masks in `sep.kron_radius()`. Formerly, using a
  boolean mask resulted in nonsense results.

* Fix segfault in `Background.back()` when box size is same as image size.

* Fix bug in creating long error messages on Python 3.

v0.5.2 (4 January 2016)
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
