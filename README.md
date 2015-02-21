SEP
===

Python and C library for Source Extraction and Photometry

[![Build Status](http://img.shields.io/travis/kbarbary/sep.svg?style=flat-square)](https://travis-ci.org/kbarbary/sep)
[![PyPI](https://img.shields.io/pypi/v/sep.svg?style=flat-square)](https://pypi.python.org/pypi/sep)


*"... [it's] an SEP: Somebody Else's Problem."  
"Oh, good. I can relax then."*

[Source Extractor](http://www.astromatic.net/software/sextractor) is
great, but sometimes you want to use a few of the pieces from it
without running the entire executable. SEP makes available some of the
algorithms in SExtractor as stand-alone functions and classes. These
operate directly on in-memory arrays (no FITS files, configuration
files, etc). The code is derived directly from the Source Extractor
code base.

SEP can be used from either Python or directly from C. See below for
language-specific build and usage instructions.

Python
------

**Documentation:** http://sep.readthedocs.org/

**Requirements:**

- Tested on Python 2.6, 2.7, 3.3, 3.4
- numpy

**Install release version:**

With pip:
```
pip install sep
```

With conda (currently 64-bit linux only):
```
conda install -c https://conda.binstar.org/kbarbary sep
```

**Install development version:**

Bulding the development verion (from github) requires Cython (v0.16 or
higher).  Build and install in the usual place:

```
./setup.py install
```

**Run tests:** To run the tests, execute `./test.py` in the top-level
directory.  Requires the `pytest` package. Some tests require a FITS
reader (either fitsio or astropy) and will be skipped if neither is
present.

C Library
---------

_Note: The C library should not yet be considered stable. (To my knowledge,
no one is using it directly.)_

**Build:** To build the C library from source, you must have
[scons](http://scons.org/) available on your system. In the top level
directory,

```
scons          # build the library
scons --clean  # clean the built library
```

**Run tests:** The test program requires that the `cfitsio` library
and development header be installed. On Ubuntu `sudo apt-get install
libcfitsio3-dev` should do it. In the top level directory:

```
scons ctest             # build the test executable
cd ctest && ./runtests  # run tests
scons ctest --clean     # clean the built test executable
```

Note: before *running* the tests, ensure that the built shared library
in `src` can be found. On linux, you can do this by putting the path
to the library in the `LD_LIBRARY_PATH` environment variable.


**Install or link:** The static library and header can be installed with

```
scons install --prefix=/path/to/prefix
```

This will install the library in `/path/to/prefix/lib` and header file
in `/path/to/prefix/include`. If you wish to link against the static
library without installing, it will be found in the `src` subdirectory
after building.

The shared library cannot yet be automatically installed.

**API:** The C library API is documented in the header file
[sep.h](src/sep.h).

License
-------

The license for all parts of the code derived from SExtractor is
LGPLv3. The license for code not derived from SExtractor is MIT. The
license for the library as a whole is therefore LGPLv3. The license
for each file is explicitly stated at the top of each file and the
full text of the licenses can be found in `licenses`.

FAQ
---

**Why isn't the C library part of Source Extractor?**

Source Extractor is *not* designed as a library with an
executable built on top of the library. In Source Extractor, background
estimation, object detection and photometry are deeply integrated into the
Source Extractor executable. Many changes to the code were necessary in
order to put the functionality in stand-alone C functions. It's too much
to ask of the Source Extractor developer to rewrite large parts of the 
core of the Source Extractor program with little gain for the executable.

**What sort of changes?**

- Source Extractor reads in only a small portion of each image at a
  time.  This allows it to keep its memory footprint extremely low and
  to operate on images that are much larger than the system's physical
  memory. It also means that a FITS reader is deeply integrated into
  the code.  SEP operates on images in memory, so all the FITS I/O
  machinery in Source Extractor is not used here.

- Error handling: When it encounters a problem, Source Extractor
  immediately exits with an error message. This is fine for an
  executable, but a library function doesn't have that luxury. Instead
  it must ensure that allocated memory is freed and return an error
  code.

- Options: Source Extractor has many options that affect its behavior. These
  are stored in a global structure used throughout the executable. In SEP,
  options for a particular function are passed as function parameters.

- Array types: Source Extractor can operate on FITS images containing various
  types of data (float, double, int, etc). Internally, it does this by
  converting all data to `float` immediately when reading from disk.
  SEP does something similar, but in memory: SEP functions typically convert
  input arrays to float on the fly within each function, then perform
  all operations as floating point.

**Is SEP as fast as Source Extractor?**

It's fast. It should be similar to Source Extractor as a lot of the code
is identical. Source Extractor has the advantage of doing all the
operations (detection and analysis) simultaneously on each image
section, which may confer CPU cache advantages, but this hasn't been
tested at all. On the other hand, depending on your usage SEP might
let you avoid writing files to disk, which is likely to be a bigger
win.

**What happens when Source Extractor is updated in the future?**

SEP can be considered a fork of the Source Extractor codebase: it's
development will not track that of Source Extractor in any automated
way. However, the algorithms implemented so far in SEP are stable in
Source Extractor: the SEP code was forked from v2.18.11, yet it is tested
against the results of v2.8.6. This indicates that the algorithms have
not changed in SExtractor over the last few years.

**Is it "sep" or "ess eee pee"?**

I don't really care, just use it. :)