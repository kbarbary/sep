SEP
===

Python and C library for Source Extraction and Photometry

[![Build Status](http://img.shields.io/travis/kbarbary/sep.svg?style=flat-square&label=linux)](https://travis-ci.org/kbarbary/sep)
[![Build status](https://img.shields.io/appveyor/ci/kbarbary/sep.svg?style=flat-square&label=windows)](https://ci.appveyor.com/project/kbarbary/sep/branch/master)
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

- Tested on Python 2.6, 2.7, 3.3, 3.4, 3.5
- numpy

**Install release version:**

```
pip install sep
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

_Note: The build process only works on Linux and OS X._

**Build:** To build the C library from source, you must have

```
make
```

**Run tests:**

```
make test
```

**Install** The static library and header can be installed with

```
make install
make PREFIX=/path/to/prefix install
```

This will install the shared and static library in `/path/to/prefix/lib`
and header file in `/path/to/prefix/include`. The default prefix is
`/usr/local`.

**API:** The C library API is documented in the header file
[sep.h](src/sep.h).

License
-------

The license for all parts of the code derived from SExtractor is
LGPLv3. The license for code not derived from SExtractor is MIT. The
license for the library as a whole is therefore LGPLv3. The license
for each file is explicitly stated at the top of each file and the
full text of the licenses can be found in `licenses`.

Citation
--------

If you use SEP in a publication, please cite the following DOI. The link provides a variety of citation styles and BibTeX export.

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.15669.svg)](http://dx.doi.org/10.5281/zenodo.15669)

You may also wish to cite the original SourceExtractor paper (Bertin 1996).

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

**In the Python interface, why do I have to byteswap data when using
astropy.io.fits?**

This occurs because FITS files have big-endian [byte
order](http://en.wikipedia.org/wiki/Endianness), whereas most widely
used CPUs have little-endian byte order. In order for the CPU to
operate on the data, it must be byte swapped at *some point.* Some
FITS readers such as [fitsio](http://github.com/esheldon/fitsio) do
the byte swap immediately when reading the data from disk to memory,
returning numpy arrays in native (little-endian) byte order. However,
astropy.io.fits does not (for reasons having to do with memory
mapping). Most of the time you never notice this because when you do
any numpy operations on such arrays, numpy uses an intermediate buffer
to byteswap the array behind the scenes and returns the result as a
native byte order array. Internally, SEP is not using numpy
operations; it's just getting a pointer to the data in the array and
passing it to C code. As the C code does not include functionality to
do buffered byte swapping, the input array must already be in native
byte order.

It would be possible to add buffered byteswapping capability to the
SEP code, but it would increase the code complexity. A simpler
alternative would be to make a byteswapped copy of the entire input
array, whenever necessary. However, this would significantly increase
memory use, and would have to be done repeatedly in multiple SEP
functions: `Background`, `extract`, `sum_circle`, etc. Each would make
a copy of the entire data array. Given these considerations, it seemed
best to just explicitly tell the user to do the byteswap operation
themselves so they could just do it once, immediately after reading in
the data.
