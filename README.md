SEP
===

[![Build Status](https://api.travis-ci.org/kbarbary/sep.svg?branch=master)](https://travis-ci.org/kbarbary/sep)

C library for Source Extraction and Photometry

*"... [it's] an SEP: Somebody Else's Problem."  
"Oh, good. I can relax then."*

[Source Extractor](http://www.astromatic.net/software/sextractor) is
great, but sometimes you want to use a few of the pieces from it
without running the entire executable. This library implements a few
of the algorithms used in SExtractor as stand-alone pieces. So far
this includes:

* global background estimation
* source detection
* circular aperture photometry

In the future, the library may also include other functions
related to photometry that are not in Source Extractor.

SEP is designed both to be used in C programs and to be wrapped in
higher-level languages such as Python or Julia. To make the latter
easier, SEP has minimal dependencies.

Build
-----

To build, you must have [scons](http://scons.org/) available on your
system. In the top level directory,

```
scons          # build the library
scons --clean  # clean the built library
```

Test
----

The test program requires that the `cfitsio` library and
development header be installed. On Ubuntu:
```
sudo apt-get install libcfitsio3-dev
```

In the top level directory:

```
scons test          # build the test executable
scons test --clean  # clean the built test executable
```

To run the test executable, first ensure that the built shared library
(in `src`) can be found. For example, on linux, put the path to the
library in the `LD_LIBRARY_PATH` environment variable. Then do

```
./runtests
```

API
---

The API is documented in the header file [sep.h](src/sep.h).

License
-------

The license for all parts of the code derived from SExtractor is
LGPLv3. The license for code not derived from SExtractor is MIT. The
license for the library as a whole is therefore LGPLv3. The license
for each file is explicitly stated at the top of each file and the
full text of the licenses can be found in `licenses`.

FAQ
---

**Why isn't this part of Source Extractor?**

Source Extractor isn't designed as a library implementation with an
executable built on top of the library. In Source Extractor, background
estimation, object detection and photometry are deeply integrated into the
Source Extractor executable. Many changes to the code were necessary in
order to put the functionality in stand-alone C functions.

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

It should be in the same ballpark, as a lot of the core code is the
same.  Source Extractor has the advantage of doing all the operations
(detection and analysis) simultaneously on each image section, which
may confer CPU cache advantages, but this hasn't been tested.

**What happens when Source Extractor is updated in the future?**

SEP can be considered a fork of the Source Extractor codebase: it's
development will not track that of Source Extractor in any automated
way. However, the algorithms implemented so far in SEP are stable in
Source Extractor: the SEP code was forked from v2.18.11, yet tested
against the results of v2.8.6.