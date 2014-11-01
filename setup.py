#!/usr/bin/env python
import os
from glob import glob
from distutils.core import setup
from distutils.extension import Extension
import numpy

if os.path.exists("sep.pyx"):
    USE_CYTHON = True
    fname = "sep.pyx"
else:
    USE_CYTHON = False
    fname = "sep.c"

sourcefiles = [fname] + glob(os.path.join("src", "*.c"))
include_dirs=[numpy.get_include(), "src"]
extensions = [Extension("sep", sourcefiles, include_dirs=include_dirs)]
if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

description = "Source extraction and photometry library"

classifiers = [
    "Development Status :: 5 - Production/Stable",
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: GNU Lesser General Public License v3 "
    "or later (LGPLv3+)",
    "Topic :: Scientific/Engineering :: Astronomy",
    "Intended Audience :: Science/Research"]

setup(name="sep", 
      version="0.1.0",
      description=description,
      long_description=description,
      license = "LGPLv3+",
      classifiers=classifiers,
      url="https://github.com/kbarbary/sep-python",
      author="Kyle Barbary",
      author_email="kylebarbary@gmail.com",
      ext_modules=extensions)
