#!/usr/bin/env python
import os
from glob import glob
from distutils.core import setup
from distutils.extension import Extension
import re

import numpy

if os.path.exists("sep.pyx"):
    USE_CYTHON = True
    fname = "sep.pyx"
else:
    USE_CYTHON = False
    fname = "sep.c"

sourcefiles = [fname] + glob(os.path.join("src", "*.c"))
headerfiles = glob(os.path.join("src", "*.h"))
include_dirs=[numpy.get_include(), "src"]
extensions = [Extension("sep", sourcefiles, include_dirs=include_dirs,
                        depends=headerfiles)]
if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

# Synchronize version from code.
version = re.findall(r"__version__ = \"(.*?)\"", open(fname).read())[0]

description = "Astronomical source extraction and photometry library"
long_description = "http://sep.readthedocs.org"

classifiers = [
    "Development Status :: 5 - Production/Stable",
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: GNU Lesser General Public License v3 "
    "or later (LGPLv3+)",
    "Topic :: Scientific/Engineering :: Astronomy",
    "Intended Audience :: Science/Research"]

setup(name="sep", 
      version=version,
      description=description,
      long_description=long_description,
      license="LGPLv3+",
      classifiers=classifiers,
      url="https://github.com/kbarbary/sep",
      author="Kyle Barbary",
      author_email="kylebarbary@gmail.com",
      ext_modules=extensions)
