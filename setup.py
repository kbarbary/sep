#!/usr/bin/env python

import os
import re
import sys
from glob import glob

from setuptools import setup
from setuptools.dist import Distribution
from setuptools.extension import Extension

# Synchronize version from code.
fname = "sep.pyx"
version = re.findall(r"__version__ = \"(.*?)\"", open(fname).read())[0]

# Detect if setup.py is being run with an argument that doesn't require
# building anything. (egg info, help commands, etc)
options = Distribution.display_option_names + ['help-commands', 'help']
is_non_build = (
    any('--' + opt in sys.argv for opt in options)
    or len(sys.argv) == 1
    or sys.argv[1] in ('egg_info', 'clean', 'help')
)

# extension module(s): only add if we actually need to build, because we need
# to import numpy and cython to build, and we'd rather non-build commands
# work when those dependencies are not installed.
if is_non_build:
    extensions = None
else:
    import numpy
    from Cython.Build import cythonize
    sourcefiles = [fname] + glob(os.path.join("src", "*.c"))
    headerfiles = glob(os.path.join("src", "*.h"))
    include_dirs = [numpy.get_include(), "src"]
    extensions = [Extension("sep", sourcefiles, include_dirs=include_dirs,
                            depends=headerfiles)]
    extensions = cythonize(extensions)



description = "Astronomical source extraction and photometry library"
long_description = "http://sep.readthedocs.org"

classifiers = [
    "Development Status :: 5 - Production/Stable",
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
      python_requires='>=3.5',
      install_requires=['numpy'],
      ext_modules=extensions)
