#!/usr/bin/env python
from __future__ import print_function

from os.path import join
import time
import numpy as np
import sep

# try to import photutils for comparison timing
try:
    import photutils
    HAVE_PHOTUTILS = True
except ImportError:
    HAVE_PHOTUTILS = False

# Try to import any FITS reader
try:
    from fitsio import read as getdata
    HAVE_FITS = True
    NEED_BYTESWAP = False
except ImportError:
    try:
        from astropy.io.fits import getdata
        HAVE_FITS = True
    except ImportError:
        HAVE_FITS = False

CONDENSED = True

if HAVE_FITS:
    rawdata = getdata(join("data", "image.fits"))  # original is 256 x 256
    data = np.tile(rawdata, (4, 4))

    print("test image shape:", data.shape)
    print("test image dtype:", data.dtype)

    t0 = time.time()
    bkg = sep.Background(data)  # estimate background
    t1 = time.time()
    print("measure background: {0:6.2f} ms".format((t1-t0) * 1.e3))

    t0 = time.time()
    bkg.subfrom(data)  # subtract it
    t1 = time.time()
    print("subtract background: {0:6.2f} ms".format((t1-t0) * 1.e3))

    t0 = time.time()
    backarr = bkg.back(dtype=np.float64)  # background
    t1 = time.time()
    print("background array: {0:6.2f} ms".format((t1-t0) * 1.e3))

    t0 = time.time()
    rmsarr = bkg.rms()
    t1 = time.time()
    print("rms array: {0:6.2f} ms".format((t1-t0) * 1.e3))

    t0 = time.time()
    objects = sep.extract(data, 1.5 * bkg.globalrms)
    t1 = time.time()
    print("extract: {0:6.2f} ms  [{1:d} objects]"
          .format((t1-t0) * 1.e3, len(objects)))

#--------------------------------------------------------------------------
# Background subtraction

print("")
if HAVE_PHOTUTILS:
    print("sep version:      ", sep.__version__)
    print("photutils version:", photutils.__version__)
    print("""
| test                    | sep             | photutils       | ratio  |
|-------------------------|-----------------|-----------------|--------|""")
    blankline = \
"|                         |                 |                 |        |"

else:
    print("sep version: ", sep.__version__)
    print("""
| test                    | sep             |
|-------------------------|-----------------|""")
    blankline = "|                         |                 |"

nloop = 50

for ntile in [4]:
    data = np.tile(rawdata, (ntile, ntile))
    line = "| {0:4d}^2 image background |".format(data.shape[0])

    t0 = time.time()
    for _ in range(0, nloop):
        bkg = sep.Background(data)
    t1 = time.time()
    t_sep = (t1-t0) * 1.e3 / nloop
    line += "      {0:7.2f} ms |".format(t_sep)

    if HAVE_PHOTUTILS:
        t0 = time.time()
        for _ in range(0, nloop):
            try:
                bkg = photutils.Background(data, (64, 64))  # estimate background
            except AttributeError:
                bkg = photutils.Background2D(data, (64, 64))  # estimate background
        t1 = time.time()
        t_pu = (t1-t0) * 1.e3 / nloop
        line += "      {0:7.2f} ms | {1:6.2f} |".format(t_pu, t_pu/t_sep)

    print(line)

#------------------------------------------------------------------------------
# Circular aperture photometry benchmarks

if not CONDENSED:
    print(blankline)
    line = "| **aperture photometry** |                 |"
    if HAVE_PHOTUTILS:
        line += "                 |        |"
    print(line)

naper = 1000
data = np.ones((2000, 2000), dtype=np.float32)
x = np.random.uniform(200., 1800., naper)
y = np.random.uniform(200., 1800., naper)

if CONDENSED:
    r_list = [5.]
    subpix_list = [(5, "subpixel", "subpix=5"), (0, "exact", "exact")]
else:
    r_list = [3., 5., 10., 20.]
    subpix_list = [(1, "center", "subpix=1"), (5, "subpixel", "subpix=5"),
                   (0, "exact", "exact")]

for r in r_list:
    for subpix, method, label in subpix_list:

        line = "| circles  r={0:2d}  {1:8s} |".format(int(r), label)

        t0 = time.time()
        for _ in range(0, nloop):
            flux, fluxerr, flag = sep.sum_circle(data, x, y, r, subpix=subpix)
        t1 = time.time()
        t_sep = (t1-t0) * 1.e6 / naper / nloop
        line += " {0:7.2f} us/aper |".format(t_sep)

        if HAVE_PHOTUTILS:
            apertures = photutils.CircularAperture((x, y), r)
            t0 = time.time()
            for _ in range(0, nloop):
                res = photutils.aperture_photometry(
                    data, apertures, method=method, subpixels=subpix)
            t1 = time.time()
            t_pu = (t1-t0) * 1.e6 / naper / nloop
            line += " {0:7.2f} us/aper | {1:6.2f} |".format(t_pu, t_pu/t_sep)

        print(line)

if not CONDENSED:
    print(blankline)

a = 1.
b = 1.
theta = np.pi/4.

for r in r_list:
    for subpix, method, label in subpix_list:
        line = "| ellipses r={0:2d}  {1:8s} |".format(int(r), label)

        t0 = time.time()
        for _ in range(0, nloop):
            flux, fluxerr, flag = sep.sum_ellipse(data, x, y, a, b, theta, r,
                                                subpix=subpix)
        t1 = time.time()
        t_sep = (t1-t0) * 1.e6 / naper / nloop
        line += " {0:7.2f} us/aper |".format(t_sep)

        if HAVE_PHOTUTILS:
            apertures = photutils.EllipticalAperture((x, y), a*r, b*r, theta)
            t0 = time.time()
            for _ in range(0, nloop):
                res = photutils.aperture_photometry(
                    data, apertures, method=method, subpixels=subpix)
            t1 = time.time()
            t_pu = (t1-t0) * 1.e6 / naper / nloop
            line += " {0:7.2f} us/aper | {1:6.2f} |".format(t_pu, t_pu/t_sep)

        print(line)
