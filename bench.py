#!/usr/bin/env python

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
except:
    try:
        from astropy.io.fits import getdata
        HAVE_FITS = True
    except:
        HAVE_FITS = False


if HAVE_FITS:
    data = getdata(join("data", "image.fits"))  # original is 256 x 256
    data = np.tile(data, (4, 4))
    
    print "test image shape:", data.shape
    print "test image dtype:", data.dtype

    t0 = time.time()
    bkg = sep.Background(data) # estimate background
    t1 = time.time()
    print "measure background: {0:6.2f} ms".format((t1-t0) * 1.e3)

    t0 = time.time()
    bkg.subfrom(data)  # subtract it
    t1 = time.time()
    print "subtract background: {0:6.2f} ms".format((t1-t0) * 1.e3)

    t0 = time.time()
    backarr = bkg.back(dtype=np.float64)  # background
    print backarr.dtype
    t1 = time.time()
    print "background array: {0:6.2f} ms".format((t1-t0) * 1.e3)

    t0 = time.time()
    rmsarr = bkg.rms()
    t1 = time.time()
    print "rms array: {0:6.2f} ms".format((t1-t0) * 1.e3)

    t0 = time.time()
    objects = sep.extract(data, 1.5 * bkg.globalrms)
    t1 = time.time()
    print ("extract: {0:6.2f} ms  [{1:d} objects]"
           .format((t1-t0) * 1.e3, len(objects)))

#--------------------------------------------------------------------
# Aperture photometry benchmarks

datasize = (2000,2000)
naper = 1000

data = np.ones(datasize, dtype=np.float32)
x = np.random.uniform(200., 1800., naper)
y = np.random.uniform(200., 1800., naper)

for r in [3., 5., 10., 20., 100.]:
    print "r = {0:.1f}".format(r)

    if HAVE_PHOTUTILS:
        apertures = [photutils.CircularAperture(r)] * len(x)
        t0 = time.time()
        flux = photutils.aperture_photometry(data, x, y, apertures,
                                             method='subpixel',
                                             subpixels=5)
        t1 = time.time()
        print ("  photutils: {0:6.2f} us/aperture"
               .format((t1-t0) * 1.e6 / naper))

    rs = r * np.ones(naper, dtype=np.float)
    t0 = time.time()
    flux, fluxerr, flag = sep.apercirc(data, x, y, rs)
    t1 = time.time()
    print ("  sep:     {0:6.2f} us/aperture"
           .format((t1-t0) * 1.e6 / naper))
