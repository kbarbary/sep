#!/usr/bin/env py.test
from __future__ import print_function, division
# unicode_literals doesn't play well with numpy dtype field names

import os
import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_approx_equal
import sep

# Try to import any FITS reader
try:
    from fitsio import read as getdata
    NO_FITS = False
except:
    try:
        from astropy.io.fits import getdata
        NO_FITS = False
    except:
        NO_FITS = True

IMAGE_FNAME = os.path.join("data", "image.fits")
BACKIMAGE_FNAME = os.path.join("data", "back.fits")
IMAGECAT_FNAME = os.path.join("data", "image.cat")
IMAGECAT_DTYPE = [('number', np.int64),
                  ('x', np.float64),
                  ('y', np.float64),
                  ('flux', np.float64),
                  ('fluxerr', np.float64),
                  ('flags', np.int64)]
IMAGE_DTYPES = [np.float64, np.float32, np.int32]  # supported image dtypes

@pytest.mark.skipif(NO_FITS, reason="no FITS reader") 
def test_vs_sextractor():
    data = getdata(IMAGE_FNAME)

    bkg = sep.Background(data, bw=64, bh=64, fw=3, fh=3)

    # Test that SExtractor background is same as SEP:
    bkgarr = bkg.back(dtype=np.float32)
    refbackarr = getdata(BACKIMAGE_FNAME)
    assert_allclose(bkgarr, refbackarr, rtol=1.e-5)

    # Extract objects
    bkg.subfrom(data)
    objects = sep.extract(data, 1.5*bkg.globalrms)
    objects = np.sort(objects, order=['y'])

    # Read SExtractor result
    refobjects = np.loadtxt(IMAGECAT_FNAME, dtype=IMAGECAT_DTYPE)
    refobjects = np.sort(refobjects, order=['y'])

    # Found correct number of sources at the right locations?
    assert_allclose(objects['x'], refobjects['x'] - 1., atol=1.e-3)
    assert_allclose(objects['y'], refobjects['y'] - 1., atol=1.e-3)

    # Test flux
    flux, fluxerr, flag = sep.apercirc(data, objects['x'], objects['y'], 5.,
                                       err=bkg.globalrms)
    assert_allclose(flux, refobjects['flux'], rtol=2.e-4)
    assert_allclose(fluxerr, refobjects['fluxerr'], rtol=1.0e-5)

    assert sep.istruncated(flag).sum() == 4
    assert sep.hasmasked(flag).sum() == 0

def test_extract_noise_array():

    # Get some background-subtracted test data:
    data = getdata(IMAGE_FNAME)
    bkg = sep.Background(data, bw=64, bh=64, fw=3, fh=3)
    bkg.subfrom(data)

    # Ensure that extraction with constant noise array gives the expected
    # result.
    objects = sep.extract(data, 1.5*bkg.globalrms)
    noise = np.ones_like(data)
    objects2 = sep.extract(data, 1.5*bkg.globalrms, noise=noise)
    ndiff = 0
    for i in range(len(objects)):
        if objects[i] != objects2[i]:
            print(objects[i]['x'], objects[i]['y'])
            #print(i)
            #print(objects[i])
            #print(objects2[i])
    exit()


def test_apercirc_dtypes():
    naper = 100
    x = np.random.uniform(200., 800., naper)
    y = np.random.uniform(200., 800., naper)
    r = 3.
    fluxes = []
    for dt in IMAGE_DTYPES:
        data = np.ones((1000, 1000), dtype=dt)
        flux, fluxerr, flag = sep.apercirc(data, x, y, r)
        fluxes.append(flux)

    for i in range(1, len(fluxes)):
        assert_allclose(fluxes[0], fluxes[i])

def test_mask_ellipse():
    arr = np.zeros((20, 20), dtype=np.bool)

    # should mask 5 pixels:
    sep.mask_ellipse(arr, 10., 10., cxx=1.0, cyy=1.0, cxy=0.0, scale=1.001)
    assert arr.sum() == 5

    # should mask 13 pixels:
    sep.mask_ellipse(arr, 10., 10., cxx=1.0, cyy=1.0, cxy=0.0, scale=2.001)
    assert arr.sum() == 13


def test_masked_background():
    data = 0.1 * np.ones((6,6))
    data[1,1] = 1.
    data[4,1] = 1.
    data[1,4] = 1.
    data[4,4] = 1.

    mask = np.zeros((6,6), dtype=np.bool)

    # Background array without mask
    sky = sep.Background(data, bw=3, bh=3, fw=1, fh=1)
    bkg1 = sky.back()

    # Background array with all False mask
    sky = sep.Background(data, mask=mask, bw=3, bh=3, fw=1, fh=1)
    bkg2 = sky.back()

    # All False mask should be the same
    assert_allclose(bkg1, bkg2)

    # Masking high pixels should give a flat background
    mask[1, 1] = True
    mask[4, 1] = True
    mask[1, 4] = True
    mask[4, 4] = True
    sky = sep.Background(data, mask=mask, bw=3, bh=3, fw=1, fh=1)
    assert_approx_equal(sky.globalback, 0.1)
    assert_allclose(sky.back(), 0.1 * np.ones((6, 6)))
