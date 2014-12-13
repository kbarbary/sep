#!/usr/bin/env py.test
from __future__ import print_function, division
# unicode_literals doesn't play well with numpy dtype field names

import os
import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal, assert_approx_equal
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
                  ('flux_aper', np.float64),
                  ('fluxerr_aper', np.float64),
                  ('kron_radius', np.float64),
                  ('flux_auto', np.float64),
                  ('fluxerr_auto', np.float64),
                  ('flags', np.int64)]
SUPPORTED_IMAGE_DTYPES = [np.float64, np.float32, np.int32]

def assert_allclose_structured(x, y):
    """Assert that two structured arrays are close.

    Compares floats relatively and everything else exactly.
    """
    assert x.dtype == y.dtype
    for name in x.dtype.names:
        if np.issubdtype(x.dtype[name], float):
            assert_allclose(x[name], y[name])
        else:
            assert_equal(x[name], y[name])

# If we have a FITS reader, read in the necessary test images
if not NO_FITS:
    image_data = getdata(IMAGE_FNAME)
    image_refback = getdata(BACKIMAGE_FNAME)

# -----------------------------------------------------------------------------
# Test versus Source Extractor results

@pytest.mark.skipif(NO_FITS, reason="no FITS reader") 
def test_vs_sextractor():
    data = np.copy(image_data)  # make an explicit copy so we can 'subfrom'
    bkg = sep.Background(data, bw=64, bh=64, fw=3, fh=3)

    # Test that SExtractor background is same as SEP:
    bkgarr = bkg.back(dtype=np.float32)
    assert_allclose(bkgarr, image_refback, rtol=1.e-5)

    # Extract objects
    bkg.subfrom(data)
    objs = sep.extract(data, 1.5*bkg.globalrms)
    objs = np.sort(objs, order=['y'])

    # Read SExtractor result
    refobjs = np.loadtxt(IMAGECAT_FNAME, dtype=IMAGECAT_DTYPE)
    refobjs = np.sort(refobjs, order=['y'])

    # Found correct number of sources at the right locations?
    assert_allclose(objs['x'], refobjs['x'] - 1., atol=1.e-3)
    assert_allclose(objs['y'], refobjs['y'] - 1., atol=1.e-3)

    # Test aperture flux
    flux, fluxerr, flag = sep.sum_circle(data, objs['x'], objs['y'], 5.,
                                         err=bkg.globalrms)
    assert_allclose(flux, refobjs['flux_aper'], rtol=2.e-4)
    assert_allclose(fluxerr, refobjs['fluxerr_aper'], rtol=1.0e-5)

    # check if the flags work at all (comparison values 
    assert ((flag & sep.APER_TRUNC) != 0).sum() == 4
    assert ((flag & sep.APER_HASMASKED) != 0).sum() == 0

    # Test "flux_auto"
    kr, flag = sep.kron_radius(data, objs['x'], objs['y'], objs['a'],
                               objs['b'], objs['theta'], 6.0)

    flux, fluxerr, flag = sep.sum_ellipse(data, objs['x'], objs['y'],
                                          objs['a'], objs['b'],
                                          objs['theta'], r=2.5 * kr,
                                          err=bkg.globalrms, subpix=1)

    # For some reason, object at index 59 doesn't match. It's very small
    # and kron_radius is set to 0.0 in SExtractor, but 0.08 in sep.
    # Most of the other values are within 1e-4 except one which is only
    # within 0.01. This might be due to a change in SExtractor between
    # v2.8.6 (used to generate "truth" catalog) and v2.18.11.
    kr[59] = 0.0
    flux[59] = 0.0
    fluxerr[59] = 0.0
    assert_allclose(2.5*kr, refobjs['kron_radius'], rtol=0.01)
    assert_allclose(flux, refobjs['flux_auto'], rtol=0.01)
    assert_allclose(fluxerr, refobjs['fluxerr_auto'], rtol=0.01)

    # Test ellipse representation conversion
    cxx, cyy, cxy = sep.ellipse_coeffs(objs['a'], objs['b'], objs['theta'])
    assert_allclose(cxx, objs['cxx'], rtol=1.e-4)
    assert_allclose(cyy, objs['cyy'], rtol=1.e-4)
    assert_allclose(cxy, objs['cxy'], rtol=1.e-4)

    a, b, theta = sep.ellipse_axes(objs['cxx'], objs['cyy'], objs['cxy'])
    assert_allclose(a, objs['a'], rtol=1.e-4)
    assert_allclose(b, objs['b'], rtol=1.e-4)
    assert_allclose(theta, objs['theta'], rtol=1.e-4)

    #test round trip
    cxx, cyy, cxy = sep.ellipse_coeffs(a, b, theta)
    assert_allclose(cxx, objs['cxx'], rtol=1.e-4)
    assert_allclose(cyy, objs['cyy'], rtol=1.e-4)
    assert_allclose(cxy, objs['cxy'], rtol=1.e-4)


# -----------------------------------------------------------------------------
# Background

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


# -----------------------------------------------------------------------------
# Extract

@pytest.mark.skipif(NO_FITS, reason="no FITS reader") 
def test_extract_with_noise_array():

    # Get some background-subtracted test data:
    data = np.copy(image_data)
    bkg = sep.Background(data, bw=64, bh=64, fw=3, fh=3)
    bkg.subfrom(data)

    # Ensure that extraction with constant noise array gives the expected
    # result. We have to use conv=None here because the results are *not*
    # the same when convolution is on! This is because the noise map is
    # convolved. Near edges, the convolution doesn't adjust for pixels
    # off edge boundaries. As a result, the convolved noise map is not
    # all ones.
    objects = sep.extract(data, 1.5*bkg.globalrms, conv=None)
    objects2 = sep.extract(data, 1.5*bkg.globalrms, err=np.ones_like(data),
                           conv=None)
    assert_equal(objects, objects2)

    # Less trivial test where thresh is realistic. Still a flat noise map.
    noise = bkg.globalrms * np.ones_like(data)
    objects2 = sep.extract(data, 1.5, err=noise, conv=None)
    assert_equal(objects, objects2)


# -----------------------------------------------------------------------------
# aperture tests

naper = 1000
x = np.random.uniform(200., 800., naper)
y = np.random.uniform(200., 800., naper)
data_shape = (1000, 1000)

def test_aperture_dtypes():
    """Ensure that all supported image dtypes work in sum_circle() and
    give the same answer"""

    r = 3.

    fluxes = []
    for dt in SUPPORTED_IMAGE_DTYPES:
        data = np.ones(data_shape, dtype=dt)
        flux, fluxerr, flag = sep.sum_circle(data, x, y, r)
        fluxes.append(flux)

    for i in range(1, len(fluxes)):
        assert_allclose(fluxes[0], fluxes[i])


def test_apertures_small_ellipse_exact():
    """Regression test for a bug that manifested primarily when x == y."""

    data = np.ones(data_shape)
    r = 0.3
    rtol=1.e-10
    flux, fluxerr, flag = sep.sum_ellipse(data, x, x, r, r, 0., subpix=0)
    assert_allclose(flux, np.pi*r**2, rtol=rtol)

def test_apertures_all():
    """Test that aperture subpixel sampling works"""

    data = np.random.rand(*data_shape)
    r = 3.
    rtol=1.e-8

    for subpix in [0, 1, 5]:
        flux_ref, fluxerr_ref, flag_ref = sep.sum_circle(data, x, y, r,
                                                         subpix=subpix)

        flux, fluxerr, flag = sep.sum_circann(data, x, y, 0., r,
                                              subpix=subpix)
        assert_allclose(flux, flux_ref, rtol=rtol)

        flux, fluxerr, flag = sep.sum_ellipse(data, x, y, r, r, 0.,
                                              subpix=subpix)
        assert_allclose(flux, flux_ref, rtol=rtol)

        flux, fluxerr, flag = sep.sum_ellipse(data, x, y, 1., 1., 0., r=r,
                                              subpix=subpix)
        assert_allclose(flux, flux_ref, rtol=rtol)



def test_apertures_exact():
    """Test area as measured by exact aperture modes on array of ones"""

    theta = np.random.uniform(-np.pi/2., np.pi/2., naper)
    ratio = np.random.uniform(0.2, 1.0, naper)
    r = 3.

    for dt in SUPPORTED_IMAGE_DTYPES:
        data = np.ones(data_shape, dtype=dt)
        for r in [0.5, 1., 3.]:
            flux, fluxerr, flag = sep.sum_circle(data, x, y, r, subpix=0)
            assert_allclose(flux, np.pi*r**2)

            rout = r*1.1
            flux, fluxerr, flag = sep.sum_circann(data, x, y, r, rout,
                                                  subpix=0)
            assert_allclose(flux, np.pi*(rout**2 - r**2))

            flux, fluxerr, flag = sep.sum_ellipse(data, x, y, 1., ratio,
                                                  theta, r=r, subpix=0)
            assert_allclose(flux, np.pi*ratio*r**2)

            rout = r*1.1
            flux, fluxerr, flag = sep.sum_ellipann(data, x, y, 1., ratio,
                                                   theta, r, rout, subpix=0)
            assert_allclose(flux, np.pi*ratio*(rout**2 - r**2))


def test_aperture_bkgann_overlapping():
    """Test bkgann functionality in circular & elliptical apertures."""

    # If bkgann overlaps aperture exactly, result should be zero
    # (with subpix=1)
    data = np.random.rand(*data_shape)
    r = 5.
    f, _, _ = sep.sum_circle(data, x, y, r, bkgann=(0., r), subpix=1)
    assert_allclose(f, 0., rtol=0., atol=1.e-13)

    f, _, _ = sep.sum_ellipse(data, x, y, 2., 1., np.pi/4., r=r,
                              bkgann=(0., r), subpix=1)
    assert_allclose(f, 0., rtol=0., atol=1.e-13)
    

def test_aperture_bkgann_ones():
    """Test bkgann functionality with flat data"""

    data = np.ones(data_shape)
    r=5.
    bkgann=(6., 8.)

    # On flat data, result should be zero for any bkgann and subpix
    f, _, _ = sep.sum_circle(data, x, y, r, bkgann=bkgann)
    assert_allclose(f, 0., rtol=0., atol=1.e-13)
   
    f, _, _ = sep.sum_ellipse(data, x, y, 2., 1., np.pi/4., r, bkgann=bkgann)
    assert_allclose(f, 0., rtol=0., atol=1.e-13)


def test_mask_ellipse():
    arr = np.zeros((20, 20), dtype=np.bool)

    # should mask 5 pixels:
    sep.mask_ellipse(arr, 10., 10., 1.0, 1.0, 0.0, r=1.001)
    assert arr.sum() == 5

    # should mask 13 pixels:
    sep.mask_ellipse(arr, 10., 10., 1.0, 1.0, 0.0, r=2.001)
    assert arr.sum() == 13

def test_mask_ellipse_dep():
    """Deprecated version of mask_ellipse"""
    arr = np.zeros((20, 20), dtype=np.bool)

    # should mask 5 pixels:
    sep.mask_ellipse(arr, 10., 10., cxx=1.0, cyy=1.0, cxy=0.0, scale=1.001)
    assert arr.sum() == 5

    # should mask 13 pixels:
    sep.mask_ellipse(arr, 10., 10., cxx=1.0, cyy=1.0, cxy=0.0, scale=2.001)
    assert arr.sum() == 13


# -----------------------------------------------------------------------------
# General behavior and utilities

def test_byte_order_exception():
    """Test that error about byte order is raised with non-native
    byte order input array. This should happen for Background, extract,
    and aperture functions."""

    data = np.ones((100, 100), dtype=np.float64)
    data = data.byteswap(True).newbyteorder()
    with pytest.raises(ValueError) as excinfo:
        bkg = sep.Background(data)
    assert 'byte order' in str(excinfo.value)


def test_set_pixstack():
    """Ensure that setting the pixel stack size works."""
    old = sep.get_extract_pixstack()
    new = old * 2
    sep.set_extract_pixstack(new)
    assert new == sep.get_extract_pixstack()
    sep.set_extract_pixstack(old)
