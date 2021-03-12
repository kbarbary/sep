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
RMSIMAGE_FNAME = os.path.join("data", "rms.fits")
IMAGECAT_FNAME = os.path.join("data", "image.cat")
IMAGECAT_DTYPE = [('number', np.int64),
                  ('x', np.float64),
                  ('y', np.float64),
                  ('xwin', np.float64),
                  ('ywin', np.float64),
                  ('x2', np.float64),
                  ('y2', np.float64),
                  ('xy', np.float64),
                  ('errx2', np.float64),
                  ('erry2', np.float64),
                  ('errxy', np.float64),
                  ('a', np.float64),
                  ('flux_aper', np.float64),
                  ('fluxerr_aper', np.float64),
                  ('kron_radius', np.float64),
                  ('flux_auto', np.float64),
                  ('fluxerr_auto', np.float64),
                  ('flux_radius', np.float64, (3,)),
                  ('flags', np.int64)]
SUPPORTED_IMAGE_DTYPES = [np.float64, np.float32, np.int32]

# If we have a FITS reader, read in the necessary test images
if not NO_FITS:
    image_data = getdata(IMAGE_FNAME)
    image_refback = getdata(BACKIMAGE_FNAME)
    image_refrms = getdata(RMSIMAGE_FNAME)


# -----------------------------------------------------------------------------
# Helpers

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


def matched_filter_snr(data, noise, kernel):
    """Super slow implementation of matched filter SNR for testing.

    At each pixel, value is

        sum(data[i] * kernel[i] / noise[i]^2)
        -------------------------------------
        sqrt(sum(kernel[i]^2 / noise[i]^2))
    """
    ctr = kernel.shape[0] // 2, kernel.shape[1] // 2
    kslice = ((0 - ctr[0], kernel.shape[0] - ctr[0]),  # range in axis 0
              (0 - ctr[1], kernel.shape[1] - ctr[1]))  # range in axis 1
    out = np.empty_like(data)

    for y in range(data.shape[0]):
        jmin = y + kslice[0][0] # min and max indicies to sum over
        jmax = y + kslice[0][1]
        kjmin = 0               # min and max kernel indicies to sum over
        kjmax = kernel.shape[0]

        # if we're over the edge of the image, limit extent
        if jmin < 0:
            offset = -jmin
            jmin += offset
            kjmin += offset
        if jmax > data.shape[0]:
            offset = data.shape[0] - jmax
            jmax += offset
            kjmax += offset

        for x in range(data.shape[1]):
            imin = x + kslice[1][0] # min and max indicies to sum over
            imax = x + kslice[1][1]
            kimin = 0               # min and max kernel indicies to sum over
            kimax = kernel.shape[1]

            # if we're over the edge of the image, limit extent
            if imin < 0:
                offset = -imin
                imin += offset
                kimin += offset
            if imax > data.shape[1]:
                offset = data.shape[1] - imax
                imax += offset
                kimax += offset

            d = data[jmin:jmax, imin:imax]
            n = noise[jmin:jmax, imin:imax]
            w = 1. / n**2
            k = kernel[kjmin:kjmax, kimin:kimax]
            out[y, x] = np.sum(d * k * w) / np.sqrt(np.sum(k**2 * w))

    return out


# -----------------------------------------------------------------------------
# Test versus Source Extractor results

@pytest.mark.skipif(NO_FITS, reason="no FITS reader")
def test_vs_sextractor():
    """Test behavior of sep versus sextractor.

    Note: we turn deblending off for this test. This is because the
    deblending algorithm uses a random number generator. Since the sequence
    of random numbers is not the same between sextractor and sep or between
    different platforms, object member pixels (and even the number of objects)
    can differ when deblending is on.

    Deblending is turned off by setting DEBLEND_MINCONT=1.0 in the sextractor
    configuration file and by setting deblend_cont=1.0 in sep.extract().
    """

    data = np.copy(image_data)  # make an explicit copy so we can 'subfrom'
    bkg = sep.Background(data, bw=64, bh=64, fw=3, fh=3)

    # Test that SExtractor background is same as SEP:
    bkgarr = bkg.back(dtype=np.float32)
    assert_allclose(bkgarr, image_refback, rtol=1.e-5)

        # Test that SExtractor background rms is same as SEP:
    rmsarr = bkg.rms(dtype=np.float32)
    assert_allclose(rmsarr, image_refrms, rtol=1.e-4)

    # Extract objects (use deblend_cont=1.0 to disable deblending).
    bkg.subfrom(data)
    objs = sep.extract(data, 1.5, err=bkg.globalrms, deblend_cont=1.0)
    objs = np.sort(objs, order=['y'])

    # Read SExtractor result
    refobjs = np.loadtxt(IMAGECAT_FNAME, dtype=IMAGECAT_DTYPE)
    refobjs = np.sort(refobjs, order=['y'])

    # Found correct number of sources at the right locations?
    assert_allclose(objs['x'], refobjs['x'] - 1., atol=1.e-3)
    assert_allclose(objs['y'], refobjs['y'] - 1., atol=1.e-3)

    # Correct Variance and Variance Errors?
    assert_allclose(objs['x2'], refobjs['x2'], atol=1.e-4)
    assert_allclose(objs['y2'], refobjs['y2'], atol=1.e-4)
    assert_allclose(objs['xy'], refobjs['xy'], atol=1.e-4)
    assert_allclose(objs['errx2'], refobjs['errx2'], rtol=1.e-4)
    assert_allclose(objs['erry2'], refobjs['erry2'], rtol=1.e-4)
    assert_allclose(objs['errxy'], refobjs['errxy'], rtol=1.e-3)

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

    # For some reason, one object doesn't match. It's very small
    # and kron_radius is set to 0.0 in SExtractor, but 0.08 in sep.
    # Could be due to a change in SExtractor between v2.8.6 (used to
    # generate "truth" catalog) and v2.18.11 (from which sep was forked).
    i = 56  # index is 59 when deblending is on.
    kr[i] = 0.0
    flux[i] = 0.0
    fluxerr[i] = 0.0

    # We use atol for radius because it is reported to nearest 0.01 in
    # reference objects.
    assert_allclose(2.5*kr, refobjs['kron_radius'], atol=0.01, rtol=0.)
    assert_allclose(flux, refobjs['flux_auto'], rtol=0.0005)
    assert_allclose(fluxerr, refobjs['fluxerr_auto'], rtol=0.0005)

    # Test using a mask in kron_radius and sum_ellipse.
    for dtype in [np.bool_, np.int32, np.float32, np.float64]:
        mask = np.zeros_like(data, dtype=dtype)
        kr2, flag = sep.kron_radius(data, objs['x'], objs['y'],
                                    objs['a'], objs['b'], objs['theta'],
                                    6.0, mask=mask)
        kr2[i] = 0.
        assert np.all(kr == kr2)

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

    # test flux_radius
    fr, flags = sep.flux_radius(data, objs['x'], objs['y'], 6.*refobjs['a'],
                                [0.1, 0.5, 0.6], normflux=refobjs['flux_auto'],
                                subpix=5)
    assert_allclose(fr, refobjs["flux_radius"], rtol=0.04, atol=0.01)

    # test winpos
    sig = 2. / 2.35 * fr[:, 1]  # flux_radius = 0.5
    xwin, ywin, flag = sep.winpos(data, objs['x'], objs['y'], sig)
    assert_allclose(xwin, refobjs["xwin"] - 1., rtol=0., atol=0.0015)
    assert_allclose(ywin, refobjs["ywin"] - 1., rtol=0., atol=0.0015)

# -----------------------------------------------------------------------------
# Background

def test_masked_background():
    data = 0.1 * np.ones((6,6))
    data[1,1] = 1.
    data[4,1] = 1.
    data[1,4] = 1.
    data[4,4] = 1.

    mask = np.zeros((6,6), dtype=np.bool_)

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


@pytest.mark.skipif(NO_FITS, reason="no FITS reader")
def test_background_special():
    """Test special methods of Background"""

    bkg = sep.Background(image_data, bw=64, bh=64, fw=3, fh=3)

    # test __array__ method
    assert np.all(np.array(bkg) == bkg.back())

    # test __rsub__ method
    d1 = image_data - bkg

    d2 = np.copy(image_data)
    bkg.subfrom(d2)
    assert np.all(d1 == d2)


def test_background_boxsize():
    """Test that background works when boxsize is same as image"""

    ny, nx = 100, 100
    data = np.ones((ny, nx), dtype=np.float64)
    bkg = sep.Background(data, bh=ny, bw=nx, fh=1, fw=1)
    bkg.back()


def test_background_rms():
    """Test that Background.rms() at least works"""

    ny, nx = 1024, 1024
    data = np.random.randn(ny, nx)
    bkg = sep.Background(data)
    rms = bkg.rms()
    assert rms.dtype == np.float64
    assert rms.shape == (ny, nx)

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
    objects = sep.extract(data, 1.5*bkg.globalrms, filter_kernel=None)
    objects2 = sep.extract(data, 1.5*bkg.globalrms, err=np.ones_like(data),
                           filter_kernel=None)
    names_to_remove = ['errx2', 'erry2', 'errxy']
    names_to_keep = [i for i in objects.dtype.names if i not in names_to_remove]
    objects = objects[names_to_keep]
    objects2 = objects2[names_to_keep]

    assert_equal(objects, objects2)

    # Less trivial test where thresh is realistic. Still a flat noise map.
    noise = bkg.globalrms * np.ones_like(data)
    objects2 = sep.extract(data, 1.5, err=noise, filter_kernel=None)

    names_to_remove = ['errx2', 'erry2', 'errxy']
    names_to_keep = [i for i in objects.dtype.names if i not in names_to_remove]
    objects = objects[names_to_keep]
    objects2 = objects2[names_to_keep]

    assert_equal(objects, objects2)


def test_extract_with_noise_convolution():
    """Test extraction when there is both noise and convolution.

    This will use the matched filter implementation, and will handle bad pixels
    and edge effects gracefully.
    """

    # Start with an empty image where we label the noise as 1 sigma everywhere.
    image = np.zeros((20, 20))
    error = np.ones((20, 20))

    # Add some noise representing bad pixels. We do not want to detect these.
    image[17, 3] = 100.
    error[17, 3] = 100.
    image[10, 0] = 100.
    error[10, 0] = 100.
    image[17, 17] = 100.
    error[17, 17] = 100.

    # Add some real point sources that we should find.
    image[3, 17] = 10.

    image[6, 6] = 2.0
    image[7, 6] = 1.0
    image[5, 6] = 1.0
    image[6, 5] = 1.0
    image[6, 7] = 1.0

    objects = sep.extract(image, 2.0, minarea=1, err=error)
    objects.sort(order=['x', 'y'])

    # Check that we recovered the two correct objects and not the others.
    assert len(objects) == 2

    assert_approx_equal(objects[0]['x'], 6.)
    assert_approx_equal(objects[0]['y'], 6.)

    assert_approx_equal(objects[1]['x'], 17.)
    assert_approx_equal(objects[1]['y'], 3.)


def test_extract_matched_filter_at_edge():
    """Exercise bug where bright star at end of image not detected
    with noise array and matched filter on."""

    data = np.zeros((20, 20))
    err = np.ones_like(data)
    kernel = np.array([[1., 2., 1.],
                       [2., 4., 2.],
                       [1., 2., 1.]])

    data[18:20, 9:12] = kernel[0:2, :]

    objects, pix = sep.extract(data, 2.0, err=err, filter_kernel=kernel,
                               filter_type="matched", segmentation_map=True)
    assert len(objects) == 1
    assert objects["npix"][0] == 6


@pytest.mark.skipif(NO_FITS, reason="no FITS reader")
def test_extract_with_mask():

    # Get some background-subtracted test data:
    data = np.copy(image_data)
    bkg = sep.Background(data, bw=64, bh=64, fw=3, fh=3)
    bkg.subfrom(data)

    # mask half the image
    ylim = data.shape[0] // 2
    mask = np.zeros(data.shape, dtype=np.bool_)
    mask[ylim:,:] = True

    objects = sep.extract(data, 1.5*bkg.globalrms, mask=mask)

    # check that we found some objects and that they are all in the unmasked
    # region.
    assert len(objects) > 0
    assert np.all(objects['y'] < ylim)


@pytest.mark.skipif(NO_FITS, reason="no FITS reader")
def test_extract_segmentation_map():

    # Get some background-subtracted test data:
    data = np.copy(image_data)
    bkg = sep.Background(data, bw=64, bh=64, fw=3, fh=3)
    bkg.subfrom(data)

    objects, segmap = sep.extract(data, 1.5*bkg.globalrms,
                                  segmentation_map=True)

    assert type(segmap) is np.ndarray
    assert segmap.shape == data.shape
    for i in range(len(objects)):
        assert objects["npix"][i] == (segmap == i+1).sum()


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
    f, fe, _ = sep.sum_circle(data, x, y, r, bkgann=bkgann, gain=1.)
    assert_allclose(f, 0., rtol=0., atol=1.e-13)

    # for all ones data and no error array, error should be close to
    # sqrt(Npix_aper + Npix_ann * (Npix_aper**2 / Npix_ann**2))
    aper_area = np.pi * r**2
    bkg_area = np.pi * (bkgann[1]**2 - bkgann[0]**2)
    expected_error = np.sqrt(aper_area + bkg_area * (aper_area/bkg_area)**2)
    assert_allclose(fe, expected_error, rtol=0.1)

    f, _, _ = sep.sum_ellipse(data, x, y, 2., 1., np.pi/4., r, bkgann=bkgann)
    assert_allclose(f, 0., rtol=0., atol=1.e-13)

def test_masked_segmentation_measurements():
    """Test measurements with segmentation masking"""

    NX = 100
    data = np.zeros((NX*2,NX*2))
    yp, xp = np.indices(data.shape)

    ####
    # Make two 2D gaussians that slightly overlap

    # width of the 2D objects
    gsigma = 10.

    # offset between two gaussians in sigmas
    off = 4

    for xy in [[NX,NX], [NX+off*gsigma, NX+off*gsigma]]:
        R = np.sqrt((xp-xy[0])**2+(yp-xy[1])**2)
        g_i = np.exp(-R**2/2/gsigma**2)
        data += g_i

    # Absolute total
    total_exact = g_i.sum()

    # Add some noise
    rms = 0.02
    np.random.seed(1)
    data += np.random.normal(size=data.shape)*rms

    # Run source detection
    objs, segmap = sep.extract(data, thresh=1.2, err=rms, mask=None,
                               segmentation_map=True)

    seg_id = np.arange(1, len(objs)+1, dtype=np.int32)

    # Compute Kron/Auto parameters
    x, y, a, b = objs['x'], objs['y'], objs['a'], objs['b']
    theta = objs['theta']

    kronrad, krflag = sep.kron_radius(data, x, y, a, b, theta, 6.0)

    flux_auto, fluxerr, flag = sep.sum_ellipse(data, x, y, a, b, theta,
                                               2.5*kronrad,
                                               segmap=segmap, seg_id=seg_id,
                                               subpix=1)

    # Test total flux
    assert_allclose(flux_auto, total_exact, rtol=5.e-2)

    # Flux_radius
    for flux_fraction in [0.2, 0.5]:

        # Exact solution
        rhalf_exact = np.sqrt(-np.log(1-flux_fraction)*gsigma**2*2)

        # Masked measurement
        flux_radius, flag = sep.flux_radius(data, x, y, 6.*a, flux_fraction,
                                        seg_id=seg_id, segmap=segmap,
                                        normflux=flux_auto, subpix=5)

        # Test flux fraction
        assert_allclose(flux_radius, rhalf_exact, rtol=5.e-2)

    if False:
        print('test_masked_flux_radius')
        print(total_exact, flux_auto)
        print(rhalf_exact, flux_radius)

def test_mask_ellipse():
    arr = np.zeros((20, 20), dtype=np.bool_)

    # should mask 5 pixels:
    sep.mask_ellipse(arr, 10., 10., 1.0, 1.0, 0.0, r=1.001)
    assert arr.sum() == 5

    # should mask 13 pixels:
    sep.mask_ellipse(arr, 10., 10., 1.0, 1.0, 0.0, r=2.001)
    assert arr.sum() == 13


def test_flux_radius():
    data = np.ones(data_shape)
    fluxfrac = [0.2**2, 0.3**2, 0.7**2, 1.]
    true_r = [2., 3., 7., 10.]
    r, _ = sep.flux_radius(data, x, y, 10.*np.ones_like(x),
                           [0.2**2, 0.3**2, 0.7**2, 1.], subpix=5)
    for i in range(len(fluxfrac)):
        assert_allclose(r[:, i], true_r[i], rtol=0.01)


def test_mask_ellipse_alt():
    """mask_ellipse with cxx, cyy, cxy parameters."""
    arr = np.zeros((20, 20), dtype=np.bool_)

    # should mask 5 pixels:
    sep.mask_ellipse(arr, 10., 10., cxx=1.0, cyy=1.0, cxy=0.0, r=1.001)
    assert arr.sum() == 5

    # should mask 13 pixels:
    sep.mask_ellipse(arr, 10., 10., cxx=1.0, cyy=1.0, cxy=0.0, r=2.001)
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
    assert 'byte order' in excinfo.value.args[0]


def test_set_pixstack():
    """Ensure that setting the pixel stack size works."""
    old = sep.get_extract_pixstack()
    new = old * 2
    sep.set_extract_pixstack(new)
    assert new == sep.get_extract_pixstack()
    sep.set_extract_pixstack(old)


def test_set_sub_object_limit():
    """Ensure that setting the sub-object deblending limit works."""
    old = sep.get_sub_object_limit()
    new = old * 2
    sep.set_sub_object_limit(new)
    assert new == sep.get_sub_object_limit()
    sep.set_sub_object_limit(old)


def test_long_error_msg():
    """Ensure that the error message is created successfully when
    there is an error detail."""

    # set extract pixstack to an insanely small value; this will trigger
    # a detailed error message when running sep.extract()
    old = sep.get_extract_pixstack()
    sep.set_extract_pixstack(5)

    data = np.ones((10, 10), dtype=np.float64)
    with pytest.raises(Exception) as excinfo:
        sep.extract(data, 0.1)
    msg = excinfo.value.args[0]
    assert type(msg) == str  # check that message is the native string type
    assert msg.startswith("internal pixel buffer full: The limit")

    # restore
    sep.set_extract_pixstack(old)
