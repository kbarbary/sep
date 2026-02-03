#!/usr/bin/env py.test

"""Test the python functionality of SEP."""

from __future__ import division, print_function

import os

import numpy as np
import pytest
from numpy.lib import recfunctions as rfn
from numpy.testing import assert_allclose, assert_approx_equal, assert_equal

import sep

# unicode_literals doesn't play well with numpy dtype field names


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
IMAGECAT_DTYPE = [
    ("number", np.int64),
    ("x", np.float64),
    ("y", np.float64),
    ("xwin", np.float64),
    ("ywin", np.float64),
    ("x2", np.float64),
    ("y2", np.float64),
    ("xy", np.float64),
    ("errx2", np.float64),
    ("erry2", np.float64),
    ("errxy", np.float64),
    ("a", np.float64),
    ("flux_aper", np.float64),
    ("fluxerr_aper", np.float64),
    ("kron_radius", np.float64),
    ("flux_auto", np.float64),
    ("fluxerr_auto", np.float64),
    ("flux_radius", np.float64, (3,)),
    ("flags", np.int64),
]
SUPPORTED_IMAGE_DTYPES = [np.float64, np.float32, np.int32]

# If we have a FITS reader, read in the necessary test images
if not NO_FITS:
    image_data = getdata(IMAGE_FNAME)
    image_refback = getdata(BACKIMAGE_FNAME)
    image_refrms = getdata(RMSIMAGE_FNAME)


# -----------------------------------------------------------------------------
# Helpers


def assert_allclose_structured(x, y):
    """
    Assert that two structured arrays are close.

    Compares floats relatively and everything else exactly.

    Parameters
    ----------
    x, y : array-like
        Structured arrays to be compared.
    """
    assert x.dtype == y.dtype
    for name in x.dtype.names:
        if np.issubdtype(x.dtype[name], float):
            assert_allclose(x[name], y[name])
        else:
            assert_equal(x[name], y[name])


def matched_filter_snr(data, noise, kernel):
    r"""
    Super slow implementation of matched filter SNR for testing.

    At each output pixel :math:`i`, the value is:

    .. math::

        \frac{\sum(\text{data}[i] * \text{kernel}[i] / \text{noise}[i]^2)}
            {\sqrt\sum(\text{kernel}[i]^2 / \text{noise}[i]^2)}

    Parameters
    ----------
    data : array-like
        The 2D data to be tested.
    noise : array-like
        The noise corresponding to the input ``data``.
    kernel : array-like
        The kernel used for filtering.

    Returns
    -------
    array-like
        The output SNR array, the same size as ``data``.
    """
    ctr = kernel.shape[0] // 2, kernel.shape[1] // 2
    kslice = (
        (0 - ctr[0], kernel.shape[0] - ctr[0]),  # range in axis 0
        (0 - ctr[1], kernel.shape[1] - ctr[1]),
    )  # range in axis 1
    out = np.empty_like(data)

    for y in range(data.shape[0]):
        jmin = y + kslice[0][0]  # min and max indicies to sum over
        jmax = y + kslice[0][1]
        kjmin = 0  # min and max kernel indicies to sum over
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
            imin = x + kslice[1][0]  # min and max indicies to sum over
            imax = x + kslice[1][1]
            kimin = 0  # min and max kernel indicies to sum over
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
            w = 1.0 / n**2
            k = kernel[kjmin:kjmax, kimin:kimax]
            out[y, x] = np.sum(d * k * w) / np.sqrt(np.sum(k**2 * w))

    return out


# -----------------------------------------------------------------------------
# Test versus Source Extractor results


@pytest.mark.skipif(NO_FITS, reason="no FITS reader")
def test_vs_sextractor():
    """
    Test behavior of sep versus sextractor.

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
    assert_allclose(bkgarr, image_refback, rtol=1.0e-5)

    # Test that SExtractor background rms is same as SEP:
    rmsarr = bkg.rms(dtype=np.float32)
    assert_allclose(rmsarr, image_refrms, rtol=1.0e-4)

    # Extract objects (use deblend_cont=1.0 to disable deblending).
    bkg.subfrom(data)
    objs = sep.extract(data, 1.5, err=bkg.globalrms, deblend_cont=1.0)
    objs = np.sort(objs, order=["y"])

    # Read SExtractor result
    refobjs = np.loadtxt(IMAGECAT_FNAME, dtype=IMAGECAT_DTYPE)
    refobjs = np.sort(refobjs, order=["y"])

    # Found correct number of sources at the right locations?
    assert_allclose(objs["x"], refobjs["x"] - 1.0, atol=1.0e-3)
    assert_allclose(objs["y"], refobjs["y"] - 1.0, atol=1.0e-3)

    # Correct Variance and Variance Errors?
    assert_allclose(objs["x2"], refobjs["x2"], atol=1.0e-4)
    assert_allclose(objs["y2"], refobjs["y2"], atol=1.0e-4)
    assert_allclose(objs["xy"], refobjs["xy"], atol=1.0e-4)
    assert_allclose(objs["errx2"], refobjs["errx2"], rtol=1.0e-4)
    assert_allclose(objs["erry2"], refobjs["erry2"], rtol=1.0e-4)
    assert_allclose(objs["errxy"], refobjs["errxy"], rtol=1.0e-3)

    # Test aperture flux
    flux, fluxerr, flag = sep.sum_circle(
        data, objs["x"], objs["y"], 5.0, err=bkg.globalrms
    )
    assert_allclose(flux, refobjs["flux_aper"], rtol=2.0e-4)
    assert_allclose(fluxerr, refobjs["fluxerr_aper"], rtol=1.0e-5)

    # check if the flags work at all (comparison values
    assert ((flag & sep.APER_TRUNC) != 0).sum() == 4
    assert ((flag & sep.APER_HASMASKED) != 0).sum() == 0

    # Test "flux_auto"
    kr, flag = sep.kron_radius(
        data, objs["x"], objs["y"], objs["a"], objs["b"], objs["theta"], 6.0
    )

    flux, fluxerr, flag = sep.sum_ellipse(
        data,
        objs["x"],
        objs["y"],
        objs["a"],
        objs["b"],
        objs["theta"],
        r=2.5 * kr,
        err=bkg.globalrms,
        subpix=1,
    )

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
    assert_allclose(2.5 * kr, refobjs["kron_radius"], atol=0.01, rtol=0.0)
    assert_allclose(flux, refobjs["flux_auto"], rtol=0.0005)
    assert_allclose(fluxerr, refobjs["fluxerr_auto"], rtol=0.0005)

    # Test using a mask in kron_radius and sum_ellipse.
    for dtype in [np.bool_, np.int32, np.float32, np.float64]:
        mask = np.zeros_like(data, dtype=dtype)
        kr2, flag = sep.kron_radius(
            data,
            objs["x"],
            objs["y"],
            objs["a"],
            objs["b"],
            objs["theta"],
            6.0,
            mask=mask,
        )
        kr2[i] = 0.0
        assert np.all(kr == kr2)

    # Test ellipse representation conversion
    cxx, cyy, cxy = sep.ellipse_coeffs(objs["a"], objs["b"], objs["theta"])
    assert_allclose(cxx, objs["cxx"], rtol=1.0e-4)
    assert_allclose(cyy, objs["cyy"], rtol=1.0e-4)
    assert_allclose(cxy, objs["cxy"], rtol=1.0e-4)

    a, b, theta = sep.ellipse_axes(objs["cxx"], objs["cyy"], objs["cxy"])
    assert_allclose(a, objs["a"], rtol=1.0e-4)
    assert_allclose(b, objs["b"], rtol=1.0e-4)
    assert_allclose(theta, objs["theta"], rtol=1.0e-4)

    # test round trip
    cxx, cyy, cxy = sep.ellipse_coeffs(a, b, theta)
    assert_allclose(cxx, objs["cxx"], rtol=1.0e-4)
    assert_allclose(cyy, objs["cyy"], rtol=1.0e-4)
    assert_allclose(cxy, objs["cxy"], rtol=1.0e-4)

    # test flux_radius
    fr, flags = sep.flux_radius(
        data,
        objs["x"],
        objs["y"],
        6.0 * refobjs["a"],
        [0.1, 0.5, 0.6],
        normflux=refobjs["flux_auto"],
        subpix=5,
    )
    assert_allclose(fr, refobjs["flux_radius"], rtol=0.04, atol=0.01)

    # test winpos
    sig = 2.0 / 2.35 * fr[:, 1]  # flux_radius = 0.5
    xwin, ywin, flag = sep.winpos(data, objs["x"], objs["y"], sig)
    assert_allclose(xwin, refobjs["xwin"] - 1.0, rtol=0.0, atol=0.0015)
    assert_allclose(ywin, refobjs["ywin"] - 1.0, rtol=0.0, atol=0.0015)


# -----------------------------------------------------------------------------
# Background


def test_masked_background():
    """
    Check the background filtering.

    Check that the derived background is consistent with an explicit
    mask, masking no pixels. Also check that the expected result is
    returned if certain pixels are masked.
    """

    data = 0.1 * np.ones((6, 6))
    data[1, 1] = 1.0
    data[4, 1] = 1.0
    data[1, 4] = 1.0
    data[4, 4] = 1.0

    mask = np.zeros((6, 6), dtype=np.bool_)

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
    """
    Test the special methods of `sep.Background`.
    """

    bkg = sep.Background(image_data, bw=64, bh=64, fw=3, fh=3)

    # test __array__ method
    assert np.all(np.array(bkg) == bkg.back())

    # test __rsub__ method
    d1 = image_data - bkg

    d2 = np.copy(image_data)
    bkg.subfrom(d2)
    assert np.all(d1 == d2)


def test_background_boxsize():
    """
    Test that `sep.Background` works when boxsize is same as image.
    """

    ny, nx = 100, 100
    data = np.ones((ny, nx), dtype=np.float64)
    bkg = sep.Background(data, bh=ny, bw=nx, fh=1, fw=1)
    bkg.back()


def test_background_rms():
    """
    Test that `sep.Background.rms` at least works.
    """

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
    """
    Test extraction with a flat noise array.

    This checks that a constant noise array gives the same result as
    extracting without a noise array, for a given threshold.
    """

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
    # Deblending is also turned off, as this appears to differ slightly
    # across platforms - see `test_vs_sextractor()`.
    objects = sep.extract(
        data, 1.5 * bkg.globalrms, filter_kernel=None, deblend_cont=1.0
    )
    objects2 = sep.extract(
        data,
        1.5 * bkg.globalrms,
        err=np.ones_like(data),
        filter_kernel=None,
        deblend_cont=1.0,
    )

    names_to_remove = ["errx2", "erry2", "errxy"]
    names_to_keep = [i for i in objects.dtype.names if i not in names_to_remove]
    objects = objects[names_to_keep]
    objects2 = objects2[names_to_keep]

    assert_allclose_structured(objects, objects2)

    # Less trivial test where thresh is realistic. Still a flat noise map.
    noise = bkg.globalrms * np.ones_like(data)
    objects2 = sep.extract(data, 1.5, err=noise, filter_kernel=None, deblend_cont=1.0)

    names_to_remove = ["errx2", "erry2", "errxy"]
    names_to_keep = [i for i in objects.dtype.names if i not in names_to_remove]
    objects = objects[names_to_keep]
    objects2 = objects2[names_to_keep]

    assert_allclose_structured(objects, objects2)


def test_extract_with_noise_convolution():
    """
    Test extraction when there is both noise and convolution.

    This will use the matched filter implementation, and will handle bad pixels
    and edge effects gracefully.
    """

    # Start with an empty image where we label the noise as 1 sigma everywhere.
    image = np.zeros((20, 20))
    error = np.ones((20, 20))

    # Add some noise representing bad pixels. We do not want to detect these.
    image[17, 3] = 100.0
    error[17, 3] = 100.0
    image[10, 0] = 100.0
    error[10, 0] = 100.0
    image[17, 17] = 100.0
    error[17, 17] = 100.0

    # Add some real point sources that we should find.
    image[3, 17] = 10.0

    image[6, 6] = 2.0
    image[7, 6] = 1.0
    image[5, 6] = 1.0
    image[6, 5] = 1.0
    image[6, 7] = 1.0

    objects = sep.extract(image, 2.0, minarea=1, err=error)
    objects.sort(order=["x", "y"])

    # Check that we recovered the two correct objects and not the others.
    assert len(objects) == 2

    assert_approx_equal(objects[0]["x"], 6.0)
    assert_approx_equal(objects[0]["y"], 6.0)

    assert_approx_equal(objects[1]["x"], 17.0)
    assert_approx_equal(objects[1]["y"], 3.0)


def test_extract_matched_filter_at_edge():
    """
    Test bright source detection at the edge of an image.

    Exercise bug where bright star at end of image not detected
    with noise array and matched filter on.
    """

    data = np.zeros((20, 20))
    err = np.ones_like(data)
    kernel = np.array([[1.0, 2.0, 1.0], [2.0, 4.0, 2.0], [1.0, 2.0, 1.0]])

    data[18:20, 9:12] = kernel[0:2, :]

    objects, pix = sep.extract(
        data,
        2.0,
        err=err,
        filter_kernel=kernel,
        filter_type="matched",
        segmentation_map=True,
    )
    assert len(objects) == 1
    assert objects["npix"][0] == 6


@pytest.mark.skipif(NO_FITS, reason="no FITS reader")
def test_extract_with_mask():
    """
    Test that object detection only occurs in unmasked regions.
    """

    # Get some background-subtracted test data:
    data = np.copy(image_data)
    bkg = sep.Background(data, bw=64, bh=64, fw=3, fh=3)
    bkg.subfrom(data)

    # mask half the image
    ylim = data.shape[0] // 2
    mask = np.zeros(data.shape, dtype=np.bool_)
    mask[ylim:, :] = True

    objects = sep.extract(data, 1.5 * bkg.globalrms, mask=mask)

    # check that we found some objects and that they are all in the unmasked
    # region.
    assert len(objects) > 0
    assert np.all(objects["y"] < ylim)


@pytest.mark.skipif(NO_FITS, reason="no FITS reader")
def test_extract_with_maskthresh():
    """
    Test that object detection only occurs in unmasked regions.
    """

    # Get some background-subtracted test data:
    data = np.copy(image_data)
    bkg = sep.Background(data, bw=64, bh=64, fw=3, fh=3)
    bkg.subfrom(data)

    # mask half the image
    ylim = data.shape[0] // 2
    mask = np.zeros(data.shape, dtype=float)
    mask[ylim:, :] = 1.0

    objects_unmasked = sep.extract(data, 1.5 * bkg.globalrms, deblend_cont=1.0)
    objects_unmasked_w_thresh = sep.extract(
        data, 1.5 * bkg.globalrms, maskthresh=1.0, deblend_cont=1.0
    )

    # Check that changing the mask threshold does not change anything,
    # if no mask is provided
    assert_allclose_structured(objects_unmasked, objects_unmasked_w_thresh)

    objects_masked = sep.extract(data, 1.5 * bkg.globalrms, mask=mask, deblend_cont=1.0)
    objects_masked_w_hthresh = sep.extract(
        data, 1.5 * bkg.globalrms, mask=mask, maskthresh=1.0, deblend_cont=1.0
    )
    objects_masked_w_lthresh = sep.extract(
        data, 1.5 * bkg.globalrms, mask=mask, maskthresh=0.5, deblend_cont=1.0
    )

    # Applying a mask should return a different number of objects
    assert len(objects_unmasked) != len(objects_masked)

    # As long as the mask is above the threshold, the results should not change
    assert_allclose_structured(objects_masked, objects_masked_w_lthresh)

    # Object detection where the maskthresh >= max(mask) should be the same
    # as if no mask were provided
    assert_allclose_structured(objects_unmasked, objects_masked_w_hthresh)


@pytest.mark.skipif(NO_FITS, reason="no FITS reader")
def test_extract_segmentation_map():
    """
    Test the returned segmentation map.

    Check that the segmentation map has the same dimensions as the input
    image, and that the number of object pixels match the catalogue field.
    """

    # Get some background-subtracted test data:
    data = np.copy(image_data)
    bkg = sep.Background(data, bw=64, bh=64, fw=3, fh=3)
    bkg.subfrom(data)

    objects, segmap = sep.extract(data, 1.5 * bkg.globalrms, segmentation_map=True)

    assert type(segmap) is np.ndarray
    assert segmap.shape == data.shape
    for i in range(len(objects)):
        assert objects["npix"][i] == (segmap == i + 1).sum()


@pytest.mark.skipif(NO_FITS, reason="no FITS reader")
def test_extract_seg_map_array():
    """
    Test the extraction when an existing segmentation map is supplied.

    Test that the returned catalogue is equal with and without a variable
    noise array, and that the majority of fields match even when
    deblending is performed on the original extraction.
    """

    # Get some background-subtracted test data:
    data = np.copy(image_data)
    bkg = sep.Background(data, bw=64, bh=64, fw=3, fh=3)
    bkg.subfrom(data)

    noise = bkg.globalrms * np.ones_like(data)

    for err in [None, noise]:
        # err=None
        # err=noise

        objects, segmap = sep.extract(data, 1.5, err, segmentation_map=True)

        assert type(segmap) is np.ndarray
        assert segmap.shape == data.shape
        for i in range(len(objects)):
            assert objects["npix"][i] == (segmap == i + 1).sum()

        objects2, segmap2 = sep.extract(data, 1.5, err, segmentation_map=segmap)

        # Test the values for which we expect an exact match
        names_exact_match = [
            "thresh",
            "npix",
            "tnpix",
            "xmin",
            "xmax",
            "ymin",
            "ymax",
            "cflux",
            "flux",
            "cpeak",
            "peak",
            "xcpeak",
            "ycpeak",
            "xpeak",
            "ypeak",
        ]

        # The position depends on the object being deblended. As no deblending
        # is performed when a segmentation map is supplied, all derived
        # parameters may vary slightly. We test those for which we have a
        # measurement of the uncertainty
        names_close = ["x", "y"]
        names_close_var = ["x2", "y2"]

        assert segmap2.shape == data.shape
        for o_i, o_ii in zip(objects, objects2):
            o_i_exact = o_i[names_exact_match]
            o_ii_exact = o_ii[names_exact_match]
            assert_equal(o_i_exact, o_ii_exact)

            o_i_close = o_i[names_close]
            o_ii_close = o_ii[names_close]
            for n, v in zip(names_close, names_close_var):
                if o_i["flag"] == 0:
                    assert_equal(o_i[n], o_ii[n])
                else:
                    assert_allclose(o_i[n], o_ii[n], atol=np.sqrt(o_i[v]))

        # Perform a second test with deblending disabled.
        objects3, segmap3 = sep.extract(
            data, 1.5, err, segmentation_map=True, deblend_cont=1.0
        )

        objects4, segmap4 = sep.extract(
            data, 1.5, err, segmentation_map=segmap3, deblend_cont=1.0
        )

        # The flag will not be the same, as the second extraction does not test
        # for deblended objects.
        objects3 = rfn.drop_fields(objects3, "flag")
        objects4 = rfn.drop_fields(objects4, "flag")
        assert_allclose_structured(objects3, objects4)


# -----------------------------------------------------------------------------
# aperture tests

naper = 1000
x = np.random.uniform(200.0, 800.0, naper)
y = np.random.uniform(200.0, 800.0, naper)
data_shape = (1000, 1000)


def test_aperture_dtypes():
    """
    Test the aperture extraction of multiple data types.

    Ensure that all supported image dtypes work in sum_circle() and
    give the same answer.
    """

    r = 3.0

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
    rtol = 1.0e-10
    flux, fluxerr, flag = sep.sum_ellipse(data, x, x, r, r, 0.0, subpix=0)
    assert_allclose(flux, np.pi * r**2, rtol=rtol)


def test_apertures_all():
    """
    Test that aperture subpixel sampling works.
    """

    data = np.random.rand(*data_shape)
    r = 3.0
    rtol = 1.0e-8

    for subpix in [0, 1, 5]:
        flux_ref, fluxerr_ref, flag_ref = sep.sum_circle(data, x, y, r, subpix=subpix)

        flux, fluxerr, flag = sep.sum_circann(data, x, y, 0.0, r, subpix=subpix)
        assert_allclose(flux, flux_ref, rtol=rtol)

        flux, fluxerr, flag = sep.sum_ellipse(data, x, y, r, r, 0.0, subpix=subpix)
        assert_allclose(flux, flux_ref, rtol=rtol)

        flux, fluxerr, flag = sep.sum_ellipse(
            data, x, y, 1.0, 1.0, 0.0, r=r, subpix=subpix
        )
        assert_allclose(flux, flux_ref, rtol=rtol)


def test_apertures_exact():
    """
    Test area as measured by exact aperture modes on array of ones.
    """

    theta = np.random.uniform(-np.pi / 2.0, np.pi / 2.0, naper)
    ratio = np.random.uniform(0.2, 1.0, naper)
    r = 3.0

    for dt in SUPPORTED_IMAGE_DTYPES:
        data = np.ones(data_shape, dtype=dt)
        for r in [0.5, 1.0, 3.0]:
            flux, fluxerr, flag = sep.sum_circle(data, x, y, r, subpix=0)
            assert_allclose(flux, np.pi * r**2)

            rout = r * 1.1
            flux, fluxerr, flag = sep.sum_circann(data, x, y, r, rout, subpix=0)
            assert_allclose(flux, np.pi * (rout**2 - r**2))

            flux, fluxerr, flag = sep.sum_ellipse(
                data, x, y, 1.0, ratio, theta, r=r, subpix=0
            )
            assert_allclose(flux, np.pi * ratio * r**2)

            rout = r * 1.1
            flux, fluxerr, flag = sep.sum_ellipann(
                data, x, y, 1.0, ratio, theta, r, rout, subpix=0
            )
            assert_allclose(flux, np.pi * ratio * (rout**2 - r**2))


def _sigma_clip_mean(values, sigma=3.0, maxiters=5):
    mask = np.ones(values.shape, dtype=bool)
    for _ in range(maxiters):
        vals = values[mask]
        if vals.size == 0:
            return np.nan
        med = np.median(vals)
        mad = np.median(np.abs(vals - med))
        std = 1.4826 * mad
        if std == 0.0:
            return vals.mean()
        lo = med - sigma * std
        hi = med + sigma * std
        newmask = mask & (values >= lo) & (values <= hi)
        if newmask.sum() == mask.sum():
            return vals.mean()
        mask = newmask
    vals = values[mask]
    return vals.mean() if vals.size else np.nan


def test_stats_circann_flat():
    data = np.ones(data_shape)
    rin, rout = 3.0, 6.0
    mean, std, med, mad_std, mean_clip, flag = sep.stats_circann(
        data, x, y, rin, rout, subpix=1
    )
    assert_allclose(mean, 1.0)
    assert_allclose(std, 0.0)
    assert_allclose(med, 1.0)
    assert_allclose(mad_std, 0.0)
    assert_allclose(mean_clip, 1.0)


def test_stats_circann_matches_numpy_subpix1():
    rng = np.random.default_rng(12345)
    data = rng.normal(size=(64, 64))
    x0 = np.array([30.4])
    y0 = np.array([31.2])
    rin, rout = 5.0, 8.0

    mean, std, med, mad_std, mean_clip, flag = sep.stats_circann(
        data, x0, y0, rin, rout, subpix=1
    )

    yy, xx = np.indices(data.shape)
    rpix2 = (xx - x0[0]) ** 2 + (yy - y0[0]) ** 2
    mask = (rpix2 >= rin**2) & (rpix2 < rout**2)
    vals = data[mask]

    exp_mean = vals.mean()
    exp_std = vals.std()
    exp_med = np.median(vals)
    exp_mad = np.median(np.abs(vals - exp_med))
    exp_mad_std = 1.4826 * exp_mad  # robust std estimate for Gaussian data
    exp_clip = _sigma_clip_mean(vals)

    assert_allclose(mean[0], exp_mean, rtol=1.0e-10, atol=1.0e-8)
    assert_allclose(std[0], exp_std, rtol=1.0e-10, atol=1.0e-8)
    assert_allclose(med[0], exp_med, rtol=1.0e-10, atol=1.0e-8)
    assert_allclose(mad_std[0], exp_mad_std, rtol=1.0e-10, atol=1.0e-7)
    assert_allclose(mean_clip[0], exp_clip, rtol=1.0e-10, atol=1.0e-8)
    assert flag[0] == 0


def test_aperture_bkgann_overlapping():
    """
    Test bkgann functionality in circular & elliptical apertures.
    """

    # If bkgann overlaps aperture exactly, result should be zero
    # (with subpix=1)
    data = np.random.rand(*data_shape)
    r = 5.0
    f, _, _ = sep.sum_circle(data, x, y, r, bkgann=(0.0, r), subpix=1)
    assert_allclose(f, 0.0, rtol=0.0, atol=1.0e-13)

    f, _, _ = sep.sum_ellipse(
        data, x, y, 2.0, 1.0, np.pi / 4.0, r=r, bkgann=(0.0, r), subpix=1
    )
    assert_allclose(f, 0.0, rtol=0.0, atol=1.0e-13)


def test_aperture_bkgann_ones():
    """
    Test bkgann functionality with flat data.
    """

    data = np.ones(data_shape)
    r = 5.0
    bkgann = (6.0, 8.0)

    # On flat data, result should be zero for any bkgann and subpix
    f, fe, _ = sep.sum_circle(data, x, y, r, bkgann=bkgann, gain=1.0)
    assert_allclose(f, 0.0, rtol=0.0, atol=1.0e-13)

    # for all ones data and no error array, error should be close to
    # sqrt(Npix_aper + Npix_ann * (Npix_aper**2 / Npix_ann**2))
    aper_area = np.pi * r**2
    bkg_area = np.pi * (bkgann[1] ** 2 - bkgann[0] ** 2)
    expected_error = np.sqrt(aper_area + bkg_area * (aper_area / bkg_area) ** 2)
    assert_allclose(fe, expected_error, rtol=0.1)

    f, _, _ = sep.sum_ellipse(data, x, y, 2.0, 1.0, np.pi / 4.0, r, bkgann=bkgann)
    assert_allclose(f, 0.0, rtol=0.0, atol=1.0e-13)


def test_masked_segmentation_measurements():
    """
    Test measurements with segmentation masking.
    """

    NX = 100
    data = np.zeros((NX * 2, NX * 2))
    yp, xp = np.indices(data.shape)

    ####
    # Make two 2D gaussians that slightly overlap

    # width of the 2D objects
    gsigma = 10.0

    # offset between two gaussians in sigmas
    off = 4

    for xy in [[NX, NX], [NX + off * gsigma, NX + off * gsigma]]:
        R = np.sqrt((xp - xy[0]) ** 2 + (yp - xy[1]) ** 2)
        g_i = np.exp(-(R**2) / 2 / gsigma**2)
        data += g_i

    # Absolute total
    total_exact = g_i.sum()

    # Add some noise
    rms = 0.02
    np.random.seed(1)
    data += np.random.normal(size=data.shape) * rms

    # Run source detection
    objs, segmap = sep.extract(
        data, thresh=1.2, err=rms, mask=None, segmentation_map=True
    )

    seg_id = np.arange(1, len(objs) + 1, dtype=np.int32)

    # Compute Kron/Auto parameters
    x, y, a, b = objs["x"], objs["y"], objs["a"], objs["b"]
    theta = objs["theta"]

    kronrad, krflag = sep.kron_radius(data, x, y, a, b, theta, 6.0)

    flux_auto, fluxerr, flag = sep.sum_ellipse(
        data, x, y, a, b, theta, 2.5 * kronrad, segmap=segmap, seg_id=seg_id, subpix=1
    )

    # Test total flux
    assert_allclose(flux_auto, total_exact, rtol=5.0e-2)

    # Flux_radius
    for flux_fraction in [0.2, 0.5]:

        # Exact solution
        rhalf_exact = np.sqrt(-np.log(1 - flux_fraction) * gsigma**2 * 2)

        # Masked measurement
        flux_radius, flag = sep.flux_radius(
            data,
            x,
            y,
            6.0 * a,
            flux_fraction,
            seg_id=seg_id,
            segmap=segmap,
            normflux=flux_auto,
            subpix=5,
        )

        # Test flux fraction
        assert_allclose(flux_radius, rhalf_exact, rtol=5.0e-2)

    if False:
        print("test_masked_flux_radius")
        print(total_exact, flux_auto)
        print(rhalf_exact, flux_radius)


def test_mask_ellipse():
    """
    Test that the correct number of elements are masked with an ellipse.
    """
    arr = np.zeros((20, 20), dtype=np.bool_)

    # should mask 5 pixels:
    sep.mask_ellipse(arr, 10.0, 10.0, 1.0, 1.0, 0.0, r=1.001)
    assert arr.sum() == 5

    # should mask 13 pixels:
    sep.mask_ellipse(arr, 10.0, 10.0, 1.0, 1.0, 0.0, r=2.001)
    assert arr.sum() == 13


def test_flux_radius():
    """
    Test that the correct radius is returned for varying flux fractions.
    """
    data = np.ones(data_shape)
    fluxfrac = [0.2**2, 0.3**2, 0.7**2, 1.0]
    true_r = [2.0, 3.0, 7.0, 10.0]
    r, _ = sep.flux_radius(
        data, x, y, 10.0 * np.ones_like(x), [0.2**2, 0.3**2, 0.7**2, 1.0], subpix=5
    )
    for i in range(len(fluxfrac)):
        assert_allclose(r[:, i], true_r[i], rtol=0.01)


def test_mask_ellipse_alt():
    """
    Mask_ellipse with cxx, cyy, cxy parameters.
    """
    arr = np.zeros((20, 20), dtype=np.bool_)

    # should mask 5 pixels:
    sep.mask_ellipse(arr, 10.0, 10.0, cxx=1.0, cyy=1.0, cxy=0.0, r=1.001)
    assert arr.sum() == 5

    # should mask 13 pixels:
    sep.mask_ellipse(arr, 10.0, 10.0, cxx=1.0, cyy=1.0, cxy=0.0, r=2.001)
    assert arr.sum() == 13


# -----------------------------------------------------------------------------
# General behavior and utilities


def test_byte_order_exception():
    """
    Test that SEP will not run with non-native byte order.

    Test that error about byte order is raised with non-native
    byte order input array. This should happen for Background, extract,
    and aperture functions.
    """

    data = np.ones((100, 100), dtype=np.float64)
    data = data.view(data.dtype.newbyteorder("S"))
    with pytest.raises(ValueError) as excinfo:
        bkg = sep.Background(data)
    assert "byte order" in excinfo.value.args[0]


def test_set_pixstack():
    """
    Ensure that setting the pixel stack size works.
    """
    old = sep.get_extract_pixstack()
    new = old * 2
    sep.set_extract_pixstack(new)
    assert new == sep.get_extract_pixstack()
    sep.set_extract_pixstack(old)


def test_set_sub_object_limit():
    """
    Ensure that setting the sub-object deblending limit works.
    """
    old = sep.get_sub_object_limit()
    new = old * 2
    sep.set_sub_object_limit(new)
    assert new == sep.get_sub_object_limit()
    sep.set_sub_object_limit(old)


def test_long_error_msg():
    """
    Test the error handling in SEP.

    Ensure that the error message is created successfully when
    there is an error detail.
    """

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
