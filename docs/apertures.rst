Aperture photometry
===================

There are four aperture functions available:

==================  =========================
Function            Sums data within...
==================  =========================
`sep.sum_circle`    circle(s)
`sep.sum_circann`   circular annulus/annuli
`sep.sum_ellipse`   ellipse(s)
`sep.sum_ellipann`  elliptical annulus/annuli
==================  =========================

The follow examples demonstrate options for circular aperture
photometry. The other functions behave similarly. ::

   # sum flux in circles of radius=3.0
   flux, fluxerr, flag = sep.sum_circle(data, objs['x'], objs['y'], 3.0)

   # x, y and r can be arrays and obey numpy broadcasting rules.
   # Here, r is an array instead of a single number:
   flux, fluxerr, flag = sep.sum_circle(data, objs['x'], objs['y'],
                                        3.0 * np.ones(len(objs)))

   # use a different subpixel sampling (default is 5; 0 means "exact")
   flux, fluxerr, flag = sep.sum_circle(data, objs['x'], objs['y'], 3.0,
                                        subpix=0)

**Error calculation**

In the default modes illustrated above, the uncertainty ``fluxerr`` is
not calculated and values of 0 are simply returned. The uncertainty can be
flexibly and efficiently calculated, depending on the characteristics
of your data. The presence of an ``err`` or ``var`` keyword indicates
a per-pixel noise, while the presense of a ``gain`` keyword indicates
that the Poisson uncertainty on the total sum should be included. Some
illustrative examples::

   # Specify a per-pixel "background" error and a gain. This is suitable
   # when the data have been background subtracted. 
   flux, fluxerr, flag = sep.sum_circle(data, objs['x'], objs['y'], 3.0,
                                        err=bkg.globalrms, gain=1.0)

   # Variance can be passed instead of error, with identical results.
   flux, fluxerr, flag = sep.sum_circle(data, objs['x'], objs['y'], 3.0,
                                        var=bkg.globalrms**2, gain=1.0)

   # Error or variance can be arrays, indicating that the background error
   # is variable.
   bkgrms = bkg.rms()  # array, same shape as data
   flux, fluxerr, flag = sep.sum_circle(data, objs['x'], objs['y'], 3.0,
                                        err=bkgrms, gain=1.0)

   # If your uncertainty array already includes Poisson noise from the object,
   # leave gain as None (default):
   flux, fluxerr, flag = sep.sum_circle(data, objs['x'], objs['y'], 3.0,
                                        err=error_array)

   # If your data represent raw counts (it is not background-subtracted),
   # set only gain to get the poisson error:
   flux, fluxerr, flag = sep.sum_circle(data, objs['x'], objs['y'], 3.0,
                                        gain=1.0)

The error is calculated as

.. math::

   \sigma_F^2 = \sum_i \sigma_i^2 + F/g

where the sum is over pixels in the aperture, :math:`\sigma_i` is the
noise in each pixel, :math:`F` is the sum in the aperture and
:math:`g` is the gain. The last term is not added if ``gain`` is
`None`.

**Masking** 

Apply a mask (same shape as data). Pixels where the mask is True are
"corrected" to the average value within the aperture. ::

   flux, fluxerr, flag = sep.sum_circle(data, objs['x'], objs['y'], 3.0,
                                        mask=mask)

**Local background subtraction**

The `~sep.sum_circle` and `~sep.sum_ellipse` functions have options
for performing local background subtraction. For example, to subtract the background calculated in an annulus between 6 and 8 pixel radius::

   flux, fluxerr, flag = sep.sum_circle(data, objs['x'], objs['y'], 3.0,
                                        mask=mask, bkgann=(6., 8.))

Pixels in the background annulus are not subsampled and any masked
pixels in the annulus are completely igored rather than corrected.
The inner and outer radii can also be arrays. The error in the background
is included in the reported error.

Equivalent of FLUX_AUTO (e.g., MAG_AUTO) in Source Extractor
------------------------------------------------------------

This is a two-step process. First we calculate the Kron radius for each
object, then we perform elliptical aperture photometry within that radius::

   kronrad, krflag = sep.kron_radius(data, x, y, a, b, theta, 6.0)
   flux, fluxerr, flag = sep.sum_ellipse(data, x, y, a, b, theta, 2.5*kronrad,
                                         subpix=1)
   flag |= krflag  # combine flags into 'flag'

This specific example is the equilvalent of setting ``PHOT_AUTOPARAMS
2.5, 0.0`` in Source Extractor (note the 2.5 in the argument to
``sep.sum_ellipse``). The second Source Extractor parameter is a
minimum diameter. To replicate Source Extractor behavior for non-zero
values of the minimum diameter, one would put in logic to use circular
aperture photometry if the Kron radius is too small. For example::

   r_min = 1.75  # minimum diameter = 3.5
   use_circle = kronrad * np.sqrt(a * b) < r_min
   cflux, cfluxerr, cflag = sep.sum_circle(data, x[use_circle], y[use_circle],
                                           r_min, subpix=1)
   flux[use_circle] = cflux
   fluxerr[use_circle] = cfluxerr
   flag[use_circle] = cflag

Equivalent of FLUX_RADIUS in Source Extractor
---------------------------------------------

In Source Extractor, the FLUX_RADIUS parameter gives the radius of a
circle enclosing a desired fraction of the total flux. For example,
with the setting ``PHOT_FLUXFRAC 0.5``, FLUX_RADIUS will give the
radius of a circle containing half the "total flux" of the object. For
the definition of "total flux", Source Extractor uses its measurement
of FLUX_AUTO, which is taken through an elliptical aperture (see
above). Thus, with the setting ``PHOT_FLUXFRAC 1.0``, you would find
the circle containing the same flux as whatever ellipse Source
Extractor used for ``FLUX_AUTO``.

Given a previous calculation of ``flux`` as above, calculate the
radius for a flux fraction of 0.5::

    r, flag = sep.flux_radius(data, objs['x'], objs['y'], 6.*objs['a'], 0.5,
                              normflux=flux, subpix=5)

And for multiple flux fractions::

    r, flag = sep.flux_radius(data, objs['x'], objs['y'], 6.*objs['a'],
                              [0.5, 0.6], normflux=flux, subpix=5)


Equivalent of XWIN_IMAGE, YWIN_IMAGE in Source Extractor
--------------------------------------------------------

Source Extractor's XWIN_IMAGE, YWIN_IMAGE parameters can be used for
more accurate object centroids than the default X_IMAGE, Y_IMAGE.
Here, the ``winpos`` function provides this behavior.  To match Source
Extractor exactly, the right ``sig`` parameter (giving a description
of the effective width) must be used for each object.  Source
Extractor uses ``2.  / 2.35 * (half-light radius)`` where the
half-light radius is calculated using ``flux_radius`` with a fraction
of 0.5 and a normalizing flux of ``FLUX_AUTO``. The equivalent here is::

    sig = 2. / 2.35 * r  # r from sep.flux_radius() above, with fluxfrac = 0.5
    xwin, ywin, flag = sep.winpos(data, objs['x'], objs['y'], sig)


Masking image regions
---------------------

Create a boolean array with elliptical regions set to True::

   mask = np.zeros(data.shape, dtype=np.bool)
   sep.mask_ellipse(mask, objs['x'], objs['y'], obs['a'], objs['b'],
                    objs['theta'], r=3.)
