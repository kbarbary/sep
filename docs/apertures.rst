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

**Masking** 

Apply a mask (same shape as data). Pixels where the mask is True are
"corrected" to the average value within the aperture. ::

   flux, fluxerr, flag = sep.sum_circle(data, objs['x'], objs['y'], 3.0,
                                        mask=mask)

**Local background subtraction**

The `~sep.sum_circle` and `~sep.sum_ellipse` functions have options
for performing local background subtraction in an annulus between 6
and 8 pixel radius::

   flux, fluxerr, flag = sep.sum_circle(data, objs['x'], objs['y'], 3.0,
                                        mask=mask, bkgann=(6., 8.))

Pixels in the background annulus are not subsampled and any masked
pixels in the annulus are completely igored rather than corrected.
The inner and outer radii can also be arrays.

**Equivalent of FLUX_AUTO (e.g., MAG_AUTO) in Source Extractor**

This is a two-step process. First we calculate the Kron radius for each
object, then we perform elliptical aperture photometry within that radius::

   kronrad, krflag = sep.kron_radius(data, objs['x'], objs['y'], objs['a'],
                                     objs['b'], objs['theta'], 6.0)
   flux, fluxerr, flag = sep.sum_ellipse(data, objs['x'], objs['y'],
                                         objs['a'], objs['b'], objs['theta'],
                                         2.5 * kronrad, subpix=1)
   flag |= krflag  # combine flags into 'flag'

Specifically this is the equilvalent of setting ``PHOT_AUTOPARAMS 2.5,
0.0`` in Source Extractor (note the ``2.5`` in the argument to
``sum_ellipse``). The second Source Extractor
parameter is a minimum radius. To replicate Source Extractor behavior
for non-zero values, one would put in logic to use circular aperture
photometry if the kron radius is too small. For example::

   rmin = 1.75  # minimum diameter = 3.5
   use_circle = kronrad * np.sqrt(obj['a'] * obj['b']) < rmin
   cflux, cfluxerr, cflag = sep.sum_circle(data, objs['x'][use_circle],
                                           objs['y'][use_circle], rmin,
                                           subpix=1)
   flux[use_circle] = cflux
   fluxerr[use_circle] = cfluxerr
   flag[use_circle] = cflag

**Mask image regions**

Create a boolean array with elliptical regions set to True::

   mask = np.zeros(data.shape, dtype=np.bool)
   sep.mask_ellipse(mask, objs['x'], objs['y'], obs['a'], objs['b'],
                    objs['theta'], r=3.)
