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
photometry. The other functions behave similarly::

   # sum flux in circles of radius=3.0
   flux, fluxerr, flag = sep.sum_circle(data, objects['x'], objects['y'], 3.0)

   # x, y and r can be arrays and obey numpy broadcasting rules
   flux, fluxerr, flag = sep.sum_circle(data, objects['x'], objects['y'],
                                        3.0 * np.ones(len(objects)))

   # use a different subpixel sampling (default is 5; 0 means "exact")
   flux, fluxerr, flag = sep.sum_circle(data, objects['x'], objects['y'], 3.0,
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
   flux, fluxerr, flag = sep.sum_circle(data, objects['x'], objects['y'], 3.0,
                                        err=bkg.globalrms, gain=1.0)

   # Variance can be passed instead of error, with identical results.
   flux, fluxerr, flag = sep.sum_circle(data, objects['x'], objects['y'], 3.0,
                                        var=bkg.globalrms**2, gain=1.0)

   # Error or variance can be arrays, indicating that the background error
   # is variable.
   bkgrms = bkg.rms()  # array, same shape as data
   flux, fluxerr, flag = sep.sum_circle(data, objects['x'], objects['y'], 3.0,
                                        err=bkgrms, gain=1.0)

   # If your uncertainty array already includes Poisson noise from the object,
   # leave gain as None (default):
   flux, fluxerr, flag = sep.sum_circle(data, objects['x'], objects['y'], 3.0,
                                        err=error_array)

   # If your data represent raw counts (it is not background-subtracted),
   # set only gain to get the poisson error:
   flux, fluxerr, flag = sep.sum_circle(data, objects['x'], objects['y'], 3.0,
                                        gain=1.0)

**Masking** 

Apply a mask (same shape as data). Pixels are igored where mask is True::

   flux, fluxerr, flag = sep.sum_circle(data, objects['x'], objects['y'], 3.0,
                                        mask=mask)

**Local background subtraction**

Perform local background subtraction in an annulus between 6 and 8 pixel
radius::

   flux, fluxerr, flag = sep.sum_circle(data, objects['x'], objects['y'], 3.0,
                                        mask=mask, bkgann=(6., 8.))

Pixels in the background annulus are not subsampled and any masked
pixels in the annulus are completely igored rather than corrected.
The inner and outer radii can also be arrays.

**Equivalent of MAG_AUTO or FLUX_AUTO in Source Extractor**

This is a two-step process. First we calculate the Kron radius for each
object::

    kronrad, flag = sep.kron_radius(data, objs['x'], objs['y'], objs['a'],
                                    objs['b'], objs['theta'], 6.0)



**Mask image regions**

Create a boolean array with elliptical regions set to True::

   mask = np.zeros(data.shape, dtype=np.bool)
   sep.mask_ellipse(mask, objs['x'], objs['y'], obs['a'], objs['b'],
                    objs['theta'], r=3.)
