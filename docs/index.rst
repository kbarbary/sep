SEP
===

Source Extraction and Photometry in Python

Usage
-----

**Background estimation**

::

   import sep

   # Measure a spatially variable background of some image data
   # (a numpy array)
   bkg = sep.Background(data)
    
   # ... or with some optional parameters
   bkg = sep.Background(data, mask=mask, bw=64, bh=64, fw=3, fh=3)
    
   # The above creates a Background object which you can then use in 
   # several ways:

   # Evaluate the spatially variable background:
   back = bkg.back()  # creates an array, same shape and type as data

   # Evaluate the spatially variable RMS of the background:
   rms = bkg.rms()  # creates an array, same shape and type as data

   # Directly subtract the background from data without allocating any
   # new memory:
   bkg.subfrom(data)

   bkg.globalback  # Global "average" background level
   bkg.globalrms  # Global "average" RMS of background

**Object detection**

::

   # extract objects from data (second argument is threshold)
   objects = sep.extract(data, 1.5 * bkg.globalrms)

   # ... or with some optional parameters:
   kernel = np.array([[1., 2., 1.],
                      [2., 4., 2.],
		      [1., 2., 1.]])
   objects = sep.extract(data, 1.5 * bkg.globalrms, minarea=5, conv=kernel,
                         clean=True, clean_param=1.0, deblend_nthresh=32,
                         deblend_cont=0.005)

   # objects is a numpy structured array:
   len(objects)  # number of objects
   objects['x'][i]  # flux-weighted x-center coordinate of i-th object
   ...              # ... and many other fields. See objects.dtype.names

**Aperture photometry**

The follow examples demonstrate options for circular aperture photometry::

   # sum flux in circles of radius=3.0
   flux, fluxerr, flag = sep.apercirc(data, objects['x'], objects['y'], 3.0)

   # x, y and r can be arrays and obey numpy broadcasting rules
   flux, fluxerr, flag = sep.apercirc(data, objects['x'], objects['y'],
                                      3.0 * np.ones(len(objects)))

   # use a different subpixel sampling (default is 5)
   flux, fluxerr, flag = sep.apercirc(data, objects['x'], objects['y'], 3.0,
                                      subpix=10)

   # Specify a per-pixel "background" error (default is zero):
   flux, fluxerr, flag = sep.apercirc(data, objects['x'], objects['y'], 3.0,
                                      err=bkg.globalrms, gain=1.0)

   # Use "var" for variance rather than error:
   flux, fluxerr, flag = sep.apercirc(data, objects['x'], objects['y'], 3.0,
                                      var=bkg.globalrms**2, gain=1.0)

   # err/var can also be an array:
   bkgrms = bkg.rms()  # array, same shape as data
   flux, fluxerr, flag = sep.apercirc(data, objects['x'], objects['y'], 3.0,
                                      err=bkgrms, gain=1.0)

   # If your uncertainty array includes Poisson noise from the object,
   # leave gain as None (default):
   flux, fluxerr, flag = sep.apercirc(data, objects['x'], objects['y'], 3.0,
                                      err=error_array)

   # If your data represent raw counts (not background-subtracted), set only
   # gain to get the poisson error:
   flux, fluxerr, flag = sep.apercirc(data, objects['x'], objects['y'], 3.0,
                                      gain=1.0)

   # Apply a mask (same shape as data). Pixels are igored where mask is True.
   flux, fluxerr, flag = sep.apercirc(data, objects['x'], objects['y'], 3.0,
                                      mask=mask)

   # Perform local background subtraction in an annulus between 6 and 8 pixel
   # radius. Pixels in the background annulus are not subsampled and any masked
   # pixels in the annulus are completely igored rather than corrected.
   # Inner and outer radii can also be arrays: 
   flux, fluxerr, flag = sep.apercirc(data, objects['x'], objects['y'], 3.0,
                                      mask=mask, bkgann=(6., 8.))

   # Convert flag array to boolean for specific flags:
   sep.istruncated(flag)  # True where aperture was truncated by image edge.
   sep.hasmasked(flag)    # True where aperture includes masked pixels.

**Mask image regions**

::

   # Create a boolean array with elliptical regions set to True:
   mask = np.zeros(data.shape, dtype=np.bool)
   sep.mask_ellipse(mask, objects['x'], objects['y'],
                    cxx=objects['cxx'], cyy=objects['cyy'], cxy=objects['cxy'],
                    scale=3.)

API
---

**Global Background Estimation**

.. autosummary::
   :toctree: api
   
   sep.Background
   sep.Background.back
   sep.Background.rms
   sep.Background.subfrom

**Object detection**

.. autosummary::
   :toctree: api

   sep.extract

**Aperture photometry & masking**

.. autosummary::
   :toctree: api
   
   sep.apercirc
   sep.aperellip
   sep.mask_ellipse
   sep.kronrad


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

