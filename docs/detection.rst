========================================
Background estimation & source detection
========================================

Most optical/IR data must be background subtracted before sources can be detected. In SEP, background estimation and source detection are two separate steps.

Background estimation
=====================

Given a numpy array ``data`` containing the image data,

::

   import numpy as np
   import sep

   # Measure a spatially variable background of some image data
   # (a numpy array)
   bkg = sep.Background(data)
    
   # ... or with some optional parameters
   bkg = sep.Background(data, mask=mask, bw=64, bh=64, fw=3, fh=3)
    
This creates a `~sep.Background` object which you can then use in
several ways::

   # Evaluate the spatially variable background:
   back = bkg.back()  # creates an array, same shape and type as data
   back = np.array(bkg)  # equivalent to the above

   # Evaluate the spatially variable RMS of the background:
   rms = bkg.rms()  # creates an array, same shape and type as data

   # Subtract the background from the data
   bkg_subtraced_data = data - bkg

   # Directly subtract the background from the data in place
   bkg.subfrom(data)  # data is modified in-place

   # get global statistics
   bkg.globalback  # Global "average" background level
   bkg.globalrms  # Global "average" RMS of background

.. note::

   If you are using SEP to analyze data read from FITS files with
   `astropy.io.fits <http://astropy.readthedocs.org/en/stable/io/fits/>`_,
   you may see an error message such as::

        ValueError: Input array with dtype '>f4' has non-native byte order.
        Only native byte order arrays are supported. To change the byte
        order of the array 'data', do 'data = data.byteswap().newbyteorder()'

   It is usually easiest to do this byte-swap operation directly after
   reading the array from the FITS file. You can even perform the byte
   swap in-place by doing ::

       >>> data = data.byteswap(inplace=True).newbyteorder()

   If you do this in-place operation, ensure that there are no other
   references to ``data``, as they will be rendered nonsensical.

   For the interested reader, this byteswap operation is necessary because
   astropy.io.fits always returns big-endian byte order arrays, even on
   little-endian machines. For more on this, see
   `this FAQ <https://github.com/kbarbary/sep#faq>`_.


Source detection
================

Once the data has been background subtracted, detect objects in the 
data, given some threshold::

   thresh = 1.5 * bkg.globalrms
   objects = sep.extract(data, thresh)

   # objects is a numpy structured array:
   len(objects)  # number of objects
   objects['x'][i]  # flux-weighted x-center coordinate of i-th object
   ...              # ... and many other fields.

See the reference section for all the fields available in the returned
structured array.

.. warning::

   If the data array is not background-subtracted or the threshold is
   too low, you will tend to get one giant object or, more likely, an
   exception will be raised due to exceeding the internal memory
   constraints of the function.


