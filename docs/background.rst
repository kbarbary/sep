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
   np.array(bkg)  # equivalent to the above

   # Evaluate the spatially variable RMS of the background:
   rms = bkg.rms()  # creates an array, same shape and type as data

   # Subtract the background from the data
   bkg_subtraced_data = data - bkg

   # Directly subtract the background from the data in place
   # (modifies the data array)
   bkg.subfrom(data)

   bkg.globalback  # Global "average" background level
   bkg.globalrms  # Global "average" RMS of background
