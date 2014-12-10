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

   Be careful! If the data array is not background-subtracted or the
   threshold is too low, you will tend to get one giant object or,
   more likely, an exception will be raised due to exceeding the
   internal memory constraints of the function.
