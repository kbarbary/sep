Changes to the C API
====================

This page details cumulative changes to the C API between v1.2.1 and
v1.4.0, or since this package was forked from
`kbarbary/sep <https://github.com/kbarbary/sep>`_. Almost all changes to
the public API are to integer parameters, many of which have been changed
from ``int`` to ``int64_t``, to fix problems which may arise with
extremely large arrays. The other change of note is to the ``sep_image``
struct, which enables the user to pass an existing segmentation map, and
re-extract all morphological and photometric quantities on a different
image.

All changes here are transparent to users of the Python interface.

.. c:struct::  sep_image

 - Three new parameters have been added:

   .. c:var:: int64_t * segids

      The unique ids in an existing segmentation map.

   .. c:var:: int64_t * idcounts

      The number of occurrences of each unique id in an existing
      segmentation map.

   .. c:var:: int64_t numids

      The total number of unique ids in an existing segmentation map.

 - The type of ``w`` and ``h`` has changed from ``int`` to ``int64_t``.

.. c:struct:: sep_bkg

 - The type of the following parameters has changed from ``int`` to
   ``int64_t``: ``w``, ``h``, ``bw``, ``bh``, ``nx``, ``ny``, and ``n``.

.. c:struct:: sep_catalog

 - The type of the following parameters has changed from ``int`` to
   ``int64_t``: ``npix``, ``tnpix``, ``xmin``, ``xmax``, ``ymin``,
   ``ymax``, ``xcpeak``, ``ycpeak``, ``xpeak``, ``ypeak``, ``pix``, and
   ``objectspix``.

.. c:function:: int sep_background()

 - The type of the following parameters has changed from ``int`` to
   ``int64_t``: ``bw``, ``bh``, ``fw``, and ``fh``.

.. c:function:: float sep_bkg_pix()

 - The type of ``x`` and ``y`` has changed from ``int`` to ``int64_t``.

.. c:function:: int sep_bkg_line()

 - The type of ``y`` has changed from ``int`` to ``int64_t``.

.. c:function:: int sep_bkg_subline()

 - The type of ``y`` has changed from ``int`` to ``int64_t``.

.. c:function:: int sep_bkg_rmsline()

 - The type of ``y`` has changed from ``int`` to ``int64_t``.

.. c:function:: int sep_extract()

 - The type of ``convw`` and ``convh`` has changed from ``int`` to
   ``int64_t``.

.. c:function:: int sep_sum_circann_multi()

 - The type of ``n`` has changed from ``int`` to ``int64_t``.


.. c:function:: int sep_stats_circann()

 - Added ``area`` and ``sumerr`` outputs.

.. c:function:: int sep_flux_radius()

 - The type of ``n`` has changed from ``int`` to ``int64_t``.

.. c:function:: void sep_set_ellipse()

 - The type of ``w`` and ``h`` has changed from ``int`` to ``int64_t``.
