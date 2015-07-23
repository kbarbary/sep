Derivation of the matched filter formula
----------------------------------------

Assume that we have an image containing a single point source. This produces a
signal with PSF :math:`S_i` and noise :math:`N_i` at each pixel indexed by
:math:`i`. Then the measured image data :math:`D_i` (i.e. our pixel values) is
given by:

.. math::
    D_i = S_i + N_i

Then we want to apply a linear transformation :math:`T_i` which gives an
output :math:`Y`:

.. math::
    Y = \sum_i T_i D_i = T^T D

We use matrix notation from here on and drop the explicit sums. Our objective
is to find the transformation :math:`T_i` which maximizes the signal-to-noise
ratio :math:`SNR`.

.. math::
    SNR = \frac{(T^T S)^2}{E[(T^T N)^2]}

We can expand the denominator as:

.. math::
    E[(T^T N)^2] &= E[(T^T N)(N^T T)] = T^T \cdot E[N N^T] \cdot T = T^T C T

Where :math:`C_{ik}` is the covariance of the noise between pixels :math:`i`
and :math:`k`. Now using the Cauchy-Schwarz inequality on the numerator:

.. math::
    (T^T S)^2 = (T^T C^{1/2} C^{-1/2} S)^2 \le (T^T C^{1/2})^2 (C^{-1/2} S)^2 =
    (T^T C T) (S^T C^{-1} S)

since :math:`C^T = C`. The signal-to-noise ratio is therefore bounded by:

.. math::
    &SNR \le \frac{(T^T C T)(S^T C^{-1} S)}{(T^T C T)} \\
    &SNR \le S^T C^{-1} S

Choosing :math:`T = \alpha C^{-1} S` where :math:`\alpha` is an arbitrary
normalization constant, we get equality. Hence this choise of :math:`T` is the
optimal linear tranformation. We normalize this linear transformation so that
if there is no signal and only noise, we get an expected signal-to-noise ratio
of 1. With this definition, the output :math:`SNR` represents the number of
standard deviations above the background. This gives:

.. math::
    &E[(T^T N)^2] = T^T C T = \alpha^2 S^T C^{-1} C C^{-1} S = \alpha^2 S^T
    C^{-1} S = 1 \\ 
    &\alpha = \frac{1}{\sqrt{S^T C^{-1} S}}

Putting everything together, our normalized linear transformation is:

.. math::
    T = \frac{C^{-1} S}{\sqrt{S^T C^{-1} S}}

And the optimal signal-to-noise is given in terms of the known variables as:

.. math::
    SNR = \frac{S^T C^{-1} D}{\sqrt{S^T C^{-1} S}}


Implementation
--------------

The above derivation gives the optimal signal-to-noise for a single point
source. For object detection, we assume that for each output pixel the data
consists of a single point source at the location of that output pixel. We then
apply the matched filter algorithm for each output pixel individually.
In Source Extractor, the matched filter is implemented in the case where there
is no noise, or equal noise across all pixels. The above formula then
simplifies to a convolution of the data with the PSF. In SEP, we support this
same behaviour in `~sep.extract` when no error array is specified and a
convolution kernel is passed.

SEP also supports having independent errors on each of the input pixels in
`~sep.extract`, which will run the matched filter algorithm with a diagonal
covariance matrix. This is a new addition that is not available in Source
Extractor. Some benefits of this method are that flagged bad pixels are
ignored, detector sensitivity can be taken into account, and edge effects are
handled gracefully. See the SEP tests for examples of situations where you
might want to use the matched filter with a noise array.

Ideally, the convolution kernel should be set to shape of the PSF in the data
when detecting stars. For galaxy detection, a larger kernel could be used. In
practice, anything that is roughly the shape of the desired object works well
since the main goal is to negate the effects of the background, and a
reasonable estimate is good enough.
