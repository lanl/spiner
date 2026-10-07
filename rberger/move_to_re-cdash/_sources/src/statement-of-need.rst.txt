.. _statement-of-need:

Why Develop Spiner?
====================

As Moore's law comes to an end, more and more performance comes from
specialized hardware, such as GPUs. A key tool in the toolbox for many
scientific codes is tabulated data. Fluid and continuum dynamics codes
often encapsulate the equation of state as data tabulated in density
and temperature. Radiation transport uses emissivity and absorption
opacity on tables. As continuum dynamics is required for a variety of
applications, such as astrophysics, geophysics, climate science,
vehicle engineering, and national security, utilizing a very large
number of supercomputer cycles, providing interpolation on tabulated
data for these applications has the potential for significant impact.

These capabilities must be supported on all hardware a code may be run
on, whether this is an NVIDIA GPU, an Intel CPU, or a next generation
accelerator manufactured by one of any number of hardware vendors. To
our knowledge there is no performance portable interpolation library
on which these codes can rely, and there is a clear need, which we
have developed ``Spiner`` to meet.

To see some examples of software projects that leverage ``Spiner`` see
`singularity-EOS`_, `singularity-opac`_, and `Phoebus`_.

.. _singularity-eos: https://github.com/lanl/singularity-eos

.. _singularity-opac: https://github.com/lanl/singularity-eos

.. _Phoebus: https://github.com/lanl/singularity-opac

State of the Field
^^^^^^^^^^^^^^^^^^^

Interpolation is a common problem, implemented countless times across
software projects, and a core part of any introductory text on
scientific computing. In graphics applications interpolation is so
ubiquitous that hardware primitives are provided by GPUs. These
hardware intrinsics are, however, severely limited for scientific
application. For example, on NVIDIA GPUs, the values to be
interpolated must be single precision floating point, and the
interpolation coefficients themselves are only half-precision, which
is often insufficient to capture the high precision required for
scientific applications. As GPUs are inherently vector devices,
hardware interpoaltion is also vectorized in nature. However,
downstream applications may be easier to reason about if scalar
operations are available. For example, equation of state lookups often
require root finds on interpolated data, and this can be easier to
implement as a scalar operation, even if the final operation is
vectorized over warps. Texture interpolation also does not support
multi-dimensional mixed indexing/interpoaltion operations where, say,
three indices of a four-dimensional array are interpolated and one is
merely indexed into.

Moreover, relying on hardware intrinsics is not a portable solution. A
software interpolation library can, if written with care, work on not
only the current generation of accelerators, but also on general
purpose CPUs and the next generation of hardware as well.

Unfortunately, a performance-portable implementation not tuned to a
specific use-case or embedded in a larger project is (to our
knowledge) not available in the literature. A common problem in
performance-portable computing is the management of
performance-portable data structures.

Interpolation is far more ubiquitous than its application in continuum
dynamics and radiation transport, and we expect Spiner will find
applications in the broader space of applications, such as image
resampling. However, the team built Spiner with simulations in mind.
