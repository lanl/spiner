.. ======================================================================
..  © (or copyright) 2019-2026. Triad National Security, LLC. All rights
..  reserved.  This program was produced under U.S. Government contract
..  89233218CNA000001 for Los Alamos National Laboratory (LANL), which is
..  operated by Triad National Security, LLC for the U.S.  Department of
..  Energy/National Nuclear Security Administration. All rights in the
..  program are reserved by Triad National Security, LLC, and the
..  U.S. Department of Energy/National Nuclear Security
..  Administration. The Government is granted for itself and others acting
..  on its behalf a nonexclusive, paid-up, irrevocable worldwide license
..  in this material to reproduce, prepare derivative works, distribute
..  copies to the public, perform publicly and display publicly, and to
..  permit others to do so.
.. ======================================================================

.. Spiner Documentation master file, created by
   sphinx-quickstart on Tue Nov 2 16:56:44 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Spiner: Performance portable routines for generic, tabulated, multi-dimensional data
=====================================================================================

`Spiner`_ is a library for storing, indexing, and interpolating
multidimensional data in a performance-portable way. It's intended to
run on CPUs, GPUs and everything in-between. You can create a table on
a CPU, copy it to a GPU, and interpolate on it in a GPU kernel, for
example.

.. _Spiner: https://github.com/lanl/spiner

Spiner also defines (via hdf5) a file format that bundles data
together with instructions for interpolating it. This means you don't
have to specify anything to start interpolating, simple load the file
and evaluate where you want.

Interpolation is linear by default. For values and analytic derivatives in
1D and 2D, ``interpToRealWithGrads`` also supports the
``Spiner::InterpolationStencils::Rational`` policy. See
:ref:`rational-interpolation` for policy selection, the four-point minimum
per axis, and boundary behavior. Value-only ``interpToReal`` calls remain
linear.

Here's an example of 3D linear interpolation (2D
slice shown) on a GPU, with second-order convergence:

.. image:: ../../figs/convergence.png

Interpolation is fast and portable. Here's a benchmark showing
performance on CPU and GPU for several architectures and problem
sizes:

.. image:: ../../figs/spiner_interpolation_benchmark.png

See below for details of how to use spiner in your project and how to
develop for it.

Spiner also relies on `Ports of Call`_ as a simple performance
portability layer. Ports of Call is included as a submodule, and
automatically integrated into the build system.

.. _Ports of Call: https://lanl.github.io/ports-of-call/main/index.html

Getting Help and Contributing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you use Spiner and need help, submit an issue to the Spiner
repository. If you'd like to contribute, just fork and submit a pull
request. See our contributing page for more details on the
contribution guidance.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   src/statement-of-need
   src/building
   src/getting-started
   src/databox
   src/interpolation
   src/contributing
   src/sphinx-howto

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

This documentation is approved for unlimited release,
LA-UR-22-20363. Generative AI was used to assist in modifying this
file.
