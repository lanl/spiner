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

Interpolation is linear. Here's an example of 3D interpolation (2D
slice shown) on a GPU, with second-order convergence:

.. image:: ../../figs/convergence.png

See below for details of how to use spiner in your project and how to
develop for it.

Spiner also relies on `Ports of Call`_ as a simple performance
portability layer. Ports of Call is included as a submodule, and
automatically integrated into the build system.

.. _Ports of Call: https://lanl.github.io/ports-of-call/main/index.html

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   src/building
   src/getting-started
   src/databox
   src/interpolation
   src/sphinx-howto

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

This documentation is approved for unlimited release, LA-UR-22-20363.
