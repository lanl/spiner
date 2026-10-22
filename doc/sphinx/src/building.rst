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

.. _building:

.. Generative AI was used to assist with modifications to this file.

Building and Installation
==========================

``Spiner`` is header-only. It requires C++20, CMake 3.24 or newer, and
Ports-of-Call 3.1.1 or newer. Clone it as:

.. code-block:: bash

  git clone --recursive git@github.com:lanl/spiner.git


Building from source
^^^^^^^^^^^^^^^^^^^^^

To build tests and install:

.. code-block:: bash

  cd spiner  # cd to repo root
  cmake -S . -B bin -DSPINER_BUILD_TESTS=ON -DSPINER_TEST_USE_KOKKOS=OFF
  cmake --build bin --parallel 6
  ctest --test-dir bin --output-on-failure -j 6
  cmake --install bin --prefix /path/to/install

Spiner supports a few ``cmake`` configuration options:

* ``SPINER_BUILD_TESTS`` builds tests (default: ``OFF``).
* ``SPINER_USE_HDF`` enables saving and loading tables as `hdf5`_ (default: ``OFF``).
* ``SPINER_USE_PORTABLE_NQT`` selects portable NQT functions (default: ``OFF``).
* ``SPINER_TEST_USE_KOKKOS`` enables `Kokkos`_ for tests (default: ``ON`` when tests are enabled).
* ``SPINER_TEST_USE_KOKKOS_CUDA`` enables CUDA for fetched test Kokkos (default: ``ON``).
  Set this to ``OFF`` for CPU Kokkos tests.
* ``CMAKE_PREFIX_PATH`` supplies installed dependency prefixes.
* ``FETCHCONTENT_SOURCE_DIR_<NAME>`` supplies dependency sources, for example
  ``FETCHCONTENT_SOURCE_DIR_KOKKOS``.
* ``CMAKE_INSTALL_PREFIX`` sets the install location.
* ``CMAKE_BUILD_TYPE`` sets the build type.

CMake searches installed packages first, then fetches missing dependencies
where supported. Ports-of-Call must be at least 3.1.1. HDF5 must be installed
when enabled. CUDA tests need an appropriate compiler such as ``nvcc_wrapper``
and Kokkos architecture settings for the target GPU.

.. _`Spack`: https://spack.io

.. _`hdf5`: https://www.hdfgroup.org/solutions/hdf5

.. _`Kokkos`: https://github.com/kokkos/kokkos

.. _`ports-of-call`: https://lanl.github.io/ports-of-call/main/index.html

HDF5 is searched for and configured via the usual `cmake`_ machinery.

.. _`cmake`: https://cmake.org/

A ``format_spiner`` target is also added if ``clang-format`` is found, so
that ``make format_spiner`` will auto-format the repository.

Testing is enabled via `Catch2`_, which is automatically downloaded
during the cmake configure phase if needed.

.. _`Catch2`: https://github.com/catchorg/Catch2

Spack
^^^^^^

Spiner is available through `Spack`_. If you have spack installed, simply call:

.. code-block:: bash

  spack install spiner

Including Spiner in your Project
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Spiner can be included into a cmake project, either in-tree as a
submodule or after installation via ``find_package``.
The cmake system provides the ``spiner::spiner`` cmake target.
