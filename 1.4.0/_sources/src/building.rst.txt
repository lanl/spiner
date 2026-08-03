.. _building:

Building and Installation
==========================

``Spiner`` is self-contained and header-only. Clone it as:

.. code-block:: bash

  git clone --recursive git@github.com:lanl/spiner.git


Building from source
^^^^^^^^^^^^^^^^^^^^^

To build tests and install:

.. code-block:: bash

  mkdir -p spiner/bin
  cd sppiner/bin
  cmake -DBUILD_TESTING=ON
  make -j
  make test
  make install

Spiner supports a few ``cmake`` configuration options:

* ``BUILD_TESTING`` enables tests
* ``SPINER_USE_HDF5`` enables support for saving and loading tables as `hdf5`_.
* ``SPINER_HDF5_INSTALL_DIR`` tells the build system where `hdf5`_ is located.
* ``SPINER_USE_KOKKOS`` enables `Kokkos`_ as a backend
* ``SPINER_USE_KOKKOS_SRC`` tells the build system to build `Kokkos`_ from source, and where the source directory is located. Note that if you use this option, you cannot install Spiner, only build the tests.
* ``SPINER_KOKKOS_INSTALL_DIR`` tells the build system where to find pre-compiled `Kokkos`_
* ``SPINER_USE_CUDA`` enables the Kokkos cuda backend
* ``CMAKE_INSTALL_PREFIX`` sets the install location
* ``CMAKE_BUILD_TYPE`` sets the build type
* ``SPINER_FORCE_INTERNAL_PORTS`` forces use of a `ports-of-call`_ submodule rather than a system install

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

.. warning::
  The spack build is currently experimental. 
  Please report problems you have as github issues.

Although the spackage has not yet made it to the main `Spack`_
repositories, we provide a spackage for ``Spiner`` within the
the source repository. If you have spack installed,
simply call

.. _Spack: https://spack.io/

.. code-block:: bash

  spack repo add spiner/spack-repo
  spack install spiner

The spack repo supports a few variants:

* ``+kokkos`` enables the Kokkos backend
* ``+cuda`` enables the cuda backend. A ``cuda_arch`` must be specified.
* ``+hdf5`` enables HDF5 file support.
* ``+mpi`` enables parallel hdf5 support
* ``+python`` installs python, numpy, and matplotlib support
* ``+doc`` adds tooling for building the docs
* ``+format`` adds support for clang-format

Including Spiner in your Project
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Spiner can be included into a cmake project, either in-tree as a
submodule or after installation via ``find_package``.
The cmake system provides the ``spiner::spiner`` cmake target.
