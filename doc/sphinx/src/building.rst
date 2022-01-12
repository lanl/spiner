.. _building:

Building and Installation
==========================

``spiner`` is self-contained and header-only. Clone it as:

.. code-block:: bash

  git clone git@github.com:lanl/spiner.git


Building from source
^^^^^^^^^^^^^^^^^^^^^

To build tests and install:

.. code-block:: bash

  mkdir -p spiner/bin
  cd spiner/bin
  cmake -DBUILD_TESTING=ON
  make -j
  make test
  make install

Spiner supports a few ``cmake`` configuration options:

* ``BUILD_TESTING`` enables tests
* ``SPINER_USE_HDF5`` enables support for saving and loading tables as `hdf5`_.
* ``SPINER_USE_KOKKOS`` enables `Kokkos`_ as a backend
* ``SPINER_USE_CUDA`` enables the Kokkos cuda backend
* ``CMAKE_INSTALL_PREFIX`` sets the install location
* ``CMAKE_BUILD_TYPE`` sets the build type

.. _`hdf5`: https://www.hdfgroup.org/solutions/hdf5

.. _`Kokkos`: https://github.com/kokkos/kokkos

HDF5 is searched for and configured via the usual `cmake`_ machinery.

.. _`cmake`: https://cmake.org/

A ``format`` target is also added if ``clang-format`` is found, so
that ``make format`` will auto-format the repository.

Spack
^^^^^^

.. warning::
  The spack build is currently experimental. 
  Please report problems you have as github issues.

Although the spackage has not yet made it to the main `Spack`_
repositories, we provide a spackage for ``spiner`` within the
the source repository. If you have spack installed,
simply call

.. _Spack: https://spack.io/

.. code-block:: bash

  spack repo add spiner/spack-repo
  spack install spiner

The spack repo supports a few variants:

* ``+kokkos`` enables the Kokkos backend
* ``+cuda`` enables the cuda backend. A ``cuda_arch`` must be specified.
* ``+python`` installs python, numpy, and matplotlib support
* ``+doc`` adds tooling for building the docs
* ``+format`` adds support for clang-format

Including Spiner in your Project
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Spiner can be included into a cmake project, either in-tree as a
submodule or after installation. The cmake system provides
``spiner::flags`` and ``spiner::libs`` cmake targets. The former adds
appropriate compilation flags, the latter adds link flags for
dependencies such as hdf5.
