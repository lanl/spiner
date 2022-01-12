.. _ports-of-call:

Ports of Call
==============

Ports of call is a header-only library that provides a bit of
flexibility for performance portability. At the moment it mainly
provides a one-header abstraction to enable or disable `Kokkos`_ in a
code. However other backends can be added. (If you're interested in
adding a backend, please let us know!)

.. _Kokkos: https://github.com/kokkos/kokkos

We define a few portability macros which are useful:

1. ``PORTABLE_FUNCTION``: decorators necessary for compiling a kernel function
2. ``PORTABLE_INLINE_FUNCTION``: ditto, but for when functions ought to be inlined
3. ``PORTABLE_FORCEINLINE_FUNCTION``: forces the compiler to inline
4. ``PORTABLE_LAMBDA``: Resolves to a ``KOKKOS_LAMBDA`` or to ``[=]`` depending on context
5. ``_WITH_KOKKOS_``: Defined if Kokkos is enabled.
6. ``_WITH_CUDA_``: Defined when Cuda is enabled
7. ``Real``: a typedef to double (default) or float (if you define ``SINGLE_PRECISION_ENABLED``)
8. ``PORTABLE_MALLOC()``, ``PORTABLE_FREE()``: A macro or wrapper for kokkos_malloc or cudaMalloc, or raw malloc.

At compile time, you define
``PORTABILITY_STRATEGY_{KOKKOS,CUDA,NONE}`` (if you don't define it,
it defaults to NONE). The above macros then behave as expected. In
particular, ``PORTABLE_FUNCTION`` and friends resolve to ``__host__
__device__`` decorators as appropriate.

There are to be two headers in this library, for different use cases.

portability.hpp
^^^^^^^^^^^^^^^^

``portability.hpp`` provides the above-mentioned macros for decorating
functions. Also provides loop abstractions that can be leveraged by a
code. These loop abstractions are of the form:

.. cpp:function:: void portableFor(const char *name, int start, int stop, Function Function)

where ``Function`` is a template parameter and should be set to a
functor that takes one index, e.g., an index in an array. For example:

.. code-block:: cpp

  portableFor("Example", 0, 5,
    PORTABLE_LAMBDA(int i) {
      printf("hello from thread %d\n", i);
  });

``start`` is inclusive, ``stop`` is exclusive. Up to five-dimensional
``portableFor`` loops are available. For example:

.. code-block:: cpp

  template <typename Function>
  void portableFor(const char *name, int startb, int stopb, int starta, int stopa,
    int startz, int stopz, int starty, int stopy, int startx,
    int stopx, Function function) {

We also provide ``portableReduce``, however the functionality is very
limited. The syntax is:

.. code-block::

  template <typename Function, typename T>
  void portableReduce(const char *name, int starta, int stopa, int startz,
    int stopz, int starty, int stopy, int startx, int stopx,
    Function function, T &reduced) {

where ``Function`` now takes as many indices are required and
``reduced`` as arguments.

portable_arrays.hpp
^^^^^^^^^^^^^^^^^^^

``portable_arrays.hpp`` provides a wrapper class, ``PortableMDArray``,
around a contiguous block of host or device memory that knows stride
and layout, enabling one to mock up multidimensional arrays from a
pointer to memory. The design is heavily inspired by the
``AthenaArray`` class from `Athena++`_.

.. _`Athena++`: https://www.athena-astro.app

One constructs a ``PortableMDArray`` by passing it a pointer to
underlying data and a shape. For example:

.. code-block:: cpp

  #include <portability.hpp>
  #include <portable_arrays.hpp>
  constexpr int NX = 2;
  constexpr int NY = 3;
  constexpr int NZ = 4;
  Real *data = (Real*)PORTABLE_MALLOC(NX*NY*NZ*sizeof(Real));
  PortableMDArray<Real> my_3d_array(data, NZ, NY, NX);

Note that ``PortableMDArray`` is templated on underlying data
type. Note also that ``PortableMDArray`` is column-major-ordered. The
slowest moving index is ``z`` and the fastest is ``x``. You can then
set or access an element by reference as:

.. code-block:: cpp

  // z = 3, y = 2, x = 1
  my_3d_array(3,2,1) = 5.0;

You can always access the "flat" array by simply using the 1D accessor:

.. code-block:: cpp

  my_3d_array(6) = 2.0;

By default ``PortableMDArray`` has reference-semantics. In
other words, copies are shallow.

You can assign new data and a new shape to a ``PortableMDArray`` with
the ``NewPortableMDArray`` function. For example:

.. code-block:: cpp

  my_3d_array.NewPortableArray(new_data, 9, 8, 7);

would reshape ``my_3d_array`` to be of shape 7x8x9 and point it at the
``new_data`` pointer.

``PortableMDArray`` also provides a few useful methods:

.. cpp:function:: size_t PortableMDArray::GetRank()

provides the number of dimensions of the array.

.. cpp:function:: int PortableMDArray::GetDim(size_t i)

returns the size of a given dimension (indexed from 1, not 0).

.. cpp:function:: int PortableMDArray::GetSize()

returns the size of the flattened array.

.. cpp:function:: size_t PortableMDArray::GetSizeInBytes()

returns the size of the flattened array in bytes.

.. cpp:function:: bool PortableMDArray::IsEmpty()

returns true if the array is empty and false otherwise.

.. cpp:function:: T* PortableMDArray::data()

returns the underlying pointer. The ``begin()`` and ``end()``
functions return pointers to the beginning and end of the array.

.. cpp:function:: void PortableMDArray::Reshape(int nx3, int nx2, int nx1)

resets the shape of the array without pointing to a new underlying
data pointer. It accepts anywhere between 1 and 6 sizes.

``PortableMDArray`` also supports some simple boolean comparitors,
such as ``==`` and arithmetic such as ``+``, and ``-``.
