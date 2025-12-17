.. _databox:

The DataBox
===========

The fundamental data type in ``spiner`` is the ``DataBox``. A
``DataBox`` packages a multi-dimensional (up to six dimensions) array
with routines for interpolating on the array and for saving the data
to and loading the data from file.

To use databox, simply include the relevant header:

.. code-block:: cpp

  #include <databox.hpp>

``DatBox`` is templated on underyling data type, which defaults to the
``Real`` type provided by ``ports-of-call``. (This is usually a
``double``.)

.. note::
  The default type can be set to type ``float`` if the preprocessor
  macro ``SINGLE_PRECISION_ENABLED`` is defined.

Any arithmetic type is supported, although the code has only been
tested carefully with floating point numbers. To set ``DataBox`` to a
single type, you may wish to declare a type alias such as:

.. code-block:: cpp

   using DataBox = Spiner::DataBox<double>

Spiner is also templated on how the interpolation gridding works. This
template parameter is called ``Grid_t``. The available options at this time are:

* ``Spiner::RegularGrid1D<T>``
* ``Spiner::PiecewiseGrid1D<T>``

where here ``T`` is the arithmetic type as discussed above. The
default type is ``RegularGrid1D``. You can further alias ``DataBox``
as, for example:

.. code-block:: cpp

   using DataBox = Spiner::DataBox<double, Spiner::RegularGrid1D<double>>;

More detail on the interpolation gridding is available below and in
the interpolation section.

.. note::
   In C++17 and later, you can also get the default type specialization
   by simply omitting the template arguments.

.. note::
  In the function signatures below, GPU/performance portability
  decorators have been excluded for brevity. However they are present
  in the actual code.

.. note::
   In the function signatures below, we will often refer to the type
   ``Real`` and the type ``T``. These are both references to the
   underlying templated arithmetic type.
   
Creating a ``DataBox``
^^^^^^^^^^^^^^^^^^^^^^

You can create a ``DataBox`` of a given shape via the constructor:

.. code-block:: cpp

  int nx1 = 2;
  int nx2 = 3;
  int nx3 = 4;
  Spiner::DataBox<double> db(nx3, nx2, nx1);

The constructor takes any number of shape values (e.g., ``nx*``) up to
six (or ``Spiner::MAXRANK``) values. Zero shape values initializes an
empty, size-zero array.

.. note::
  ``DataBox`` is row-major ordered. By convention, ``x3`` is the
  slowest moving index and ``x1`` is the fastest.

If GPU support is enabled, a ``DataBox`` can be allocated on either
host or device, depending on the ``AllocationTarget``. For example, to
explicitly allocate one array on the host and one on the device, you
might call:

.. code-block:: cpp

  // Allocates on the host (CPU)
  Spiner::DataBox<double> db_host(Spiner::AllocationTarget::Host, nx2, nx1);
  // Allocates on the device (GPU)
  Spiner::DataBox<double> db_dev(Spiner::AllocationTarget::Device, nx2, nx1);

.. note::
  If GPU support is not enabled, these both allocate on host.

You can also wrap a ``DataBox`` around a pointer you allocated
yourself. For example:

.. code-block:: cpp

  std::vector<double> mydata(nx1*nx2);
  Spiner::DataBox<double> db(mydata.data(), nx2, nx1);

You can also resize a ``DataBox``, which you can use to modify a
``DataBox`` in-place. For example:

.. code-block:: cpp

  Spiner::DataBox<double> db; // empty
  // clears old memory, resizes the underlying array,
  // and resets strides
  db.resize(nx3, nx2, nx1);
 
Just like the constructor, ``resize`` takes an optional (first)
argument for the ``AllocationTarget``.

.. warning::
  ``DataBox::resize`` is destructive. The underlying data is not preserved.

If you want to change the stride without changing the underlying data,
you can use ``reshape``, which modifies the dimensions of the
array, without modifying the underlying memory. For example:

.. code-block:: cpp

  // allocate a 1D databox
  Spiner::DataBox<double> db(nx3*nx2*nx1);
  // interpret it as a 3D object
  db.reshape(nx3, nx2, nx1);

.. warning::

  Make sure not to change the underlying size of the array
  when using ``reshape``. This is checked with an ``assert``
  statement, so you will get errors when compiling without
  the ``NDEBUG`` preprocessor macro.

The method

.. cpp:function:: void DataBox::reset();

sets the ``DataBox`` to be empty with zero rank.

Copying a ``DataBox`` to device
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If GPU support is enabled, you can deep-copy a ``DataBox`` and any
data contained in it from host to device with the function

.. cpp:function:: DataBox getOnDeviceDataBox(DataBox &db_host);

which returns a new databox with the data in ``db_host`` copied to
GPU. An object-oriented method

.. cpp:function:: DataBox Databox::getOnDevice() const;

exists as well, which returns a new object with the underlying data
copied to GPU.

.. note::
  If GPU support is not enabled, ``getOnDevice`` and friends are
  no-ops.

Semantics and Memory Management
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``DataBox`` has reference semantics---meaning that copying a
``DataBox`` does not copy the underlying data. In other words,

.. code-block:: cpp

  Spiner::DataBox<double> db1(size);
  Spiner::DataBox<double> db2 = db1;

shallow-copies ``db1`` into ``db2``. Especially for `Kokkos`_ like
workflows, this is very useful.

.. _`Kokkos`: https://github.com/kokkos/kokkos

.. warning::
  ``DataBox`` is neither reference-counted nor garbage-collected.
  If you create a ``DataBox`` you must clear the memory allocated
  just like you would for a pointer.

Two functions are provided for freeing memory in ``DataBox``:

.. cpp:function:: void free(DataBox& db);

and

.. cpp:function:: DataBox::finalize();

both will do the same thing and free the memory in a ``DataBox`` in a
context-dependent way. I.e., no matter what the ``AllocationTarget``
was, the appropriate memory will be freed.

.. warning::
  Do not free a ``DataBox`` if its memory is managed externally, e.g.,
  via a ``std::vector``. ``DataBox`` checks for this use-case
  via an ``assert`` statement.

You can check whether a given ``DataBox`` is empty, unmanaged, or
allocated on host or device with the

.. cpp:function:: DataBox::dataStatus() const;

method. It returns an ``enum class``, ``Spiner::DataStatus``, which
can take on the values ``Empty``, ``Unmanaged``, ``AllocatedHost``, or
``AllocatedDevice``. You can also check whether or not ``free`` should
be called with the method

.. cpp:function:: bool DataBox::ownsAllocatedMemory();

which returns ``true`` if a given databox is managing memory and
``false`` otherwise. The method

.. cpp:function:: bool DataBox::isReference();

returns ``false`` if the databox is managing memory and ``true``
otherwise.

Using ``DataBox`` with smart pointers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Smart pointers can be used to manage a ``DataBox`` and automatically
call ``free`` for you, so long as you use them with a custom
deleter. Spiner provides the following deleter for use in this
scenario:

.. code-block:: cpp

  struct DBDeleter {
    template <typename T>
    void operator()(T *ptr) {
      ptr->finalize();
      delete ptr;
    }
  };

It can be used, for example, with a ``std::unique_ptr`` via:

.. code-block:: cpp

  // needed for smart pointers
  #include <memory>

  // Creates a unique pointer pointing to a DataBox
  // with memory allocated on device
  std::unique_ptr<DataBox, Spiner::DBDeleter> pdb(
    new DataBox(Spiner::AllocationTarget::Device, N));
  
  // Before using the databox in, e.g., a GPU or Kokkos kernel, get a
  // shallow copy:
  auto db = *pdb;
  // some kokkos code...
  
  // when you leave scope, the data box will be freed.

Serialization and de-serialization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Shared memory models, such as `MPI Windows`_, require allocation of
memory through an external API call (e.g.,
``MPI_Win_allocate_shared``), which tabulated data must be written
to. ``Spiner`` supports this model through **serialization** and
**de-serialization**. The relevant methods are as follows. The
function

.. cpp:function:: std::size_t DataBox::serializedSizeInBytes() const;

reports how much memory a ``DataBox`` object requires to be externally
allocated. The function

.. cpp:function:: std::size_t serialize(char *dst) const;

takes a ``char*`` pointer, assumed to contain enough space for a
``DataBox``, and stores all information needed for the ``DataBox`` to
reconstruct itself. The return value is the amount of memory in bytes
used in the array by the serialized ``DataBox`` object. This method is
non-destructive; the original ``DataBox`` is unchanged. The function

.. cpp:function:: std::size_t DataBox::setPointer(T *src);

with the overload

.. cpp:function:: std::size_t DataBox::setPointer(char *src);

sets the underlying tabulated data from the src pointer, which is
assumed to be the right size and shape. This is useful for the
deSerialize function (described below) and for building your own
serialization/de-serialization routines in composite objects. The
function

.. cpp:function:: std::size_t DataBox::deSerialize(char *src);

initializes a ``DataBox`` to match the serialized ``DataBox``
contained in the ``src`` pointer.

.. note::

  Note that the de-serialized ``DataBox`` has **unmanaged** memory, as
  it is assumed that the ``src`` pointer manages its memory for
  it. Therefore, one **cannot** ``free`` the ``src`` pointer until
  everything you want to do with the de-serialized ``DataBox`` is
  over.

Putting this all together, an application of
serialization/de-serialization probably looks like this:

.. code-block:: cpp

  // load a databox from, e.g., file
  Spiner::DataBox<double> db;
  db.loadHDF(filename);
  
  // get size of databox
  std::size_t allocate_size = db.serialSizeInBytes();
  
  // Allocate the memory for the new databox.
  // In practice this would be an API call for, e.g., shared memory
  char *memory = (char*)malloc(allocate_size);
  
  // serialize the old databox
  std::size_t write_size = db.serialize(memory);
  
  // make a new databox and de-serialize it
  Spiner::DataBox<double> db2;
  std::size_t read_size = db2.deSerialize(memory);

  // read_size, write_size, and allocate_size should all be the same.
  assert((read_size == write_size) && (write_size == allocate_size));

.. warning::

  The serialization routines described here are **not** architecture
  aware. Serializing and de-serializing on a single architecture
  inside a single executable will work fine. However, do not use
  serialization as a file I/O strategy, as there is no guarantee that
  the serialized format for a ``DataBox`` on one architecture will be
  the same as on another. This is due to, for example,
  architecture-specific differences in endianness and padding.

.. _`MPI Windows`: https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node311.htm

Accessing Elements of a ``DataBox``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Elements of a ``DataBox`` can be accessed and set via the ``()``
operator. For example:

.. code-block:: cpp

  Spiner::DataBox<double> db(nx3, nx2, nx1);
  db(2,1,0) = 5.0;

The ``()`` operator accepts between one and six indexes. If you pass
in more indexes than the rank of the array, the excess indices are
ignored. If you pass in fewer, the unset indices are assumed to be
zero. The exception is the one-dimensional operator. You can always
stride through the "flattened" array by using the one-dimensional
accessor. For example:

.. code-block:: cpp

  for (int i = 0; i < nx3*nx2*nx1; ++i) {
    db(i) = static_cast<double>(i);
  }

fills the three-dimensional array above with the flat index of each
element.

Slicing
^^^^^^^^

A new ``DataBox`` containing a shallow slice of another ``DataBox``
can be constructed with the ``slice`` method:

.. cpp:function:: DataBox DataBox::slice(const int dim, const int indx, const int nvar) const;

this is fairly limited functionality. It returns a new ``DataBox``
containing only elements from ``indx`` to ``indx + nvar - 1`` in the
``dim`` direction. All other directions are unchanged. The slowest
moving dimension can be sliced to a single index with

.. cpp:function:: DataBox DataBox::slice(const int indx) const;
   
and the slowst-moving two dimensions can be sliced to a single pair of
indicies with

.. cpp:function:: DataBox DataBox::slice(const int i2, int i1) const;

Index Types and Interpolation Ranges
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Often-times an array mixes "continuous" and "discrete" variables. In
other words, some indices of an array are discretizations of a
continuous quantity, and we want to interpolate in those directions,
but other indices are discrete---they may index a particle species,
for example. A common example is in neutrino transport, where an array
of emissivities may depend on fluid density, fluid temperature,
electron fraction, neutrino energy, and neutrino species. The species
can only take three discrete values, but the density, temperature, and
electron fraction are all continuous.

``Spiner`` accounts for this by assigning each dimension in the array
a "type," represented as an ``enum class``, ``IndexType``. Currently
the type can be either ``Interpolated`` or ``Indexed``. When a new
``DataBox`` is created, all dimensions are set to
``IndexType::Indexed``. A dimension can be set to ``Interpolated`` via
the ``setRange`` method.

.. cpp:function:: void DataBox::setRange(int i, Grid_t g);

where here ``i`` is the dimension and ``g`` is the gridding object for
this index. In the default setup, where grids are uniformly spaced
(i.e., you use a ``RegularGrid1D``), this is:

.. cpp:function:: void DataBox::setRange(int i, T min, T max, int N);
   
where here ``i`` is the dimension, ``min`` is the minimum value of the
independent variable, ``max`` is the maximum value of the indpendent
variable, and ``N`` is the number of points in the ``i``
dimension. (Here ``T`` is the underlying templated data type.)

.. warning::
  In these routines, the dimension is indexed from zero. Moreover, the
  indexing possibly the opposite direction you're expecting. The zeroth
  dimension is the right most index when ``interpToReal`` is called.

.. note::
  There is a set of lower-level objects, ``RegularGrid1D``, and
  ``PiecewiseGrid1D``, which represent these interpolation ranges
  internally. There is a getter method ``range`` that works
  with the underlying ``Grid_t`` class directly. For
  more details, see the relevant documentation.

It's often desirable to have multiple databoxes with the exact same
shape and interpolation structure (i.e., independent variable
ranges). In this case, the method

.. cpp:function:: void DataBox::copyMetadata(const DataBox &src);

can assist. This method resets and re-allocates the data in a
``DataBox`` to the exact same size and shape as ``src``. More
importantly, it also copies the relevant ``IndexType`` and independent
variable range for each dimension.

One can also manually set the ``IndexType`` in a given dimension with

.. cpp:function:: void DataBox::setIndexType(int i, IndexType t);
   
and retrieve the ``IndexType`` with

.. cpp:function:: IndexType &DataBox::indexType(const int i);

to see if a dimension is interpolatable.

Interpolation to a real number
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The family of ``DataBox::interpToReal`` methods interpolate the
"entire" ``DataBox`` to a real number. Up to four-dimensional
interpolation is supported:

.. cpp:function:: T DataBox::interpToReal(const T x) const;

.. cpp:function:: T DataBox::interpToReal(const T x2, const T x1) const;

.. cpp:function:: T DataBox::interpToReal(const T x3, const T x2, const T x1) const;

.. cpp:function:: T DataBox::interpToReal(const T x4, const T x3, const T x2, const T x1) const;

where ``x1`` is the fastest moving direction, ``x2`` is less fast, and
so on. These interpolation routines are hand-tuned for performance.

.. warning::
  Do not call ``interpToReal`` with a ``DataBox`` that is the wrong shape
  or try to interpolate on indices that are not interpolatable.
  This is checked with an ``assert`` statement.

Mixed interpolation and indexing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the case where an array has some dimensions that are discrete and
some that are interpolatable, one can fuse interpolation and indexing
into a single operation. These operations are still named
``DataBox::interpToReal``, but one of the input arguments is an
integer instead of a floating point number. The location of the
integer in the function signature indicates which dimension in the
``DataBox`` is indexed. For example:

.. cpp:function:: T DataBox::interpToReal(const T x3, const T x2, const T x1, const int idx) const;

interpolates the three slower-moving indices and indexes the fastest
moving index. On the other hand,

.. cpp:function:: T DataBox::interpToReal(const T x4, const T x3, const T x2, const int idx, const T x1) const;

interpolates the fastest moving index, then indexes the
second-fastest, then interpolates the remaining three slower. The
above fused operations are the only ones currently supported.

Interpolating into another ``DataBox``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There is limited functionality for filling a ``DataBox`` with the
interpolated values of another ``DataBox``. For example, the method

.. cpp:function:: void DataBox::interpFromDB(const DataBox &src, const T x);

allocates the ``DataBox`` to have a rank one lower than ``src`` and
fill it with the faster moving elements of ``src`` interpolated to
``x`` in the slowest-moving direction. Similarly for

.. cpp:function:: void DataBox::interpFromDB(const DataBox &src, const T x2, const T x1);

The methods

.. cpp:function:: DataBox Databox::InterpToDB(const T x) const;

and

.. cpp:function:: DataBox Databox::InterpToDB(const T x2, const T x1);

return a new ``DataBox`` object, rather than setting it from a source ``DataBox``.

File I/O
^^^^^^^^^

If `hdf5`_ is enabled, ``Spiner`` can save an array to or load an
array from disk. Each array so-saved is also saved with the
``IndexType`` and independent variable ranges bundled with it, so that
knowledge of how to interpolate the data is automatically
available.

.. _`hdf5`: https://www.hdfgroup.org/solutions/hdf5/

The following methods are supported:

.. cpp:function:: herr_t DataBox::saveHDF(const std::string &filename) const;
   
saves the ``DataBox`` to a file with ``filename``.

.. cpp:function:: herr_t DataBox::saveHDF(hid_t loc, const std::string &groupname) const;

saves the ``DataBox`` as an hdf5 group at the location ``loc`` in an hdf5 file.

.. cpp:function:: DataBox::loadHDF(const std::string &filename);

fills the ``DataBox`` from information in the root of a file with ``filename``.

.. cpp:function:: DataBox::loadHDF(hid_t loc, const std::string &groupname);

fills the ``DataBox`` from information in the group with ``groupname``
based at location ``loc`` in the file.

.. warning::
  HDF5 I/O is only supported for single- and double-precision types at this time.

Miscellany
^^^^^^^^^^^

Here we list a few convenience functions available that were not
covered elsewhere.

.. cpp:function:: T DataBox::min() const;

and

.. cpp:function:: T DataBox::max() const;

compute and return the minimum and maximum values (respectively) in the array.

.. cpp:function:: int rank() const;

returns the rank (number of dimensions) of the array.

.. cpp:function:: int size() const;

returns the total number of elements in the underlying array.

.. cpp:function:: int sizeBytes() const;

returns the total size of the underlying array in bytes.

.. cpp:function:: int dim(int i) const;

returns the size in a given dimension/direction, indexed from zero.
