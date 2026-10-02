.. _interpolation:

Gridding for Interpolation
===========================

nSpiner performs interpolation on Cartesian-product
grids. There are two lower-level objects:

* ``RegularGrid1D``
* ``PiecewiseGrid1D``

These objects contain the metadata required for interpolation
operations and have a few useful userspace functions, which are
described here.

Like ``DataBox``, these grid objects are templated on
underlying data type, the default type being a ``Real`` as provided by
``ports-of-call``. You may wish to specialize to a specific type with
a type alias such as:

.. code-block:: cpp

   using RegularGrid1D = Spiner::RegularGrid1D<double>;
   using PiecewiseGrid1D = Spiner::PiecewiseGrid1D<double>;

.. note::
   In the function signature below we refer to ``T`` and ``Real`` as
   the underlying arithmetic data type.

When constructing a ``DataBox``, you may wish to specify which
interpolation object you are using. It is a template parameter.

``RegularGrid1D``
------------------

We begin by discussing ``RegularGrid1D``, as the ``PiecewiseGrid1D``
object is built on top of it.

Construction
^^^^^^^^^^^^^

A ``RegularGrid1D`` requires three values to specify an interpolation
grid: the minimum value of the independent variable, the maximum value
of the independent variable, and the number of points on the
grid. These are passed into the constructor:

.. cpp:function:: RegularGrid1D::RegularGrid1D(T min, T max, size_t N);

Default constructors and copy constructors are also provided.

Mapping an index to a real number and vice-versa
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The function

.. cpp:function:: T RegularGrid1D::x(const int i) const;

returns a "physical" position on the grid given an index ``i``.

The function

.. cpp:function:: int RegularGrid1D::index(const T x) const;

returns the index on the grid of a "physical" value ``x``.

The function

.. cpp:function:: T RegularGrid1D::min() const;

returns the minimum value on the independent variable grid.

The function

.. cpp:function:: T RegularGrid1D::max() const;

returns the maximum value on the independent variable grid.

The function

.. cpp:function:: T RegularGrid1D::dx() const;

returns the grid spacing for the independent variable.

The function

.. cpp:function:: int RegularGrid1D::nPoints() const;

returns the number of points in the independent variable grid.

The ``PiecewiseGrid1D``
------------------------

A ``PiecewiseGrid1D`` is a non-intersecting, contiguous, ordered
collection ``RegularGrid1D`` s. It can be used to construct grids with
non-uniform spacing, so long as the grid spacing is piecewise
constant.

The maximum number of ``RegularGrid1D``s that can be used to construct
a ``PiecewiseGrid1D`` is a compile-time parameter (default is 5). You
can specify a different value with, e.g.,

.. code-block:: cpp

   // Maximum number of "pieces" in a grid = 10
   using PiecewiseGrid1D = Spiner::PiecewiseGrid1D<double, 10>;

Constructiong a ``PiecewiseGrid1D``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A ``PiecewiseGrid1D`` is constructed from either a ``std::vector`` or
a ``std::initializer_list`` of ``RegularGrid1D`` s. For example:

.. code-block:: cpp

   // Initialize the regular grids
   // Note that the start and end points match
   // for each consecutive pair of grids.
   // g1 ends when g2 starts, etc.
   Spiner::RegularGrid1D<double> g1(0, 0.25, 3);
   Spiner::RegularGrid1D<double> g2(0.25, 0.75, 11);
   Spiner::RegularGrid1D<double> g3(0.75, 1, 7);

   // Build the piecewise grid. The double bracket notation
   // is an "initalizer list" and is very convenient,
   // as it is a C++ language feature.
   Spiner::PiecewiseGrid1D<double> h = {{g1, g2, g3}};

Default constructors and copy constructors are also provided.

Index Mapping with ``PiecewiseGrid1D``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A ``PiecewiseGrid1D`` has all the same functionality as
``RegularGrid1D``, but it automatically uses the relevant underlying
grid spacing.

The function

.. cpp:function:: T PiecewiseGrid1D::x(const int i) const;

returns a "physical" position on the grid given an index ``i``.

The function

.. cpp:function:: int PiecewiseGrid1D::index(const T x) const;

returns the index on the grid of a "physical" value ``x``.

The function

.. cpp:function:: T PiecewiseGrid1D::min() const;

returns the minimum value on the independent variable grid.

The function

.. cpp:function:: T PiecewiseGrid1D::max() const;

returns the maximum value on the independent variable grid.

The function

.. cpp:function:: T PiecewiseGrid1D::dx() const;

returns the grid spacing for the independent variable.

The function

.. cpp:function:: int PiecewiseGrid1D::nPoints() const;

returns the number of points in the independent variable grid.


Developer functionality
------------------------

For developers, additional functionality is available. Please consult
the code.
