.. _interpolation:

Gridding for Interpolation
===========================

Spiner performs interpolation on uniform, Cartesian-product
grids. There is a lower-level object, ``RegularGrid1D`` which contains
the metadata required for these operations. ``RegularGrid1D`` has a
few useful userspace functions, which are described here.

Like ``DataBox``, the ``RegularGrid1D`` object is templated on
underlying data type, the default type being a ``Real`` as provided by
``ports-of-call``. You may wish to specialize to a specific type with
a type alias such as:

.. code-block:: cpp

   using RegularGrid1D = Spiner::RegularGrid1D<double>;

.. note::
   In the function signature below we refer to ``T`` and ``Real`` as
   the underlying arithmetic data type.

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

.. cpp:function:: int index(const T x) const;

returns the index on the grid of a "physical" value ``x``.

The function

.. cpp:function:: T min() const;

returns the minimum value on the independent variable grid.

The function

.. cpp:function:: T max() const;

returns the maximum value on the independent variable grid.

The function

.. cpp:function:: T dx() const;

returns the grid spacing for the independent variable.

The function

.. cpp:function:: int nPoints() const;

returns the number of points in the independent variable grid.

Developer functionality
^^^^^^^^^^^^^^^^^^^^^^^^

For developers, additional functionality is available. Please consult
the code.
