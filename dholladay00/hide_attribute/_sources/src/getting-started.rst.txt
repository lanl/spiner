.. _getting-started:

Getting Started
================

The following provides a simple example of utilizing a ``DataBox``.

.. code-block:: cpp

  #include <iostream>
  #include <databox.hpp>
  using namespace Spiner;

  int main() {
    // create a databox
    constexpr int NX1 = 2;
    constexpr int NX2 = 3;
    constexpr int NX3 = 4;
    DataBox db(NX3, NX2, NX1);
    
    // fill the databox with the flat index of each element
    for (int i = 0; i < db.size(); ++i) {
      db(i) = static_cast<double>(i);
    }
    
    // set the interpolation ranges to [0,1] or each dimension
    for (int d = 0; d < db.rank(); ++d) {
      db.setRange(d, 0, 1, db.dim(d));
    }
    
    // interpolate
    double val = db.interpToReal(0.2, 0.3, 0.4);
    
    // save to file
    db.saveHDF("my_data.sp5");
    
    // load a new databox from file
    DataBox db2;
    db2.loadHDF("my_data.sp5");
    
    // interpolate new databox to the same location
    double val2 = db2.itnerpToReal(0.2, 0.3, 0.4);
    
    // print the interpolated values and see they're the same
    std::cout << val1 << ", " val2 << ": " << (val1 - val2) << std::endl;
    
    // free the databoxes
    free(db);
    free(db2);

    return 0;
  }

For more examples, please consult the test directory.
