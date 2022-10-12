//------------------------------------------------------------------------------
// © 2021-2022. Triad National Security, LLC. All rights reserved.  This
// program was produced under U.S. Government contract 89233218CNA000001
// for Los Alamos National Laboratory (LANL), which is operated by Triad
// National Security, LLC for the U.S.  Department of Energy/National
// Nuclear Security Administration. All rights in the program are
// reserved by Triad National Security, LLC, and the U.S. Department of
// Energy/National Nuclear Security Administration. The Government is
// granted for itself and others acting on its behalf a nonexclusive,
// paid-up, irrevocable worldwide license in this material to reproduce,
// prepare derivative works, distribute copies to the public, perform
// publicly and display publicly, and to permit others to do so.
//------------------------------------------------------------------------------
#include "module.hpp"

template <typename T>
py::class_<T> shifted_eos_class(py::module_ &m) {
  // define Shifted utility function
  m.def(
      "Shifted",
      [](typename T::BaseType eos, Real shift) { return T(std::move(eos), shift); },
      py::arg("eos"), py::arg("shift"));

  // define shifted class
  return eos_class<T>(m, T::EosPyType())
      .def(py::init<typename T::BaseType, Real>(), py::arg("eos"), py::arg("shift"));
}

template <typename... Ts>
void create_shifted(py::module_ &m, tl<Ts...>) {
  // C++14 workaround, since we don't have C++17 fold expressions
  auto l = {(shifted_eos_class<Ts>(m), 0)...};
}

void create_shifted_eos_classes(py::module_ &m) {
  create_shifted(m, singularity::shifted);
}
