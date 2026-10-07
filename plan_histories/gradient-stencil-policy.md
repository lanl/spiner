# DataBox gradient stencil policy

Copyright 2026 Triad National Security, LLC. All rights reserved.
Generative AI was used to assist with this document.

Completed 2026-10-06.

## Plan and results

1. Inspect the user-added interpolation_stencils.hpp and in-progress result
   aliases/include changes. Preserve those changes and correct the stencil's
   leftover DataBox member references and missing result template arguments.
   Completed; the two static Linear methods use their supplied grids/data and
   return Result1D<T>/Result2D<T> with portable noexcept annotations.
2. Add GradStencil as the third DataBox template parameter, defaulting to
   InterpolationStencils::Linear. Propagate it through all member definitions,
   slices, copies, and the free getOnDeviceDataBox helper. Completed, including
   a GradStencilType alias. Existing one-/two-parameter instantiations retain
   their behavior; explicitly supplied Concept arguments move to position four.
3. Delegate both interpToRealWithGrads overloads after their existing asserts.
   Completed without changing value-only interpolation or ownership semantics.
4. Add a distinguishable custom policy regression test. Completed for 1D/2D
   dispatch, value-only behavior, slice types, shallow/deep copies, assignment,
   and portable execution on device copies. Existing all-grid tests remain.
5. Update the DataBox reference with the policy interface and parameter order.
   Completed. Format changed C++ regions and check whitespace. Completed.

## Validation

- cmake --build build-debug --parallel 6: passed.
- ctest --test-dir build-debug --output-on-failure -j 6: 22/22 passed.
- git diff --check: passed.
- The configured portability backend is Kokkos Serial. CUDA compilation and
  GPU execution were not tested.
