# Interpolation gradient documentation

Copyright 2026 Triad National Security, LLC. All rights reserved.
Generative AI was used to assist with this document.

## Plan and completion

1. Locate the existing DataBox and grid interpolation references. Complete.
2. Document both interpToRealWithGrads overloads, aggregate fields, usage,
   derivative coordinate order, knot/extrapolation semantics, and execution
   space requirements. Complete in doc/sphinx/src/databox.rst.
3. Document weightsWithGrad for all built-in grids and the custom-grid
   interface. Complete in doc/sphinx/src/interpolation.rst.
4. Verify documentation formatting and the worked example. Complete:
   git diff --check passed; the extracted 2D example compiled with C++20
   and the serial portability strategy and ran with assertions checking all
   three documented results. An initial C++17 compile failed because the
   current library/dependencies require C++20.

Sphinx and docutils are unavailable in the environment, so rendered HTML
was not built. No library code changed.
