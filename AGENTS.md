# Agent Instructions

## Build and test

- Run builds with parallelism enabled: use `make -j` (or the equivalent parallel build option for the selected build system). Use 6 threads.

## Portability and memory spaces

- Spiner is a performance-portable library. Take care to preserve the distinction between host and device memory spaces in all code changes.
- Ensure allocations, accesses, copies, and execution contexts are valid for the memory space in which they occur. Do not assume that host and device pointers are interchangeable.

## Databox ownership and copy semantics

- Spiner `Databox` objects use pointer semantics: they behave like references and are shallow copies unless a deep copy is explicitly stated or requested.
- `Databox` objects are not reference counted and do not obey RAII. They must be explicitly freed when they are no longer needed, using the appropriate library deallocation mechanism.
- Track only the supported allocation states: trivial, allocated on device, allocated on host, empty, or unmanaged.
- Do not add RAII copy constructors, move constructors, or move-assignment semantics to `Databox` types or `Grid` types.

## Prior agentic work

- Machine-readable plans describing previous agentic efforts are available in the `plan_histories/` folder. Consult them when prior implementation context is relevant.
