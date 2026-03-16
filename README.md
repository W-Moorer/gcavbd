# gcavbd

Minimal C++20 + CMake scaffold for a dynamics project.

## What is included

- Static library target: `gcavbd`
- Test target: `test_math`
- Dependencies via CMake FetchContent:
  - Eigen
  - GoogleTest

## Quick start

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
```
