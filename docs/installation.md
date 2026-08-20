# Installation

The library is a set of plain C sources and headers. There is nothing to
"install" — you add the files to your project. No third-party libraries are
required.

## What to include

| Purpose | Files |
|---------|-------|
| All three filters | `kalman/src/kalman_matrix.c`, `kalman_kf.c`, `kalman_ekf.c`, `kalman_ukf.c`, `kalman_common.c` |
| Linear KF only | `kalman/src/kalman_matrix.c`, `kalman_kf.c`, `kalman_common.c` (and compile with `-DKF_ENABLE_EKF=0 -DKF_ENABLE_UKF=0`) |

Add `kalman/include` to the include path and `#include "kalman.h"` where
needed.

## STM32CubeIDE / CubeMX

1. Copy the `kalman/` directory into your project (e.g. `Core/`).
2. Project → Properties → C/C++ Build → Settings → MCU GCC Compiler →
   Include paths → add `../Core/kalman/include`.
3. Add the `.c` files under `kalman/src/` to the build (right-click →
   "Add to build", or add the folder to the source path).
4. `#include "kalman.h"` in `main.c`.

The core never includes STM32 headers, so it compiles identically in a plain C
project and in a CubeIDE project.

## Make / CMake

A top-level `Makefile` and `CMakeLists.txt` are provided for host development:

```sh
make            # static lib + examples + tests + benchmark
make test       # run tests
```

To use the library from your own Makefile, compile the sources into your build
and add `-I kalman/include`:

```makefile
CFLAGS += -I kalman/include
OBJS += kalman/src/kalman_matrix.o kalman/src/kalman_kf.o ...
```

## Compiler requirements

- C11 (or later) compiler. `kalman_config.h` asserts this at compile time.
- No C standard-library *runtime* requirements beyond `math.h` for `sqrt`
  (used by the matrix engine). There is no `printf`, no `malloc`, and no
  floating-point I/O in the core.

## First build check

To confirm the library builds in a minimal configuration, run:

```sh
make CFLAGS="-std=c11 -O2 -Wall -Wextra -DKF_ENABLE_EKF=0 -DKF_ENABLE_UKF=0" test
```

A clean build with `-Wall -Wextra -Werror` should produce no warnings.
