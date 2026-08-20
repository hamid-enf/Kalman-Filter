#!/usr/bin/env bash
#
# Run the test suite under AddressSanitizer + UndefinedBehaviorSanitizer
# (plus -fsanitize=bounds and pointer-overflow) to catch buffer overruns,
# memory leaks, and undefined behaviour.
#
# Usage: ./tools/sanitize.sh [gcc|clang]
#   (defaults to gcc; clang is also supported)

set -eu

CC="${1:-gcc}"

FLAGS="-std=c11 -O1 -g -Wall -Wextra \
       -fsanitize=address,undefined,bounds,pointer-overflow \
       -fno-sanitize-recover=all"

echo "=== Building test suite with $CC under ASan+UBSan ==="

# Prefer the compiler's ASan flags; fall back if clang is not present.
if [ "$CC" = "clang" ] && ! command -v clang >/dev/null 2>&1; then
    echo "clang not found; using gcc"
    CC=gcc
fi

"$CC" $FLAGS -Ikalman/include -Itests \
    tests/test_framework.c \
    tests/test_matrix.c \
    tests/test_kf.c \
    tests/test_ekf.c \
    tests/test_ukf.c \
    tests/test_extensions.c \
    tests/test_reference.c \
    tests/test_boundary.c \
    tests/test_stress.c \
    tests/main.c \
    kalman/src/*.c \
    -lm -o /tmp/kalman_test_sanitized

echo "=== Running ==="
/tmp/kalman_test_sanitized
