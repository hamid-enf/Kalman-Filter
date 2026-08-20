#!/usr/bin/env bash
#
# MISRA-C:2012 check (cppcheck + official misra.py addon).
#
# Usage:
#   CPPCHECK=/path/to/cppcheck ./tools/misra_check.sh
#
# The script needs the cppcheck source tree (for the addon/misra.py file).
# Optionally set CPPCHECK_ADDON to point at a local copy of addon/misra.py.
#
# NOTE: cppcheck's MISRA support is *partial* and can report false positives;
# it is not a substitute for a certified MISRA checker (LDRA, Helix QAC,
# PC-lint Plus, Parasoft). Treat its output as guidance, and see
# docs/misra_compliance.md for the authoritative, reviewed statement.

set -u

CPPCHECK="${CPPCHECK:-cppcheck}"
ADDON="${CPPCHECK_ADDON:-/tmp/cppcheck-src/addons/misra.py}"
STD="${STD:-c11}"
INCDIR="kalman/include"

if ! command -v "$CPPCHECK" >/dev/null 2>&1; then
    echo "error: cppcheck not found (set CPPCHECK=...)" >&2
    exit 1
fi
if [ ! -f "$ADDON" ]; then
    echo "error: misra.py not found (set CPPCHECK_ADDON=...)" >&2
    exit 1
fi

srcs="kalman/src/kalman_matrix.c kalman/src/kalman_kf.c \
      kalman/src/kalman_ekf.c kalman/src/kalman_ukf.c kalman/src/kalman_common.c"

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

for f in $srcs; do
    "$CPPCHECK" --std="$STD" --dump -I"$INCDIR" "$f" >/dev/null 2>&1
    python3 "$ADDON" "${f}.dump" > "$tmpdir/$(basename "$f").misra" 2>/dev/null
    rm -f "${f}.dump"
done

echo "=== Per-file summary ==="
for f in "$tmpdir"/*.misra; do
    echo "--- $(basename "$f" .misra) ---"
    grep -A10 'MISRA rules violated' "$f" | grep -E 'misra-c2012-'
done
