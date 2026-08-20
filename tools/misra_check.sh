#!/usr/bin/env bash
#
# MISRA-C:2012 check (cppcheck + official misra.py addon).
#
# Usage:
#   ./tools/misra_check.sh
#
# Optionally override:
#   CPPCHECK=<path to cppcheck binary>
#   CPPCHECK_ADDON=<path to addon/misra.py>
#
# The script tries to auto-locate both (cppcheck on PATH, and misra.py under
# the common cppcheck data directories or a source checkout).
#
# NOTE: cppcheck's MISRA support is *partial* and can report false positives;
# it is not a substitute for a certified MISRA checker (LDRA, Helix QAC,
# PC-lint Plus, Parasoft). Treat its output as guidance, and see
# docs/misra_compliance.md for the authoritative, reviewed statement.

set -u

CPPCHECK="${CPPCHECK:-cppcheck}"
ADDON="${CPPCHECK_ADDON:-}"
STD="${STD:-c11}"
INCDIR="kalman/include"

if ! command -v "$CPPCHECK" >/dev/null 2>&1; then
    echo "error: cppcheck not found (set CPPCHECK=...)" >&2
    exit 1
fi

# Auto-locate misra.py if not given.
if [ -z "$ADDON" ] || [ ! -f "$ADDON" ]; then
    for cand in \
        /usr/share/cppcheck/addons/misra.py \
        /usr/local/share/cppcheck/addons/misra.py \
        /tmp/cppcheck-src/addons/misra.py \
        "$(dirname "$0")/../cppcheck/addons/misra.py"; do
        if [ -f "$cand" ]; then ADDON="$cand"; break; fi
    done
fi
if [ -z "$ADDON" ] || [ ! -f "$ADDON" ]; then
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
