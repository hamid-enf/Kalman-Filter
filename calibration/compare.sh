#!/usr/bin/env bash
# Compare the C and Python calibration outputs and check they agree to 5
# decimals (the C port runs in float, so a 1e-5 tolerance is expected).
#
# Usage: ./compare.sh
set -eu

cd "$(dirname "$0")"

gcc -O2 -w fusion.c \
    ../kalman/src/kalman_matrix.c ../kalman/src/kalman_kf.c \
    ../kalman/src/kalman_ekf.c ../kalman/src/kalman_ukf.c \
    ../kalman/src/kalman_common.c -I../kalman/include -lm -o /tmp/fusion_c

/tmp/fusion_c | tail -n +2 | awk '{print $1, $3, $5}' > /tmp/fusion_c.txt
PYTHONPATH=../python python3 fusion.py | tail -n +2 | awk '{print $1, $3, $5}' > /tmp/fusion_py.txt

echo "=== C ==="; cat /tmp/fusion_c.txt
echo "=== Python ==="; cat /tmp/fusion_py.txt

python3 - <<'PY'
c = [l.split() for l in open('/tmp/fusion_c.txt')]
p = [l.split() for l in open('/tmp/fusion_py.txt')]
assert len(c) == len(p) == 4, "row count mismatch"
ok = True
for (tc, c1, c2), (tp, p1, p2) in zip(c, p):
    for a, b in ((c1, p1), (c2, p2)):
        if abs(float(a) - float(b)) > 1e-4:
            ok = False
            print(f"MISMATCH at t={tc}: {a} vs {b}")
print("MATCH (within 1e-4)" if ok else "MISMATCH")
raise SystemExit(0 if ok else 1)
PY
