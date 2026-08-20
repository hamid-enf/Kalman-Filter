# Cross-language calibration

One **identical scenario**, implemented in the C library and its three ports
(Python, MATLAB, C#), so the numerical outputs can be compared side by side.

The scenario is the **sensor-fusion** example (odometry velocity + GNSS
position), because it exercises predict, update, two different `H`/`R` matrices
and sequential updates — the most code paths per line.

## The shared scenario (spec)

- **Model**: 2-state constant velocity, `x = [pos, vel]`,
  `F = [[1, dt],[0, 1]]`, `dt = 0.01`.
- **Truth**: cart moves at `2.0 m/s`; `true_pos += 2.0*dt` each step.
- **Noise** (identical in every language — this is what makes outputs
  comparable):
  - PRNG: **xorshift32**, seed `0xC0FFEE`.
  - uniform `u ∈ [-1, 1]` from `(int32)(state % 2000001) - 1000000` / 1e6.
  - gaussian ≈ `(u1 + u2 + u3) / 3`.
- **Sensors**:
  - *Odometry* (100 Hz): `vel + 0.2*gauss`, updates **velocity** with
    `H = [0, 1]`, `R = 0.04`.
  - *GNSS* (1 Hz, every 100 steps): `true_pos + 2.0*gauss`, updates
    **position** with `H = [1, 0]`, `R = 4.0`.
- **Process noise**: `Q = [[0,0],[0,0.05]]`, `P0 = diag(10, 10)`.
- **Run**: 400 fast steps (4 s). Print the fused position and velocity at
  `t = 1, 2, 3, 4` seconds.

## Expected output

Every implementation must print the **same four rows** (to within float
rounding). The verified reference (C and Python agree to 5 decimals):

```
  t(s) | fused pos | fused vel
  0.99 |   1.20936 |   1.92865
  1.99 |   3.66873 |   2.05031
  2.99 |   6.22010 |   2.05000
  3.99 |   8.30493 |   1.99125
```

(`t` is printed as `step*dt`, i.e. the GNSS fixes land just before the whole
seconds.) The velocity converges to ~2.0 m/s and the position tracks the true
trajectory, both drift-free thanks to the GNSS correction.

## Files

| File | Language | Run |
|------|----------|-----|
| `fusion.c` | C (the library) | `gcc fusion.c ../kalman/src/kalman_*.c -I../kalman/include -lm && ./a.out` |
| `fusion.py` | Python | `PYTHONPATH=../python python3 fusion.py` |
| `fusion.m` | MATLAB | `addpath('../matlab'); fusion;` |
| `fusion.cs` | C# (.NET) | `dotnet run --project ../csharp/Calibration` |

Run any two and diff their output — the numbers must agree to ~5 decimal
places. `compare.sh` automates the C↔Python comparison (tolerance 1e-4: the C
library runs in `float`, so a ~1e-5 accumulation over 400 steps is expected
and correct).
