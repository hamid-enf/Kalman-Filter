# MISRA-C:2012 compliance

## What this means here — read this first

The library is **written to MISRA-C:2012**, and the status below is **verified
with a real static-analysis pass** (cppcheck 2.22 + the official `misra.py`
addon), not just a claim. However, two honest caveats apply:

1. **cppcheck's MISRA support is partial.** It covers a useful subset of rules
   and occasionally reports false positives. It is *not* a certified MISRA
   checker (LDRA, Helix QAC, PC-lint Plus, Parasoft C/C++test). Full formal
   certification requires one of those tools plus a compliance matrix agreed
   with the project.
2. **Advisory rules are handled as documented deviations**, not all suppressed
   to zero. This is standard practice: most MISRA projects deviate from a
   documented set of advisory rules with rationale.

The actionable, reproducible result: **every Mandatory and Required rule that
cppcheck checks is clean; the remaining findings are Advisory and are either
false positives or justified deviations, listed below.**

## How to reproduce

```sh
CPPCHECK=/path/to/cppcheck ./tools/misra_check.sh
```

The script dumps each source with `cppcheck --dump` and runs the `misra.py`
addon. (The MISRA rule texts are copyrighted and not distributed with the
addon, so messages show rule numbers only — which is all that is needed here.)

## Results (cppcheck 2.22, `misra.py`, C11)

| Rule | Category | Status | Notes |
|------|----------|--------|-------|
| **10.8** | Required | ✅ Compliant | Composite expressions are not cast to a different essential type |
| **11.8** | Required | ✅ Compliant | No cast removes `const`/`volatile` |
| **14.2** | Required | ✅ Compliant | All `for` loops well-formed |
| **17.7** | Required | ✅ Compliant | Every non-void return value is used or explicitly `(void)`-cast |
| 12.1 | Advisory | ⚠️ Deviation | See below |
| 12.3 | Advisory | ⚠️ False positive | No comma operators exist; addon flags declarator lists |
| 15.5 | Advisory | ⚠️ Deviation | See below |
| 20.10 | Advisory | ✅ Compliant | No `#`/`##` in shipped macros |

## Documented deviations

### Rule 15.5 — "A function should have a single point of exit" (Advisory)

The library returns early on error throughout:

```c
kf_status_t kf_kf_predict(kf_kf_t *kf, ...)
{
    if (kf == NULL) { return KF_ERROR_NULL_POINTER; }   /* early return */
    ...
    return KF_OK;
}
```

**Rationale:** early-return is idiomatic, reduces nesting, and makes the error
path explicit at the point of failure. Refactoring ~300 sites to a single exit
would increase code size and indentation and add risk, for no functional
benefit. This deviation is deliberate and applied uniformly. (It is one of the
most commonly deviated advisory rules in production embedded code.)

### Rule 12.1 — "Operator precedence should be made explicit" (Advisory)

The config macros were given explicit parentheses; a handful of arithmetic
expressions in the matrix engine (e.g. `sum += a * b;`) rely on standard
operator precedence.

**Rationale:** these expressions are short, use only `*`/`+`/`-` on
`kf_real_t`, and are unambiguous. Teams requiring zero advisory findings can
add parentheses mechanically; no correctness issue is masked.

### Rule 12.3 — "The comma operator shall not be used" (Advisory)

The code contains **no comma operators**. cppcheck's addon flags comma-
separated declarator lists (e.g. `uint16_t a, b;`); block-scope declarations
were split into single-identifier form (also satisfying advisory 8.9), and the
remaining reports are function-parameter lists, which are not comma operators.

## Rules not checkable by cppcheck

Several MISRA-C:2012 rules require information cppcheck does not track (e.g.
Directive 4.x documentation rules, Rule 21.x standard-library requirements).
These are addressed by design and documented here:

- **Dir 4.7 / 4.10** (error handling, header guards): every fallible function
  returns `kf_status_t`; every header has an include guard.
- **Rule 21.2** (reserved identifiers): no identifiers begin with `__` or `_` +
  uppercase; the `kf_` prefix is used consistently.
- **Rule 21.3/21.5/21.6** (library functions): the only `<math.h>` function used
  in the core is `sqrt`/`sqrtf`/`isfinite` (via the `kf_sqrt`/`kf_isfinite`
  wrappers); `<stdint.h>`/`<stddef.h>` for explicit types. No `printf`, no
  allocation, no `string.h` in the core.
- **Rule 21.7/21.8** (stdlib functions): not used in the core.
- **Rule 8.4** (external linkage declarations): every externally-linked
  function is declared in a header included by the defining translation unit.

## Compliance summary

| Category | Result |
|----------|--------|
| Mandatory rules | ✅ No findings (design satisfies them; cppcheck reports none) |
| Required rules checkable by cppcheck | ✅ 0 findings |
| Advisory rules | ⚠️ 2 documented deviations (15.5, 12.1) + 1 false positive (12.3) |

**Bottom line:** the core is free of every Mandatory/Required-rule violation
that the available tool can detect, and the advisory deviations are explicit,
uniform, and justified. For a safety-certified product, run a certified MISRA
checker against this same code and adopt this document as the starting
compliance matrix.
