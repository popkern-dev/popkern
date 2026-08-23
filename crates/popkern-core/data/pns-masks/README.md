# PnS mask tables (popkern, AGPL-3.0)

These 19 CSV tables map PnS input-net stencils to bi-3 BB-coefficients for
the polyhedral-spline cap families consumed by
`src/geometry/polyhedral_spline/tables.rs`. They are **popkern's own
tables**, derived clean-room from the published papers (PK-317); nothing in
this directory is copied from the ACM ToMS 1032 package, which remains
deliberately untracked (`docs/code/ACM-ToMS-PnS-1032`, see PK-315) and is
used only as a numerical verification oracle and for OBJ test fixtures.

Regenerate with:

```
~/venvs/ocp/bin/python tools/reference/pns-tables/emit_all.py
```

Generators (all under `tools/reference/pns-tables/`, with per-rule paper
citations inline):

| Tables | Generator | Source | Status vs reference package |
|---|---|---|---|
| polarSct3..8 | `polar.py` | 20polar.pdf closed forms | reproduces ≤ 7.1e-10 |
| T0, T1, T2 | `t_tables.py` | 20quaT.pdf | reproduces ≤ 8.1e-10 |
| ngonSct3, ngonSct5 | `eop_ngon.py` (`ngon_small`) | 14bi234.pdf §4.2 | reproduces ≤ 7.8e-10 |
| ngonSct6..8 | `eop_ngon.py` (`ngon_macro`) | 14bi234.pdf App. C | self-derived (see below) |
| eopSct3..8 | `eop_ngon.py` (`eop_small`/`eop_macro`) | composition | follows the ngon tables |

Two structural facts keep the star (`eop`) family from being independent
data: every eop table equals the corresponding ngon table composed with the
bi-2 virtual-refinement stencil (9,3,3,1)/16 (14bi234 Fig. 5/6a, App. A) —
verified against the reference package at ≤ 7.2e-10 for **all** valences.

## The macro family (n ≥ 6) is self-derived, not package-identical

The reference package's 2×2 bi-3 macro caps set their interior coefficients
by an unpublished guided-sampling pipeline (19bi34-style; TOMS 1032
describes the residual degrees of freedom as line-searched tuned weights).
Our tables instead follow the *published* recipe of 14bi234 Appendix C /
Fig. 23c: centre by eq. (5), tensor-border and G1 machinery as printed, and
F3 minimization for the remaining coefficients. Everything that carries a
continuity contract is identical to the reference (verified numerically):

- split borders + transversal layers under β = (3^{2×2}) with γ·c,
  γ = (n−2.3)/(2n−7.72) — 4.7e-11;
- sector G1 with weights (2c, γc) inner / (γc, 0) outer in per-sub-patch
  derivatives — ~2e-9;
- internal C1 junctions; centre by eq. (5), e = (0.02n+1.5)/(n+1.75);
- G2 at the centre via the quadratic-jet-through-σ model, with σ_ss pinned
  by the paper's assignment (2̄) and σ_st the single measured tuned scalar
  per valence (`SIGMA_Z` in `eop_ngon.py`).

Only shape-interior cells differ from the reference (≤ 4.3e-2 in mask
weight, ≤ 3.5e-3 after the eop composition). All popkern validation gates
(C0 seams, G1 sector continuity, finiteness, regularity) hold identically.
