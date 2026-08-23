# Clean-room PnS mask table generators (PK-317)

Regenerates the 19 polyhedral-net-spline mask tables that
`popkern-core`'s PnS machinery consumes, from the PUBLISHED
constructions — so the repository can ship self-derived tables with
no ACM-package content.

Sources of truth (papers in `docs/papers/`):

- `20polar.pdf`  — Karciauskas & Peters, "Smooth polar caps for
  locally quad-dominant meshes", CAGD 81 (2020): `polar.py`.
- `20quaT.pdf`   — K&P, "Low degree splines for locally quad-dominant
  meshes", CAGD 83 (2020): T0/T1/T2 (GT-splines).
- `14bi234.pdf`  — K&P, "Smooth multi-sided blending of biquadratic
  splines", C&G 46 (2015): star (eop) + n-gon caps.

Clean-room protocol: the generators implement the papers' formulas;
the ACM ToMS 1032 package tables (deliberately untracked,
`docs/code/ACM-ToMS-PnS-1032/`) are used ONLY as verification
oracles. The package tables carry ~5e-10 of their own numeric noise
(numeric solves printed at 10 digits; e.g. polarSct3 prints
0.1666666661 where the closed form is exactly 1/6), so the
verification gate is 2e-9 max-abs deviation (5e-9 for the composed
eop tables, which stack two tables' noise).

Reproduction status: polar (6 tables), T0/T1/T2 and the n in {3,5}
star/ngon caps reproduce the reference package within the gate. The
2x2 bi-3 macro caps (star/ngon, n in {6,7,8}) are SELF-DERIVED per
14bi234 Appendix C / Fig. 23c: the reference package fills their
interior coefficients from an unpublished guided-sampling pipeline,
so interiors differ by up to ~4e-2 while every continuity-bearing
quantity (borders, G1 machinery, centre, central G2 jet) matches the
reference exactly. See `eop_ngon.py` docstrings and
`crates/popkern-core/data/pns-masks/README.md`.

`emit_all.py` regenerates and writes all 19 CSVs into
`crates/popkern-core/data/pns-masks/` (the tables `tables.rs` embeds).

Run with the ocp venv python (`~/venvs/ocp/bin/python`, has numpy).
