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
verification gate is 2e-9 max-abs deviation.

Run with the ocp venv python (`~/venvs/ocp/bin/python`, has numpy).
