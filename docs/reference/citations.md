# Research Paper Citation Trail

A running index of which published papers our code implements content from, which modules contain that content, and the commits that landed it. This is the authoritative map from paper → code for all non-classical work in the fillet pipeline.

Keep this in sync when adding new paper-derived modules. One paper = one section. Group sections by workstream.

## Clean-room policy

All paper-derived implementations in this codebase are **clean-room from papers only**. We do not read or copy reference implementations (e.g. the SurfLab Bitbucket repo for PnS, Peters' unpublished Maple scripts, or any other author-released source code). The papers themselves and their cited prior work are the sole inputs. This policy enables AGPL licensing of the output while respecting any non-commercial / share-alike restrictions on authors' own code.

---

## Karčiauskas–Peters bidegree-9 rectangular vertex blend (N ∈ {3, 4})

### K-P 2004 — Guided surface construction

**Citation:** Karčiauskas, K. & Peters, J. (2004). "Guided Subdivision and Geometric Constructions." (The 2004 guided-surface framework.)

**Files:** `/home/graeme/code/popkern/docs/papers/` — multiple 2004 talk / paper files; primary is the guided-subdivision framework that the bi-nonic construction extends.

**What's implemented:**
- Bi-nonic guided-surface cap (K-P 2004 §4) — 9×9 Gregory patches per sub-patch at N ∈ {3, 4} trihedral / square-pyramid corners
- Iterative C-twist optimization at central vertex
- Bi-quintic guide substrate (Option C extension)
- L-M solver for free guide CPs

**Modules:** `crates/popkern-fillet/src/construct/kp_binonic_optimizer.rs`, `kp_biquintic_optimizer.rs`, `kp_lm_solver.rs`, `kp_twist_optimizer.rs`; `crates/popkern-core/src/geometry/gregory/rectangular_binonic.rs`, `kp_guide_patch.rs`, `kp_guide_patch_biquintic.rs`, `kp_ring_build.rs`, `kp_row1.rs`, `kp_row2.rs`, `kp_interior_fit.rs`

**Key commits:** `657dc98` (H.6 dispatch), `b8466c3` (H.7 gate), `98a6148` (Option B twist-opt), `002bc27` (L-M phase 3), numerous H.4 sub-phase commits.

### Chiyokura-Kimura 1983 — Rectangular Gregory patches

**Citation:** Chiyokura, H. & Kimura, F. (1983). "Design of Solids with Free-Form Surfaces." *SIGGRAPH*.

**What's implemented:**
- Exact rational Gregory rectangular patch evaluation
- Twist-extraction from incident-surface normals at hex corners (commercial-grade #2)
- Parallelogram-rule inner CPs + tangent-plane projection

**Modules:** `crates/popkern-core/src/geometry/gregory/rectangular.rs`; `crates/popkern-fillet/src/blend/gregory_parallelogram.rs`, `gregory_inner_cp_shrink.rs`, `gregory_three_quad.rs`

### Shirman-Séquin 1987 — Vertex enclosure

**Citation:** Shirman, L. A. & Séquin, C. H. (1987). "Local surface interpolation with Bézier patches." *Computer Aided Geometric Design* 4(4).

**What's implemented:**
- Shirman-Séquin vertex enclosure at central point of three-quad decomposition
- Full G² vertex enclosure extension (commercial-grade #3)

**Modules:** `crates/popkern-fillet/src/blend/gregory_vertex_enclosure.rs`; central-point twist construction inside `kp_guide_ring.rs`

### Longhi 1987 — Triangular Gregory patches

**Citation:** Longhi, L. (1987). "Triangular Bézier patches for interpolation." (Interior twist-free triangular Gregory evaluator.)

**What's implemented:**
- Triangular Gregory evaluator with analytic partials

**Modules:** `crates/popkern-core/src/geometry/gregory/triangular.rs`

---

## Polar Karčiauskas–Peters 2013 / 2015 caps (N ≥ 5)

### KP15 — Bi-6 polar cap

**Citation:** Karčiauskas, K. & Peters, J. (2015). "Biquintic G² surfaces" (and the 2015 bi-6 polar construction for N ≥ 5).

**Local file:** `docs/papers/15bi6g2.pdf`

**What's implemented:**
- Full KP15 bi-6 polar pipeline: tensor-border → bi-6 cap → F_5 SPD minimisation → bi-9 degree elevation
- `ρ`, `β`, `γ` reparameterisations; `μ` coefficient tables (`µ̇, µ̈, ⃛µ, µ̃, µ̄, w₁₂₃`)
- BB rules S1–S4; weights for F_5 functional

**Modules:** `crates/popkern-core/src/geometry/gregory/kp15_bi6_reparam.rs`, `kp15_bi6_mu.rs`, `kp15_bi6_rules.rs`, `kp15_bi6_rule_s3.rs`, `kp15_bi6_cap.rs`, `kp15_bi6_cap_input.rs`, `kp15_bi6_polar_bridge.rs`, `kp15_bi6_tensor_border_extract.rs`; fillet wrapper in `crates/popkern-fillet/src/construct/kp_polar_bi6_apply.rs`

**Key commits:** I.1–I.10f sequence; ADR `docs/decisions/2026-04-20-phase-i-polar-kp-formula-choice.md`

**Research gap:** `n6` Maple script not available; deep-interior CP closed forms in the F_5 optimisation were worked around with the F_κ ladder (options A–F). Citation-chasing path noted in the polar-cap WIP notes.

### KP13 — Bi-5 macro-patch construction

**Citation:** Karčiauskas, K. & Peters, J. (2013). "G² Constructions from Biquintic Macro-patches." (Referenced as KP13 throughout.)

**Local file:** `docs/papers/kk12bi5.pdf` (preprint 2012, published 2013)

**What's implemented:**
- KP13 bi-5 macro-patch pipeline: guide g△ → g → auxiliary a → final p → bi-9 elevation
- `σ`, `ρ`, `β` reparameterisations (bi-5 variants)
- LSQ g△ fitting + post-fit eq (6) match at `m_k`
- I.14.a–e cap assembly (F_1 LSQ, eq 17 O piece, C² internal, eq 20 chain rule)
- 4N BRep subdivision (`poke_4n`) + env-gated dispatch

**Modules:** `crates/popkern-core/src/geometry/gregory/kp13_bi5_*.rs` (22 files); fillet wrappers in `crates/popkern-fillet/src/construct/kp_bi5_macro_apply.rs`, `kp13_bi5_promote.rs`, `poke_4n.rs`

**Key commits:** I.11–I.16 + deferred #1–#4 sequence; session-tail commits `d18f788..1688995`

**Research gaps:**
- β reparameterization §4.1 left at β=identity (eq 11/12/13 indexing ambiguous without Fig 8(b) visual confirmation).
- Upstream two-stripe promotion rejects N=5 folded seeds — pipeline reaches polar but needs `BlendBoundary → DirectBoundaryInput` extractor (polar-cap WIP).

### KP09a / KP 2009 — n ∈ {3, 4} special cases

**Citation:** Karčiauskas, K. & Peters, J. (2009). "Generalizing Bicubic Splines for Arbitrary Topology" and related 2009 polar-spline work.

**Local file:** `docs/papers/09tri8c2_smi.pdf`, `docs/papers/09guiSpline_cagd.pdf`

**What's implemented:**
- §6 pseudo-inverse path for circulant rank-deficient systems at n ∈ {3, 4}
- Polar boundary input validation for low valence

**Modules:** `crates/popkern-core/src/geometry/gregory/kp_polar_boundary_input.rs` (low-valence handling); Deferred #4 logic across KP13 circulant solvers

**Key commits:** `1031cb1` (Deferred #4)

### Stam 1998 — Catmull-Clark characteristic map

**Citation:** Stam, J. (1998). "Exact Evaluation of Catmull-Clark Subdivision Surfaces at Arbitrary Parameter Values." *SIGGRAPH*.

**What's implemented:**
- Full Stam-1998 CC χ (characteristic map) evaluation
- Used by both KP13 bi-5 and KP15 bi-6 in their LSQ / F_κ setups

**Modules:** `crates/popkern-core/src/geometry/gregory/kp13_bi5_stam_chi.rs`

**Key commits:** Deferred #2 (`6897e4b`)

### KP16 — Refinability reference

**Citation:** Karčiauskas, K. & Peters, J. (2021). "Least Degree G¹-refinable Multi-sided Surfaces Suitable for Inclusion into C¹ bi-2 Splines." *Computer-Aided Design* 130. (KP16 preprint label.)

**Local file:** `docs/papers/16gref.pdf`

**What's implemented:** Referenced for degree-bound theory only; no direct code.

---

## PnS (Polyhedral-Net Splines) n-gon caps

### KK-Peters 2015 — Smooth multi-sided blending of biquadratic splines

**Citation:** Karčiauskas, K. & Peters, J. (2015). "Smooth multi-sided blending of biquadratic splines." *Computers & Graphics* 46, 172–185. (Published 2015; preprint 2014 — local file uses preprint dating.)

**Local file:** `docs/papers/14bi234.pdf`

**What's implemented:**
- β reparameterization (eq 3 bi-4 variant, eq 3³ bi-3 n ∈ {3, 5} variant) + Appendix B defaults
- G¹ cap-cap constraints (eq 2 bi-4, eq 2³ bi-3)
- Central p_22 symmetric average (eq 4)
- Central-point weighted average `eop` (eq 5) + three `e_def` variants (bi-4, bi-3 n ∈ {3,5}, bi-3 2×2 Appendix C)
- σ central-reparameterization jet (Appendix B) — order-≤2 Taylor jet at pole
- Bi-3 vertex-enclosure twist solver (eq 7) — closed-form pair-sum circulant for n ∈ {3, 5}
- Bi-2 tensor-border type + b∘β composition producing cap's outer two layers

**Modules:** `crates/popkern-core/src/geometry/gregory/bernstein_bi3.rs`, `pns_beta.rs`, `pns_g1_constraints.rs`, `pns_center.rs`, `pns_sigma.rs`, `pns_bi3_twists.rs`, `pns_boundary.rs`, `pns_f_kappa.rs`, `pns_gram.rs`, `pns_central_quad.rs`

**Key commits:** `f224e9b..545f8c1` (12 commits, 67 tests)

**Paper coverage (14bi234.pdf):**
- ✅ Eq 1 (G¹ formal identity)
- ✅ Eq 2 / Eq 2³ (G¹ cap-cap constraints, with polynomial-identity verification)
- ✅ Eq 3 / Eq 3³ (β reparameterisation, bi-4 + bi-3 variants)
- ✅ Eq 4 (p_22 symmetric average)
- ✅ Eq 5 (central point eop) + all 3 `e_def` variants
- ✅ Eq 6 (F_κ minimisation framework — Gram matrix + `(n+5)×(n+5)` solver)
- ✅ Eq 7 (bi-3 vertex-enclosure twist solver, closed-form circulant)
- ✅ Appendix B (σ center jet + `d_1^def`, `d_2^def`)
- ✅ b ∘ β tensor-border composition
- ✅ `a ∘ σ` chain-rule central quadratic expansion
- ✅ F_κ functional evaluator (F_1 through F_4) via Gauss-Legendre quadrature

**Research gap:** The derived `σ_{01}, σ_{02}` BB coefs (from flip-symmetry-about-bisector) are a geometric interpretation of "flip-symmetry of each sector about its diagonal" rather than a direct paper transcription. Verified internally consistent by the σ rotational-symmetry tests (|ds_u|=|ds_v|, ds_uv along bisector), but paper doesn't explicitly state the image-space rotation. First thing to revisit if output shape looks off on a concrete fixture.

**Scope-explicit limitations (WIP-tracked future work):**
- Full cap assembly (integrates the primitives into `Vec<[[Point3; 5]; 5]>`).
- Bi-3 variant for n ∈ {3, 5} (parallel assembly using `pns_bi3_twists`).
- 2×2 bi-3 Appendix C variant for n > 5.
- BlendBoundary → Bi2TensorBorderRing CAD adapter.
- Dispatch wiring.
- σ's full BB-net beyond order-≤2 jet (only needed if F_κ optimisation over σ's free params is added later).

### ACM TOMS Algorithm 1032 — Polyhedral splines

**Citation:** Peters, J., Lo, K. S.-H., & Karčiauskas, K. (2023). "Algorithm 1032: Bi-cubic Splines for Polyhedral Control Nets." *ACM Transactions on Mathematical Software* 49(1), Article 7.

**Local file:** `docs/papers/3570158.pdf`

**Reference package:** `docs/code/ACM-ToMS-PnS-1032`

**Permission/provenance note:** The checked-in reference package remains documentation and validation material. The popkern implementation is a Rust-native architecture under `crates/popkern-core/src/geometry/polyhedral_spline`, using the paper's published polyhedral spline families and the ACM 1032 mask tables with explicit permission from the author to reimplement the code/tables. Do not copy the C++ architecture into production code; keep topology, diagnostics, and validation popkern-native.

**What's implemented:** Standalone Bernstein-Bezier mask backend for regular quads, extraordinary stars, n-gons, polar fans, T0, T1, T2, and two-triangle/two-quad layouts; degree elevation to bicubic; Bezier-NURBS conversion; finite-control validation; patch-atlas construction; sampled C0/tangent-plane seam validation including collapsed polar singular seams; interior Jacobian regularity gates; BRep-readiness metadata for patch domains, ordinary seams, boundary edges, and collapsed singular edge roles; open-shell BRep export with exact NURBS patch faces, shared seam topology, boundary-edge marking, linear pcurves on patch loops, face provenance, transactional no-mutation failure behavior, work estimates, and OBJ/debug summaries. `ControlNet` adapts the same local layout families. The commercial validation matrix covers the reference polar/eop/n-gon valence fixtures, T0/T1/T2 fixtures, scaled adversarial fixtures, BRep export, and debug artifact gates.

---

## Karčiauskas–Peters 2020 G² macro cap (setback vertex blends)

### KP20 — G² vertex-blend macro cap

**Citation:** Karčiauskas, K. & Peters, J. (2020). "High quality refinable G-splines for locally quad-dominant meshes with T-gons" / the 2020 G² refinable surface construction (referenced as KP20 throughout).

**Local file:** `docs/papers/20g2Ref.pdf`

**What's implemented:**
- Full §4.4 macro-cap assembly: 4n bi-9 pieces per n-gon cap on the setback border — the production cap for KP20-class vertex blends (`CapFamily` history: the setback surgery's watertight tiling)
- §4.3 planar τ family: appendix A fig 6c τ solves (F₃ energy, paper row sets, rotational variants, active-set solver)
- Characteristic cap net / Stam-eval ring (PR08 §5–6.1 + Stam 1998 inner rings)
- Piece maps with C²-stitch, orientation and domain validity reports
- eq (23) µ₀ scaling (κn/10⁵ — the misread that cost a campaign is recorded in the journal), S₀ solve, G² strip weights

**Modules:** `crates/popkern-core/src/geometry/gregory/kp20_bi6_*.rs` (assemble, cap, cnet, cc_ring, math, tau, tau_f3_solve, tau_solve_assemble, piece_maps, s0_solve, guide_build/sample/lsq, stam_eval, validate); fillet-side setback integration in `crates/popkern-fillet/src/construct/kp20_*.rs` and `setback_vertex_topology.rs`

### KP18 — Rapid degree-6 guide construction

**Citation:** Karčiauskas, K. & Peters, J. (2018). "Rapidly contracting subdivision yields finite, effectively C² surfaces" / the 2018 rapid guide construction KP20 delegates its guide to (referenced as KP18; KP20's [23]).

**Local file:** `docs/papers/18rapid.pdf`

**What's implemented:**
- §3–§4.2 total-degree-6 triangular guide: eq (1) join weights, guide fit consumed through `solve_kp20_guide_lsq` (the authoritative spoke-identifying encoding; the eq (2)–(5) side transcription with a mismatched index convention was deleted after audit)
- Guide sampling to bi-6/bi-9 corner blocks (`kp20_bi6_guide_sample`, `kp20_bi6_math` Hermite-to-BB tables)

**Modules:** `crates/popkern-core/src/geometry/gregory/kp20_bi6_guide_*.rs`, `kp20_bi6_math.rs`

---

## Optional / deferred PnS follow-ons

The following papers cover other PnS configurations. We haven't implemented these because our current CAD vertex-blend use case only needs the n-gon configuration. Adding them later would be additive extensions.

### KK-Peters 2020 — Polar caps

**Citation:** Karčiauskas, K. & Peters, J. (2020). "Smooth polar caps for locally quad-dominant meshes." *Computer Aided Geometric Design* 81, 101908.

**Local file:** `docs/papers/20polar.pdf`

**Would add:** Polar configuration (n triangles around a pole, surrounded by a quad ring) for cylinder/sphere cap blending.

### KK-Peters 2020 — Low-degree T-junction splines

**Citation:** Karčiauskas, K. & Peters, J. (2020). "Low degree splines for locally quad-dominant meshes." *Computer Aided Geometric Design* 83, 101934.

**Local file:** `docs/papers/20quaT.pdf`

**Would add:** T0 / T1 / T2 junction configurations for adaptive quad-strip merging.

### Gupta-Karčiauskas-Peters — PnS23 localized upgrade

**Citation:** Gupta, P., Karčiauskas, K., & Peters, J. "Localized '2-4-3' conversion of bi-2 C¹ splines to bi-3 C² splines." (SSRN preprint 6346314.)

**Local file:** `docs/papers/ssrn-6346314-1.pdf`

**Would add:** C² upgrade of PnS2's regular-grid bi-2 parts while preserving non-regular configuration patches.

---

## Classical foundations (no paper trail needed)

These are standard textbook / widely-known constructions; no specific paper implementation mapping.

- **de Boor** — B-spline theory, NURBS evaluation (`popkern-core::geometry::nurbs_*`).
- **de Casteljau** — Bezier subdivision (`popkern-core::geometry::bezier_*`).
- **Catmull-Clark 1978** — Subdivision framework (referenced by Stam 1998; our CC use is via Stam's exact-evaluation formula).
- **Bernstein / Bernstein-Bézier basis** — All `bernstein_*` modules.

---

## Maintenance

When a new paper-derived module lands:
1. Add a section under the appropriate workstream above.
2. Cite the paper with full bibliographic info, link to local PDF if present.
3. List the modules and key commits.
4. Note any research gaps / deferred work.

When citing prior work from within a module's header comment, point back to this file for the full citation.
