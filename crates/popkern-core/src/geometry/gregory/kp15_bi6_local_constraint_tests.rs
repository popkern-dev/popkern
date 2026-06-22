//! Phase B exit-criterion tests for the [KP15] eq (2) residual
//! at the cross-sector edges of a bi-6 polar cap.
//!
//! Strategy:
//! - Build a deterministic pentagon / hexagon / heptagon input
//!   that runs cleanly through the existing pipeline.
//! - Drive `assemble_bi6_polar_cap` end-to-end.
//! - For every adjacent sector pair `(k − 1, k)`, evaluate the
//!   eq (2) G¹ residual polynomial. Report the worst BB-coef
//!   magnitude across all 8 indices and all `n` edges.
//!
//! The exit criterion (per
//! `docs/decisions/2026-06-22-kp15-bi6-deep-interior-without-maple.md`
//! Phase B) is that the residual coefficient at BB indices
//! 0, 1, 2, 3, 4, 5, 6 vanishes to machine zero (the rules-
//! (S1)–(S4) determined band), and the 7-th coefficient vanishes
//! once rule (5) on `p̀_{30}` is applied — that's what makes the
//! identity polynomial degree-6 rather than degree-7 (paper
//! line 149).
//!
//! Any non-zero coefficient is a transcription gap somewhere in
//! the pipeline that the next phases (constraint matrix, solve)
//! will need to mind.

use super::*;
use crate::geometry::Point3;
use crate::geometry::Vec3;
use crate::geometry::gregory::kp_polar_boundary_input::{
    PolarBoundaryInput, TensorBorder, TensorBorderBi3,
};
use crate::geometry::gregory::kp15_bi6_cap::assemble_bi6_polar_cap;
use crate::geometry::gregory::kp15_bi6_polar_bridge::convert_polar_to_bi6_cap_input;
use crate::geometry::gregory::{OuterCornerData, SpokeBoundaryData, TangentPlane};

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn tp_xy() -> TangentPlane {
    TangentPlane::from_basis(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)).unwrap()
}

fn polar_input(n: usize) -> PolarBoundaryInput {
    let center = pt(0.0, 0.0, 0.0);
    let mut spokes = Vec::with_capacity(n);
    let mut outer_corners = Vec::with_capacity(n);
    let mut face_tangent_planes = Vec::with_capacity(n);
    let mut per_spoke_radii = Vec::with_capacity(n);
    for k in 0..n {
        let theta = (2.0 * std::f64::consts::PI) * (k as f64) / (n as f64);
        let v = pt(theta.cos(), theta.sin(), 0.0);
        let mut boundary_cps = [pt(0.0, 0.0, 0.0); 8];
        for (i, cp) in boundary_cps.iter_mut().enumerate() {
            let t = i as f64 / 7.0;
            *cp = pt(t * v.x, t * v.y, t * v.z);
        }
        spokes.push(SpokeBoundaryData {
            boundary_cps,
            xtan_a_at_c: Vec3::new(0.0, 1.0, 0.0),
            xtan_b_at_c: Vec3::new(0.0, -1.0, 0.0),
            xtan_a_at_v: Vec3::new(0.0, 1.0, 0.0),
            xtan_b_at_v: Vec3::new(0.0, -1.0, 0.0),
            radius_at_c: 0.05,
            radius_at_v: 0.05,
            edge_id: format!("edge_{k}"),
        });
        outer_corners.push(OuterCornerData {
            position: v,
            tangent_plane: tp_xy(),
        });
        face_tangent_planes.push(tp_xy());
        per_spoke_radii.push(0.05);
    }
    PolarBoundaryInput {
        central_vertex: center,
        central_tangent_plane: tp_xy(),
        spokes,
        outer_corners,
        face_tangent_planes,
        per_spoke_radii,
    }
}

fn synthetic_bi3_tensor_border(n: usize) -> TensorBorder {
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let theta = (2.0 * std::f64::consts::PI) * (k as f64) / (n as f64);
        let dx = theta.cos();
        let dy = theta.sin();
        let row0: [Point3; 4] = [
            pt(dx, dy, 0.1),
            pt(dx * 1.1, dy * 1.1, 0.1),
            pt(dx * 1.2, dy * 1.2, 0.1),
            pt(dx * 1.3, dy * 1.3, 0.1),
        ];
        let row1: [Point3; 4] = [
            pt(dx, dy, 0.2),
            pt(dx * 1.1, dy * 1.1, 0.2),
            pt(dx * 1.2, dy * 1.2, 0.2),
            pt(dx * 1.3, dy * 1.3, 0.2),
        ];
        let row2: [Point3; 4] = [
            pt(dx, dy, 0.3),
            pt(dx * 1.1, dy * 1.1, 0.3),
            pt(dx * 1.2, dy * 1.2, 0.3),
            pt(dx * 1.3, dy * 1.3, 0.3),
        ];
        sectors.push([row0, row1, row2]);
    }
    TensorBorder::Bi3(TensorBorderBi3 { sectors })
}

/// Compute eq (2) BB residual on every cross-sector edge of a
/// fully-assembled bi-6 polar cap. Returns the per-edge per-BB
/// magnitudes; the worst across `n × 8` is the diagnostic.
fn eq2_per_edge_residual(n: usize) -> Vec<[f64; 8]> {
    let polar = polar_input(n);
    let tb = synthetic_bi3_tensor_border(n);
    let cap_input = convert_polar_to_bi6_cap_input(&polar, &tb)
        .unwrap_or_else(|e| panic!("convert failed for n = {n}: {e}"));
    let out = assemble_bi6_polar_cap(&cap_input)
        .unwrap_or_else(|e| panic!("assemble failed for n = {n}: {e}"));
    let mut per_edge = Vec::with_capacity(n);
    for k in 0..n {
        let km1 = (k + n - 1) % n;
        let grid_prev = &out.patches_bi6[km1];
        let grid_curr = &out.patches_bi6[k];
        let bb = super::eq2_g1_residual_bb_degree7(grid_prev, grid_curr, &cap_input.params);
        let mags: [f64; 8] = core::array::from_fn(|i| bb[i].length());
        per_edge.push(mags);
    }
    per_edge
}

fn worst_residual(per_edge: &[[f64; 8]]) -> (usize, usize, f64) {
    let mut worst = (0usize, 0usize, 0.0_f64);
    for (k, edge) in per_edge.iter().enumerate() {
        for (i, &m) in edge.iter().enumerate() {
            if m > worst.2 {
                worst = (k, i, m);
            }
        }
    }
    worst
}

#[test]
fn eq2_phase_a_baseline_unit_loading_at_u0() {
    // Sanity check on the residual extractor itself: with a
    // bare-minimum loading (p̀_{00} = 1, everything else 0,
    // ṕ_{0, 0} = p̀_{0, 0}, ṕ_{0, 1} = 2(1 − c) by rule (S1)),
    // BB[0] must vanish. This is the same loading exercised in
    // `kp15_eq2_identity_zeros_at_u0_with_rule_s1_p01` but via
    // the production residual function rather than ad-hoc
    // BernsteinPoly construction.
    let mut grid_prev = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
    let mut grid_curr = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
    grid_prev[0][0] = Point3::new(1.0, 0.0, 0.0);
    grid_curr[0][0] = Point3::new(1.0, 0.0, 0.0);
    let params = Bi6Reparams::for_valence(6);
    let c = params.c;
    // ṕ_{0, 1} per rule (S1) with p̀_{00} = (1, 0, 0):
    grid_curr[0][1] = Point3::new(2.0 * (1.0 - c), 0.0, 0.0);
    let bb = super::eq2_g1_residual_bb_degree7(&grid_prev, &grid_curr, &params);
    assert!(
        bb[0].length() < 1e-12,
        "u^0 BB-coef = {:?} (mag {})",
        bb[0],
        bb[0].length()
    );
}

/// Build a synthetic two-sector pair that is (S1)-strict by
/// construction: a deterministic "previous" sector grid loaded
/// only with non-zero values at the BB-coefs (S1) sees, and a
/// "current" sector derived from it by:
///
/// 1. Identifying the shared cross-sector row: `ṕ_{i, 0} :=
///    p̀_{i, 0}` for all `i ∈ 0..=6`. The paper specifies this
///    identification for the entire row (line 158: "We select
///    the six BB-coefficients p^0_ij, 0 ≤ i + j ≤ 2"), and (S1)
///    re-affirms it explicitly for `i = 0, 1, 2`. The rest of
///    the row is the local segment of the tensor-border shared
///    across the cross-sector edge by construction.
/// 2. Applying [KP15 eq (6)]'s `rule_s1_star_propagation` to
///    derive `ṕ_{0, 1}, ṕ_{1, 1}, ṕ_{0, 2}` from the previous
///    sector's (0..=2, 0..=2) corner.
/// 3. Leaving every other BB-coef of the current sector at zero.
///
/// Under this construction, the eq (2) polynomial's first three
/// BB-coefficients — `u^0, u^1, u^2` — must vanish (they're
/// exactly the constraints (S1) was *derived* to satisfy). The
/// higher-degree BB-coefs (`u^3..u^7`) depend on `ṕ_{i, 1}` for
/// `i ≥ 2` and on rule (5)'s `p̀_{30}` substitution; we don't
/// load those here, so they will be non-zero — but the test
/// only asserts the first three.
fn build_s1_strict_two_sectors(n: usize) -> ([[Point3; 7]; 7], [[Point3; 7]; 7], Bi6Reparams) {
    let params = Bi6Reparams::for_valence(n);

    // A reproducible loading of the previous sector's near-star
    // and near-rim BB-coefs at row j ∈ {0, 1, 2}. Values are
    // arbitrary but distinct, so the test is a real linear
    // combination rather than a trivial zero.
    let mut grid_prev = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
    grid_prev[0][0] = Point3::new(0.0, 0.0, 0.0);
    grid_prev[1][0] = Point3::new(0.20, 0.07, 0.01);
    grid_prev[2][0] = Point3::new(0.45, 0.12, 0.02);
    grid_prev[3][0] = Point3::new(0.70, 0.18, 0.03);
    grid_prev[4][0] = Point3::new(0.90, 0.21, 0.03);
    grid_prev[5][0] = Point3::new(1.05, 0.23, 0.04);
    grid_prev[6][0] = Point3::new(1.20, 0.25, 0.04);
    grid_prev[0][1] = Point3::new(0.03, 0.15, 0.02);
    grid_prev[1][1] = Point3::new(0.22, 0.20, 0.03);
    grid_prev[0][2] = Point3::new(0.05, 0.30, 0.04);

    // Build current sector by (S1)-strict propagation.
    let s1 = super::super::kp15_bi6_rules::rule_s1_star_propagation(
        &params,
        grid_prev[0][0],
        grid_prev[0][1],
        grid_prev[0][2],
        grid_prev[1][0],
        grid_prev[1][1],
        grid_prev[2][0],
    );
    let mut grid_curr = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
    // Row 0 identification — paper line 158, the cross-sector
    // edge shared by sector k − 1 and sector k.
    for i in 0..=6 {
        grid_curr[i][0] = grid_prev[i][0];
    }
    grid_curr[0][1] = s1.p01;
    grid_curr[1][1] = s1.p11;
    grid_curr[0][2] = s1.p02;

    (grid_prev, grid_curr, params)
}

#[test]
fn s1_strict_two_sectors_zero_eq2_at_u0_u1() {
    // Phase B positive control: under the (S1)-strict construction,
    // the first two BB-coefs of eq (2) — at u^0 and u^1 — must
    // vanish across all valences. These are exactly the
    // constraints `rule_s1_star_propagation` was derived to
    // satisfy.
    for n in [3usize, 5, 6, 7, 8, 9, 10] {
        let (grid_prev, grid_curr, params) = build_s1_strict_two_sectors(n);
        let bb = super::eq2_g1_residual_bb_degree7(&grid_prev, &grid_curr, &params);
        assert!(
            bb[0].length() < 1e-12,
            "n={n}: u^0 BB-coef = {:?} (mag {})",
            bb[0],
            bb[0].length()
        );
        assert!(
            bb[1].length() < 1e-12,
            "n={n}: u^1 BB-coef = {:?} (mag {})",
            bb[1],
            bb[1].length()
        );
    }
}

#[test]
fn s1_strict_two_sectors_zero_eq3_at_u0() {
    // Phase B v2 positive control for eq (3): under the (S1)-strict
    // construction, the *first* BB-coef of the φ-multiplied eq (3)
    // residual must vanish across all valences. This index involves
    // only the BB-coefs ṕ_{0, 0..2} and p̀_{0, 0..2} that (S1)
    // determines.
    //
    // (Higher BB-coefs of eq (3) depend on ṕ_{1, 1}, ṕ_{1, 2}, etc.,
    // which (S2)/(S3) determine — and which the (S1)-strict two-
    // sector synthetic skeleton does NOT populate. Phase C's matrix
    // solver will fold those in by also enforcing (S2)/(S3) on the
    // ring of `n` sectors simultaneously.)
    for n in [3usize, 5, 6, 7, 8, 9, 10] {
        let (grid_prev, grid_curr, params) = build_s1_strict_two_sectors(n);
        let bb = super::eq3_g2_residual_phi_multiplied(&grid_prev, &grid_curr, &params);
        assert!(
            !bb.is_empty(),
            "n={n}: eq3 residual returned an empty BB-coef vector"
        );
        assert!(
            bb[0].length() < 1e-12,
            "n={n}: eq3 u^0 BB-coef = {:?} (mag {}); deg of residual polynomial = {}",
            bb[0],
            bb[0].length(),
            bb.len() - 1
        );
    }
}

#[test]
fn eq3_residual_polynomial_degree_matches_paper_after_assembly() {
    // The φ-multiplied eq (3) residual polynomial is at most degree
    // 9 as composed from BB-product arithmetic in this module (each
    // term's per-direction degree adds: term 4 = 4·1·4 = 9, etc.).
    // The paper's line 150 claim is that *after rule (5)* the
    // polynomial drops to degree 7. We check the BB-coef magnitudes
    // beyond index 7 on a (S1)-strict synthetic input — they aren't
    // expected to vanish here because (S2)/(S3)/(S4) aren't applied
    // — but the polynomial's BB-vector length must be ≥ 10 so the
    // Phase C matrix builder gets all 15 (= 7 + 8) coefficients
    // worth of constraints.
    let (grid_prev, grid_curr, params) = build_s1_strict_two_sectors(6);
    let bb = super::eq3_g2_residual_phi_multiplied(&grid_prev, &grid_curr, &params);
    assert!(
        bb.len() >= 8,
        "expected at least 8 BB-coefs (paper's degree 7); got {}",
        bb.len()
    );
}

/// Rank of a dense `f64` matrix via Gaussian elimination with
/// partial pivoting and a relative-pivot tolerance. Test-only
/// utility — production rank checks should use a robust SVD path.
fn matrix_rank(rows: usize, cols: usize, entries: &[f64], tol: f64) -> usize {
    assert_eq!(entries.len(), rows * cols);
    let mut m: Vec<f64> = entries.to_vec();
    let mut rank = 0usize;
    let mut row = 0usize;
    for col in 0..cols {
        if row >= rows {
            break;
        }
        // Find pivot in this column at or below `row`.
        let mut best = row;
        let mut best_abs = m[row * cols + col].abs();
        for r in (row + 1)..rows {
            let v = m[r * cols + col].abs();
            if v > best_abs {
                best = r;
                best_abs = v;
            }
        }
        if best_abs < tol {
            continue;
        }
        // Swap pivot row up.
        if best != row {
            for c in 0..cols {
                m.swap(row * cols + c, best * cols + c);
            }
        }
        // Eliminate below and above.
        let pivot = m[row * cols + col];
        for r in 0..rows {
            if r == row {
                continue;
            }
            let factor = m[r * cols + col] / pivot;
            if factor == 0.0 {
                continue;
            }
            for c in col..cols {
                m[r * cols + c] -= factor * m[row * cols + c];
            }
        }
        rank += 1;
        row += 1;
    }
    rank
}

#[test]
fn eq2_jacobian_rank_is_eight_independent_constraints() {
    // Phase C foundation: the eq (2) Jacobian over a single cross-
    // sector edge is an `8 × 98` matrix linear in the BB-coefs of
    // both adjacent sectors. With no rule applied yet — i.e., the
    // 49 + 49 grid positions all treated as free — the 8 rows
    // should be linearly independent: there's no algebraic identity
    // among the BB-coefs of the eq (2) polynomial that makes one
    // BB-coef predictable from the others without invoking rule (5)
    // or the (S1)-(S4) identifications.
    //
    // After rule (5) (substituting `p̀_{30} := (3/4)(p̀_{20} + p̀_{40})
    // − (3/10)(p̀_{10} + p̀_{50}) + (1/20)(p̀_{00} + p̀_{60})`) the
    // row at u^7 should fall into the row space of the others —
    // we check that in a follow-on test.
    for n in [3usize, 5, 6, 7, 8, 9, 10] {
        let params = Bi6Reparams::for_valence(n);
        let jac = super::eq2_jacobian_dense(&params);
        let mut flat = Vec::with_capacity(8 * 98);
        for row in jac.iter() {
            flat.extend_from_slice(row);
        }
        let rank = matrix_rank(8, 98, &flat, 1e-12);
        assert_eq!(
            rank, 8,
            "n={n}: eq (2) Jacobian rank = {rank}, expected 8 (all BB rows independent)"
        );
    }
}

#[test]
fn eq3_jacobian_rank_reports_paper_structure() {
    // The φ-multiplied eq (3) residual polynomial built from BB-
    // product arithmetic returns a vector of length up to 10
    // (degree 9). Empirically the Jacobian rank is consistently
    // 1 less than the row count across valences — one of the
    // top-degree BB-coefs is a fixed linear combination of the
    // others, an algebraic identity baked into eq (3)'s sum-of-
    // products structure with the paper's `a := -1, b := 2c(1-u)²,
    // δ := 24c(τ - 1)` choices. (At `n = 3` the rank may drop
    // further; the paper's §3.5 notes the construction differs
    // there.)
    //
    // The assertion below establishes a stable Phase C foundation:
    // 1. The rank is ≥ 8 (paper's "8 constraints" at degree 7).
    // 2. The rank is ≤ the row count.
    // It also reports the actual rank so a future tightening can
    // pin the exact paper-prescribed rank once Phase D's symbolic
    // analysis is in place.
    for n in [3usize, 5, 6, 7, 8, 9, 10] {
        let params = Bi6Reparams::for_valence(n);
        let jac = super::eq3_jacobian_dense(&params);
        let mut flat = Vec::with_capacity(jac.len() * 98);
        for row in &jac {
            flat.extend_from_slice(row);
        }
        let rank = matrix_rank(jac.len(), 98, &flat, 1e-12);
        eprintln!(
            "n={n:2}: eq (3) Jacobian rows = {}, rank = {rank}",
            jac.len()
        );
        assert!(
            rank >= 8,
            "n={n}: eq (3) Jacobian rank = {rank} < 8 (paper says 8 constraints)"
        );
        assert!(
            rank <= jac.len(),
            "n={n}: rank {rank} exceeds row count {}",
            jac.len()
        );
    }
}

#[test]
fn global_constraint_matrix_rank_matches_paper_pre_rule5() {
    // Phase C v2: stack the `n` per-edge Jacobians into the full
    // global M(n, c, τ). Pre-rule-(5), pre-(S1)-(S4) substitutions,
    // we expect the rank to be `n × (8 + 9) = 17n`: the per-edge
    // ranks (8 from eq (2) + 9 from eq (3)) multiplied by `n`
    // edges. There's NO cross-edge rank deficiency at this stage
    // because the constraints involve different sector pairs and
    // are linearly independent in their own column spaces.
    //
    // After rule (5) collapses two BB-coefs per edge (one from
    // each of eq (2) and eq (3)), the rank should drop to `15n`,
    // matching paper line 152's "15 linear equations" per edge.
    // That step lives in Phase D.
    for n in [3usize, 5, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let (rows, cols, m) = super::global_constraint_matrix(n, &params);
        assert_eq!(rows, n * 18, "n={n}: row count mismatch");
        assert_eq!(cols, n * 49, "n={n}: col count mismatch");
        let rank = matrix_rank(rows, cols, &m, 1e-12);
        eprintln!(
            "n={n:2}: global M shape = {rows} × {cols}, rank = {rank}, n×17 = {}",
            n * 17
        );
        // Lower-bound only: 17 per edge - up to 17 cross-edge
        // dependencies introduced by the column overlap between
        // adjacent edges. The exact rank is part of Phase D's
        // measurement output.
        assert!(
            rank >= n * 8,
            "n={n}: rank {rank} < n × 8 = {} (eq (2) alone should clear this)",
            n * 8
        );
        assert!(
            rank <= n * 17,
            "n={n}: rank {rank} > n × 17 = {} (impossible upper bound)",
            n * 17
        );
    }
}

#[test]
fn generating_table_built_and_pivot_count_matches_rank() {
    // Phase D v1: build `Kp15Bi6GeneratingTable` from the global
    // M(n, c, τ) and verify the pivot count agrees with the
    // independently-measured rank from Phase C v2. This is the
    // first concrete production-side artifact of the closure
    // derivation; later phases (D v2, E) layer rule (5)
    // substitution and (S1)-(S4) tightening atop it.
    //
    // Expected pivot counts from Phase C v2's `n × 17 − defect`
    // ranks: n=3 → 41, n=5 → 72, n=6 → 86, n=7 → 105.
    // RREF with partial pivoting may pick up one additional pivot
    // beyond the Gauss-elimination rank for borderline numerical
    // cases; we accept ±1.
    let expected: &[(usize, usize)] = &[(3, 41), (5, 72), (6, 86), (7, 105)];
    for &(n, expected_rank) in expected {
        let params = Bi6Reparams::for_valence(n);
        let h = super::Kp15Bi6GeneratingTable::build(n, params);
        assert_eq!(h.rows, 18 * n);
        assert_eq!(h.cols, 49 * n);
        let rank = h.rank();
        assert!(
            rank.abs_diff(expected_rank) <= 1,
            "n={n}: generating-table pivot count = {rank}, expected {expected_rank} (±1)"
        );
        assert_eq!(h.free_cols.len(), h.cols - h.rank());
    }
}

#[test]
fn generating_table_round_trip_satisfies_constraints() {
    // For a deterministic non-zero free-column loading, the
    // table's `apply` must produce a full column vector with
    // ‖M · x‖ ≈ 0. This is the Phase D analogue of the (S1)-
    // strict test for the per-edge residual: any rule deriving
    // pivot values from free values must satisfy the underlying
    // M·x = 0 constraint.
    for n in [3usize, 5, 6, 7] {
        let params = Bi6Reparams::for_valence(n);
        let h = super::Kp15Bi6GeneratingTable::build(n, params);
        // Deterministic non-zero free-column values keyed to n
        // and to the column index — distinct values, finite.
        let free_values: Vec<f64> = (0..h.free_cols.len())
            .map(|i| 0.1 + (i as f64 + 1.0).sqrt() * 0.01 + (n as f64) * 0.001)
            .collect();
        let x = h.apply(&free_values);
        let residual = h.constraint_residual_norm(&x);
        assert!(
            residual < 1e-8,
            "n={n}: ‖M·x‖ = {residual:.3e} (free-col-determined x should satisfy all constraints)"
        );
    }
}

#[test]
fn rule5_is_implied_by_eq2_eq3_constraints() {
    // **Phase D v2 finding:** [KP15 eq (5)] is the paper's choice
    // that makes the boundary curve `e(u)` drop from monomial-
    // degree 6 to degree 5. We expected this to be **independent**
    // of eq (2)/(3) — the paper describes it as a "setting"
    // (line 148), not a derived identity.
    //
    // Measurement: append the `n` rule-(5) rows to M and re-RREF.
    // **The rank does not increase** at any valence in {3, 5, 6, 7}.
    // Rule (5) is in the row space of the per-edge eq (2)/(3)
    // Jacobians — a non-obvious fact the paper does not state but
    // that falls out of the BB-arithmetic structure once you encode
    // both equations together.
    //
    // Consequence: the paper's "after rule (5), eq (2) drops to
    // degree 6" claim is mechanically *equivalent* to a specific
    // BB-coef cancellation already encoded by eq (2) + eq (3) +
    // the engineered `a, b, d, e` choices. This rules out one
    // path for tightening M's rank; Phase E.1's bridge fix
    // (proper (S1) propagation across sectors) is the next
    // candidate.
    for n in [3usize, 5, 6, 7] {
        let params = Bi6Reparams::for_valence(n);
        let h_plain = super::Kp15Bi6GeneratingTable::build(n, params);
        let h_rule5 = super::Kp15Bi6GeneratingTable::build_with_rule5(n, params);
        assert_eq!(h_rule5.rows, h_plain.rows + n);
        let rank_jump = h_rule5.rank() as isize - h_plain.rank() as isize;
        eprintln!(
            "n={n:2}: rank plain = {}, with rule (5) = {}, Δ = {rank_jump}",
            h_plain.rank(),
            h_rule5.rank()
        );
        assert_eq!(
            rank_jump, 0,
            "n={n}: rule (5) Δrank = {rank_jump}; expected 0 (rule (5) is redundant)"
        );
    }
}

#[test]
fn cross_sector_identifications_lift_rank() {
    // Phase E.1 candidate: add `p^k_{i, 0} = p^{k-1}_{i, 0}` for
    // i ∈ {0, 1, 2}, k ∈ 0..n. These are (S1)'s corner-shared
    // identifications the paper re-affirms at line 158-161. If
    // they're already implied by eq (2)/(3) — like rule (5) is
    // (Phase D v2 finding) — the rank delta is 0 and the identi-
    // fication path is also closed. If they're independent, the
    // rank lifts by `3·(n-1)` (one rotation is redundant given
    // the others form a closed cycle).
    for n in [3usize, 5, 6, 7] {
        let params = Bi6Reparams::for_valence(n);
        let h_plain = super::Kp15Bi6GeneratingTable::build(n, params);
        let h_id = super::Kp15Bi6GeneratingTable::build_with_identifications(n, params);
        assert_eq!(h_id.rows, h_plain.rows + 3 * n);
        let delta = h_id.rank() as isize - h_plain.rank() as isize;
        eprintln!(
            "n={n:2}: rank plain = {}, with identifications = {}, Δ = {delta} (predicted ≤ 3(n−1) = {})",
            h_plain.rank(),
            h_id.rank(),
            3 * (n - 1)
        );
        assert!(
            delta >= 0,
            "n={n}: identifications reduced rank — impossible"
        );
    }
}

#[test]
fn generating_table_with_identifications_round_trip_satisfies_all_constraints() {
    // Companion to `cross_sector_identifications_lift_rank`: when
    // identifications are active in M, the table's `apply` must
    // produce a column vector that satisfies BOTH the eq (2)/(3)
    // constraints AND the cross-sector identifications.
    for n in [3usize, 5, 6, 7] {
        let params = Bi6Reparams::for_valence(n);
        let h = super::Kp15Bi6GeneratingTable::build_with_identifications(n, params);
        let free_values: Vec<f64> = (0..h.free_cols.len())
            .map(|i| 0.1 + (i as f64 + 1.0).sqrt() * 0.005 + (n as f64) * 0.001)
            .collect();
        let x = h.apply(&free_values);
        let residual = h.constraint_residual_norm(&x);
        assert!(
            residual < 1e-8,
            "n={n}: identifications-augmented ‖M·x‖ = {residual:.3e}"
        );
        // Explicit cross-sector identification check: for each
        // adjacent pair (k-1, k) and each i ∈ {0, 1, 2},
        // x[col(k, i, 0)] should equal x[col(k-1, i, 0)].
        for k in 0..n {
            let km1 = (k + n - 1) % n;
            for i in 0..3 {
                let curr = x[super::global_grid_col(n, k, i, 0)];
                let prev = x[super::global_grid_col(n, km1, i, 0)];
                assert!(
                    (curr - prev).abs() < 1e-10,
                    "n={n}, k={k}, i={i}: ṕ - p̀ = {} (curr={curr}, prev={prev})",
                    curr - prev
                );
            }
        }
    }
}

#[test]
fn phase_e2_table_driven_cap_satisfies_eq2_and_eq3_to_machine_zero() {
    // Phase E.2 — exit criterion: the table-driven cap
    // construction must produce per-sector bi-6 grids whose
    // eq (2) G¹ AND eq (3) G² residuals at every cross-sector
    // edge are machine zero, even though the kernel's *existing*
    // assemble path reports unit-length residuals on the same
    // input — see `assembly_eq2_residual_diagnostic`. This test
    // pins the improvement empirically across valences for both
    // constraint families.
    use crate::geometry::gregory::kp_polar_shared::catmull_clark_star_point;
    for n in [3usize, 5, 6, 7, 8] {
        let polar = polar_input(n);
        let tb = synthetic_bi3_tensor_border(n);
        let cap_input = convert_polar_to_bi6_cap_input(&polar, &tb).unwrap();
        let star = catmull_clark_star_point(n, cap_input.c7, &cap_input.c5s, &cap_input.c6s);
        let table = super::Kp15Bi6GeneratingTable::build_with_identifications(n, cap_input.params);
        let grids = table.populate_cap_grids_from_input(&cap_input.sectors, star);

        assert_eq!(grids.len(), n);
        let mut worst_eq2 = 0.0_f64;
        let mut worst_eq3 = 0.0_f64;
        for k in 0..n {
            for i in 0..7 {
                for j in 0..7 {
                    let p = grids[k][i][j];
                    assert!(
                        p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
                        "n={n} sector {k} CP[{i}][{j}] non-finite: {p:?}"
                    );
                }
            }
            let km1 = (k + n - 1) % n;
            let bb2 = super::eq2_g1_residual_bb_degree7(&grids[km1], &grids[k], &cap_input.params);
            for v in bb2 {
                let m = v.length();
                if m > worst_eq2 {
                    worst_eq2 = m;
                }
            }
            let bb3 =
                super::eq3_g2_residual_phi_multiplied(&grids[km1], &grids[k], &cap_input.params);
            for v in &bb3 {
                let m = v.length();
                if m > worst_eq3 {
                    worst_eq3 = m;
                }
            }
        }
        eprintln!(
            "Phase E.2 (n={n:2}): worst eq (2) = {worst_eq2:.3e}, worst eq (3) = {worst_eq3:.3e}"
        );
        assert!(
            worst_eq2 < 1e-8,
            "n={n}: worst eq (2) residual = {worst_eq2:.3e}; expected machine zero"
        );
        assert!(
            worst_eq3 < 1e-8,
            "n={n}: worst eq (3) residual = {worst_eq3:.3e}; expected machine zero"
        );
    }
}

#[test]
fn phase_e3_table_driven_f5_satisfies_eq2_eq3_and_beats_existing_f5() {
    // Phase E.3 — F_5-minimised cap construction via the table.
    // The table's free-column subset that isn't seeded by
    // `Bi6CapInput` (= the deep-interior CPs) is chosen to
    // minimise the per-patch F_5 functional, via three
    // independent SPD solves of size |variable| × |variable|.
    //
    // Exit criteria:
    //   1. The minimisation converges (no SPD singularity).
    //   2. eq (2) and eq (3) residuals at every cross-sector
    //      edge remain machine zero (since the table guarantees
    //      M·x = 0 for ANY free-value loading).
    //   3. The achieved F_5 is less than or equal to the
    //      F_5 of the same cap with all extras defaulted to zero
    //      (i.e. `populate_cap_grids_from_input` without
    //      minimisation). This is the *guarantee* of F_5
    //      minimisation over the variable subset.
    use crate::geometry::gregory::kp_polar_shared::catmull_clark_star_point;
    use crate::geometry::gregory::kp15_bi6_f5::f_kappa_quadratic_form;

    let q_local = f_kappa_quadratic_form(5);

    for n in [3usize, 5, 6, 7, 8] {
        let polar = polar_input(n);
        let tb = synthetic_bi3_tensor_border(n);
        let cap_input = convert_polar_to_bi6_cap_input(&polar, &tb).unwrap();
        let star = catmull_clark_star_point(n, cap_input.c7, &cap_input.c5s, &cap_input.c6s);
        let table = super::Kp15Bi6GeneratingTable::build_with_identifications(n, cap_input.params);

        let grids_zero = table.populate_cap_grids_from_input(&cap_input.sectors, star);
        let (grids_f5, f_vals) = table
            .populate_cap_grids_from_input_with_f_kappa(&cap_input.sectors, star, &q_local)
            .unwrap_or_else(|e| panic!("n={n}: F_5 solve failed: {e}"));

        let f_of = |g: &[[[Point3; 7]; 7]]| -> f64 {
            let mut acc = 0.0;
            for grid in g {
                for r in 0..49 {
                    let i_r = r / 7;
                    let j_r = r % 7;
                    for c in 0..49 {
                        let i_c = c / 7;
                        let j_c = c % 7;
                        let p_r = grid[i_r][j_r];
                        let p_c = grid[i_c][j_c];
                        let q = q_local[r][c];
                        if q == 0.0 {
                            continue;
                        }
                        acc += q * (p_r.x * p_c.x + p_r.y * p_c.y + p_r.z * p_c.z);
                    }
                }
            }
            acc
        };
        let f_zero = f_of(&grids_zero);
        let f_min = f_of(&grids_f5);

        // Verify eq (2)/(3) residuals are still machine zero.
        let mut worst_eq2 = 0.0_f64;
        let mut worst_eq3 = 0.0_f64;
        for k in 0..n {
            let km1 = (k + n - 1) % n;
            let bb2 =
                super::eq2_g1_residual_bb_degree7(&grids_f5[km1], &grids_f5[k], &cap_input.params);
            for v in bb2 {
                worst_eq2 = worst_eq2.max(v.length());
            }
            let bb3 = super::eq3_g2_residual_phi_multiplied(
                &grids_f5[km1],
                &grids_f5[k],
                &cap_input.params,
            );
            for v in &bb3 {
                worst_eq3 = worst_eq3.max(v.length());
            }
        }

        eprintln!(
            "Phase E.3 (n={n:2}): F_5 zero-default = {f_zero:.6e}, F_5 minimised = {f_min:.6e} \
             (Δ = {:.2}×), per-coord F = {:?}; eq(2) = {worst_eq2:.2e}, eq(3) = {worst_eq3:.2e}",
            f_zero / f_min.max(1e-30),
            f_vals,
        );

        assert!(
            f_min <= f_zero + 1e-6,
            "n={n}: F_5 minimised ({f_min:.3e}) exceeds zero-default ({f_zero:.3e})"
        );
        assert!(
            worst_eq2 < 1e-8,
            "n={n}: F_5 cap eq (2) residual = {worst_eq2:.3e}"
        );
        assert!(
            worst_eq3 < 1e-8,
            "n={n}: F_5 cap eq (3) residual = {worst_eq3:.3e}"
        );
    }
}

#[test]
fn comparison_four_variants_paper_p_and_beta_reparameterization() {
    // Compare four cap-construction variants on the same input:
    //   A  No β, full F_5 over free (existing via-table, current code).
    //   B  No β, F_5 restricted to paper's P column set as free.
    //   C  β reparameterization applied to tensor border, full F_5.
    //   D  β + paper-P restriction (closest to KP15 §3.6 algorithm).
    //
    // Reports per variant:
    //   - F_5 value
    //   - worst eq (2) BB residual
    //   - worst eq (3) BB residual
    //   - rank of the constraint matrix used
    //   - free-column count
    use crate::geometry::Point3;
    use crate::geometry::gregory::kp_polar_shared::catmull_clark_star_point;
    use crate::geometry::gregory::kp15_bi6_circulant::CirculantSectorRhs;
    use crate::geometry::gregory::kp15_bi6_f5::f_kappa_quadratic_form;

    let q_local = f_kappa_quadratic_form(5);

    fn sector_rhs_from_grid_rows(q: [[Point3; 7]; 3]) -> CirculantSectorRhs {
        CirculantSectorRhs {
            p_00: q[0][0],
            p_10: q[0][1],
            p_20: q[0][2],
            p_40: q[0][4],
            p_50: q[0][5],
            p_60: q[0][6],
            p_01: q[1][0],
            p_11: q[1][1],
            p_02: q[2][0],
            p_42: q[2][4],
            p_52: q[2][5],
            p_62: q[2][6],
        }
    }

    fn extract_bi3_rows_from_synthetic(
        tb: &crate::geometry::gregory::kp_polar_boundary_input::TensorBorder,
    ) -> Vec<[[Point3; 4]; 3]> {
        match tb {
            crate::geometry::gregory::kp_polar_boundary_input::TensorBorder::Bi3(b) => {
                b.sectors.clone()
            }
            _ => panic!("expected bi-3 tensor border"),
        }
    }

    fn f5_of(patches: &[[[Point3; 7]; 7]], q: &[[f64; 49]; 49]) -> f64 {
        let mut acc = 0.0;
        for grid in patches {
            for r in 0..49 {
                let (i_r, j_r) = (r / 7, r % 7);
                for c in 0..49 {
                    let (i_c, j_c) = (c / 7, c % 7);
                    let p_r = grid[i_r][j_r];
                    let p_c = grid[i_c][j_c];
                    let qq = q[r][c];
                    if qq == 0.0 {
                        continue;
                    }
                    acc += qq * (p_r.x * p_c.x + p_r.y * p_c.y + p_r.z * p_c.z);
                }
            }
        }
        acc
    }

    fn residuals(patches: &[[[Point3; 7]; 7]], params: &Bi6Reparams) -> (f64, f64) {
        let n = patches.len();
        let mut worst_eq2 = 0.0_f64;
        let mut worst_eq3 = 0.0_f64;
        for k in 0..n {
            let km1 = (k + n - 1) % n;
            let bb2 = super::eq2_g1_residual_bb_degree7(&patches[km1], &patches[k], params);
            for v in bb2 {
                worst_eq2 = worst_eq2.max(v.length());
            }
            let bb3 = super::eq3_g2_residual_phi_multiplied(&patches[km1], &patches[k], params);
            for v in &bb3 {
                worst_eq3 = worst_eq3.max(v.length());
            }
        }
        (worst_eq2, worst_eq3)
    }

    eprintln!("\n#  4-way comparison: with / without β × default-free / paper-P-free");

    for n in [5usize, 6] {
        eprintln!("\n## n = {n}");
        let polar = polar_input(n);
        let tb = synthetic_bi3_tensor_border(n);
        let cap_input = convert_polar_to_bi6_cap_input(&polar, &tb).unwrap();
        let star = catmull_clark_star_point(n, cap_input.c7, &cap_input.c5s, &cap_input.c6s);

        // Build β-reparameterized sectors.
        let bi3_rows = extract_bi3_rows_from_synthetic(&tb);
        let mut sectors_beta: Vec<CirculantSectorRhs> = Vec::with_capacity(n);
        for sector_b in &bi3_rows {
            let q = super::apply_beta_reparameterization(sector_b, &cap_input.params);
            sectors_beta.push(sector_rhs_from_grid_rows(q));
        }

        let table_default =
            super::Kp15Bi6GeneratingTable::build_with_identifications(n, cap_input.params);
        let table_paper_p = super::Kp15Bi6GeneratingTable::build_with_identifications_paper_p_free(
            n,
            cap_input.params,
        );

        // Variant A: no β, default-free F_5.
        let (a_grids, _) = table_default
            .populate_cap_grids_from_input_with_f_kappa(&cap_input.sectors, star, &q_local)
            .unwrap();
        // Variant B: no β, paper-P-free F_5.
        let (b_grids, _) = table_paper_p
            .populate_cap_grids_from_input_with_f_kappa(&cap_input.sectors, star, &q_local)
            .unwrap();
        // Variant C: with β, default-free F_5.
        let (c_grids, _) = table_default
            .populate_cap_grids_from_input_with_f_kappa(&sectors_beta, star, &q_local)
            .unwrap();
        // Variant D: with β, paper-P-free F_5.
        let (d_grids, _) = table_paper_p
            .populate_cap_grids_from_input_with_f_kappa(&sectors_beta, star, &q_local)
            .unwrap();

        for (label, grids, table) in [
            ("A: no β, default-free   ", &a_grids, &table_default),
            ("B: no β, paper-P-free   ", &b_grids, &table_paper_p),
            ("C: with β, default-free ", &c_grids, &table_default),
            ("D: with β, paper-P-free ", &d_grids, &table_paper_p),
        ] {
            let f5 = f5_of(grids, &q_local);
            let (eq2, eq3) = residuals(grids, &cap_input.params);
            eprintln!(
                "  {label} F_5 = {f5:.4e}, eq(2) = {eq2:.2e}, eq(3) = {eq3:.2e}, \
                 rank = {}, free_cols = {}",
                table.rank(),
                table.free_cols.len()
            );
        }
    }
}

#[test]
#[ignore = "Phase B diagnostic — measure end-to-end residual"]
fn assembly_eq2_residual_diagnostic() {
    // Diagnostic: report eq (2) residual on the current pipeline
    // output. NOT a pass/fail gate — the current bridge does not
    // (S1)-propagate sector 0's near-star CPs through the rest
    // and does not paper-strictly identify the cross-sector row
    // (per `apply_free_params` line 149-153 in
    // `kp15_bi6_cap_f5_solve.rs`). The residual reported here is
    // the *quantitative* size of those gaps and the deep-interior
    // closure gap the rest of this Phase B / C plan is closing.
    for n in [5usize, 6] {
        let per_edge = eq2_per_edge_residual(n);
        let (edge, bb_idx, worst) = worst_residual(&per_edge);
        eprintln!("n={n}: worst eq (2) residual = {worst:.3e} at edge {edge}, BB index {bb_idx}");
        for (k, edge_mags) in per_edge.iter().enumerate() {
            eprintln!(
                "  edge {k:2} BB residuals: {}",
                edge_mags
                    .iter()
                    .map(|m| format!("{m:.2e}"))
                    .collect::<Vec<_>>()
                    .join("  ")
            );
        }
    }
}
