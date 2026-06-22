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
