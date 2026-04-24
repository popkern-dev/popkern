//! Tests for [`pns_cap_assemble`] — end-to-end bi-4 cap assembly.

#![allow(clippy::needless_range_loop)]

use super::super::pns_beta::BetaReparam;
use super::super::pns_boundary::{Bi2BorderStrip, Bi2TensorBorderRing};
use super::super::pns_sigma::sigma_center_jet_default;
use super::*;
use crate::geometry::Point3;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn approx_pt(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

/// Build a trivial boundary ring where every sector's strip is
/// a planar horizontal curve (z = 0) — useful as a first-cut
/// smoke test for the scaffolding.
fn planar_ring(n: usize, radius: f64) -> Bi2TensorBorderRing {
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let a0 = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
        let a1 = 2.0 * std::f64::consts::PI * ((k as f64) + 0.5) / (n as f64);
        let a2 = 2.0 * std::f64::consts::PI * ((k + 1) as f64) / (n as f64);
        sectors.push(Bi2BorderStrip {
            pos: [
                pt(radius * a0.cos(), radius * a0.sin(), 0.0),
                pt(radius * a1.cos(), radius * a1.sin(), 0.0),
                pt(radius * a2.cos(), radius * a2.sin(), 0.0),
            ],
            deriv: [pt(0.0, 0.0, -1.0); 3], // inward-pointing (toward pole)
        });
    }
    Bi2TensorBorderRing { sectors }
}

#[test]
fn scaffold_even_valence_is_unsupported() {
    for n in [4_usize, 6, 8] {
        let input = CapAssembleInput {
            boundary: planar_ring(n, 1.0),
            beta: BetaReparam::bi4(n),
            sigma: sigma_center_jet_default(n),
            pole: pt(0.0, 0.0, 0.5),
            kappa: 3,
        };
        let result = assemble_pns_bi4_cap(&input);
        assert!(
            matches!(result, Err(CapAssembleError::EvenValenceUnsupported(_))),
            "n={n}: expected EvenValenceUnsupported, got {:?}",
            result.err()
        );
    }
}

#[test]
fn scaffold_invalid_valence_errors() {
    for n in [1_usize, 2, 13, 20] {
        let input = CapAssembleInput {
            // Even though some valences are "invalid", we still
            // need a well-formed boundary ring to instantiate.
            // Use 5 sectors but tell the function to expect these.
            boundary: planar_ring(n.max(3), 1.0),
            beta: BetaReparam::bi4(n.max(3)),
            sigma: sigma_center_jet_default(n.max(3)),
            pole: pt(0.0, 0.0, 0.5),
            kappa: 3,
        };
        // Actually, we pass planar_ring(n) so n must be ≥ 3.
        // Skip n < 3 here; they're caught elsewhere.
        if n < 3 {
            continue;
        }
        let result = assemble_pns_bi4_cap(&input);
        // n in {13, 20} is > 12 and should error.
        if n > 12 {
            assert!(
                matches!(result, Err(CapAssembleError::InvalidValence(_))),
                "n={n}: expected InvalidValence, got {:?}",
                result.err()
            );
        }
    }
}

#[test]
fn scaffold_pole_pinned_at_all_sector_origins() {
    let n = 5;
    let pole = pt(0.1, 0.2, 0.3);
    let input = CapAssembleInput {
        boundary: planar_ring(n, 1.0),
        beta: BetaReparam::bi4(n),
        sigma: sigma_center_jet_default(n),
        pole,
        kappa: 3,
    };
    let result = assemble_pns_bi4_cap(&input).expect("n=5 should succeed");
    assert_eq!(result.patches.len(), n);
    for k in 0..n {
        assert!(
            approx_pt(result.patches[k][0][0], pole, 1e-12),
            "sector {k} pole at (0,0) is {:?}, expected {pole:?}",
            result.patches[k][0][0]
        );
    }
}

#[test]
fn scaffold_boundary_row_i4_matches_b_o_beta_row0() {
    use super::super::pns_boundary::compose_bi4_outer_layers;
    let n = 5;
    let pole = pt(0.0, 0.0, 1.0);
    let boundary = planar_ring(n, 2.0);
    let beta = BetaReparam::bi4(n);
    let input = CapAssembleInput {
        boundary: boundary.clone(),
        beta,
        sigma: sigma_center_jet_default(n),
        pole,
        kappa: 3,
    };
    let result = assemble_pns_bi4_cap(&input).expect("n=5 should succeed");

    // Independently compute the expected row0 for each sector and
    // compare.
    for k in 0..n {
        let (row0, _row1) = compose_bi4_outer_layers(&boundary.sectors[k], &beta);
        for j in 0..5 {
            assert!(
                approx_pt(result.patches[k][4][j], row0[j], 1e-12),
                "sector {k} patch[4][{j}] = {:?}, expected row0 {:?}",
                result.patches[k][4][j],
                row0[j]
            );
        }
    }
}

#[test]
fn scaffold_returns_n_patches() {
    for n in [3_usize, 5, 7, 9, 11] {
        let input = CapAssembleInput {
            boundary: planar_ring(n, 1.0),
            beta: BetaReparam::bi4(n),
            sigma: sigma_center_jet_default(n),
            pole: pt(0.0, 0.0, 0.0),
            kappa: 3,
        };
        let result = assemble_pns_bi4_cap(&input).expect("odd n should succeed");
        assert_eq!(result.patches.len(), n, "expected {n} patches");
    }
}

#[test]
fn j4_row_equals_next_sector_j0_row() {
    // Cap topology: P[k][i, 4] = P[(k+1) mod N][i, 0] for i=0..4.
    let n = 5;
    let input = CapAssembleInput {
        boundary: planar_ring(n, 1.0),
        beta: BetaReparam::bi4(n),
        sigma: sigma_center_jet_default(n),
        pole: pt(0.0, 0.0, 0.0),
        kappa: 3,
    };
    let result = assemble_pns_bi4_cap(&input).expect("n=5 should succeed");
    for k in 0..n {
        let nxt = (k + 1) % n;
        for i in 0..5 {
            assert!(
                approx_pt(result.patches[k][i][4], result.patches[nxt][i][0], 1e-10),
                "P[{k}][{i}, 4] = {:?}, P[{nxt}][{i}, 0] = {:?}",
                result.patches[k][i][4],
                result.patches[nxt][i][0]
            );
        }
    }
}

#[test]
fn a_0_pinned_to_pole() {
    let pole = pt(0.3, -0.7, 0.2);
    let n = 5;
    let input = CapAssembleInput {
        boundary: planar_ring(n, 1.0),
        beta: BetaReparam::bi4(n),
        sigma: sigma_center_jet_default(n),
        pole,
        kappa: 3,
    };
    let result = assemble_pns_bi4_cap(&input).expect("n=5 should succeed");
    assert!(
        approx_pt(result.a.a_0, pole, 1e-12),
        "a.a_0 = {:?}, expected {pole:?}",
        result.a.a_0
    );
}

/// Diagnostic probe: report the L¹ residual between paper's
/// `cur_{10}` from eq (2) and our construction's `p^k` canonical
/// value at the matching location, given the axis swap.
///
/// Per the KK-Peters reference implementation author's
/// clarification (2026-04-24):
/// - σ^{k+1}_{ij} = R_α · σ^k_{ij} between adjacent sectors
///   (coefficient-by-coefficient image-space rotation)
/// - Paper's ``p_{ij}` (cur side) = canonical `p^k_{j, i}` (swapped)
/// - Paper's `´p_{ij}` (nxt side) = canonical `p^{k+1}_{i, j}` (no swap)
///
/// Under these two corrections and the paper's σ_{10} = (d_1, 0)
/// and σ_{01} = (d_1·cos α, d_1·sin α), eq (2)'s cur_{10} matches
/// `p^k_{0, 1}` canonical (paper (1, 0) = canonical (0, 1)).
#[test]
#[ignore = "diagnostic probe — verifies eq (2) cur_{10} ≈ 0 residual"]
fn probe_eq2_cur10_residual() {
    use super::super::pns_g1_constraints::bi4_g1_row_from_row1;
    let n = 5;
    let input = CapAssembleInput {
        boundary: planar_ring(n, 1.5),
        beta: BetaReparam::bi4(n),
        sigma: sigma_center_jet_default(n),
        pole: pt(0.0, 0.0, 0.3),
        kappa: 3,
    };
    let result = assemble_pns_bi4_cap(&input).expect("n=5 should succeed");
    for k in 0..n {
        let nxt = (k + 1) % n;
        // Paper's cur's j=1 row = canonical `p^k_{j, i}` for paper
        // i=0..4 → canonical (1, 0..4), i.e., fixed u_c=1 with
        // v_c varying.
        let cur_row1 = [
            result.patches[k][1][0], // paper cur_{0, 1} = canonical (1, 0)
            result.patches[k][1][1], // paper cur_{1, 1} = canonical (1, 1)
            result.patches[k][1][2], // paper cur_{2, 1} = canonical (1, 2)
            result.patches[k][1][3], // paper cur_{3, 1} = canonical (1, 3)
            result.patches[k][1][4], // paper cur_{4, 1} = canonical (1, 4)
        ];
        // Paper's nxt's j=1 row = canonical (i, 1) for paper i=0..4
        // (no swap).
        let nxt_row1 = [
            result.patches[nxt][0][1],
            result.patches[nxt][1][1],
            result.patches[nxt][2][1],
            result.patches[nxt][3][1],
            result.patches[nxt][4][1],
        ];
        let g1 = bi4_g1_row_from_row1(result.patches[k][0][0], cur_row1, nxt_row1, n);
        // Paper's cur_{1, 0} = canonical p^k_{0, 1}.
        let expected = result.patches[k][0][1];
        let diff = (g1.cur10.x - expected.x).abs()
            + (g1.cur10.y - expected.y).abs()
            + (g1.cur10.z - expected.z).abs();
        println!("sector {k}: eq2 cur_{{10}} - P[{k}][0, 1] L¹ = {diff:.3e}");
    }
}

#[test]
fn f_kappa_reduced_relative_to_zero_fill() {
    // Verify the F_κ solve actually REDUCES the functional relative
    // to a trivial baseline (free CPs set to zero, a=0). If the
    // joint-solve routine is broken, F_κ would be >= baseline.
    use super::super::pns_f_kappa::f_kappa_bi4;

    let n = 5;
    let pole = pt(0.0, 0.0, 0.3);
    let kappa = 3;
    let input = CapAssembleInput {
        boundary: planar_ring(n, 1.5),
        beta: BetaReparam::bi4(n),
        sigma: sigma_center_jet_default(n),
        pole,
        kappa,
    };
    let result = assemble_pns_bi4_cap(&input).expect("n=5 should succeed");

    // F_κ sum over all sectors' patches.
    let f_solved: f64 = (0..n).map(|k| f_kappa_bi4(&result.patches[k], kappa)).sum();

    // Baseline: place free CPs at zero, a=0 — reconstruct by
    // calling assemble with a proxy via scaffold-style. Instead,
    // just use a patch with pole + boundary + zeros elsewhere.
    // Compute baseline F_κ by reconstructing the cap WITHOUT the
    // joint solve (free CPs = 0, a = 0).
    use super::super::pns_boundary::compose_bi4_outer_layers;
    let zero = pt(0.0, 0.0, 0.0);
    let mut baseline_patches = vec![[[zero; 5]; 5]; n];
    for k in 0..n {
        let (row0, row1) = compose_bi4_outer_layers(&input.boundary.sectors[k], &input.beta);
        for j in 0..5 {
            baseline_patches[k][4][j] = row0[j];
            baseline_patches[k][3][j] = pt(
                row0[j].x + 0.25 * row1[j].x,
                row0[j].y + 0.25 * row1[j].y,
                row0[j].z + 0.25 * row1[j].z,
            );
        }
        baseline_patches[k][0][0] = pole;
    }
    let f_baseline: f64 = (0..n)
        .map(|k| f_kappa_bi4(&baseline_patches[k], kappa))
        .sum();

    println!("F_κ solved = {f_solved:.6e}, baseline = {f_baseline:.6e}");
    assert!(
        f_solved < f_baseline,
        "F_κ solve did not reduce objective (solved = {f_solved:.3e}, baseline = {f_baseline:.3e})"
    );
}
