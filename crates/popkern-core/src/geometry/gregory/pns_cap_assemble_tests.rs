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
fn even_valence_assembles_successfully() {
    // n=4 is rejected (c=0 σ singularity); n ∈ {6, 8, 10, 12} succeed.
    for n in [6_usize, 8, 10, 12] {
        let input = CapAssembleInput {
            boundary: planar_ring(n, 1.0),
            beta: BetaReparam::bi4(n),
            sigma: sigma_center_jet_default(n),
            pole: pt(0.0, 0.0, 0.5),
            kappa: 3,
        };
        let result = assemble_pns_bi4_cap(&input);
        assert!(
            result.is_ok(),
            "n={n} should succeed, got {:?}",
            result.err()
        );
        let cap = result.unwrap();
        assert_eq!(cap.patches.len(), n);
    }
}

#[test]
fn pole_from_boundary_gives_centroid_for_symmetric_ring() {
    // For the planar_ring test fixture with deriv = (0, 0, -1),
    // b^k_{i, 1} = pos[i] + (0, 0, -0.5). Symmetric around the
    // z-axis, so the x, y components average to 0 and the z
    // component averages to -0.5 for both d_1 and (d_2+d_4)/2.
    // Therefore pole ≈ (0, 0, -0.5) for any e.
    let n = 5;
    let ring = planar_ring(n, 1.0);
    let pole = pole_from_boundary(&ring);
    assert!(pole.x.abs() < 1e-12, "pole.x = {}", pole.x);
    assert!(pole.y.abs() < 1e-12, "pole.y = {}", pole.y);
    assert!(
        (pole.z - (-0.5)).abs() < 1e-12,
        "pole.z = {} (expected -0.5)",
        pole.z
    );
}

#[test]
fn pole_from_boundary_feeds_assembly() {
    // End-to-end: auto-pole + assembly produces a cap consistent
    // with the derived pole.
    let n = 5;
    let ring = planar_ring(n, 1.5);
    let pole = pole_from_boundary(&ring);
    let input = CapAssembleInput {
        boundary: ring,
        beta: BetaReparam::bi4(n),
        sigma: sigma_center_jet_default(n),
        pole,
        kappa: 3,
    };
    let result = assemble_pns_bi4_cap(&input).expect("auto-pole assembly should succeed");
    for k in 0..n {
        assert!(approx_pt(result.patches[k][0][0], pole, 1e-12));
    }
}

#[test]
fn n4_rejected_with_informative_error() {
    let input = CapAssembleInput {
        boundary: planar_ring(4, 1.0),
        beta: BetaReparam::bi4(4),
        sigma: sigma_center_jet_default(4),
        pole: pt(0.0, 0.0, 0.5),
        kappa: 3,
    };
    let result = assemble_pns_bi4_cap(&input);
    assert!(
        matches!(result, Err(CapAssembleError::N4NotSupported)),
        "n=4 should return N4NotSupported, got {:?}",
        result.err()
    );
}

#[test]
fn eq2_g1_holds_for_even_n() {
    // Verify all 5 eq (2) outputs hold to floating-point precision
    // for even n ∈ {6, 8} (n=4 is singular; other even n work).
    use super::super::pns_g1_constraints::bi4_g1_row_from_row1;
    for n in [6_usize, 8] {
        let input = CapAssembleInput {
            boundary: planar_ring(n, 1.3),
            beta: BetaReparam::bi4(n),
            sigma: sigma_center_jet_default(n),
            pole: pt(0.0, 0.0, 0.2),
            kappa: 3,
        };
        let result = assemble_pns_bi4_cap(&input).expect("even n should succeed");
        let tol = 1e-10;
        for k in 0..n {
            let nxt = (k + 1) % n;
            let cur_row1 = [
                result.patches[k][1][0],
                result.patches[k][1][1],
                result.patches[k][1][2],
                result.patches[k][1][3],
                result.patches[k][1][4],
            ];
            let nxt_row1 = [
                result.patches[nxt][0][1],
                result.patches[nxt][1][1],
                result.patches[nxt][2][1],
                result.patches[nxt][3][1],
                result.patches[nxt][4][1],
            ];
            let g1 = bi4_g1_row_from_row1(result.patches[k][0][0], cur_row1, nxt_row1, n);
            assert!(
                approx_pt(g1.cur10, result.patches[k][0][1], tol),
                "n={n} sector {k} cur_{{10}} mismatch"
            );
            assert!(
                approx_pt(g1.cur20, result.patches[k][0][2], tol),
                "n={n} sector {k} cur_{{20}} mismatch"
            );
            assert!(
                approx_pt(g1.cur30, result.patches[k][0][3], tol),
                "n={n} sector {k} cur_{{30}} mismatch"
            );
            assert!(
                approx_pt(g1.cur40, result.patches[k][0][4], tol),
                "n={n} sector {k} cur_{{40}} mismatch"
            );
            assert!(
                approx_pt(g1.nxt21, result.patches[nxt][2][1], tol),
                "n={n} sector {k} nxt_{{21}} mismatch"
            );
        }
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
fn canonical_v4_edge_matches_b_o_beta_row0() {
    // Cap's outer boundary = canonical v_c=4 edge (per paper's
    // axis-swap convention). Check patches[k][i][4] == row0[i].
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

    for k in 0..n {
        let (row0, _row1) = compose_bi4_outer_layers(&boundary.sectors[k], &beta);
        for i in 0..5 {
            assert!(
                approx_pt(result.patches[k][i][4], row0[i], 1e-12),
                "sector {k} patch[{i}][4] = {:?}, expected row0 {:?}",
                result.patches[k][i][4],
                row0[i]
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
fn cur_nxt_seam_sharing_holds() {
    // Cap topology (canonical): cur's (0, v) = nxt's (v, 0) for
    // v=0..4. This encodes the shared seam between sectors k and
    // k+1 along their respective canonical u_c=0 / v_c=0 edges.
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
        for v in 0..5 {
            assert!(
                approx_pt(result.patches[k][0][v], result.patches[nxt][v][0], 1e-10),
                "P[{k}][0, {v}] = {:?}, P[{nxt}][{v}, 0] = {:?}",
                result.patches[k][0][v],
                result.patches[nxt][v][0]
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

/// Verify eq (2) G¹ identities hold across every seam of the
/// assembled cap. Paper's eq (2) pins 5 outputs
/// (cur_{10}, cur_{20}, cur_{30}, cur_{40}, nxt_{21}); under the
/// correct σ convention + axis swap, the assembled cap should
/// satisfy all 5 to floating-point precision.
///
/// Paper-to-canonical mapping (cur side has swap, nxt doesn't):
/// - paper cur_{i, 0} = canonical p^k_{0, i}
/// - paper cur_{i, 1} = canonical p^k_{1, i}
/// - paper nxt_{i, 1} = canonical p^{k+1}_{i, 1}
#[test]
fn eq2_g1_all_outputs_hold_after_solve() {
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
    let tol = 1e-10;
    for k in 0..n {
        let nxt = (k + 1) % n;
        // Paper's cur_{i, 1} = canonical (1, i) — swap on cur.
        let cur_row1 = [
            result.patches[k][1][0],
            result.patches[k][1][1],
            result.patches[k][1][2],
            result.patches[k][1][3],
            result.patches[k][1][4],
        ];
        // Paper's nxt_{i, 1} = canonical (i, 1) — no swap on nxt.
        let nxt_row1 = [
            result.patches[nxt][0][1],
            result.patches[nxt][1][1],
            result.patches[nxt][2][1],
            result.patches[nxt][3][1],
            result.patches[nxt][4][1],
        ];
        let g1 = bi4_g1_row_from_row1(result.patches[k][0][0], cur_row1, nxt_row1, n);
        // Check all 5 eq (2) outputs against the assembled cap's
        // canonical values (with paper→canonical swap on cur side).
        let expected_cur10 = result.patches[k][0][1]; // paper (1, 0) = canonical (0, 1)
        let expected_cur20 = result.patches[k][0][2]; // paper (2, 0) = canonical (0, 2)
        let expected_cur30 = result.patches[k][0][3]; // paper (3, 0) = canonical (0, 3)
        let expected_cur40 = result.patches[k][0][4]; // paper (4, 0) = canonical (0, 4)
        let expected_nxt21 = result.patches[nxt][2][1]; // paper (2, 1) no swap
        assert!(
            approx_pt(g1.cur10, expected_cur10, tol),
            "sector {k} cur_{{10}}: eq2={:?} vs canonical={:?}",
            g1.cur10,
            expected_cur10
        );
        assert!(
            approx_pt(g1.cur20, expected_cur20, tol),
            "sector {k} cur_{{20}}: eq2={:?} vs canonical={:?}",
            g1.cur20,
            expected_cur20
        );
        assert!(
            approx_pt(g1.cur30, expected_cur30, tol),
            "sector {k} cur_{{30}}: eq2={:?} vs canonical={:?}",
            g1.cur30,
            expected_cur30
        );
        assert!(
            approx_pt(g1.cur40, expected_cur40, tol),
            "sector {k} cur_{{40}}: eq2={:?} vs canonical={:?}",
            g1.cur40,
            expected_cur40
        );
        assert!(
            approx_pt(g1.nxt21, expected_nxt21, tol),
            "sector {k} nxt_{{21}}: eq2={:?} vs canonical={:?}",
            g1.nxt21,
            expected_nxt21
        );
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
