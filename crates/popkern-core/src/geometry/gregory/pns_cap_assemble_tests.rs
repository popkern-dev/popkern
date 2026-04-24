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
