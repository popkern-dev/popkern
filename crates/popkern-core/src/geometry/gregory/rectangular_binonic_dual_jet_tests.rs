//! Tests for the dual cross-jet bi-nonic patch: bit-identity lifts,
//! per-boundary jet exactness of the squared-weight blend, analytic
//! partials against finite differences, and the exact rational NURBS
//! conversion.

use super::*;
use crate::geometry::gregory::RectangularGregoryBiNonic;

fn sample_plain() -> RectangularGregoryBiNonic {
    // A curved but well-behaved grid: z = sin-ish bump over a sheared
    // lattice, deterministic.
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, p) in row.iter_mut().enumerate() {
            let (x, y) = (i as f64 / 9.0, j as f64 / 9.0);
            *p = Point3::new(
                x + 0.03 * y,
                y - 0.02 * x,
                0.2 * (2.1 * x + 0.3).sin() * (1.7 * y + 0.5).cos(),
            );
        }
    }
    RectangularGregoryBiNonic { grid }
}

fn perturbed_dual_jet(corner: usize) -> RectangularGregoryBiNonicDualJet {
    let mut patch = RectangularGregoryBiNonicDualJet::from_plain(&sample_plain());
    patch.activate_corner(corner);
    let block = patch.corners[corner].as_mut().unwrap();
    // Distinct claims per cell and per side.
    let bump = |p: Point3, d: f64| Point3::new(p.x + d, p.y - 0.5 * d, p.z + 2.0 * d);
    block.cell_11 = (bump(block.cell_11.0, 0.008), bump(block.cell_11.1, -0.010));
    block.cell_12 = (bump(block.cell_12.0, 0.011), bump(block.cell_12.1, -0.007));
    block.cell_21 = (bump(block.cell_21.0, -0.009), bump(block.cell_21.1, 0.013));
    block.cell_22 = (bump(block.cell_22.0, 0.005), bump(block.cell_22.1, -0.012));
    patch
}

/// Plain tensor patch equal to `patch` with a corner's blended cells
/// replaced by ONE boundary's claims.
fn tensor_with_claims(
    patch: &RectangularGregoryBiNonicDualJet,
    corner: usize,
    u_boundary: bool,
) -> RectangularGregoryBiNonic {
    let mut grid = patch.grid;
    let block = patch.corners[corner].unwrap();
    for (cell, pair) in [
        (DualJetCell::Cell11, block.cell_11),
        (DualJetCell::Cell12, block.cell_12),
        (DualJetCell::Cell21, block.cell_21),
        (DualJetCell::Cell22, block.cell_22),
    ] {
        let (i, j) = dual_jet_cell_position(corner, cell);
        grid[i][j] = if u_boundary { pair.0 } else { pair.1 };
    }
    RectangularGregoryBiNonic { grid }
}

#[test]
fn from_plain_is_bit_identical() {
    let plain = sample_plain();
    let mut lifted = RectangularGregoryBiNonicDualJet::from_plain(&plain);
    for corner in 0..4 {
        lifted.activate_corner(corner);
    }
    for i in 0..=16 {
        for j in 0..=16 {
            let (u, v) = (i as f64 / 16.0, j as f64 / 16.0);
            let a = plain.point_at(u, v);
            let b = lifted.point_at(u, v);
            assert_eq!((a.x, a.y, a.z), (b.x, b.y, b.z), "at ({u}, {v})");
            let pu_a = plain.partial_u(u, v);
            let pu_b = lifted.partial_u(u, v);
            assert!(
                (pu_a - pu_b).length() < 1e-12,
                "partial_u drift at ({u}, {v})"
            );
        }
    }
}

/// Each boundary through an activated corner must realize its own
/// claims exactly through second derivatives: the patch's d0/d1/d2
/// jets ON the boundary equal the plain tensor patch built with that
/// boundary's claims.
#[test]
fn boundary_jets_realize_own_claims() {
    for corner in 0..4 {
        let patch = perturbed_dual_jet(corner);
        let (ic, jc) = DUAL_JET_CORNER_POSITIONS[corner];
        let u_fixed = if ic == 0 { 0.0 } else { 1.0 };
        let v_fixed = if jc == 0 { 0.0 } else { 1.0 };
        // u-constant boundary (u = u_fixed): claims .0.
        let ref_u = tensor_with_claims(&patch, corner, true);
        // v-constant boundary: claims .1.
        let ref_v = tensor_with_claims(&patch, corner, false);
        for s in 0..=32 {
            let t = s as f64 / 32.0;
            // u-constant boundary jets.
            for (label, got, want) in [
                (
                    "p",
                    patch.point_at(u_fixed, t) - ref_u.point_at(u_fixed, t),
                    0.0,
                ),
                (
                    "du",
                    Point3::new(0.0, 0.0, 0.0)
                        + (patch.partial_u(u_fixed, t) - ref_u.partial_u(u_fixed, t))
                        - Point3::new(0.0, 0.0, 0.0),
                    0.0,
                ),
                (
                    "duu",
                    Point3::new(0.0, 0.0, 0.0)
                        + (patch.second_partial_uu(u_fixed, t)
                            - ref_u.second_partial_uu(u_fixed, t))
                        - Point3::new(0.0, 0.0, 0.0),
                    0.0,
                ),
            ] {
                let err = got.length();
                // The FD second partial carries ~1e-8·|claim delta|
                // truncation; the analytic terms are exact.
                let bar = if label == "duu" { 1e-6 } else { 1e-12 };
                assert!(
                    err < bar + want,
                    "corner {corner} u-boundary {label} err {err:.3e} at t={t}"
                );
            }
            // v-constant boundary jets.
            let p_err = (patch.point_at(t, v_fixed) - ref_v.point_at(t, v_fixed)).length();
            let dv_err = (patch.partial_v(t, v_fixed) - ref_v.partial_v(t, v_fixed)).length();
            let dvv_err = (patch.second_partial_vv(t, v_fixed)
                - ref_v.second_partial_vv(t, v_fixed))
            .length();
            assert!(p_err < 1e-12, "corner {corner} v-boundary p {p_err:.3e}");
            assert!(dv_err < 1e-12, "corner {corner} v-boundary dv {dv_err:.3e}");
            assert!(
                dvv_err < 1e-6,
                "corner {corner} v-boundary dvv {dvv_err:.3e} at t={t}"
            );
        }
    }
}

#[test]
fn analytic_partials_match_finite_differences() {
    let patch = perturbed_dual_jet(2);
    let h = 1e-6;
    for &(u, v) in &[
        (0.3, 0.4),
        (0.85, 0.9),
        (0.95, 0.7),
        (0.7, 0.95),
        (0.98, 0.98),
        (0.05, 0.05),
    ] {
        let fd_u = (patch.point_at(u + h, v) - patch.point_at(u - h, v)) * (0.5 / h);
        let fd_v = (patch.point_at(u, v + h) - patch.point_at(u, v - h)) * (0.5 / h);
        let an_u = patch.partial_u(u, v);
        let an_v = patch.partial_v(u, v);
        assert!(
            (fd_u - an_u).length() < 1e-6,
            "partial_u mismatch at ({u}, {v}): {:.3e}",
            (fd_u - an_u).length()
        );
        assert!(
            (fd_v - an_v).length() < 1e-6,
            "partial_v mismatch at ({u}, {v}): {:.3e}",
            (fd_v - an_v).length()
        );
    }
}

#[test]
fn corner_limit_convention_is_midpoint() {
    let patch = perturbed_dual_jet(0);
    // At the exact corner the surface value equals the grid corner CP
    // (all blended cells' Bernstein factors vanish there).
    let p = patch.point_at(0.0, 0.0);
    let expected = patch.grid[0][0];
    assert!((p.x - expected.x).abs() < 1e-12);
    assert!((p.y - expected.y).abs() < 1e-12);
    assert!((p.z - expected.z).abs() < 1e-12);
}

#[test]
fn to_nurbs_exact_agrees_with_evaluation() {
    use crate::geometry::surface::SurfaceEval;
    for corner in 0..4 {
        let patch = perturbed_dual_jet(corner);
        let nurbs = patch.to_nurbs_exact();
        assert_eq!(nurbs.degree_u, 11);
        assert_eq!(nurbs.degree_v, 11);
        for (cp, w) in nurbs.control_points.iter().zip(nurbs.weights.iter()) {
            assert!(
                cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite() && w.is_finite(),
                "corner {corner}: non-finite CP/weight in conversion"
            );
            assert!(*w > 0.0, "corner {corner}: non-positive weight {w}");
        }
        let mut worst = 0.0_f64;
        for i in 0..=24 {
            for j in 0..=24 {
                let (u, v) = (i as f64 / 24.0, j as f64 / 24.0);
                let a = patch.point_at(u, v);
                let b = nurbs.point_at(u, v);
                let d = a.distance_to(&b);
                assert!(d.is_finite(), "corner {corner}: NaN eval at ({u}, {v})");
                worst = worst.max(d);
            }
        }
        assert!(
            worst < 1e-10,
            "corner {corner}: NURBS conversion worst deviation {worst:.3e}"
        );
    }
}

#[test]
fn to_nurbs_plain_delegates_to_binonic() {
    let patch = RectangularGregoryBiNonicDualJet::from_plain(&sample_plain());
    let nurbs = patch.to_nurbs_exact();
    assert_eq!(nurbs.degree_u, 9);
    assert_eq!(nurbs.degree_v, 9);
}
