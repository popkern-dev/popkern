//! PK-275 acceptance: the analytic ring Jacobian (inner-spoke block)
//! and the outer-fillet Jacobian block match central finite
//! differences of the actual residual functions.
//!
//! The chain under test is `g ↦ ring CPs (linear, via M) ↦ residual`,
//! exactly the structure `BinonicEval::jacobian` relies on: any base
//! ring and any dense linear map `M` exercise the analytic math in
//! full, so the fixture uses a smooth synthetic ring (no heavyweight
//! K-P builder) and a deterministic pseudo-random `M`.

use super::kp_binonic_jacobian::{analytic_ring_jacobian, outer_fillet_jacobian};
use super::{outer_fillet_residual_vector, residual_vector};
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};

const N_PATCHES: usize = 2;
const D_GUIDE: usize = 12;

/// Smooth, well-conditioned synthetic ring: each patch is a gently
/// curved height-field sheet, offset per patch so spokes are
/// non-degenerate. Normals are far from parallel/antiparallel so no
/// residual sample sits on the |·| kink.
fn base_ring() -> Vec<RectangularGregoryBiNonic> {
    (0..N_PATCHES)
        .map(|p| {
            let mut patch = RectangularGregoryBiNonic::zeros();
            for i in 0..10 {
                for j in 0..10 {
                    let u = i as f64 / 9.0;
                    let v = j as f64 / 9.0;
                    patch.grid[i][j] = Point3::new(
                        u + 0.3 * p as f64,
                        v - 0.2 * p as f64,
                        0.15 * (u * 2.1 + 0.7 * p as f64).sin() + 0.1 * (v * 1.7).cos(),
                    );
                }
            }
            patch
        })
        .collect()
}

/// Deterministic dense pseudo-random linear map `M`:
/// `D_ring × D_GUIDE`, entries in roughly [-0.5, 0.5].
fn synthetic_m() -> Vec<Vec<f64>> {
    let d_ring = N_PATCHES * 300;
    (0..d_ring)
        .map(|i| {
            (0..D_GUIDE)
                .map(|k| {
                    let h = (i * 31 + k * 17 + 7) % 97;
                    (h as f64 - 48.0) / 97.0
                })
                .collect()
        })
        .collect()
}

/// Apply `base + M · g` and rebuild the ring structs.
fn ring_at(
    base: &[RectangularGregoryBiNonic],
    m: &[Vec<f64>],
    g: &[f64],
) -> Vec<RectangularGregoryBiNonic> {
    let mut ring = base.to_vec();
    for (p, patch) in ring.iter_mut().enumerate() {
        for i in 0..10 {
            for j in 0..10 {
                let mut delta = [0.0_f64; 3];
                for (a, d) in delta.iter_mut().enumerate() {
                    let row = &m[p * 300 + (i * 10 + j) * 3 + a];
                    for (k, gk) in g.iter().enumerate() {
                        *d += row[k] * gk;
                    }
                }
                patch.grid[i][j].x += delta[0];
                patch.grid[i][j].y += delta[1];
                patch.grid[i][j].z += delta[2];
            }
        }
    }
    ring
}

/// The cylinder fillet the outer block projects onto: axis Z through
/// the origin, radius large enough that every u=0 boundary point
/// projects successfully and the footpoint normal varies smoothly.
fn fillet_surfaces() -> Vec<SurfaceGeometry> {
    (0..N_PATCHES)
        .map(|_| SurfaceGeometry::Cylinder {
            origin: Point3::new(0.5, 0.5, -1.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
        })
        .collect()
}

/// Central-difference Jacobian of an arbitrary residual closure.
fn fd_jacobian<F: Fn(&[f64]) -> Vec<f64>>(residual: F, g0: &[f64], h: f64) -> Vec<Vec<f64>> {
    let m_rows = residual(g0).len();
    let mut jac = vec![vec![0.0; g0.len()]; m_rows];
    let mut g = g0.to_vec();
    for k in 0..g0.len() {
        g[k] = g0[k] + h;
        let f_plus = residual(&g);
        g[k] = g0[k] - h;
        let f_minus = residual(&g);
        g[k] = g0[k];
        for (row, (fp, fm)) in jac.iter_mut().zip(f_plus.iter().zip(f_minus.iter())) {
            row[k] = (fp - fm) / (2.0 * h);
        }
    }
    jac
}

/// Compare two Jacobians entry-wise, scaled to the FD matrix's
/// largest magnitude (the contact_jacobian_fd_tests.rs convention).
fn assert_jacobians_agree(analytic: &[Vec<f64>], fd: &[Vec<f64>], tol: f64, label: &str) {
    assert_eq!(analytic.len(), fd.len(), "{label}: row count");
    let scale = fd
        .iter()
        .flatten()
        .fold(0.0_f64, |acc, v| acc.max(v.abs()))
        .max(1e-12);
    let mut worst = 0.0_f64;
    let mut worst_at = (0usize, 0usize);
    for (i, (ra, rf)) in analytic.iter().zip(fd.iter()).enumerate() {
        assert_eq!(ra.len(), rf.len(), "{label}: col count row {i}");
        for (k, (a, f)) in ra.iter().zip(rf.iter()).enumerate() {
            let err = (a - f).abs() / scale;
            if err > worst {
                worst = err;
                worst_at = (i, k);
            }
        }
    }
    assert!(
        worst <= tol,
        "{label}: worst scaled deviation {worst:.3e} at (row {}, col {}) exceeds {tol:.1e} \
         (analytic {:.6e} vs FD {:.6e}, scale {scale:.3e})",
        worst_at.0,
        worst_at.1,
        analytic[worst_at.0][worst_at.1],
        fd[worst_at.0][worst_at.1],
    );
}

#[test]
fn inner_spoke_analytic_jacobian_matches_fd() {
    let base = base_ring();
    let m = synthetic_m();
    let g0 = vec![0.0; D_GUIDE];
    let analytic = analytic_ring_jacobian(&ring_at(&base, &m, &g0), &m, D_GUIDE);
    let fd = fd_jacobian(|g| residual_vector(&ring_at(&base, &m, g)), &g0, 1e-6);
    assert_jacobians_agree(&analytic, &fd, 1e-5, "inner-spoke block");
}

#[test]
fn outer_fillet_analytic_jacobian_matches_fd() {
    let base = base_ring();
    let m = synthetic_m();
    let fillets = fillet_surfaces();
    let g0 = vec![0.0; D_GUIDE];
    let ring0 = ring_at(&base, &m, &g0);
    // Fixture sanity: every sample must project (no zero-padded
    // rows), otherwise the agreement below is vacuous.
    let r0 = outer_fillet_residual_vector(&ring0, &fillets);
    assert_eq!(r0.len(), N_PATCHES * 50);
    assert!(
        r0.iter().filter(|v| **v > 1e-6).count() > N_PATCHES * 40,
        "fixture must produce non-trivial outer residuals (got {} non-trivial of {})",
        r0.iter().filter(|v| **v > 1e-6).count(),
        r0.len()
    );
    let analytic = outer_fillet_jacobian(&ring0, &fillets, &m, D_GUIDE);
    let fd = fd_jacobian(
        |g| outer_fillet_residual_vector(&ring_at(&base, &m, g), &fillets),
        &g0,
        1e-6,
    );
    assert_jacobians_agree(&analytic, &fd, 1e-5, "outer-fillet block");
}

/// Off-origin agreement: repeat both comparisons at a non-zero guide
/// vector so the test isn't blessed by any special structure at g=0.
#[test]
fn jacobians_match_fd_at_nonzero_guide_vector() {
    let base = base_ring();
    let m = synthetic_m();
    let fillets = fillet_surfaces();
    let g0: Vec<f64> = (0..D_GUIDE)
        .map(|k| 0.02 * ((k * 13 % 7) as f64 - 3.0))
        .collect();
    let ring0 = ring_at(&base, &m, &g0);
    let inner_analytic = analytic_ring_jacobian(&ring0, &m, D_GUIDE);
    let inner_fd = fd_jacobian(|g| residual_vector(&ring_at(&base, &m, g)), &g0, 1e-6);
    assert_jacobians_agree(&inner_analytic, &inner_fd, 1e-5, "inner-spoke block (g≠0)");
    let outer_analytic = outer_fillet_jacobian(&ring0, &fillets, &m, D_GUIDE);
    let outer_fd = fd_jacobian(
        |g| outer_fillet_residual_vector(&ring_at(&base, &m, g), &fillets),
        &g0,
        1e-6,
    );
    assert_jacobians_agree(&outer_analytic, &outer_fd, 1e-5, "outer-fillet block (g≠0)");
}

/// Patches without a fillet entry must produce zero Jacobian rows,
/// matching the residual's zero padding.
#[test]
fn missing_fillet_entries_produce_zero_rows() {
    let base = base_ring();
    let m = synthetic_m();
    let fillets: Vec<SurfaceGeometry> = fillet_surfaces().into_iter().take(1).collect();
    let ring0 = ring_at(&base, &m, &[0.0; D_GUIDE]);
    let jac = outer_fillet_jacobian(&ring0, &fillets, &m, D_GUIDE);
    assert_eq!(jac.len(), N_PATCHES * 50);
    for (i, row) in jac.iter().enumerate().skip(50) {
        assert!(
            row.iter().all(|v| *v == 0.0),
            "row {i} (patch without fillet) must be zero-padded"
        );
    }
}
