//! PK-275 acceptance: the seam-polish mixed Jacobian matches finite
//! differences of `SeamPolishEval::residual` on every block —
//! coplanarity (analytic), folded-G2 (block FD), outer-fillet sine
//! (analytic footpoint chain), and tether (exact diagonal).
//!
//! The fixture is a CLOSED 4-patch ring (each patch a linear frame
//! `C + u·U + v·V` under a smooth height field, consecutive frames
//! satisfying the spoke-sharing convention `P_{k+1}(1−t, 0) =
//! P_k(1, t)` exactly), so shared boundary rows coincide and the
//! folded-G2 gate (`cos ≥ 0.9` across the spoke) passes — a fixture
//! with discontinuous normals would zero the G2 block and make its
//! agreement vacuous.

// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]

use super::{SeamPolishEval, pack_seam_rows};
use crate::construct::kp_lm_solver::LmEval;
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};

fn height(x: f64, y: f64) -> f64 {
    0.15 * (x + 0.3).sin() * (y - 0.2).cos()
}

/// Closed 4-ring of linear-frame patches under `height`.
fn ring_fixture() -> Vec<RectangularGregoryBiNonic> {
    let frames: [([f64; 2], [f64; 2], [f64; 2]); 4] = [
        ([0.0, 0.0], [1.0, 0.0], [0.0, 1.0]),
        ([1.0, 1.0], [0.0, -1.0], [-1.0, 0.0]),
        ([0.0, 0.0], [1.0, 0.0], [0.0, -1.0]),
        ([1.0, -1.0], [0.0, 1.0], [-1.0, 0.0]),
    ];
    frames
        .iter()
        .map(|(c, u, v)| {
            let mut patch = RectangularGregoryBiNonic::zeros();
            for i in 0..10 {
                for j in 0..10 {
                    let a = i as f64 / 9.0;
                    let b = j as f64 / 9.0;
                    let x = c[0] + a * u[0] + b * v[0];
                    let y = c[1] + a * u[1] + b * v[1];
                    patch.grid[i][j] = Point3::new(x, y, height(x, y));
                }
            }
            patch
        })
        .collect()
}

fn eval_fixture(with_fillets: bool) -> SeamPolishEval {
    let base = ring_fixture();
    let references: Vec<(Vec3, f64)> = base
        .iter()
        .map(|p| {
            let n = p.partial_u(0.5, 0.5).cross(&p.partial_v(0.5, 0.5));
            let len = n.length();
            (n, len)
        })
        .collect();
    let x0 = pack_seam_rows(&base);
    let senses = vec![true; base.len()];
    SeamPolishEval {
        x0,
        references,
        senses,
        max_spoke_folds: 0,
        outer_fillets: with_fillets.then(|| {
            (0..base.len())
                .map(|_| SurfaceGeometry::Cylinder {
                    origin: Point3::new(0.5, 0.0, -3.0),
                    axis: Vec3::new(0.0, 1.0, 0.0),
                    radius: 4.0,
                })
                .collect()
        }),
        base,
        char_len: 1.0,
        selector_m: None,
    }
}

/// Central-difference Jacobian of the eval's own residual.
fn fd_jacobian(eval: &mut SeamPolishEval, x0: &[f64], h: f64) -> Vec<Vec<f64>> {
    let m_rows = eval.residual(x0).expect("residual evaluates").len();
    let mut jac = vec![vec![0.0; x0.len()]; m_rows];
    let mut x = x0.to_vec();
    for k in 0..x0.len() {
        x[k] = x0[k] + h;
        let f_plus = eval.residual(&x).expect("residual evaluates");
        x[k] = x0[k] - h;
        let f_minus = eval.residual(&x).expect("residual evaluates");
        x[k] = x0[k];
        for (row, (fp, fm)) in jac.iter_mut().zip(f_plus.iter().zip(f_minus.iter())) {
            row[k] = (fp - fm) / (2.0 * h);
        }
    }
    jac
}

/// Scaled worst deviation over a row range.
fn worst_scaled(analytic: &[Vec<f64>], fd: &[Vec<f64>], rows: std::ops::Range<usize>) -> f64 {
    let scale = fd[rows.clone()]
        .iter()
        .flatten()
        .fold(0.0_f64, |acc, v| acc.max(v.abs()))
        .max(1e-12);
    let mut worst = 0.0_f64;
    for i in rows {
        for (a, f) in analytic[i].iter().zip(fd[i].iter()) {
            worst = worst.max((a - f).abs() / scale);
        }
    }
    worst
}

/// Evaluation point: the packed base rows plus a small deterministic
/// perturbation. The unperturbed ring is near-G1, which parks the
/// coplanarity samples on the `|triple|` kink where central FD is
/// undefined; the perturbation moves every sample well off the kink
/// (|f| ≫ h) while keeping the folded-G2 `cos ≥ 0.9` gate open.
fn eval_point(base: &[RectangularGregoryBiNonic]) -> Vec<f64> {
    pack_seam_rows(base)
        .iter()
        .enumerate()
        .map(|(k, v)| v + 0.02 * (k as f64 * 1.7).sin())
        .collect()
}

#[test]
fn polish_mixed_jacobian_matches_fd_per_block() {
    let mut eval = eval_fixture(true);
    let x = eval_point(&eval.base);
    let n = eval.base.len();
    let d = x.len();

    let analytic = eval.jacobian(&x).expect("mixed jacobian evaluates");
    let fd = fd_jacobian(&mut eval, &x, 1e-6);
    assert_eq!(analytic.len(), fd.len(), "row count must match residual");
    // Chord Sn, lean barrier 2 × Sn over the whole-spoke parameters
    // (PK-547), folded G2 25n, outer fillet 25n, tether d.
    let sn = n * crate::construct::kp_spoke_normal_sine::whole_spoke_parameters().len();
    assert_eq!(analytic.len(), 3 * sn + 25 * n + 25 * n + d);

    // Non-vacuity: the folded-G2 block must have live samples at
    // the evaluation point.
    let mut perturbed = ring_fixture();
    super::unpack_seam_rows(&x, &mut perturbed);
    let ring = perturbed;
    let g2_live: usize = (0..n)
        .map(|k| {
            let a = &ring[k];
            let b = &ring[(k + 1) % n];
            (0..25)
                .filter(|s| {
                    let t = 0.05 + 0.90 * (*s as f64) / 24.0;
                    super::g2_residual_at(a, b, t).is_some()
                })
                .count()
        })
        .sum();
    assert!(
        g2_live > n * 20,
        "fixture must keep the folded-G2 gate open (live {g2_live} of {})",
        n * 25
    );

    let w_coplanarity = worst_scaled(&analytic, &fd, 0..sn);
    assert!(
        w_coplanarity <= 1e-5,
        "coplanarity block deviates: {w_coplanarity:.3e}"
    );
    // The lean rows are inactive on this upright dome fixture (the
    // active gradient is FD-checked in `kp_spoke_normal_sine`); both
    // Jacobians must agree that they are.
    let w_lean = worst_scaled(&analytic, &fd, sn..3 * sn);
    assert!(w_lean <= 1e-5, "lean block deviates: {w_lean:.3e}");
    // The G2 block is forward-FD inside the mixed Jacobian vs the
    // central FD here — agreement bounded by the scheme difference.
    let w_g2 = worst_scaled(&analytic, &fd, 3 * sn..3 * sn + 25 * n);
    assert!(w_g2 <= 5e-4, "folded-G2 block deviates: {w_g2:.3e}");
    let w_outer = worst_scaled(&analytic, &fd, 3 * sn + 25 * n..3 * sn + 50 * n);
    assert!(
        w_outer <= 1e-5,
        "outer-fillet block deviates: {w_outer:.3e}"
    );
    let w_tether = worst_scaled(&analytic, &fd, 3 * sn + 50 * n..3 * sn + 50 * n + d);
    assert!(w_tether <= 1e-9, "tether block deviates: {w_tether:.3e}");
}

#[test]
fn polish_mixed_jacobian_matches_fd_without_fillets() {
    let mut eval = eval_fixture(false);
    let x = eval_point(&eval.base);
    let n = eval.base.len();
    let d = x.len();
    let analytic = eval.jacobian(&x).expect("mixed jacobian evaluates");
    let fd = fd_jacobian(&mut eval, &x, 1e-6);
    assert_eq!(analytic.len(), fd.len());
    let sn = n * crate::construct::kp_spoke_normal_sine::whole_spoke_parameters().len();
    assert_eq!(analytic.len(), 3 * sn + 25 * n + d);
    let w_coplanarity = worst_scaled(&analytic, &fd, 0..sn);
    assert!(
        w_coplanarity <= 1e-5,
        "coplanarity block deviates: {w_coplanarity:.3e}"
    );
    let w_g2 = worst_scaled(&analytic, &fd, 3 * sn..3 * sn + 25 * n);
    assert!(w_g2 <= 5e-4, "folded-G2 block deviates: {w_g2:.3e}");
}
