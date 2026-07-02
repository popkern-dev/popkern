use super::*;
use crate::geometry::gregory::kp20_paper_tau::solve_kp20_paper_tau;
use crate::geometry::gregory::tri6_flat_index;

fn identity_guide(n: usize, sector: usize) -> Kp20Tri6GuideSector {
    let alpha = core::f64::consts::TAU / n as f64;
    let a0 = sector as f64 * alpha;
    let a1 = (sector + 1) as f64 * alpha;
    let spoke0 = [a0.cos(), a0.sin()];
    let spoke1 = [a1.cos(), a1.sin()];
    let mut guide = Kp20Tri6GuideSector::zero();
    for i in 0..=6_usize {
        for j in 0..=(6 - i) {
            let u = i as f64 / 6.0;
            let v = j as f64 / 6.0;
            guide.cps[tri6_flat_index(6, i, j)] = Point3::new(
                u * spoke0[0] + v * spoke1[0],
                u * spoke0[1] + v * spoke1[1],
                0.0,
            );
        }
    }
    guide
}

/// A bi-6 polynomial is reconstructed exactly from its four (3,3)
/// corner blocks, so sampling `identity ∘ τ` must reproduce τ's
/// pieces to roundoff — the decisive exactness gate for the whole
/// composition/sampling pipeline.
#[test]
fn identity_guide_sampling_reproduces_tau_exactly() {
    for n in [3, 5, 6] {
        let tau = solve_kp20_paper_tau(n).expect("tau").tau;
        let guide = identity_guide(n, 0);
        let sampled = kp20_sample_paper_cap_sector(n, 0, &guide, &tau);
        let mut worst = 0.0_f64;
        for piece_id in Kp20PieceId::all() {
            let a = sampled.piece(piece_id);
            let b = tau.piece(piece_id);
            for i in 0..7 {
                for j in 0..7 {
                    worst = worst
                        .max((a.cps[i][j].x - b.cps[i][j].x).abs())
                        .max((a.cps[i][j].y - b.cps[i][j].y).abs())
                        .max(a.cps[i][j].z.abs());
                }
            }
        }
        eprintln!("sampling exactness n={n}: worst={worst:e}");
        assert!(
            worst < 1.0e-9,
            "n={n}: sampled cap deviates from tau by {worst:e}",
        );
    }
}

/// Linear lift: with the guide `(x, y) -> (x, y, ax + by)` the sampled
/// cap must be the lifted τ exactly.
#[test]
fn linear_lift_guide_sampling_is_exact() {
    let n = 5;
    let (alpha_lift, beta_lift) = (0.375, -0.2);
    let tau = solve_kp20_paper_tau(n).expect("tau").tau;
    let mut guide = identity_guide(n, 0);
    for cp in guide.cps.iter_mut() {
        *cp = Point3::new(cp.x, cp.y, alpha_lift * cp.x + beta_lift * cp.y);
    }
    let sampled = kp20_sample_paper_cap_sector(n, 0, &guide, &tau);
    let mut worst = 0.0_f64;
    for piece_id in Kp20PieceId::all() {
        let a = sampled.piece(piece_id);
        let b = tau.piece(piece_id);
        for i in 0..7 {
            for j in 0..7 {
                worst = worst
                    .max((a.cps[i][j].x - b.cps[i][j].x).abs())
                    .max((a.cps[i][j].y - b.cps[i][j].y).abs())
                    .max(
                        (a.cps[i][j].z - alpha_lift * b.cps[i][j].x - beta_lift * b.cps[i][j].y)
                            .abs(),
                    );
            }
        }
    }
    assert!(worst < 1.0e-9, "lifted sampling deviates by {worst:e}");
}
