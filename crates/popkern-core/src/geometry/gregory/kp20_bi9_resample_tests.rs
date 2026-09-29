use super::*;
use crate::geometry::gregory::kp20_bi6_guide_sample::Kp20Tri6GuideSector;
use crate::geometry::gregory::kp20_paper_guide::{Kp20PaperGuideBorder, Kp20PaperGuideBorderPiece};
use crate::geometry::gregory::{Kp20PieceId, assemble_kp20_paper_cap};

/// The lifted Catmull–Clark ring (genuine guide misfit, C² across sectors).
fn production_border(n: usize) -> Kp20PaperGuideBorder {
    let pieces =
        crate::geometry::gregory::cc_characteristic_tensor_border::cc_lifted_ring_border_pieces(
            n,
            |x, y| 0.4 * (x * x + y * y),
        )
        .expect("lifted ring");
    Kp20PaperGuideBorder {
        sectors: pieces
            .into_iter()
            .map(|[near_c, near_m]| {
                [
                    Kp20PaperGuideBorderPiece::from_bi3(near_c),
                    Kp20PaperGuideBorderPiece::from_bi3(near_m),
                ]
            })
            .collect(),
    }
}

fn identity_guide(n: usize, sector: usize) -> Kp20Tri6GuideSector {
    let alpha = core::f64::consts::TAU / n as f64;
    let (a0, a1) = (sector as f64 * alpha, (sector + 1) as f64 * alpha);
    let (s0, s1) = ([a0.cos(), a0.sin()], [a1.cos(), a1.sin()]);
    let mut guide = Kp20Tri6GuideSector::zero();
    for i in 0..=6_usize {
        for j in 0..=(6 - i) {
            let (u, v) = (i as f64 / 6.0, j as f64 / 6.0);
            guide.cps[crate::geometry::gregory::tri6_flat_index(6, i, j)] =
                Point3::new(u * s0[0] + v * s1[0], u * s0[1] + v * s1[1], 0.0);
        }
    }
    guide
}

/// Exactness gate: with the identity guide `g ∘ τ` is τ itself, a bi-6
/// polynomial, so the bi-6 sampling carries no artefact and the bi-9
/// resampling must be a no-op to roundoff (the same gate as the bi-6
/// sampler's). The paper τ satisfies eq (8)/(9) by construction, so it
/// is a feasible start.
#[test]
fn identity_guide_resampling_is_a_no_op() {
    for n in [3, 4, 5, 6, 8] {
        let tau = crate::geometry::gregory::kp20_paper_tau::solve_kp20_paper_tau(n)
            .expect("tau")
            .tau;
        let alpha = core::f64::consts::TAU / n as f64;
        let guides: Vec<Kp20Tri6GuideSector> = (0..n).map(|s| identity_guide(n, s)).collect();
        let sectors: Vec<Kp20MacroPatch> = (0..n)
            .map(|s| {
                let (c, sn) = ((s as f64 * alpha).cos(), (s as f64 * alpha).sin());
                let mut rotated = tau;
                for id in Kp20PieceId::all() {
                    for row in rotated.piece_mut(id).cps.iter_mut() {
                        for p in row.iter_mut() {
                            *p = Point3::new(p.x * c - p.y * sn, p.x * sn + p.y * c, 0.0);
                        }
                    }
                }
                crate::geometry::gregory::kp20_paper_cap_sample::kp20_sample_paper_cap_sector(
                    n, s, &guides[s], &rotated,
                )
            })
            .collect();
        let resample = resample_kp20_cap_bi9(n, &guides, &tau, &sectors, &sectors)
            .unwrap_or_else(|e| panic!("n={n}: {e}"));
        assert!(
            resample.deviation_after < 1e-10,
            "n={n}: identity resample moved {:e}",
            resample.deviation_after
        );
    }
}

/// A curved border with genuine guide misfit: every relation the cap
/// keeps (internal C², sector C⁰, eq (8)/(9)) holds to roundoff after
/// resampling.
#[test]
fn production_cap_resampling_keeps_every_relation() {
    for n in [4, 5, 6, 7] {
        let report = assemble_kp20_paper_cap(&production_border(n), Point3::new(0.0, 0.0, 0.0))
            .expect("assemble");
        let resample = report.bi9.as_ref().unwrap_or_else(|e| panic!("n={n}: {e}"));
        let [c2, c0, g1, g2] = resample.relation_residual;
        assert!(
            c2 < 1e-10 && c0 < 1e-10 && g1 < 1e-9 && g2 < 1e-8,
            "n={n}: {:?}",
            resample.relation_residual
        );
        assert_eq!(resample.pieces.len(), 4 * n);
    }
}

/// The step (iv) `SoftG` cap does not satisfy eq (8)/(9) exactly: the
/// resampling declines typed and the bi-6 cap serves.
#[test]
fn soft_enforcement_declines() {
    let report = crate::geometry::gregory::assemble_kp20_paper_cap_conditioned(
        &production_border(5),
        Point3::new(0.0, 0.0, 0.0),
        crate::geometry::gregory::Kp20StripEnforcement::SoftG,
        0.0,
    )
    .expect("assemble");
    assert_eq!(report.bi9, Err(Kp20Bi9Decline::SoftEnforcement));
}

/// The closed-form composition jets agree with central differences.
#[test]
fn exact_jets_match_differences() {
    let n = 5;
    let border = production_border(n);
    let guide = crate::geometry::gregory::kp20_paper_guide::solve_kp20_paper_guide(
        &border,
        Point3::new(0.0, 0.0, 0.0),
    )
    .expect("guide");
    let tau = crate::geometry::gregory::kp20_paper_tau::solve_kp20_paper_tau(n)
        .expect("tau")
        .tau;
    for k in 0..4 {
        let exact = ExactPiece::new(n, 2, &guide.sectors[2], &tau.pieces[k]);
        for &(u, v) in &[(0.3, 0.6), (0.7, 0.2), (0.5, 0.5)] {
            let (h, j) = (1e-5, exact.jet(u, v));
            let p = |a: f64, b: f64| exact.jet(a, b).p;
            let du = (p(u + h, v) - p(u - h, v)) * (0.5 / h);
            let dv = (p(u, v + h) - p(u, v - h)) * (0.5 / h);
            let duu = (exact.jet(u + h, v).u - exact.jet(u - h, v).u) * (0.5 / h);
            let dvv = (exact.jet(u, v + h).v - exact.jet(u, v - h).v) * (0.5 / h);
            let scale = 1.0 + j.u.length() + j.v.length();
            assert!((du - j.u).length() < 1e-6 * scale, "u");
            assert!((dv - j.v).length() < 1e-6 * scale, "v");
            assert!((duu - j.uu).length() < 1e-4 * scale, "uu");
            assert!((dvv - j.vv).length() < 1e-4 * scale, "vv");
        }
    }
}

/// PK-600: the fair interior serves only when it measures fairer, and
/// the relations hold either way. On the lifted ring at n = 5 it serves
/// (measured 7.25 → 6.44).
#[test]
fn fair_interior_serves_only_when_fairer() {
    for n in [4, 5, 7] {
        let report = assemble_kp20_paper_cap(&production_border(n), Point3::new(0.0, 0.0, 0.0))
            .expect("assemble");
        let resample = report.bi9.as_ref().unwrap_or_else(|e| panic!("n={n}: {e}"));
        let fair = resample.fair;
        if let Some(after) = fair.max_curvature_after {
            assert_eq!(
                fair.served,
                after <= fair.max_curvature_before,
                "n={n}: {fair:?}"
            );
        } else {
            assert!(!fair.served, "n={n}: {fair:?}");
        }
        if n == 5 {
            assert!(fair.served, "n=5 fair interior should serve: {fair:?}");
        }
        let [c2, c0, g1, g2] = resample.relation_residual;
        assert!(
            c2 < 1e-10 && c0 < 1e-10 && g1 < 1e-9 && g2 < 1e-8,
            "n={n}: {:?}",
            resample.relation_residual
        );
    }
}
