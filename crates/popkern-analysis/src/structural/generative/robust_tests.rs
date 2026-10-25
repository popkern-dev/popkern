//! The robust formulation (PH-124; Wang, Lazarov & Sigmund 2011 §5).

use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{box_mesh, model_of};
use crate::model::FaceId;
use crate::structural::assemble::FaceMechanics;
use crate::structural::generative::Evaluation;
use crate::structural::generative::projection::{Projection, ProjectionKind};

const E: f64 = 2000.0;
const NU: f64 = 0.3;

fn options(h: f64, tol: f64) -> AnalysisOptions {
    AnalysisOptions {
        target_cell_size: Some(h),
        cg_rel_tol: tol,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    }
}

/// A `size` cantilever clamped at x = 0, pulled down at the far end.
fn cantilever(size: [f64; 3]) -> (AnalysisModel, ElasticProblem) {
    let model = model_of(box_mesh([0.0; 3], size), 6);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    faces[1] = FaceMechanics::Traction([0.0, -1.0, 0.0]);
    (
        model,
        ElasticProblem {
            youngs: E,
            poisson: NU,
            faces,
        },
    )
}

/// GD1 for the robust formulation (η = 0.3, β = 4): the eroded
/// design's compliance and the intermediate design's volume against
/// central differences, and the eroded design the least stiff of the
/// three.
#[test]
fn robust_sensitivities_match_central_differences() {
    let (model, problem) = cantilever([4.0, 2.0, 1.0]);
    let mut design = GenerativeOptions::new(0.5, 0.75);
    design.method = UpdateMethod::Mma;
    design.projection_beta = 4.0;
    design.projection_kind = ProjectionKind::Threshold(0.5);
    design.robust = Some(0.3);
    let mut ev = Evaluator::new(&model, &problem, design, &[], &options(0.5, 1e-12))
        .unwrap_or_else(|e| panic!("{e}"));
    let n = ev.entries().len();
    let x: Vec<f64> = (0..n)
        .map(|i| 0.1 + 0.8 * ((7 * i) % 11) as f64 / 10.0)
        .collect();
    let at = |eta: f64| {
        Some(Projection {
            beta: 4.0,
            kind: ProjectionKind::Threshold(eta),
        })
    };
    // Plain evaluations at the three thresholds.
    let plain = |ev: &mut Evaluator, eta: f64| {
        ev.set_projection(at(eta));
        ev.set_robust(None);
        ev.evaluate(&x).unwrap_or_else(|e| panic!("{e}"))
    };
    let eroded = plain(&mut ev, 0.7);
    let intermediate = plain(&mut ev, 0.5);
    let dilated = plain(&mut ev, 0.3);
    assert!(
        eroded.compliance > intermediate.compliance && intermediate.compliance > dilated.compliance,
        "{} {} {}",
        eroded.compliance,
        intermediate.compliance,
        dilated.compliance
    );
    ev.set_projection(at(0.5));
    ev.set_robust(Some(0.3));
    let base = ev.evaluate(&x).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(base.compliance, eroded.compliance);
    assert_eq!(base.physical, intermediate.physical);
    let volumes: Vec<f64> = ev.entries().iter().map(|d| d.volume).collect();
    let volume =
        |e: &Evaluation| -> f64 { e.physical.iter().zip(&volumes).map(|(p, v)| p * v).sum() };
    let delta = 1e-6;
    let mut checked = 0;
    for j in (0..n).step_by(n / 8) {
        let mut up = x.clone();
        up[j] += delta;
        let mut down = x.clone();
        down[j] -= delta;
        let u = ev.evaluate(&up).unwrap_or_else(|e| panic!("{e}"));
        let d = ev.evaluate(&down).unwrap_or_else(|e| panic!("{e}"));
        let fv = (volume(&u) - volume(&d)) / (2.0 * delta);
        let rv = (fv - base.volume_gradient[j]).abs() / fv.abs().max(1e-12);
        assert!(
            rv < 1e-6,
            "volume {j}: {fv} against {}",
            base.volume_gradient[j]
        );
        let fc = (u.compliance - d.compliance) / (2.0 * delta);
        if fc.abs() > 1e-6 * base.compliance {
            checked += 1;
            let rel = (fc - base.gradient[j]).abs() / fc.abs();
            assert!(rel < 1e-5, "entry {j}: {fc} against {}", base.gradient[j]);
        }
    }
    assert!(checked >= 4, "{checked}");
}

/// A robust threshold outside (0, ½), or without the threshold
/// projection and MMA, declines.
#[test]
fn a_bad_robust_setting_declines() {
    let (model, problem) = cantilever([4.0, 2.0, 1.0]);
    let mut good = GenerativeOptions::new(0.5, 0.75);
    good.method = UpdateMethod::Mma;
    good.projection_beta = 4.0;
    good.projection_kind = ProjectionKind::Threshold(0.5);
    good.robust = Some(0.3);
    let mut bad = Vec::new();
    for eta in [0.0, 0.5, f64::NAN] {
        bad.push(GenerativeOptions {
            robust: Some(eta),
            ..good
        });
    }
    bad.push(GenerativeOptions {
        method: UpdateMethod::OptimalityCriteria,
        ..good
    });
    bad.push(GenerativeOptions {
        projection_kind: ProjectionKind::Guest,
        ..good
    });
    bad.push(GenerativeOptions {
        projection_beta: 0.0,
        ..good
    });
    for g in bad {
        let r = optimise(&model, &problem, &g, &[], &options(0.5, 1e-8));
        assert!(
            matches!(r, Err(AnalysisError::InvalidGenerative { .. })),
            "{g:?}"
        );
    }
}

/// The crisp 24 × 8 × 4 cantilever (volume 0.4, r = 1.5 mm, the loaded
/// end kept) at cell size `h`, robust at `eta` or plain.
fn crisp(h: f64, eta: Option<f64>) -> Design {
    let (model, problem) = cantilever([24.0, 8.0, 4.0]);
    let mut g = GenerativeOptions::new(0.4, 1.5);
    g.method = UpdateMethod::Mma;
    g.projection_beta = 32.0;
    g.projection_kind = ProjectionKind::Threshold(0.5);
    g.max_iterations = 400;
    g.robust = eta;
    optimise(&model, &problem, &g, &[FaceId(1)], &options(h, 1e-8))
        .unwrap_or_else(|e| panic!("{e}"))
}

/// Refinement with the filter resolved (slow, about an hour on the
/// loaded machine: run with --ignored). From 3 to 4.5 cells per filter
/// radius the crisp design moves 2.2% (27.28 → 26.68 N·mm, grey 1.9%,
/// 4.2%), the robust (η = 0.3) 7.8% (33.62 → 31.00, eroded 47.5,
/// 52.3), measured 2026-10-13; from 1.5 to 3 cells (12 × 4 × 2,
/// r = 0.75 mm) the crisp design moved 9.3% and the robust 15%.
#[test]
#[ignore = "heavy: slow robust topology optimisations (minutes to an hour); run on demand with --ignored (physics suite)"]
fn resolved_crisp_designs_hold_under_refinement() {
    let coarse = crisp(0.5, None);
    let fine = crisp(1.0 / 3.0, None);
    let change = fine.compliance / coarse.compliance - 1.0;
    assert!(
        change.abs() < 0.05,
        "{} -> {}",
        coarse.compliance,
        fine.compliance
    );
    assert!(coarse.grey_fraction < 0.06 && fine.grey_fraction < 0.06);
}

/// PH-128 (slow, about 40 minutes on the loaded machine: run with
/// --ignored): the robust design's thinnest members are wider than the
/// plain design's. Measured on the 24 × 8 × 4 cantilever at 0.5 mm
/// (r = 1.5 mm): plain 0.5 mm at the 2–10% quantiles, 1.0 at the
/// median; robust (η = 0.3) 1.0 mm from 2% to the median; both 1.5 at
/// 70–99%.
#[test]
#[ignore = "heavy: slow robust topology optimisations (minutes to an hour); run on demand with --ignored (physics suite)"]
fn robust_members_are_wider() {
    let plain = crisp(0.5, None);
    let robust = crisp(0.5, Some(0.3));
    assert!(
        robust.member_thickness > plain.member_thickness,
        "{} against {}",
        robust.member_thickness,
        plain.member_thickness
    );
}

/// PH-124 (slow, about two hours on the loaded machine: run with
/// --ignored): with the filter at 5 and 7 cells (12 × 4 × 2,
/// r = 1 mm, 0.2 and 1/7 mm cells) the robust design (η = 0.3) holds
/// its compliance and its thinnest members in millimetres, the plain
/// design's thinnest member is a cell at each (the negative control).
/// Measured 2026-10-14: robust 4.4826 → 4.4776 (−0.11%), members 0.40
/// and 0.43 mm (0.2–0.21 of 2r), eroded 11.7 and 12.7; plain 3.5749 →
/// 3.6013 (+0.74%), members 0.200 and 0.143 mm.
#[test]
#[ignore = "heavy: slow robust topology optimisations (minutes to an hour); run on demand with --ignored (physics suite)"]
fn robust_designs_hold_their_members_under_refinement() {
    let design = |h: f64, eta: Option<f64>| {
        let (model, problem) = cantilever([12.0, 4.0, 2.0]);
        let mut g = GenerativeOptions::new(0.4, 1.0);
        g.method = UpdateMethod::Mma;
        g.projection_beta = 32.0;
        g.projection_kind = ProjectionKind::Threshold(0.5);
        g.max_iterations = 400;
        g.update_tolerance = 1e-3;
        g.robust = eta;
        optimise(&model, &problem, &g, &[FaceId(1)], &options(h, 1e-8))
            .unwrap_or_else(|e| panic!("{e}"))
    };
    let (coarse, fine) = (design(0.2, Some(0.3)), design(1.0 / 7.0, Some(0.3)));
    let change = fine.compliance / coarse.compliance - 1.0;
    assert!(change.abs() < 0.05, "{change}");
    let (a, b) = (coarse.member_thickness, fine.member_thickness);
    assert!((b / a - 1.0).abs() < 0.2 && a > 1.5 * 0.2, "{a} {b}");
    let plain = design(1.0 / 7.0, None);
    assert!(
        plain.member_thickness < b,
        "{} against {b}",
        plain.member_thickness
    );
}
