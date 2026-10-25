//! MMA, the projection and the print filter through the optimisation
//! loop (PH-122, PH-118, PH-119;
//! docs/architecture/analysis-generative-design.md).

use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{box_mesh, model_of};
use crate::model::FaceId;
use crate::structural::assemble::FaceMechanics;
use crate::structural::generative::print::{BuildDirection, PrintFilter};
use crate::structural::generative::projection::ProjectionKind;

const E: f64 = 2000.0;
const NU: f64 = 0.3;

fn options(h: f64) -> AnalysisOptions {
    AnalysisOptions {
        target_cell_size: Some(h),
        cg_rel_tol: 1e-8,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    }
}

/// A 12 × 4 × 2 cantilever clamped at x = 0, pulled down at x = 12.
fn cantilever() -> (AnalysisModel, ElasticProblem) {
    let model = model_of(box_mesh([0.0; 3], [12.0, 4.0, 2.0]), 6);
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

/// MMA, the threshold projection to β = 32, volume 0.4, r = 0.75 mm,
/// the loaded end kept.
fn crisp(build: Option<BuildDirection>, h: f64, shift: [f64; 3]) -> Design {
    let (model, problem) = cantilever();
    let mut g = GenerativeOptions::new(0.4, 0.75);
    g.method = UpdateMethod::Mma;
    g.projection_beta = 32.0;
    g.projection_kind = ProjectionKind::Threshold(0.5);
    g.max_iterations = 400;
    g.build_direction = build;
    let mut o = options(h);
    o.grid_shift = shift.map(|f| f * h);
    optimise(&model, &problem, &g, &[FaceId(1)], &o).unwrap_or_else(|e| panic!("{e}"))
}

const UP_Y: BuildDirection = BuildDirection {
    axis: 1,
    negative: false,
};

/// GD2 under MMA: with p = 1 the bar in tension on symmetry rollers
/// reaches c* = F²L/(E f A), as the optimality criteria do.
#[test]
fn mma_reaches_the_convex_bar_optimum() {
    let (l, side) = (12.0, 3.0);
    let model = model_of(box_mesh([0.0; 3], [l, side, side]), 6);
    let roller = FaceMechanics::Support(PieceSupport::Roller);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = roller;
    faces[2] = roller;
    faces[4] = roller;
    faces[1] = FaceMechanics::Traction([1.0, 0.0, 0.0]);
    let problem = ElasticProblem {
        youngs: E,
        poisson: NU,
        faces,
    };
    let area = side * side;
    let f = 0.5;
    let mut g = GenerativeOptions::new(f, 0.75);
    g.penalty = 1.0;
    g.method = UpdateMethod::Mma;
    let d = optimise(&model, &problem, &g, &[], &options(0.5)).unwrap_or_else(|e| panic!("{e}"));
    assert!(d.converged, "{} updates", d.history.len());
    assert!(
        (d.volume_fraction - f).abs() < 1e-9,
        "{}",
        d.volume_fraction
    );
    let ratio = d.compliance / (area * area * l / (E * f * area));
    assert!((ratio - 1.0).abs() < 0.01, "{ratio}");
}

/// The crisp cantilever, then the same built along +y. Unconstrained,
/// the design is crisp (Mnd 4.4% against 47% unprojected, measured) at
/// its volume, but 14.4% of its solid would print unsupported along +y
/// (the negative control). With the print filter none does, for 9.3%
/// more compliance (Langelaar 2017 reports 0–11% on the MBB beam).
#[test]
fn a_crisp_design_prints_without_support_along_its_build_direction() {
    let free = crisp(None, 0.5, [0.0; 3]);
    assert!((free.volume_fraction - 0.4).abs() < 1e-9);
    assert!(free.grey_fraction < 0.06, "{}", free.grey_fraction);
    assert_eq!(free.unsupported_fraction, None);
    // PH-128: its thinnest members are at least a cell, within 2r.
    assert!(
        free.member_thickness >= free.cell_size && free.member_thickness <= 1.5,
        "{}",
        free.member_thickness
    );
    let volumes: Vec<f64> = free.entries.iter().map(|e| e.volume).collect();
    let along_y = PrintFilter::new(&free.entries, free.cell_size, UP_Y).unsupported(
        &free.densities,
        &free.floors,
        &volumes,
        0.5,
    );
    assert!(along_y > 0.1, "{along_y}");
    let printed = crisp(Some(UP_Y), 0.5, [0.0; 3]);
    assert!((printed.volume_fraction - 0.4).abs() < 1e-9);
    assert!(printed.grey_fraction < 0.06, "{}", printed.grey_fraction);
    let unsupported = printed.unsupported_fraction.unwrap_or(1.0);
    assert!(unsupported < 0.005, "{unsupported}");
    let cost = printed.compliance / free.compliance;
    assert!(cost > 1.0 && cost < 1.15, "{cost}");
}

/// GD4 for crisp designs (slow: run with --ignored). Measured
/// 2026-10-13: compliance −1.0% to +1.9% over the 5 offsets (grey
/// 0.4–4.4%); refinement to 0.25 mm in the plan.
#[test]
#[ignore = "heavy: five full crisp optimisations; run on demand with --ignored (physics suite)"]
fn crisp_designs_hold_across_grid_offsets() {
    let shifts = [
        [0.0; 3],
        [0.21, 0.0, 0.0],
        [0.0, 0.37, 0.11],
        [0.43, 0.19, 0.31],
        [0.61; 3],
    ];
    let c: Vec<f64> = shifts
        .iter()
        .map(|&s| {
            let d = crisp(None, 0.5, s);
            assert!(d.grey_fraction < 0.06, "{s:?}: {}", d.grey_fraction);
            d.compliance
        })
        .collect();
    for x in &c {
        assert!((x / c[0] - 1.0).abs() < 0.05, "{c:?}");
    }
}

/// PH-126: with a loose update tolerance the design moves on solves to
/// 1e-3, and the final design is solved again at the analysis's (its
/// last residual under 1e-8). Measured on the crisp cantilever at
/// 0.5 mm (400 updates at most): at 1e-3 3852 CG iterations against
/// 8009 at 1e-8, compliance 3.7498 against 3.7578 (−0.21%), grey 3.9%
/// against 4.7%; at 1e-4 4576, +0.22%; at 1e-6 6559, +0.11%.
#[test]
fn a_loose_update_tolerance_ends_on_an_exact_solve() {
    let (model, problem) = cantilever();
    let mut g = GenerativeOptions::new(0.4, 0.75);
    g.method = UpdateMethod::Mma;
    g.max_iterations = 15;
    g.update_tolerance = 1e-3;
    let d = optimise(&model, &problem, &g, &[FaceId(1)], &options(0.5))
        .unwrap_or_else(|e| panic!("{e}"));
    let last = d.report.residual_history.last().copied().unwrap_or(1.0);
    assert!(last <= 1e-8, "{last}");
    let loose = d
        .history
        .iter()
        .skip(1)
        .map(|i| i.cg_iterations)
        .max()
        .unwrap_or(0);
    assert!(
        loose < d.report.cg_iterations,
        "{loose} against {}",
        d.report.cg_iterations
    );
    let bad = GenerativeOptions {
        update_tolerance: 1.0,
        ..g
    };
    assert!(optimise(&model, &problem, &bad, &[], &options(0.5)).is_err());
}
