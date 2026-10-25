//! Singularities by refinement growth (PH-038): a sharp internal corner
//! whose BRep edge is not reported grows from h to h/2; a resolved
//! fillet does not.

use popkern_analysis::assess::surface_peaks;
use popkern_analysis::assess::{AssessOptions, Singularity, label_by_growth, peak_growth};
use popkern_analysis::model::{AnalysisModel, FaceId, Provenance, Surface};
use popkern_analysis::structural::assemble::FaceMechanics;
use popkern_analysis::structural::{ElasticProblem, StructuralResult, solve_elastic};
use popkern_analysis::{AnalysisOptions, element::vector_boundary::PieceSupport};
use popkern_analysis_brep::{FeatureTolerances, model_from_breps};
use popkern_analysis_cli::step_model::tessellation_options;
use popkern_analysis_fixtures::{SteppedBar, stepped_bar_with};

fn model_of(part: &popkern_part::Part) -> AnalysisModel {
    model_from_breps(
        &[part.brep().clone()],
        &tessellation_options(),
        &FeatureTolerances::default(),
        Provenance::Popkern {
            document: "singularity_fixtures".to_string(),
            revision: 0,
        },
    )
    .unwrap_or_else(|e| panic!("{e}"))
}

/// Plane faces of `m` with outward normal `n` at offset `d` (n·x = d).
fn plane_faces(m: &AnalysisModel, n: [f64; 3], d: f64) -> Vec<usize> {
    m.faces
        .iter()
        .enumerate()
        .filter(|(_, f)| match f.surface {
            Surface::Plane { normal, offset } => {
                (0..3).all(|k| (normal[k] - n[k]).abs() < 1e-9) && (offset - d).abs() < 1e-9
            }
            _ => false,
        })
        .map(|(i, _)| i)
        .collect()
}

fn clamped_faces(m: &AnalysisModel, n: [f64; 3], d: f64) -> Vec<FaceId> {
    plane_faces(m, n, d)
        .into_iter()
        .map(|f| FaceId(f as u32))
        .collect()
}

fn clamp() -> FaceMechanics {
    FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    })
}

fn problem(faces: Vec<FaceMechanics>) -> ElasticProblem {
    ElasticProblem {
        youngs: 2000.0,
        poisson: 0.3,
        faces,
    }
}

fn options(cell: f64) -> AnalysisOptions {
    AnalysisOptions {
        target_cell_size: Some(cell),
        max_active_cells: 200_000,
        ..AnalysisOptions::default()
    }
}

fn solve(m: &AnalysisModel, faces: Vec<FaceMechanics>, cell: f64) -> StructuralResult {
    solve_elastic(m, &problem(faces), &options(cell)).unwrap_or_else(|e| panic!("{e}"))
}

/// The L-bracket clamped on x = 40 and pushed along x on y = 30.
fn bracket(l: &AnalysisModel, cell: f64) -> StructuralResult {
    solve(l, bracket_faces(l), cell)
}

fn bracket_faces(l: &AnalysisModel) -> Vec<FaceMechanics> {
    let mut faces = vec![FaceMechanics::Free; l.faces.len()];
    for f in plane_faces(l, [1.0, 0.0, 0.0], 40.0) {
        faces[f] = clamp();
    }
    for f in plane_faces(l, [0.0, 1.0, 0.0], 30.0) {
        faces[f] = FaceMechanics::Traction([1.0, 0.0, 0.0]);
    }
    faces
}

/// A short C8-style shoulder (r = 3, D = 24, d = 12) in tension, half
/// thickness on a mid-plane roller: the fillet is the peak.
const SHORT_BAR: SteppedBar = SteppedBar {
    big_d: 24.0,
    step: 6.0,
    radius: 3.0,
    wide_length: 30.0,
    narrow_length: 30.0,
    thickness: 2.0,
    split_end: false,
};

fn shoulder(m: &AnalysisModel, cell: f64) -> StructuralResult {
    solve(m, shoulder_faces(m), cell)
}

fn shoulder_faces(m: &AnalysisModel) -> Vec<FaceMechanics> {
    let mut faces = vec![FaceMechanics::Free; m.faces.len()];
    for f in plane_faces(m, [-1.0, 0.0, 0.0], 0.0) {
        faces[f] = clamp();
    }
    for f in plane_faces(m, [0.0, 0.0, -1.0], 0.0) {
        faces[f] = FaceMechanics::Support(PieceSupport::Roller);
    }
    let end = SHORT_BAR.wide_length + SHORT_BAR.narrow_length;
    for f in plane_faces(m, [1.0, 0.0, 0.0], end) {
        faces[f] = FaceMechanics::Traction([10.0, 0.0, 0.0]);
    }
    faces
}

#[test]
fn an_unreported_sharp_corner_grows_under_refinement_and_is_labelled() {
    // The L-bracket with its concave edges removed: no geometric label,
    // so only growth can find the corner. Measured: 1.196 (2 -> 1 mm),
    // 1.361 (1 -> 0.5 mm).
    let o = AssessOptions::default();
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/analysis/l_bracket.step");
    let mut l = popkern_analysis_cli::step_model::model_from_step_file(&path)
        .unwrap_or_else(|e| panic!("{e}"));
    l.concave_edges.clear();
    let mut coarse = bracket(&l, 2.0);
    for cell in [1.0, 0.5] {
        let fine = bracket(&l, cell);
        let clamped = clamped_faces(&l, [1.0, 0.0, 0.0], 40.0);
        let growth =
            peak_growth(&l, &coarse, &fine, &clamped, &o).unwrap_or_else(|e| panic!("{e}"));
        let mut peaks = surface_peaks(&l, &fine, &o).unwrap_or_else(|e| panic!("{e}"));
        assert!(peaks[0].singularity.is_none(), "{cell}: geometric label");
        label_by_growth(&mut peaks, &growth, &o);
        assert!(
            matches!(peaks[0].singularity, Some(Singularity::RefinementGrowth { ratio }) if ratio > 1.15),
            // Not at the clamp: a corner label, not a support label.
            "{cell}: {growth:?}"
        );
        // The corner line is x = 5, y = 5 (both legs 5 thick).
        let p = peaks[0].von_mises.position;
        assert!(
            (p[0] - 5.0).abs() <= 1.0 && (p[1] - 5.0).abs() <= 1.0,
            "{p:?}"
        );
        coarse = fine;
    }
}

#[test]
fn a_resolved_fillet_settles_under_refinement_and_is_not_labelled() {
    // Negative control: the shoulder fillet at 4 cells (h = 0.75) grows
    // 1.020 to h/2 (3 cells: 1.078; 5 cells: 1.068). Cost of the second
    // solve, release build, one process each: 24,504 -> 166,653 dofs
    // (6.8x), 3.5 -> 18.7 s, peak RSS 248 MB -> 1.34 GB (5.4x).
    let o = AssessOptions::default();
    let m = model_of(&stepped_bar_with("short", &SHORT_BAR).unwrap_or_else(|e| panic!("{e}")));
    let (coarse, fine) = (shoulder(&m, 0.75), shoulder(&m, 0.375));
    let clamped = clamped_faces(&m, [-1.0, 0.0, 0.0], 0.0);
    let growth = peak_growth(&m, &coarse, &fine, &clamped, &o).unwrap_or_else(|e| panic!("{e}"));
    assert!(growth[0].ratio < 1.05, "{growth:?}");
    let mut peaks = surface_peaks(&m, &fine, &o).unwrap_or_else(|e| panic!("{e}"));
    label_by_growth(&mut peaks, &growth, &o);
    assert!(peaks[0].singularity.is_none(), "{:?}", peaks[0]);
}

/// PH-047, measured 2026-10-11: growth h -> h/2, local submodel (mixed
/// caps) against the whole part: the L-bracket corner (2 -> 1 mm) 1.176
/// against 1.196, with 12,474 against 24,948 dofs; the 4-cell shoulder
/// (0.75 -> 0.375) 1.015 against 1.016, with 8,409 against 166,653 dofs
/// (0.34x the first solve's 24,504). Displacement caps miss the corner
/// (0.985: the coarse part is too stiff there and they hold its
/// motion); traction caps overshoot (1.86: the coarse tractions on
/// several caps are not in balance).
#[test]
fn the_local_check_matches_the_whole_part_check_at_a_fraction_of_its_cost() {
    use popkern_analysis::assess::local_peak_growth;
    let o = AssessOptions::default();
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/analysis/l_bracket.step");
    let mut l = popkern_analysis_cli::step_model::model_from_step_file(&path)
        .unwrap_or_else(|e| panic!("{e}"));
    l.concave_edges.clear();
    let bar = model_of(&stepped_bar_with("short", &SHORT_BAR).unwrap_or_else(|e| panic!("{e}")));
    let cases = [
        (
            &l,
            bracket_faces(&l),
            2.0,
            clamped_faces(&l, [1.0, 0.0, 0.0], 40.0),
        ),
        (
            &bar,
            shoulder_faces(&bar),
            0.75,
            clamped_faces(&bar, [-1.0, 0.0, 0.0], 0.0),
        ),
    ];
    for (m, faces, cell, clamped) in cases {
        let p = problem(faces);
        let coarse = solve_elastic(m, &p, &options(cell)).unwrap_or_else(|e| panic!("{e}"));
        let fine = solve_elastic(m, &p, &options(0.5 * cell)).unwrap_or_else(|e| panic!("{e}"));
        let global = peak_growth(m, &coarse, &fine, &clamped, &o).unwrap_or_else(|e| panic!("{e}"));
        let (local, dofs) = local_peak_growth(m, &p, &coarse, &clamped, &options(cell), &o)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(
            (local[0].ratio - global[0].ratio).abs() <= 0.03 * global[0].ratio,
            "{cell}: local {local:?} vs whole part {global:?}"
        );
        // The second solve's dofs: at most the whole part's at h/2, and on
        // the shoulder (a part much larger than its peak's box) under
        // 1.5x the first solve's.
        assert!(
            dofs <= fine.report.dof,
            "{cell}: {dofs} vs {}",
            fine.report.dof
        );
        if cell < 1.0 {
            assert!(
                2 * dofs < 3 * coarse.report.dof,
                "{dofs} vs {}",
                coarse.report.dof
            );
        }
    }
}
