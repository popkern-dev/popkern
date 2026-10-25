//! Stress concentrations on popkern-built parts (PH-027): C8 and C9
//! (filleted stepped bars against Roark Table 17.1 case 5) and C10
//! (U-notch, case 3a).

use popkern_analysis::model::{AnalysisModel, FaceRole, Provenance, Surface};
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
            document: "peak_fixtures".to_string(),
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

/// Largest σ₁ over the cell-resolution surface samples on `faces`.
fn peak_on(_m: &AnalysisModel, r: &StructuralResult, faces: &[usize]) -> f64 {
    r.stresses[0]
        .surface
        .iter()
        .filter(|s| faces.contains(&(s.face.0 as usize)))
        .map(|s| s.stress.max_principal)
        .fold(f64::NEG_INFINITY, f64::max)
}

/// Roark Table 17.1 case 5 K_t; `bending` selects 5b.
fn roark_stepped(step: f64, radius: f64, big_d: f64, bending: bool) -> f64 {
    let q = step / radius;
    let s = q.sqrt();
    let x = 2.0 * step / big_d;
    let c1 = if q <= 2.0 {
        1.007 + 1.000 * s - 0.031 * q
    } else {
        1.042 + 0.982 * s - 0.036 * q
    };
    let (c2, c3, c4) = match (bending, q <= 2.0) {
        (false, true) => (
            -0.114 - 0.585 * s + 0.314 * q,
            0.241 - 0.992 * s - 0.271 * q,
            -0.134 + 0.577 * s - 0.012 * q,
        ),
        (false, false) => (
            -0.074 - 0.156 * s - 0.010 * q,
            -3.418 + 1.220 * s - 0.005 * q,
            3.450 - 2.046 * s + 0.051 * q,
        ),
        (true, true) => (
            -0.270 - 2.404 * s + 0.749 * q,
            0.677 + 1.133 * s - 0.904 * q,
            -0.414 + 0.271 * s + 0.186 * q,
        ),
        (true, false) => (
            -3.599 + 1.619 * s - 0.431 * q,
            6.084 - 5.607 * s + 1.158 * q,
            -2.527 + 3.006 * s - 0.691 * q,
        ),
    };
    c1 + c2 * x + c3 * x * x + c4 * x * x * x
}

/// The C8/C9 bar: r = 3, h = 2r, D = 8r (2h/D = 0.5, h/r = 2), long
/// enough for both validity limits; half thickness with a mid-plane
/// roller at z = 0.
const BAR: SteppedBar = SteppedBar {
    big_d: 24.0,
    step: 6.0,
    radius: 3.0,
    wide_length: 108.0,
    narrow_length: 72.0,
    thickness: 4.0,
    split_end: true,
};

/// Solves a bar clamped at x = 0, with a mid-plane roller at z = 0
/// and tractions `lower` / `upper` (along x) on the two halves of its
/// split end (faces at x = `x_end` and `x_end + 1`).
fn solve_bar(
    m: &AnalysisModel,
    cell: f64,
    shift: [f64; 3],
    x_end: f64,
    lower: f64,
    upper: f64,
) -> StructuralResult {
    solve_bar_problem(m, cell, shift, x_end, lower, upper).1
}

/// [`solve_bar`], with its problem.
fn solve_bar_problem(
    m: &AnalysisModel,
    cell: f64,
    shift: [f64; 3],
    x_end: f64,
    lower: f64,
    upper: f64,
) -> (ElasticProblem, StructuralResult) {
    let mut faces = vec![FaceMechanics::Free; m.faces.len()];
    for f in plane_faces(m, [-1.0, 0.0, 0.0], 0.0) {
        faces[f] = FaceMechanics::Support(PieceSupport::Fixed {
            u0: [0.0; 3],
            gradient: [[0.0; 3]; 3],
        });
    }
    for f in plane_faces(m, [0.0, 0.0, -1.0], 0.0) {
        faces[f] = FaceMechanics::Support(PieceSupport::Roller);
    }
    let (lo, hi) = (
        plane_faces(m, [1.0, 0.0, 0.0], x_end),
        plane_faces(m, [1.0, 0.0, 0.0], x_end + 1.0),
    );
    assert_eq!((lo.len(), hi.len()), (1, 1), "split end");
    faces[lo[0]] = FaceMechanics::Traction([lower, 0.0, 0.0]);
    faces[hi[0]] = FaceMechanics::Traction([upper, 0.0, 0.0]);
    let problem = ElasticProblem {
        youngs: 2000.0,
        poisson: 0.3,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(cell),
        grid_shift: shift,
        max_active_cells: 200_000,
        ..AnalysisOptions::default()
    };
    let r = solve_elastic(m, &problem, &options).unwrap_or_else(|e| panic!("{e}"));
    (problem, r)
}

/// Faces of `m` selected by `keep`.
fn faces_where(
    m: &AnalysisModel,
    keep: impl Fn(&popkern_analysis::FaceInfo) -> bool,
) -> Vec<usize> {
    m.faces
        .iter()
        .enumerate()
        .filter(|(_, f)| keep(f))
        .map(|(i, _)| i)
        .collect()
}

/// (fillet peak σ₁ / σ_nom, raw peak / σ_nom) for the stepped bar.
fn stepped_ratio(m: &AnalysisModel, cell: f64, shift: [f64; 3], bending: bool) -> (f64, f64) {
    let s = 10.0;
    let x_end = BAR.wide_length + BAR.narrow_length;
    let r = solve_bar(m, cell, shift, x_end, if bending { -s } else { s }, s);
    let fillets = faces_where(m, |f| matches!(f.role, Some(FaceRole::Fillet { .. })));
    assert_eq!(fillets.len(), 2);
    // Couple of ±s on the two halves of the narrow end: M = s t d²/4,
    // σ_nom = 6M/(t d²) = 1.5 s.
    let nominal = if bending { 1.5 * s } else { s };
    (
        peak_on(m, &r, &fillets) / nominal,
        r.stresses[0].raw_peak_principal / nominal,
    )
}

/// C10 bar: D = 12, one semicircular notch r = h = 3 (h/D = 0.25).
const NOTCH: (f64, f64, f64, f64) = (12.0, 3.0, 96.0, 4.0);

/// Roark Table 17.1 case 3 semicircular-notch K_t (3a tension along the
/// net-section centroid, 3b in-plane bending).
fn roark_notch(h_over_d: f64, bending: bool) -> f64 {
    let x = h_over_d;
    if bending {
        2.988 - 7.735 * x + 10.674 * x * x - 4.927 * x * x * x
    } else {
        2.988 - 7.300 * x + 9.742 * x * x - 4.429 * x * x * x
    }
}

/// (notch peak σ₁ / σ_nom, raw / σ_nom) for the U-notch bar.
fn notch_ratio(m: &AnalysisModel, cell: f64, shift: [f64; 3], bending: bool) -> (f64, f64) {
    let (big_d, h, length, _) = NOTCH;
    let d = big_d - h;
    let s = 10.0;
    let (lower, upper, nominal) = if bending {
        // M = s t D²/4 on the net section: σ_nom = 6M/(t d²).
        (-s, s, 1.5 * s * (big_d / d) * (big_d / d))
    } else {
        // Resultant on the net-section centroid y = −h/2 (the figure's
        // load line at d/2): P = s t D, σ_nom = P/(t d).
        let k = 2.0 * h / big_d;
        (s * (1.0 + k), s * (1.0 - k), s * big_d / d)
    };
    let r = solve_bar(m, cell, shift, length, lower, upper);
    let notch = faces_where(
        m,
        |f| matches!(f.surface, Surface::Cylinder { radius, .. } if (radius - h).abs() < 1e-9),
    );
    assert_eq!(notch.len(), 1);
    (
        peak_on(m, &r, &notch) / nominal,
        r.stresses[0].raw_peak_principal / nominal,
    )
}

const SHIFTS: [[f64; 3]; 4] = [
    [0.0; 3],
    [0.31, 0.17, 0.0],
    [0.53, 0.71, 0.0],
    [0.11, 0.43, 0.0],
];

/// Ratios to K_t over the four grid shifts: (recovered, raw).
fn over_shifts(f: impl Fn([f64; 3]) -> (f64, f64), cell: f64, kt: f64) -> Vec<(f64, f64)> {
    SHIFTS
        .iter()
        .map(|s| {
            let (p, raw) = f(s.map(|x| x * cell));
            (p / kt, raw / kt)
        })
        .collect()
}

fn stepped() -> AnalysisModel {
    model_of(&stepped_bar_with("c8", &BAR).unwrap_or_else(|e| panic!("{e}")))
}

fn notched() -> AnalysisModel {
    let (big_d, radius, length, thickness) = NOTCH;
    model_of(
        &popkern_analysis_fixtures::u_notch_bar("c10", big_d, radius, length, thickness)
            .unwrap_or_else(|e| panic!("{e}")),
    )
}

#[test]
fn c8_shoulder_fillet_in_tension_within_10_percent_at_4_cells() {
    let m = stepped();
    let kt = roark_stepped(BAR.step, BAR.radius, BAR.big_d, false);
    let cell = BAR.radius / 4.0;
    let r = over_shifts(|s| stepped_ratio(&m, cell, s, false), cell, kt);
    // Measured 0.933–1.033.
    assert!(r.iter().all(|x| (x.0 - 1.0).abs() <= 0.10), "{r:?}");
    // Negative controls: the raw fragment-centroid peak misses by more
    // than 10% at some shift (measured up to +26%), and one cell per
    // radius reads about 40% low.
    assert!(r.iter().any(|x| (x.1 - 1.0).abs() > 0.10), "{r:?}");
    let coarse = stepped_ratio(&m, BAR.radius, [0.0; 3], false).0 / kt;
    assert!(coarse < 0.9, "{coarse}");
}

#[test]
fn c9_shoulder_fillet_in_bending_within_10_percent_at_4_cells() {
    let m = stepped();
    let kt = roark_stepped(BAR.step, BAR.radius, BAR.big_d, true);
    let cell = BAR.radius / 4.0;
    let r = over_shifts(|s| stepped_ratio(&m, cell, s, true), cell, kt);
    // Measured 0.960–1.042.
    assert!(r.iter().all(|x| (x.0 - 1.0).abs() <= 0.10), "{r:?}");
    let coarse = stepped_ratio(&m, BAR.radius, [0.0; 3], true).0 / kt;
    assert!(coarse < 0.9, "{coarse}");
}

#[test]
fn c10_u_notch_bending_at_4_cells_and_tension_at_6() {
    let m = notched();
    let (big_d, radius, _, _) = NOTCH;
    let bending = roark_notch(radius / big_d, true);
    let cell = radius / 4.0;
    let r = over_shifts(|s| notch_ratio(&m, cell, s, true), cell, bending);
    // Measured 0.978–1.016.
    assert!(r.iter().all(|x| (x.0 - 1.0).abs() <= 0.10), "{r:?}");
    // Tension converges more slowly (0.861–0.916 at 4 cells, 0.895–0.935
    // at 5; PH-036): gated at 6, measured 0.915–0.947.
    let tension = roark_notch(radius / big_d, false);
    let cell = radius / 6.0;
    let r = over_shifts(|s| notch_ratio(&m, cell, s, false), cell, tension);
    assert!(r.iter().all(|x| (x.0 - 1.0).abs() <= 0.10), "{r:?}");
    // Negative control: raw centroid peaks read 17–30% high there.
    assert!(r.iter().any(|x| (x.1 - 1.0).abs() > 0.10), "{r:?}");
}

#[test]
fn fr1_fillet_resolution_warns_below_4_cells_only() {
    use popkern_analysis::assess::fillet::FilletStatus;
    use popkern_analysis::assess::{AssessOptions, assess};
    let m = stepped();
    let pla = popkern_analysis::material::preset("PLA (FDM)").unwrap_or_else(|e| panic!("{e}"));
    let x_end = BAR.wide_length + BAR.narrow_length;
    for cells in [1.0, 2.0, 3.0, 4.0, 5.0] {
        let r = solve_bar(&m, BAR.radius / cells, [0.0; 3], x_end, 10.0, 10.0);
        let a = assess(&m, &r, pla, &AssessOptions::default()).unwrap_or_else(|e| panic!("{e}"));
        // Both shoulder fillets carry the peak, so both are checked.
        assert_eq!(a.fillets.len(), 2, "{cells}: {:?}", a.fillets);
        for f in &a.fillets {
            assert!((f.cells - cells).abs() < 1e-9);
            match f.status {
                FilletStatus::UnderResolved { expected_low_by } => {
                    assert!(cells < 4.0, "{cells}");
                    assert!(expected_low_by > 0.0);
                }
                FilletStatus::Resolved => assert!(cells >= 4.0, "{cells}"),
                FilletStatus::Sharp => panic!("{cells}: sharp"),
            }
        }
        // A filleted shoulder has no concave sharp edge: no singularity.
        assert!(a.peaks[0].singularity.is_none(), "{cells}");
    }
}

#[test]
fn the_l_bracket_corner_peak_is_labelled_a_singularity() {
    use popkern_analysis::assess::{AssessOptions, Singularity, assess};
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/analysis/l_bracket.step");
    let mut l = popkern_analysis_cli::step_model::model_from_step_file(&path)
        .unwrap_or_else(|e| panic!("{e}"));
    let pla = popkern_analysis::material::preset("PLA (FDM)").unwrap_or_else(|e| panic!("{e}"));
    let corner = l.concave_edges[0].faces;
    let solve = |l: &AnalysisModel, cell: f64| {
        let mut faces = vec![FaceMechanics::Free; l.faces.len()];
        for f in plane_faces(l, [1.0, 0.0, 0.0], 40.0) {
            faces[f] = FaceMechanics::Support(PieceSupport::Fixed {
                u0: [0.0; 3],
                gradient: [[0.0; 3]; 3],
            });
        }
        for f in plane_faces(l, [0.0, 1.0, 0.0], 30.0) {
            faces[f] = FaceMechanics::Traction([1.0, 0.0, 0.0]);
        }
        let problem = ElasticProblem {
            youngs: 2000.0,
            poisson: 0.3,
            faces,
        };
        let options = AnalysisOptions {
            target_cell_size: Some(cell),
            ..AnalysisOptions::default()
        };
        let r = solve_elastic(l, &problem, &options).unwrap_or_else(|e| panic!("{e}"));
        assess(l, &r, pla, &AssessOptions::default()).unwrap_or_else(|e| panic!("{e}"))
    };
    for cell in [2.0, 1.0, 0.5] {
        let a = solve(&l, cell);
        assert_eq!(
            a.peaks[0].singularity,
            Some(Singularity::ConcaveEdge { faces: corner }),
            "{cell}: {:?}",
            a.peaks[0]
        );
    }
    // Negative control: without the BRep's concave edge the label is
    // gone (it comes from geometry, not from the stress statistics).
    l.concave_edges.clear();
    assert!(solve(&l, 1.0).peaks[0].singularity.is_none());
}

/// The largest refined (`--refine`) principal peak over `faces` (radius
/// `radius`) of the bar at cell size `cell` and each shift, over
/// `nominal · kt`.
fn refined_over_shifts(
    m: &AnalysisModel,
    (faces, radius): (&[usize], f64),
    (cell, x_end): (f64, f64),
    (lower, upper): (f64, f64),
    nominal_kt: f64,
) -> Vec<f64> {
    use popkern_analysis::structural::submodel::refine_raiser;
    SHIFTS
        .iter()
        .map(|s| {
            let shift = s.map(|x| x * cell);
            let (problem, r) = solve_bar_problem(m, cell, shift, x_end, lower, upper);
            let options = AnalysisOptions {
                target_cell_size: Some(cell),
                grid_shift: shift,
                max_active_cells: 200_000,
                ..AnalysisOptions::default()
            };
            faces
                .iter()
                .map(|&f| {
                    let face = popkern_analysis::FaceId(f as u32);
                    refine_raiser(m, &problem, &r, face, radius, &options)
                        .unwrap_or_else(|e| panic!("{e}"))
                        .peak_principal
                })
                .fold(0.0, f64::max)
                / nominal_kt
        })
        .collect()
}

/// PH-036: C8/C9 at 3 cells per radius through submodels. Measured
/// 2026-10-11 over Roark: tension 1.010–1.021 (coarse 0.856–1.002),
/// bending 1.027–1.051 (coarse 0.893–1.047).
#[test]
fn c8_c9_fillets_refine_to_roark_from_3_cells() {
    let m = stepped();
    let fillets = faces_where(&m, |f| matches!(f.role, Some(FaceRole::Fillet { .. })));
    let (s, x_end) = (10.0, BAR.wide_length + BAR.narrow_length);
    for bending in [false, true] {
        let kt = roark_stepped(BAR.step, BAR.radius, BAR.big_d, bending);
        let nominal = if bending { 1.5 * s } else { s };
        let r = refined_over_shifts(
            &m,
            (&fillets, BAR.radius),
            (BAR.radius / 3.0, x_end),
            (if bending { -s } else { s }, s),
            nominal * kt,
        );
        assert!(
            r.iter().all(|x| (x - 1.0).abs() <= 0.10),
            "{bending}: {r:?}"
        );
    }
}

/// PH-036: C10 tension at 4 cells per radius through a submodel.
/// Measured 2026-10-11 over Roark 3a: 0.933–0.937 (coarse 0.861–0.916).
#[test]
fn c10_tension_refines_to_roark_from_4_cells() {
    let m = notched();
    let (big_d, h, length, _) = NOTCH;
    let notch = faces_where(
        &m,
        |f| matches!(f.surface, Surface::Cylinder { radius, .. } if (radius - h).abs() < 1e-9),
    );
    let (s, k) = (10.0, 2.0 * h / big_d);
    let nominal = s * big_d / (big_d - h);
    let r = refined_over_shifts(
        &m,
        (&notch, h),
        (h / 4.0, length),
        (s * (1.0 + k), s * (1.0 - k)),
        nominal * roark_notch(h / big_d, false),
    );
    assert!(r.iter().all(|x| (x - 1.0).abs() <= 0.10), "{r:?}");
}
