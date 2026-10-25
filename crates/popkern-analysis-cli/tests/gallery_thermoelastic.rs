//! Thermoelastic corpus gate (PH-095): on every frozen gallery part, a
//! uniform temperature rise θ with the free expansion u = αθ·x held on
//! every face is stress-free: every vertex reads that displacement and
//! no stress, at the five LF1 grid offsets, unless the part is pinned
//! below to a typed decline. It checks the thermal load against the
//! stiffness on every cut cell, fragment and enhanced mode of the
//! corpus.

use std::path::{Path, PathBuf};

use popkern_analysis::AnalysisOptions;
use popkern_analysis::element::vector_boundary::PieceSupport;
use popkern_analysis::model::AnalysisModel;
use popkern_analysis::structural::assemble::FaceMechanics;
use popkern_analysis::structural::infill::Infill;
use popkern_analysis::structural::thermoelastic::ThermalLoad;
use popkern_analysis::structural::{ElasticProblem, solve_with};
use popkern_analysis_cli::run::analysis_error_kind;
use popkern_analysis_cli::step_model::model_from_step_file;

/// Offsets as fractions of the cell size (LF1's).
const SHIFTS: [[f64; 3]; 5] = [
    [0.0; 3],
    [0.21, 0.0, 0.0],
    [0.0, 0.37, 0.11],
    [0.43, 0.19, 0.31],
    [0.61; 3],
];

/// Cell budget per solve (LF1's).
const CELLS: usize = 3000;

/// Young's modulus, MPa.
const E: f64 = 2000.0;
/// Poisson's ratio.
const NU: f64 = 0.35;
/// Expansion, 1/K.
const ALPHA: f64 = 1e-4;
/// Uniform rise, K.
const THETA: f64 = 50.0;

/// Worst vertex displacement error over the largest |αθ·x| for a pass.
const EXACT: f64 = 1e-8;

/// Largest vertex stress component over Eαθ (the stress of the same
/// rise fully restrained along one axis) for a pass. Components, not
/// von Mises: a held expansion without the rise is hydrostatic, which
/// von Mises cannot see. Measured 2026-10-12: piano_hinge 6.0e-8 and
/// living_hinges 1.4e-8 (thin flexures, solver round-off), every other
/// part at most 1e-8; a wrong thermal load reads 1e-2 to 1.
const STRESS_FREE: f64 = 1e-6;

/// One part's result over the five offsets.
#[derive(Debug, Clone, PartialEq)]
enum Outcome {
    /// Every offset solved: worst relative displacement error and
    /// largest relative stress component.
    Solved { displacement: f64, stress: f64 },
    /// The first decline's kind.
    Declined(String),
}

/// One solve's relative displacement error, relative stress and cell
/// size; `rise` is the temperature rise applied, the held field is
/// always αθ·x.
fn expansion(
    model: &AnalysisModel,
    shift: [f64; 3],
    h: Option<f64>,
    rise: f64,
) -> Result<(f64, f64, f64), String> {
    expansion_with(model, shift, h, rise, None).map(|(d, s, h, _)| (d, s, h))
}

/// [`expansion`] with an infill, also returning the core's volume.
fn expansion_with(
    model: &AnalysisModel,
    shift: [f64; 3],
    h: Option<f64>,
    rise: f64,
    infill: Option<&Infill>,
) -> Result<(f64, f64, f64, f64), String> {
    let g = ALPHA * THETA;
    let held = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[g, 0.0, 0.0], [0.0, g, 0.0], [0.0, 0.0, g]],
    });
    let problem = ElasticProblem {
        youngs: E,
        poisson: NU,
        faces: vec![held; model.faces.len()],
    };
    let temperature = |_: [f64; 3]| Some(20.0 + rise);
    let load = ThermalLoad {
        alpha: ALPHA,
        reference_c: 20.0,
        temperature: &temperature,
    };
    let options = AnalysisOptions {
        target_cell_size: h,
        max_active_cells: CELLS,
        grid_shift: h.map_or([0.0; 3], |h| shift.map(|f| f * h)),
        cg_rel_tol: 1e-12,
        cg_max_iters: 1000,
        ..AnalysisOptions::default()
    };
    let r = solve_with(model, &problem, Some(&load), infill, &options)
        .map_err(|e| analysis_error_kind(&e).to_string())?;
    let (mut worst, mut scale) = (0.0f64, 0.0f64);
    for (mesh, disp) in model.solids.iter().zip(&r.vertex_displacement) {
        for (v, d) in mesh.vertices.iter().zip(disp) {
            for i in 0..3 {
                worst = worst.max((d[i] - g * v[i]).abs());
                scale = scale.max((g * v[i]).abs());
            }
        }
    }
    let stress = r
        .stresses
        .iter()
        .flat_map(|s| &s.vertices)
        .flat_map(|v| v.stress)
        .map(f64::abs)
        .fold(0.0, f64::max);
    Ok((
        worst / scale,
        stress / (E * ALPHA * THETA),
        r.report.cell_size.unwrap_or(0.0),
        r.report.infill_core_volume,
    ))
}

fn outcome(path: &Path) -> Outcome {
    let model = match model_from_step_file(path) {
        Ok(m) => m,
        Err(_) => return Outcome::Declined("StepModel".into()),
    };
    let (mut displacement, mut stress, mut h) = (0.0f64, 0.0f64, None);
    for shift in SHIFTS {
        match expansion(&model, shift, h, THETA) {
            Ok((d, s, cell)) => {
                displacement = displacement.max(d);
                stress = stress.max(s);
                h = Some(cell);
            }
            Err(kind) => return Outcome::Declined(kind),
        }
    }
    Outcome::Solved {
        displacement,
        stress,
    }
}

fn steps(dir: &Path) -> Vec<(String, PathBuf)> {
    let mut files: Vec<(String, PathBuf)> = std::fs::read_dir(dir)
        .unwrap_or_else(|e| panic!("{}: {e}", dir.display()))
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "step"))
        .map(|p| {
            let name = p
                .file_stem()
                .map_or(String::new(), |s| s.to_string_lossy().into_owned());
            (name, p)
        })
        .collect();
    files.sort();
    files
}

/// Parts pinned to a decline, with why (LF1's pins).
// None since 2026-10-25: the die opened by the STEP round trip
// (PK-783) had sphere faces exported before PK-787 fixed their sense;
// with that correction applied to the frozen file it solves and holds
// every gate (docs/architecture/analysis-integration.md fixture log).
const PINNED: [(&str, &str); 0] = [];

#[test]
fn a_uniform_rise_held_at_free_expansion_is_stress_free_on_every_gallery_part() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/analysis/gallery");
    let parts = steps(&dir);
    assert!(parts.len() >= 20, "{} gallery parts", parts.len());
    let outcomes: Vec<Outcome> = std::thread::scope(|scope| {
        let handles: Vec<_> = parts
            .iter()
            .map(|(_, path)| scope.spawn(move || outcome(path)))
            .collect();
        handles
            .into_iter()
            .map(|h| {
                h.join()
                    .unwrap_or_else(|_| Outcome::Declined("panic".into()))
            })
            .collect()
    });
    let mut failures = Vec::new();
    for ((name, _), got) in parts.iter().zip(&outcomes) {
        let pinned = PINNED.iter().find(|(n, _)| n == name).map(|(_, k)| *k);
        let ok = match (got, pinned) {
            (
                Outcome::Solved {
                    displacement,
                    stress,
                },
                None,
            ) => *displacement <= EXACT && *stress <= STRESS_FREE,
            (Outcome::Declined(kind), Some(p)) => kind == p,
            _ => false,
        };
        if !ok {
            failures.push(format!("{name}: {got:?} (pinned {pinned:?})"));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// The negative control: the same held expansion without the rise is
/// stressed, hydrostatically at E/(1 − 2ν) × αθ (the faces stretch an
/// unheated part), so the stress check is not vacuous.
#[test]
fn the_held_expansion_without_the_rise_is_stressed() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/analysis/gallery");
    let model =
        model_from_step_file(&dir.join("bored_block.step")).unwrap_or_else(|e| panic!("{e}"));
    let (_, stress, _) = expansion(&model, [0.0; 3], None, 0.0).unwrap_or_else(|e| panic!("{e}"));
    assert!((stress - 1.0 / (1.0 - 2.0 * NU)).abs() < 0.01, "{stress}");
}

/// With a gyroid core (walls and core share α), the held free expansion
/// is still stress-free: every part reads it exactly, which covers full,
/// partial and mixed-material cells. At LF1's budget with 0.3 mm walls
/// five parts have a core (measured 2026-10-12); the rest are all
/// walls, since cut cells always are.
#[test]
fn an_infilled_part_at_free_expansion_is_stress_free_on_every_gallery_part() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/analysis/gallery");
    let parts = steps(&dir);
    let infill = Infill {
        density: 0.25,
        walls_mm: 0.3,
        line_width_mm: 0.4,
    };
    let outcomes: Vec<Result<(f64, f64, f64, f64), String>> = std::thread::scope(|scope| {
        let handles: Vec<_> = parts
            .iter()
            .map(|(_, path)| {
                scope.spawn(move || {
                    let model = model_from_step_file(path).map_err(|_| "StepModel".to_string())?;
                    expansion_with(&model, [0.43, 0.19, 0.31], None, THETA, Some(&infill))
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().unwrap_or_else(|_| Err("panic".into())))
            .collect()
    });
    let mut failures = Vec::new();
    let mut cored = 0;
    for ((name, _), got) in parts.iter().zip(&outcomes) {
        let pinned = PINNED.iter().find(|(n, _)| n == name).map(|(_, k)| *k);
        let ok = match (got, pinned) {
            (Ok((displacement, stress, _, core)), None) => {
                cored += usize::from(*core > 0.0);
                *displacement <= EXACT && *stress <= STRESS_FREE
            }
            (Err(kind), Some(p)) => kind == p,
            _ => false,
        };
        if !ok {
            failures.push(format!("{name}: {got:?} (pinned {pinned:?})"));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
    assert!(cored >= 5, "only {cored} parts have a core");
}
