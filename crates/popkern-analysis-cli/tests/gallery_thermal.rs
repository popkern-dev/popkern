//! Thermal corpus gate (PH-091): on every frozen gallery part, a linear
//! temperature held on every face is reproduced at every vertex at the
//! five LF1 grid offsets, and with a uniform source and convection on
//! every face the heat leaving equals the heat generated, unless the part
//! is pinned below to a typed decline. Measured 2026-10-12 on the 22
//! solved parts: worst vertex error 1.7e-9 (centrifugal_impeller),
//! worst imbalance 1.7e-12 W per W (piano_hinge).

use std::path::{Path, PathBuf};

use popkern_analysis::AnalysisOptions;
use popkern_analysis::model::AnalysisModel;
use popkern_analysis::thermal::{ConductionProblem, FaceCondition, solve_conduction};
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

/// The held field: T = 20 + G·x, °C with x in mm.
const GRADIENT: [f64; 3] = [0.7, -0.4, 1.1];

/// Worst vertex error over the largest |T − 20| for a pass.
const EXACT: f64 = 1e-8;

/// |heat in − heat out| over the heat in for a pass.
const BALANCED: f64 = 1e-6;

/// One part's result.
#[derive(Debug, Clone, PartialEq)]
enum Outcome {
    /// Every solve ran: the worst relative vertex error over the
    /// offsets, and the relative imbalance of the convective solve.
    Solved { error: f64, imbalance: f64 },
    /// The first decline's kind.
    Declined(String),
}

fn options(h: Option<f64>, shift: [f64; 3]) -> AnalysisOptions {
    AnalysisOptions {
        target_cell_size: h,
        max_active_cells: CELLS,
        grid_shift: h.map_or([0.0; 3], |h| shift.map(|f| f * h)),
        cg_rel_tol: 1e-12,
        cg_max_iters: 1000,
        ..AnalysisOptions::default()
    }
}

/// The worst relative vertex error of the held linear field, and the
/// cell size used.
fn linear(
    model: &AnalysisModel,
    shift: [f64; 3],
    h: Option<f64>,
    gradient: [f64; 3],
) -> Result<(f64, f64), String> {
    let held = FaceCondition::Fixed {
        t0_c: 20.0,
        gradient,
    };
    let problem = ConductionProblem {
        conductivity: 0.2,
        faces: vec![held; model.faces.len()],
        source_w: 0.0,
        body_sources: Vec::new(),
        hdt_c: None,
    };
    let r = solve_conduction(model, &problem, &options(h, shift))
        .map_err(|e| analysis_error_kind(&e).to_string())?;
    let (mut worst, mut scale) = (0.0f64, 0.0f64);
    for (mesh, temps) in model.solids.iter().zip(&r.vertex_temperature) {
        for (v, t) in mesh.vertices.iter().zip(temps) {
            let exact: f64 = (0..3).map(|j| GRADIENT[j] * v[j]).sum();
            worst = worst.max((t - 20.0 - exact).abs());
            scale = scale.max(exact.abs());
        }
    }
    Ok((worst / scale, r.report.cell_size.unwrap_or(0.0)))
}

/// The relative heat imbalance of 1 W generated in the part with
/// convection on every face.
fn balance(model: &AnalysisModel, h: Option<f64>) -> Result<f64, String> {
    let problem = ConductionProblem {
        conductivity: 0.2,
        faces: vec![
            FaceCondition::Convection {
                h: 10.0,
                ambient_c: 20.0,
            };
            model.faces.len()
        ],
        source_w: 1.0,
        body_sources: Vec::new(),
        hdt_c: None,
    };
    let r = solve_conduction(model, &problem, &options(h, [0.0; 3]))
        .map_err(|e| analysis_error_kind(&e).to_string())?;
    Ok(r.heat.imbalance().abs() / r.heat.source)
}

fn outcome(path: &Path) -> Outcome {
    let model = match model_from_step_file(path) {
        Ok(m) => m,
        Err(_) => return Outcome::Declined("StepModel".into()),
    };
    let (mut error, mut h) = (0.0f64, None);
    for shift in SHIFTS {
        match linear(&model, shift, h, GRADIENT) {
            Ok((e, cell)) => {
                error = error.max(e);
                h = Some(cell);
            }
            Err(kind) => return Outcome::Declined(kind),
        }
    }
    match balance(&model, h) {
        Ok(imbalance) => Outcome::Solved { error, imbalance },
        Err(kind) => Outcome::Declined(kind),
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
fn a_linear_field_is_exact_and_heat_balances_on_every_gallery_part() {
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
            (Outcome::Solved { error, imbalance }, None) => {
                *error <= EXACT && *imbalance <= BALANCED
            }
            (Outcome::Declined(kind), Some(p)) => kind == p,
            _ => false,
        };
        if !ok {
            failures.push(format!("{name}: {got:?} (pinned {pinned:?})"));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// The negative control: a held field that is not the one checked
/// (one gradient component off by 1%) misses the gate on a gallery part.
#[test]
fn a_perturbed_field_misses_the_gate() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/analysis/gallery");
    let model =
        model_from_step_file(&dir.join("bored_block.step")).unwrap_or_else(|e| panic!("{e}"));
    let mut off = GRADIENT;
    off[2] *= 1.01;
    let (error, _) = linear(&model, [0.0; 3], None, off).unwrap_or_else(|e| panic!("{e}"));
    assert!(error > 1e3 * EXACT, "{error}");
}
