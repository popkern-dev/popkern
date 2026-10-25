//! Modal corpus gate (PH-092): every frozen gallery part, free of
//! supports, has six rigid modes per solid at (numerically) zero, a
//! positive first elastic mode, every mode converged, and the mass of
//! its rasterised volume, unless the part is pinned below to a typed
//! decline. (The rasterised volume is the tessellation's, X5; its gap to
//! the exact BRep volume is chord error, PH-094.)

use std::path::{Path, PathBuf};

use popkern_analysis::AnalysisOptions;
use popkern_analysis::element::vector_boundary::PieceSupport;
use popkern_analysis::model::AnalysisModel;
use popkern_analysis::structural::ElasticProblem;
use popkern_analysis::structural::assemble::FaceMechanics;
use popkern_analysis::structural::modal::{MODE_TOLERANCE, ModalResult, solve_modal};
use popkern_analysis_cli::run::analysis_error_kind;
use popkern_analysis_cli::step_model::model_from_step_file;

/// Cell budget per solve (LF1's).
const CELLS: usize = 3000;

/// Density, kg/m³ (PLA's order).
const RHO: f64 = 1240.0;

/// A rigid mode's |ω²| over the first elastic ω² for a pass. The rigid
/// modes are exact (locked, PH-093); their ω² is round-off on the
/// operator's scale, so a part with a very soft first mode reads
/// highest: piano_hinge at 6.1e-5 (first elastic 45 Hz), every other
/// part at or below 8.7e-9 (measured 2026-10-12).
const RIGID: f64 = 1e-4;

/// |mass − ρV| over ρV for a pass, V the rasterised volume.
const MASS: f64 = 1e-9;

/// One part's result.
#[derive(Debug, Clone, PartialEq)]
enum Outcome {
    /// Solved: the solids, the worst rigid ratio, the first elastic ω²,
    /// the worst residual and the relative mass error.
    Solved {
        solids: usize,
        rigid: f64,
        elastic: f64,
        residual: f64,
        mass: f64,
    },
    /// The decline's kind.
    Declined(String),
}

fn problem(faces: Vec<FaceMechanics>) -> ElasticProblem {
    ElasticProblem {
        youngs: 3000.0,
        poisson: 0.36,
        faces,
    }
}

fn options() -> AnalysisOptions {
    AnalysisOptions {
        max_active_cells: CELLS,
        cg_rel_tol: 1e-12,
        cg_max_iters: 2000,
        ..AnalysisOptions::default()
    }
}

fn modes(
    model: &AnalysisModel,
    faces: Vec<FaceMechanics>,
    count: usize,
) -> Result<ModalResult, String> {
    solve_modal(model, &problem(faces), RHO, count, &options())
        .map_err(|e| analysis_error_kind(&e).to_string())
}

fn outcome(path: &Path) -> Outcome {
    let model = match model_from_step_file(path) {
        Ok(m) => m,
        Err(_) => return Outcome::Declined("StepModel".into()),
    };
    let solids = model.solids.len();
    let rigid_count = 6 * solids;
    let r = match modes(
        &model,
        vec![FaceMechanics::Free; model.faces.len()],
        rigid_count + 1,
    ) {
        Ok(r) => r,
        Err(kind) => return Outcome::Declined(kind),
    };
    let elastic = r
        .modes
        .get(rigid_count)
        .map_or(f64::NAN, |m| m.omega_squared);
    let rigid = r.modes[..rigid_count.min(r.modes.len())]
        .iter()
        .map(|m| m.omega_squared.abs() / elastic)
        .fold(0.0, f64::max);
    let residual = r.modes.iter().map(|m| m.residual).fold(0.0, f64::max);
    let volume: f64 = r.report.volume_parity.iter().map(|v| v.rasterised).sum();
    let exact_mass = RHO * 1e-9 * volume;
    Outcome::Solved {
        solids,
        rigid,
        elastic,
        residual,
        mass: (r.mass_kg - exact_mass).abs() / exact_mass,
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
fn every_free_gallery_part_has_its_rigid_modes_and_its_mass() {
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
                    rigid,
                    elastic,
                    residual,
                    mass,
                    ..
                },
                None,
            ) => *rigid <= RIGID && *elastic > 0.0 && *residual <= MODE_TOLERANCE && *mass <= MASS,
            (Outcome::Declined(kind), Some(p)) => kind == p,
            _ => false,
        };
        if !ok {
            failures.push(format!("{name}: {got:?} (pinned {pinned:?})"));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// The negative control: clamping one face of a single-solid part
/// leaves no rigid mode, so the rigid check above is not vacuous.
#[test]
fn a_clamped_part_has_no_rigid_mode() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/analysis/gallery");
    let model =
        model_from_step_file(&dir.join("bored_block.step")).unwrap_or_else(|e| panic!("{e}"));
    let mut faces = vec![FaceMechanics::Free; model.faces.len()];
    faces[0] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    let r = modes(&model, faces, 2).unwrap_or_else(|e| panic!("{e}"));
    let first = r.modes[0].omega_squared;
    assert!(first > 1e-3 * r.modes[1].omega_squared, "{:?}", r.modes);
}
