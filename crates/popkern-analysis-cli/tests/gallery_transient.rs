//! Transient corpus gate (PH-098): on every frozen gallery part, a
//! uniform source heating a fully insulated part rises uniformly at
//! P/(ρcV), which backward Euler reproduces exactly: every vertex reads
//! T₀ + P·t/(ρcV) and the stored energy is P·t, unless the part is
//! pinned below to a typed decline. It checks the capacity matrix and
//! the source against each other on every cut cell and fragment of the
//! corpus.

use std::path::{Path, PathBuf};

use popkern_analysis::AnalysisOptions;
use popkern_analysis::model::AnalysisModel;
use popkern_analysis::thermal::transient::{TransientProblem, solve_transient};
use popkern_analysis::thermal::{ConductionProblem, FaceCondition};
use popkern_analysis_cli::run::analysis_error_kind;
use popkern_analysis_cli::step_model::model_from_step_file;

/// Cell budget per solve (LF1's).
const CELLS: usize = 3000;

/// Aluminium-like: conductivity W/(m·K), density kg/m³, specific heat
/// J/(kg·K).
const K: f64 = 167.0;
const RHO: f64 = 2700.0;
const C: f64 = 900.0;

/// Power, W; step, s; steps.
const POWER: f64 = 1.0;
const STEP: f64 = 0.01;
const STEPS: usize = 5;

/// Worst vertex error over the exact rise for a pass. Measured
/// 2026-10-12: piano_hinge 6.6e-8 at every CG tolerance from 1e-8 to
/// 1e-12 (so not the solver), every other part below 1e-8. Before the
/// steps were solved as increments (PH-100) the tolerance applied to
/// the 20 °C level, and hinged_box's 5e-4 K rise read 5.8e-7 at 1e-12
/// and 3.6e-3 at 1e-8.
const EXACT: f64 = 1e-7;

/// |stored − P·t| over P·t for a pass. Measured: piano_hinge 1.9e-9,
/// every other part below 1e-11. The source is spread by mesh volume
/// and the capacity integrates the moments, which agree to the volume
/// parity tolerance (1e-9, X5).
const STORED: f64 = 1e-8;

/// One part's result.
#[derive(Debug, Clone, PartialEq)]
enum Outcome {
    /// Solved: worst relative vertex error and relative stored-energy
    /// error at the last step.
    Solved { vertex: f64, stored: f64 },
    /// The decline's kind.
    Declined(String),
}

/// The relative vertex and stored-energy errors with `power` W heating
/// the part (the exact rise uses [`POWER`]).
fn heating(model: &AnalysisModel, power: f64) -> Result<(f64, f64), String> {
    let problem = TransientProblem {
        conduction: ConductionProblem {
            conductivity: K,
            faces: vec![FaceCondition::Insulated; model.faces.len()],
            source_w: power,
            body_sources: Vec::new(),
            hdt_c: None,
        },
        density: RHO,
        specific_heat: C,
        initial_c: 20.0,
        step_s: STEP,
        steps: STEPS,
    };
    let options = AnalysisOptions {
        max_active_cells: CELLS,
        cg_rel_tol: 1e-12,
        cg_max_iters: 1000,
        ..AnalysisOptions::default()
    };
    let r = solve_transient(model, &problem, &options)
        .map_err(|e| analysis_error_kind(&e).to_string())?;
    let volume: f64 = r
        .last
        .report
        .volume_parity
        .iter()
        .map(|v| v.rasterised)
        .sum();
    let t = STEP * STEPS as f64;
    let rise = POWER * t / (1e-9 * RHO * C * volume);
    let vertex = r
        .last
        .vertex_temperature
        .iter()
        .flatten()
        .map(|temp| (temp - 20.0 - rise).abs() / rise)
        .fold(0.0, f64::max);
    let stored = r.steps.last().map_or(f64::NAN, |s| s.stored_j);
    Ok((vertex, (stored - POWER * t).abs() / (POWER * t)))
}

fn outcome(path: &Path) -> Outcome {
    let model = match model_from_step_file(path) {
        Ok(m) => m,
        Err(_) => return Outcome::Declined("StepModel".into()),
    };
    match heating(&model, POWER) {
        Ok((vertex, stored)) => Outcome::Solved { vertex, stored },
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
fn an_insulated_part_heats_uniformly_at_its_capacity_on_every_gallery_part() {
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
            (Outcome::Solved { vertex, stored }, None) => *vertex <= EXACT && *stored <= STORED,
            (Outcome::Declined(kind), Some(p)) => kind == p,
            _ => false,
        };
        if !ok {
            failures.push(format!("{name}: {got:?} (pinned {pinned:?})"));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// The negative control: 1% more power than the exact rise assumes
/// misses both checks on a gallery part.
#[test]
fn a_wrong_power_misses_the_gate() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/analysis/gallery");
    let model =
        model_from_step_file(&dir.join("bored_block.step")).unwrap_or_else(|e| panic!("{e}"));
    let (vertex, stored) = heating(&model, 1.01 * POWER).unwrap_or_else(|e| panic!("{e}"));
    assert!(
        vertex > 1e3 * EXACT && stored > 1e3 * STORED,
        "{vertex} {stored}"
    );
}
