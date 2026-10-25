//! Corpus gate LF1 (PH-029): a linear displacement field prescribed on
//! every face of every frozen gallery part (`fixtures/gallery`, exported
//! with popkern's CLI at the pinned commit) is reproduced at every
//! vertex at five grid offsets, unless the part is pinned below to a
//! typed decline or to lattice-moment fallback cells (counted in the
//! report). The larger gallery parts (but six_axis_arm, under live
//! development) run from a fresh export
//! (`tools/gallery-export.sh`) in the ignored `extended_corpus`.

use std::path::{Path, PathBuf};

use popkern_analysis::AnalysisOptions;
use popkern_analysis::element::vector_boundary::PieceSupport;
use popkern_analysis::model::AnalysisModel;
use popkern_analysis::structural::assemble::FaceMechanics;
use popkern_analysis::structural::{ElasticProblem, solve_elastic};
use popkern_analysis_cli::run::analysis_error_kind;
use popkern_analysis_cli::step_model::model_from_step_file;

/// Offsets as fractions of the cell size.
const SHIFTS: [[f64; 3]; 5] = [
    [0.0; 3],
    [0.21, 0.0, 0.0],
    [0.0, 0.37, 0.11],
    [0.43, 0.19, 0.31],
    [0.61; 3],
];

/// Cell budget per solve: coarse on purpose, so thin features fall below
/// the cell and exercise fragment splitting.
const CELLS: usize = 3000;

const GRADIENT: [[f64; 3]; 3] = [[1e-3, 2e-3, 0.0], [-1e-3, 5e-4, 1e-3], [0.0, 1e-3, -2e-3]];

/// LF1 passes at this worst vertex error over the largest |G x|.
const EXACT: f64 = 1e-8;

/// One part's result over the five offsets.
#[derive(Debug, Clone, PartialEq)]
enum Outcome {
    /// Every offset solved: the worst relative error and the most
    /// split cells left on lattice moments at any offset.
    Solved { error: f64, lattice_cells: usize },
    /// The first decline's kind.
    Declined(String),
}

/// The worst relative vertex error, the lattice-moment cells and the
/// cell size of one solve.
fn lf1(
    model: &AnalysisModel,
    shift: [f64; 3],
    h: Option<f64>,
) -> Result<(f64, usize, f64), String> {
    lf1_at(model, shift, h, CELLS)
}

/// [`lf1`] under a budget of `cells`.
fn lf1_at(
    model: &AnalysisModel,
    shift: [f64; 3],
    h: Option<f64>,
    cells: usize,
) -> Result<(f64, usize, f64), String> {
    let fixed = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: GRADIENT,
    });
    let problem = ElasticProblem {
        youngs: 2000.0,
        poisson: 0.35,
        faces: vec![fixed; model.faces.len()],
    };
    let options = AnalysisOptions {
        target_cell_size: h,
        max_active_cells: cells,
        grid_shift: h.map_or([0.0; 3], |h| shift.map(|f| f * h)),
        cg_rel_tol: 1e-12,
        cg_max_iters: 1000,
        ..AnalysisOptions::default()
    };
    let r = solve_elastic(model, &problem, &options)
        .map_err(|e| analysis_error_kind(&e).to_string())?;
    let (mut worst, mut scale) = (0.0f64, 0.0f64);
    for (mesh, disp) in model.solids.iter().zip(&r.vertex_displacement) {
        for (v, d) in mesh.vertices.iter().zip(disp) {
            for i in 0..3 {
                let exact: f64 = (0..3).map(|j| GRADIENT[i][j] * v[j]).sum();
                worst = worst.max((d[i] - exact).abs());
                scale = scale.max(exact.abs());
            }
        }
    }
    Ok((
        worst / scale,
        r.report.split_cells_lattice_moments,
        r.report.cell_size.unwrap_or(0.0),
    ))
}

/// The outcome of the part at `path` (the budget sets h at the first
/// offset; the other offsets reuse it).
fn outcome(path: &Path) -> Outcome {
    let model = match model_from_step_file(path) {
        Ok(m) => m,
        Err(_) => return Outcome::Declined("StepModel".into()),
    };
    let (mut error, mut lattice_cells, mut h) = (0.0f64, 0usize, None);
    for shift in SHIFTS {
        match lf1(&model, shift, h) {
            Ok((e, l, cell)) => {
                error = error.max(e);
                lattice_cells = lattice_cells.max(l);
                h = Some(cell);
            }
            Err(kind) => return Outcome::Declined(kind),
        }
    }
    Outcome::Solved {
        error,
        lattice_cells,
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

/// Parts whose LF1 is pinned to something other than exact, with why.
/// Measured 2026-10-10 at 3,000 cells over the five offsets.
// None since 2026-10-25: the die opened by the STEP round trip
// (PK-783) had sphere faces exported before PK-787 fixed their sense;
// with that correction applied to the frozen file it solves and holds
// every gate (docs/architecture/analysis-integration.md fixture log).
const PINNED: [(&str, &str); 0] = [];

#[test]
fn lf1_holds_on_every_frozen_gallery_part_or_its_pinned_decline() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/analysis/gallery");
    let parts = steps(&dir);
    assert!(parts.len() >= 20, "{} gallery parts", parts.len());
    // Parts are independent: solve them on parallel threads, then check
    // in name order.
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
            // Exact, with every split cell on exact moments.
            (
                Outcome::Solved {
                    error,
                    lattice_cells,
                },
                None,
            ) => *error <= EXACT && *lattice_cells == 0,
            (Outcome::Solved { lattice_cells, .. }, Some("lattice")) => *lattice_cells > 0,
            (Outcome::Declined(kind), Some(p)) => kind == p,
            _ => false,
        };
        if !ok {
            failures.push(format!("{name}: {got:?} (pinned {pinned:?})"));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

#[test]
#[ignore = "env: needs tools/analysis/gallery-export.sh output in target/gallery-export"]
fn extended_corpus() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target/gallery-export");
    // six_axis_arm is under live development in popkern; skip it.
    let parts: Vec<_> = steps(&dir)
        .into_iter()
        .filter(|(name, _)| name != "six_axis_arm")
        .collect();
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
    for ((name, _), got) in parts.iter().zip(&outcomes) {
        eprintln!("{name}: {got:?}");
    }
}
