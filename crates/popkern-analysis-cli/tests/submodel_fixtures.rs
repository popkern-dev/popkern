//! PH-063: submodels around the PH-050 raisers against CalculiX
//! (boss_plate's 0.08 fillet 28.4376 MPa, bored_block's 0.18 bore
//! 109.142 MPa at their CalculiX reference meshes).

use std::path::{Path, PathBuf};

use popkern_analysis::assess::{AssessOptions, assess};
use popkern_analysis::structural::from_setup::{problem_from_setup, solve_structural};
use popkern_analysis::structural::submodel::refine_raiser_with;
use popkern_analysis_brep::resolve_faces;
use popkern_analysis_cli::setup_file::parse_setup;
use popkern_analysis_cli::step_model::model_from_step_file;

#[path = "support/report.rs"]
mod report;

fn repo(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
}

/// The part's case (optionally at cell size `h`), its unresolved
/// raiser's submodel peaks for each box reach in `radii`, over
/// `reference`; and the coarse peak over it.
fn refined(part: &str, h: Option<f64>, radii: &[f64], reference: f64) -> (f64, Vec<f64>) {
    let mut text = std::fs::read_to_string(repo(&format!(
        "tools/analysis/calculix/cases/{part}.setup.toml"
    )))
    .unwrap_or_else(|e| panic!("{e}"));
    if let Some(h) = h {
        text += &format!("\n[options]\ntarget_cell_size = {h}\n");
    }
    let setup = parse_setup(&text).unwrap_or_else(|e| panic!("{e}"));
    let model = model_from_step_file(&repo(&format!(
        "tests/fixtures/analysis/gallery/{part}.step"
    )))
    .unwrap_or_else(|e| panic!("{e}"));
    let resolution =
        resolve_faces(&model, &setup.faces, &setup.tolerance).unwrap_or_else(|e| panic!("{e:?}"));
    let (result, material) =
        solve_structural(&model, &setup, &resolution.faces).unwrap_or_else(|e| panic!("{e}"));
    let assessment = assess(&model, &result, material, &AssessOptions::default())
        .unwrap_or_else(|e| panic!("{e}"));
    let raiser = assessment
        .fillets
        .iter()
        .chain(&assessment.bores)
        .max_by(|a, b| a.peak_von_mises.total_cmp(&b.peak_von_mises))
        .unwrap_or_else(|| panic!("no raiser"));
    let (problem, _) =
        problem_from_setup(&model, &setup, &resolution.faces).unwrap_or_else(|e| panic!("{e}"));
    let peaks = radii
        .iter()
        .map(|&r| {
            refine_raiser_with(
                &model,
                &problem,
                &result,
                (raiser.face, raiser.radius),
                (r, 8.0),
                &setup.options,
            )
            .unwrap_or_else(|e| panic!("{e}"))
            .peak_von_mises
                / reference
        })
        .collect();
    (raiser.peak_von_mises / reference, peaks)
}

/// Measured 2026-10-11 (over CalculiX's peak): bored_block at 3 cells
/// per bore radius reads 0.803, its submodel at 8 cells 0.967 (the box,
/// 4 radii, already holds the whole small part, so it has no caps and is
/// the whole part at the finer cell size; 8 radii the same); boss_plate
/// at the default budget (fillet 1.5 cells) 0.650, its submodel 0.992
/// at 4 radii and 0.994 at 6, with caps on its box.
#[test]
fn submodels_recover_the_calculix_peaks_at_any_box_size() {
    let (coarse, peaks) = refined("bored_block", Some(0.06), &[4.0], 109.142);
    assert!(coarse < 0.85, "bored coarse {coarse}");
    assert!((peaks[0] - 1.0).abs() <= 0.10, "bored {peaks:?}");
    let (coarse, peaks) = refined("boss_plate", None, &[4.0, 6.0], 28.4376);
    assert!(coarse < 0.7, "boss coarse {coarse}");
    assert!((peaks[0] - 1.0).abs() <= 0.10, "boss {peaks:?}");
    assert!(
        (peaks[0] - peaks[1]).abs() <= 0.03,
        "boss box sensitivity {peaks:?}"
    );
}

#[test]
fn refine_reports_each_unresolved_raiser() {
    use popkern_analysis_cli::run::{RunOptions, run_with};
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("refine");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let setup = dir.join("bored.toml");
    let text =
        std::fs::read_to_string(repo("tools/analysis/calculix/cases/bored_block.setup.toml"))
            .unwrap_or_else(|e| panic!("{e}"));
    std::fs::write(&setup, text + "\n[options]\ntarget_cell_size = 0.06\n")
        .unwrap_or_else(|e| panic!("{e}"));
    let options = RunOptions {
        refine: true,
        ..RunOptions::default()
    };
    let out = run_with(
        &repo("tests/fixtures/analysis/gallery/bored_block.step"),
        &setup,
        &options,
    );
    report::assert_valid(&out.report);
    let subs = out.report["result"]["submodels"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    assert_eq!(subs.len(), 1, "{subs:?}");
    let peak = subs[0]["peak_von_mises"].as_f64().unwrap_or(0.0);
    let coarse = subs[0]["coarse_peak"].as_f64().unwrap_or(f64::INFINITY);
    assert!(peak > 1.15 * coarse, "{peak} vs coarse {coarse}");
    // Without the option there is nothing.
    let plain = popkern_analysis_cli::run::run(
        &repo("tests/fixtures/analysis/gallery/bored_block.step"),
        &setup,
    );
    assert!(plain.report["result"]["submodels"].is_null());
}

/// PH-064: heat flowing along bored_block past its insulated bore (the
/// thermal analogue of the pull case), at 3 cells per bore radius.
/// Measured 2026-10-11 (PH-065, over CalculiX's bore peak 6.2195e6
/// W/m², DC3D10 meshes 0.06 and 0.03 within 0.06%): the coarse solve
/// 0.951, the submodel (8 cells per radius; its box holds the whole small
/// part) 0.984, the whole part at h 0.03 0.984.
#[test]
fn refine_reports_the_bore_flux_of_a_thermal_run() {
    use popkern_analysis_cli::run::{RunOptions, run_with};
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("refine_thermal");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let setup = dir.join("bored_thermal.toml");
    let text = std::fs::read_to_string(repo(
        "tools/analysis/calculix/cases/bored_block.thermal.setup.toml",
    ))
    .unwrap_or_else(|e| panic!("{e}"));
    std::fs::write(&setup, text + "\n[options]\ntarget_cell_size = 0.06\n")
        .unwrap_or_else(|e| panic!("{e}"));
    let options = RunOptions {
        refine: true,
        ..RunOptions::default()
    };
    let step = repo("tests/fixtures/analysis/gallery/bored_block.step");
    let out = run_with(&step, &setup, &options);
    report::assert_valid(&out.report);
    let subs = out.report["result"]["submodels"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    assert_eq!(subs.len(), 1, "{subs:?}");
    let peak = subs[0]["peak_heat_flux"].as_f64().unwrap_or(0.0);
    let coarse = subs[0]["coarse_peak"].as_f64().unwrap_or(f64::INFINITY);
    let calculix = 6.2195e6;
    assert!(coarse < 0.96 * calculix, "coarse {}", coarse / calculix);
    assert!(
        (peak / calculix - 1.0).abs() <= 0.05,
        "submodel {}",
        peak / calculix
    );
    let plain = popkern_analysis_cli::run::run(&step, &setup);
    assert!(plain.report["result"]["submodels"].is_null());
}

/// PH-104: a raiser far below the coarse cell needs a box of two coarse
/// cells at its fine cell size: living_hinges' 0.4 mm raisers under
/// 2.6 mm coarse cells would take 4-5M cells each, and are declined as
/// CellBudgetExceeded rather than solved for half an hour; its 0.8 mm
/// raiser still refines (180k cells).
#[test]
fn a_raiser_far_below_the_coarse_cell_is_declined_by_budget() {
    use popkern_analysis_cli::run::run_with;
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("submodel_budget");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let setup = dir.join("heat.toml");
    let text = "analysis = \"thermal\"\nunits = \"mm\"\nmaterial = \"Aluminium 6061-T6\"\n\
        [faces]\ncold = { id = 0 }\n\
        [[fixed_temperature]]\nfaces = [\"cold\"]\ntemperature_c = 20\n\
        [[heat_source]]\npower_w = 1\n\
        [[convection]]\nfaces = \"all_other\"\nh = 10\nambient_c = 20\n\
        [options]\nmax_active_cells = 4000\n";
    std::fs::write(&setup, text).unwrap_or_else(|e| panic!("{e}"));
    let options = popkern_analysis_cli::run::RunOptions {
        refine: true,
        ..Default::default()
    };
    let out = run_with(
        &repo("tests/fixtures/analysis/gallery/living_hinges.step"),
        &setup,
        &options,
    );
    assert_eq!(out.report["status"], "solved", "{}", out.report);
    let subs = out.report["result"]["submodels"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    let declined = subs
        .iter()
        .filter(|s| s["error"]["kind"] == "CellBudgetExceeded")
        .count();
    let solved: Vec<f64> = subs
        .iter()
        .filter_map(|s| s["cells_active"].as_f64())
        .collect();
    assert_eq!((declined, solved.len()), (5, 1), "{subs:?}");
    assert!(solved[0] <= 1_000_000.0, "{solved:?}");
}
