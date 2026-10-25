//! Stage 2 runs end to end: radiation (PH-070), transients (PH-071),
//! thermoelastic loads (PH-072) and modal analysis (PH-073).

use std::path::PathBuf;

use popkern_analysis_cli::run::{RunOptions, run, run_with};

#[path = "support/report.rs"]
mod report;
use report::{assert_valid, repo};

/// PH-070: convection with an emittance radiates, solved to its fixed
/// point; the report records the passes and the emittance used.
#[test]
fn a_radiating_face_is_iterated_and_reported() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("radiating");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let setup = dir.join("radiating.toml");
    let text = std::fs::read_to_string(repo(
        "tools/analysis/calculix/cases/bored_block.thermal.setup.toml",
    ))
    .unwrap_or_else(|e| panic!("{e}"));
    let text = text
        + "\n[[convection]]\nfaces = \"all_other\"\nh = 5\nambient_c = 20\nemittance = \"material\"\n";
    std::fs::write(&setup, text).unwrap_or_else(|e| panic!("{e}"));
    let out = run(
        &repo("tests/fixtures/analysis/gallery/bored_block.step"),
        &setup,
    );
    assert_valid(&out.report);
    assert_eq!(out.report["status"], "solved", "{}", out.report);
    let passes = out.report["solve"]["radiation_passes"]
        .as_u64()
        .unwrap_or(0);
    assert!((2..20).contains(&passes), "{passes}");
    assert_eq!(
        out.report["result"]["convection"][0]["emittance"],
        "material"
    );
}

/// PH-071: a transient run reports every step, its final state as the
/// fields, and declines where it cannot be honest.
#[test]
fn a_transient_run_reports_its_steps() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("transient");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let base = std::fs::read_to_string(repo(
        "tools/analysis/calculix/cases/bored_block.thermal.setup.toml",
    ))
    .unwrap_or_else(|e| panic!("{e}"));
    let step = repo("tests/fixtures/analysis/gallery/bored_block.step");
    let run_text = |name: &str, text: String, refine: bool| {
        let setup = dir.join(name);
        std::fs::write(&setup, text).unwrap_or_else(|e| panic!("{e}"));
        let options = RunOptions {
            refine,
            ..RunOptions::default()
        };
        run_with(&step, &setup, &options)
    };
    let transient = "\n[transient]\ninitial_c = 20\nstep_s = 0.01\nduration_s = 0.1\n";
    let out = run_text("al.toml", base.clone() + transient, false);
    assert_valid(&out.report);
    assert_eq!(out.report["status"], "solved", "{}", out.report);
    let steps = out.report["result"]["transient"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    assert_eq!(steps.len(), 10);
    let last = &steps[9];
    assert!((last["time_s"].as_f64().unwrap_or(0.0) - 0.1).abs() < 1e-12);
    assert!(last["imbalance"].as_f64().unwrap_or(1.0) <= 1e-6, "{last}");
    // 1 W in for 0.1 s, minus what left through the cold end.
    let stored = last["stored_j"].as_f64().unwrap_or(0.0);
    assert!(stored > 0.0 && stored < 0.1, "{stored}");
    // ABS's conductivity is cited but its specific heat is not:
    // declined, not guessed.
    let abs = base.replace("Aluminium 6061-T6", "ABS (FDM)") + transient;
    let out = run_text("abs.toml", abs, false);
    assert_eq!(out.report["status"], "declined", "{}", out.report);
    let message = out.report["error"]["message"].as_str().unwrap_or("");
    assert!(message.contains("specific heat"), "{message}");
    let out = run_text("refine.toml", base + transient, true);
    assert_eq!(out.report["error"]["kind"], "Unsupported", "{}", out.report);
}

/// PH-072: a thermal strain on a structural run, uniform or from a
/// thermal setup solved first; a material without a cited expansion
/// declines.
#[test]
fn a_thermoelastic_run_loads_the_part_with_its_temperature() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("thermoelastic");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let base =
        std::fs::read_to_string(repo("tools/analysis/calculix/cases/bored_block.setup.toml"))
            .unwrap_or_else(|e| panic!("{e}"));
    let al = base.replace("PLA (FDM)", "Aluminium 6061-T6");
    let step = repo("tests/fixtures/analysis/gallery/bored_block.step");
    let run_text = |name: &str, text: String| {
        let setup = dir.join(name);
        std::fs::write(&setup, text).unwrap_or_else(|e| panic!("{e}"));
        run(&step, &setup)
    };
    // Uniform: the clamped end restrains the expansion.
    let uniform = al.clone() + "\n[thermal]\nreference_c = 20\nuniform_c = 80\n";
    let out = run_text("uniform.toml", uniform);
    assert_valid(&out.report);
    assert_eq!(out.report["status"], "solved", "{}", out.report);
    assert_eq!(out.report["result"]["thermal"]["alpha_per_k"], 23.0e-6);
    // From a thermal setup beside it.
    let thermal = std::fs::read_to_string(repo(
        "tools/analysis/calculix/cases/bored_block.thermal.setup.toml",
    ))
    .unwrap_or_else(|e| panic!("{e}"));
    std::fs::write(dir.join("heat.toml"), thermal).unwrap_or_else(|e| panic!("{e}"));
    let coupled = al + "\n[thermal]\nreference_c = 20\nsetup = \"heat.toml\"\n";
    let out = run_text("coupled.toml", coupled);
    assert_valid(&out.report);
    assert_eq!(out.report["status"], "solved", "{}", out.report);
    assert!(
        out.report["thermal_solve"]["max_temperature_c"]
            .as_f64()
            .unwrap_or(0.0)
            > 20.0
    );
    // PETG has no cited expansion (PLA's is cited since PH-009).
    let petg =
        base.replace("PLA (FDM)", "PETG (FDM)") + "\n[thermal]\nreference_c = 20\nuniform_c = 80\n";
    let out = run_text("petg.toml", petg);
    assert_eq!(out.report["status"], "declined", "{}", out.report);
    let message = out.report["error"]["message"].as_str().unwrap_or("");
    assert!(message.contains("thermal expansion"), "{message}");
}

/// PH-073: `[modal]` reports the lowest natural frequencies; without
/// supports the part is free and its first six modes are rigid.
#[test]
fn a_modal_run_reports_frequencies() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("modal");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let base =
        std::fs::read_to_string(repo("tools/analysis/calculix/cases/bored_block.setup.toml"))
            .unwrap_or_else(|e| panic!("{e}"));
    let step = repo("tests/fixtures/analysis/gallery/bored_block.step");
    let setup = dir.join("clamped.toml");
    // Coarse cells: the gate checks the path, not the frequencies.
    let coarse = "\n[options]\ntarget_cell_size = 0.1\n";
    std::fs::write(&setup, base + "\n[modal]\nmodes = 4\n" + coarse)
        .unwrap_or_else(|e| panic!("{e}"));
    let out = run(&step, &setup);
    assert_valid(&out.report);
    assert_eq!(out.report["status"], "solved", "{}", out.report);
    let modes = out.report["result"]["modes"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    assert_eq!(modes.len(), 4);
    let f: Vec<f64> = modes
        .iter()
        .map(|m| m["frequency_hz"].as_f64().unwrap_or(0.0))
        .collect();
    assert!(f[0] > 0.0 && f.windows(2).all(|w| w[0] <= w[1]), "{f:?}");
    assert!(out.report["result"]["mass_kg"].as_f64().unwrap_or(0.0) > 0.0);
    // Free: no supports.
    let free =
        "units = \"mm\"\nmaterial = \"PLA (FDM)\"\n[modal]\nmodes = 7\n".to_string() + coarse;
    let setup = dir.join("free.toml");
    std::fs::write(&setup, free).unwrap_or_else(|e| panic!("{e}"));
    let out = run(&step, &setup);
    assert_eq!(out.report["status"], "solved", "{}", out.report);
    let modes = out.report["result"]["modes"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    let w: Vec<f64> = modes
        .iter()
        .map(|m| m["omega_squared"].as_f64().unwrap_or(f64::NAN))
        .collect();
    assert!(w[..6].iter().all(|x| x.abs() <= 1e-6 * w[6]), "{w:?}");
}

/// PH-076: a gyroid infill softens the part and reports its core; other
/// patterns and metals decline.
#[test]
fn an_infilled_run_softens_the_part_and_declines_other_patterns() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("infill");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let base = std::fs::read_to_string(repo("tests/fixtures/analysis/block_two_bores.setup.toml"))
        .unwrap_or_else(|e| panic!("{e}"));
    let base = base + "\n[options]\ntarget_cell_size = 0.5\n";
    let step = repo("tests/fixtures/analysis/block_two_bores.step");
    let run_text = |name: &str, text: String| {
        let setup = dir.join(name);
        std::fs::write(&setup, text).unwrap_or_else(|e| panic!("{e}"));
        run(&step, &setup)
    };
    let solid = run_text("solid.toml", base.clone());
    let infill = "\n[infill]\npattern = \"gyroid\"\ndensity = 0.25\nwalls_mm = 1.0\n";
    let out = run_text("gyroid.toml", base.clone() + infill);
    assert_valid(&out.report);
    assert_eq!(out.report["status"], "solved", "{}", out.report);
    let soft = out.report["result"]["max_displacement_mm"]
        .as_f64()
        .unwrap_or(0.0);
    let hard = solid.report["result"]["max_displacement_mm"]
        .as_f64()
        .unwrap_or(f64::INFINITY);
    assert!(soft > 1.2 * hard, "{soft} vs solid {hard}");
    assert_eq!(out.report["result"]["infill"]["pattern"], "gyroid");
    // PH-101: the core's volume, none without infill.
    let core = out.report["solve"]["infill_core_volume"]
        .as_f64()
        .unwrap_or(0.0);
    assert!(core > 0.0, "{}", out.report["solve"]);
    assert_eq!(solid.report["solve"]["infill_core_volume"], 0.0);
    // PH-079: the core is checked against Bean 2025's gyroid yield.
    let safety = &out.report["result"]["assessment"]["safety"];
    assert!(
        !safety["unchecked"].to_string().contains("InfillCore"),
        "{safety}"
    );
    let checks = safety["checks"].as_array().cloned().unwrap_or_default();
    let core = checks.iter().find(|c| c["check"] == "InfillCore");
    assert!(
        core.is_some_and(|c| c["factor"].as_f64().is_some_and(|f| f > 0.0)),
        "{checks:?}"
    );
    assert_eq!(out.report["solve"]["infill_density"], 0.25);
    let grid = run_text(
        "grid.toml",
        base.clone() + &infill.replace("gyroid", "grid"),
    );
    assert_eq!(grid.report["status"], "declined");
    let message = grid.report["error"]["message"].as_str().unwrap_or("");
    assert!(message.contains("not available for `grid`"), "{message}");
    let metal = run_text(
        "metal.toml",
        base.replace("PLA (FDM)", "Aluminium 6061-T6") + infill,
    );
    assert_eq!(metal.report["status"], "declined");
}
