//! The `run` command in the gallery gate: the arguments a script's
//! exports take, the OCC oracle on a STEP, and one script through the
//! whole command (each script's test exports the document its session
//! built instead of building it a second time, PK-1223).

use std::path::Path;

use popkern_cli::{Cli, Outcome, execute, parse_args};

/// `popkern run <script> --step <step> --png <png>` at the gallery's
/// thumbnail size.
pub(crate) fn cli(script: &Path, step: &Path, png: &Path) -> Cli {
    let args: Vec<String> = [
        "run",
        script.to_str().expect("utf8"),
        "--step",
        step.to_str().expect("utf8"),
        "--png",
        png.to_str().expect("utf8"),
        "--width",
        "480",
        "--height",
        "360",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    parse_args(&args).expect("args")
}

/// The OCC oracle's verdict on `step`; `None` when the oracle is absent.
pub(crate) fn oracle_passes(step: &Path) -> Option<bool> {
    let oracle = Path::new(&std::env::var("HOME").unwrap_or_default()).join("venvs/ocp/bin/python");
    let check = super::root().join("tools/occ-oracle/occ_check.py");
    if !(oracle.exists() && check.exists()) {
        return None;
    }
    let output = std::process::Command::new(&oracle)
        .arg(&check)
        .arg(step)
        .output()
        .expect("oracle");
    Some(String::from_utf8_lossy(&output.stdout).contains("\"passed\": true"))
}

/// The command end to end on a gallery script: it builds the script
/// itself, writes a STEP the oracle accepts and a thumbnail.
#[test]
fn the_run_command_writes_step_and_thumbnail() {
    let script = super::examples_dir().join("boss_plate.rhai");
    let step = super::tmp("run_command_boss_plate.step");
    let png = super::tmp("run_command_boss_plate.png");
    for file in [&step, &png] {
        let _ = std::fs::remove_file(file);
    }
    let summary = execute(&cli(&script, &step, &png));
    assert!(
        matches!(summary.outcome, Outcome::Served),
        "{:?}",
        summary.outcome
    );
    assert_eq!(summary.outputs, vec![step.clone(), png.clone()]);
    for file in [&step, &png] {
        let size = std::fs::metadata(file).map(|m| m.len()).unwrap_or(0);
        assert!(size > 0, "{} written", file.display());
    }
    if let Some(passed) = oracle_passes(&step) {
        assert!(passed, "OCC oracle rejected the run command's STEP");
    }
}
