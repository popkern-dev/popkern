//! `popkern-analysis run <part.step> <setup.toml> [--out <dir>]
//! [--refine-check] [--refine] [--png] [--samples]`
//!
//! Writes `<dir>/<setup stem>.report.json` and, when the solve
//! succeeds, `<dir>/<setup stem>.vtu` (default dir: the setup's own);
//! `--samples` adds `<dir>/<setup stem>.samples.csv` (structural);
//! `--refine` re-solves each unresolved fillet or bore in a submodel.
//! Exit codes: 0 solved, 2 declined (the report says why), 1 usage or
//! I/O failure.

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use popkern_analysis_cli::run::{RunOptions, run_with};

const USAGE: &str = "usage: popkern-analysis run <part.step> <setup.toml> [--out <dir>] [--refine-check] [--refine] [--png] [--samples]";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut options = RunOptions::default();
    let mut out = None;
    let mut png = false;
    let mut samples = false;
    let mut positional = Vec::new();
    let mut rest = args.iter();
    while let Some(a) = rest.next() {
        match a.as_str() {
            "--out" => match rest.next() {
                Some(d) => out = Some(d),
                None => {
                    eprintln!("{USAGE}");
                    return ExitCode::from(1);
                }
            },
            "--refine-check" => options.refinement_check = true,
            "--png" => png = true,
            "--samples" => samples = true,
            "--refine" => options.refine = true,
            _ => positional.push(a),
        }
    }
    let (step, setup) = match positional.as_slice() {
        [cmd, step, setup] if *cmd == "run" => (*step, *setup),
        _ => {
            eprintln!("{USAGE}");
            return ExitCode::from(1);
        }
    };
    let setup = Path::new(setup);
    let dir: PathBuf = match out {
        Some(d) => PathBuf::from(d),
        None => setup
            .parent()
            .map_or_else(|| PathBuf::from("."), Path::to_path_buf),
    };
    let stem = setup
        .file_stem()
        .map_or_else(|| "setup".into(), |s| s.to_string_lossy().into_owned());
    if png {
        options.png = Some(dir.join(format!("{stem}.png")));
    }
    if samples {
        options.samples = Some(dir.join(format!("{stem}.samples.csv")));
    }
    // Written only by a generative run.
    options.design_mesh = Some(dir.join(format!("{stem}.design.stl")));
    let output = run_with(Path::new(step), setup, &options);
    let write = |name: String, text: &str| -> Result<(), String> {
        let path = dir.join(name);
        std::fs::write(&path, text).map_err(|e| format!("cannot write {}: {e}", path.display()))
    };
    let report = match serde_json::to_string_pretty(&output.report) {
        Ok(r) => r + "\n",
        Err(e) => {
            eprintln!("cannot serialise the report: {e}");
            return ExitCode::from(1);
        }
    };
    if let Err(e) = write(format!("{stem}.report.json"), &report) {
        eprintln!("{e}");
        return ExitCode::from(1);
    }
    if let Some(vtu) = &output.vtu
        && let Err(e) = write(format!("{stem}.vtu"), vtu)
    {
        eprintln!("{e}");
        return ExitCode::from(1);
    }
    if output.solved() {
        ExitCode::SUCCESS
    } else {
        eprintln!("declined: {}", output.report["error"]["message"]);
        ExitCode::from(2)
    }
}
