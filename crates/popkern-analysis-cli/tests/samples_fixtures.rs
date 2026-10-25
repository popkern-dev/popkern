//! `--samples` (PH-050): the cell-resolution surface samples as CSV,
//! holding the reported peak stress and largest displacement.

use std::path::PathBuf;

use popkern_analysis_cli::run::{RunOptions, run_with};

#[path = "support/report.rs"]
mod report;
use report::{assert_valid, repo};

#[test]
fn samples_cover_the_surface_and_hold_the_reported_peak() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("samples");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let step = repo("tests/fixtures/analysis/block_two_bores.step");
    let setup = repo("tests/fixtures/analysis/block_two_bores.setup.toml");
    let path = dir.join("b2b.samples.csv");
    let options = RunOptions {
        samples: Some(path.clone()),
        ..RunOptions::default()
    };
    let out = run_with(&step, &setup, &options);
    assert_valid(&out.report);
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{e}"));
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some(popkern_analysis_cli::samples::HEADER));
    let rows: Vec<Vec<f64>> = lines
        .map(|l| {
            l.split(',')
                .map(|v| v.parse::<f64>().unwrap_or(f64::NAN))
                .collect()
        })
        .collect();
    let count = out.report["result"]["samples"]["count"].as_u64();
    assert_eq!(count, Some(rows.len() as u64));
    assert!(rows.len() > 10_000, "{} samples", rows.len());
    // Every sample has a displacement, and the peak is the reported one.
    assert!(rows.iter().all(|r| r[8..11].iter().all(|v| v.is_finite())));
    // The reported largest displacement is found among the samples
    // (or the coarser display vertices).
    let largest = rows
        .iter()
        .map(|r| (r[8] * r[8] + r[9] * r[9] + r[10] * r[10]).sqrt())
        .fold(0.0, f64::max);
    let reported_u = out.report["result"]["max_displacement_mm"].as_f64();
    assert_eq!(Some(largest), reported_u);
    let peak = rows.iter().map(|r| r[6]).fold(0.0, f64::max);
    let reported = out.report["result"]["assessment"]["peaks"][0]["von_mises"]["value"].as_f64();
    assert_eq!(Some(peak), reported);
    // A thermal run refuses the option instead of ignoring it.
    let thermal = dir.join("thermal.toml");
    std::fs::write(&thermal, "analysis = \"thermal\"\nunits = \"mm\"\n")
        .unwrap_or_else(|e| panic!("{e}"));
    let refused = run_with(&step, &thermal, &options);
    assert_eq!(refused.report["error"]["kind"], "Samples");
}
