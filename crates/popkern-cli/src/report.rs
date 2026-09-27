//! What a run produced, as text or JSON (serde of the spec §10/§11
//! payloads — the same shapes scripts see).

use std::path::PathBuf;
use std::time::Duration;

use serde::Serialize;
use serde_json::Value;

use crate::args::Cli;

pub use popkern_rhai::Outcome;

/// One run's summary.
#[derive(Debug, Clone, Serialize)]
pub struct RunSummary {
    /// The script that ran.
    pub script: PathBuf,
    /// How it ended.
    #[serde(flatten)]
    pub outcome: Outcome,
    /// Wall time of the script run, when a clock exists.
    pub elapsed_ms: Option<f64>,
    /// Lines the script printed.
    pub printed: Vec<String>,
    /// The replay report (spec §11) after the run and any overrides.
    pub report: Option<Value>,
    /// Feature count of the document.
    pub features: usize,
    /// Files written.
    pub outputs: Vec<PathBuf>,
    /// Zebra render disagreement, when `--zebra` rendered the PNG.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub zebra: Option<crate::render::ZebraSummary>,
}

impl RunSummary {
    pub(crate) fn new(script: &std::path::Path) -> RunSummary {
        RunSummary {
            script: script.to_path_buf(),
            outcome: Outcome::Served,
            elapsed_ms: None,
            printed: Vec::new(),
            report: None,
            features: 0,
            outputs: Vec::new(),
            zebra: None,
        }
    }

    pub(crate) fn with_elapsed(mut self, elapsed: Option<Duration>) -> RunSummary {
        self.elapsed_ms = elapsed.map(|d| d.as_secs_f64() * 1000.0);
        self
    }
}

fn feature_lines(report: &Value) -> Vec<String> {
    let mut lines = Vec::new();
    if let Some(features) = report.get("features").and_then(Value::as_array) {
        for f in features {
            let label = f
                .pointer("/feature/label")
                .or_else(|| f.get("label"))
                .and_then(Value::as_str)
                .unwrap_or("?");
            let operation = f
                .pointer("/feature/operation")
                .or_else(|| f.get("operation"))
                .and_then(Value::as_str)
                .unwrap_or("?");
            let outcome = match f.get("outcome") {
                Some(Value::String(s)) => s.clone(),
                Some(Value::Object(map)) => map.keys().next().cloned().unwrap_or_default(),
                _ => "?".into(),
            };
            let elapsed = f
                .get("elapsed_ms")
                .and_then(Value::as_f64)
                .map(|ms| format!(" {ms:.1} ms"))
                .unwrap_or_default();
            let downgraded = f
                .pointer("/blend/downgraded")
                .and_then(Value::as_bool)
                .unwrap_or(false);
            let mut line = format!("  {label:<18} {operation:<10} {outcome}{elapsed}");
            if downgraded {
                line.push_str("  (continuity downgraded)");
            }
            lines.push(line);
        }
    }
    lines
}

/// Render the summary as text.
pub fn text(summary: &RunSummary) -> String {
    let mut out = String::new();
    match &summary.outcome {
        Outcome::Served => out.push_str(&format!("served {}", summary.script.display())),
        Outcome::ScriptError { rendered, .. } => out.push_str(&format!("script error: {rendered}")),
        Outcome::KernelDecline { rendered, .. } => {
            out.push_str(&format!("kernel decline: {rendered}"))
        }
        Outcome::GateFailure { step, detail } => out.push_str(&format!("{step} failed: {detail}")),
    }
    if let Some(ms) = summary.elapsed_ms {
        out.push_str(&format!(" ({ms:.0} ms)"));
    }
    out.push('\n');
    for line in &summary.printed {
        out.push_str(&format!("> {line}\n"));
    }
    if let Some(report) = &summary.report {
        out.push_str(&format!("{} feature(s):\n", summary.features));
        for line in feature_lines(report) {
            out.push_str(&line);
            out.push('\n');
        }
    }
    for path in &summary.outputs {
        out.push_str(&format!("wrote {}\n", path.display()));
    }
    if let Some(z) = &summary.zebra {
        out.push_str(&format!(
            "zebra ({} normals, {} stripes, {}° facets): mesh vs surface normal max {:.3}° mean {:.4}° p99 {:.3}° over {} px ({} ray hits); worst {} at {:?}; {} px change band\n",
            z.normals,
            z.stripes,
            z.facet_angle_deg,
            z.max_angle_deg,
            z.mean_angle_deg,
            z.p99_angle_deg,
            z.evaluated_pixels,
            z.ray_hit_pixels,
            z.worst_surface.unwrap_or("?"),
            z.worst_pixel.unwrap_or([0, 0]),
            z.band_mismatch_pixels,
        ));
    }
    out
}

/// Print the summary as the CLI asked (text to stdout, or JSON).
pub fn print(cli: &Cli, summary: &RunSummary) {
    if cli.json {
        match serde_json::to_string_pretty(summary) {
            Ok(json) => println!("{json}"),
            Err(e) => eprintln!("popkern: could not serialise the report: {e}"),
        }
    } else {
        print!("{}", text(summary));
    }
}
