//! One run: script under budgets → parameter overrides → report →
//! exports. Every failure lands in the summary's outcome; nothing
//! panics, nothing is half-written (an export that fails leaves the
//! previous file untouched because it fails before writing).

use popkern_rhai::{HostOptions, ScriptHost};
use popkern_script::Model;

use crate::args::Cli;
use crate::report::{Outcome, RunSummary};

/// Why the CLI could not even start a run.
#[derive(Debug, thiserror::Error)]
pub enum CliError {
    /// The script file could not be read.
    #[error("cannot read '{path}': {detail}")]
    Read {
        /// The path.
        path: String,
        /// The error.
        detail: String,
    },
}

/// The §10 classes that mean the kernel declined a feature (exit 3)
/// rather than the script being wrong (exit 2).
const KERNEL_CLASSES: &[&str] = &[
    "FeatureFailed",
    "BlendDeclined",
    "BooleanFailed",
    "SweepFailed",
    "ShellFailed",
    "DraftFailed",
    "ResolutionFailed",
    "CoreError",
];

/// The §10 classes that are gates or I/O (exit 4).
const GATE_CLASSES: &[&str] = &["NotWatertight", "ValidationFailed", "IoFailed"];

fn classify(kind: &str, rendered: String, payload: serde_json::Value) -> Outcome {
    if KERNEL_CLASSES.contains(&kind) {
        Outcome::KernelDecline {
            kind: kind.to_string(),
            rendered,
            payload,
        }
    } else if GATE_CLASSES.contains(&kind) {
        Outcome::GateFailure {
            step: kind.to_string(),
            detail: rendered,
        }
    } else {
        Outcome::ScriptError {
            kind: kind.to_string(),
            rendered,
            payload,
        }
    }
}

/// The exit code for a summary.
pub fn exit_code(summary: &RunSummary) -> u8 {
    match summary.outcome {
        Outcome::Served => 0,
        Outcome::ScriptError { .. } => 2,
        Outcome::KernelDecline { .. } => 3,
        Outcome::GateFailure { .. } => 4,
    }
}

fn host(cli: &Cli) -> ScriptHost {
    ScriptHost::with_options(HostOptions {
        max_operations: cli.max_ops,
        timeout: cli.timeout,
    })
}

fn report_value(model: &Model) -> serde_json::Value {
    let report = model.report();
    let mut value = serde_json::to_value(report).unwrap_or(serde_json::Value::Null);
    if let serde_json::Value::Object(map) = &mut value {
        map.insert("is_clean".into(), report.is_clean().into());
        map.insert("failures".into(), report.failures().count().into());
        map.insert(
            "downgraded_blends".into(),
            report.downgraded_blends().count().into(),
        );
    }
    value
}

/// Run the script and produce every requested output.
pub fn execute(cli: &Cli) -> RunSummary {
    let mut summary = RunSummary::new(&cli.script);
    let source = match std::fs::read_to_string(&cli.script) {
        Ok(text) => text,
        Err(e) => {
            summary.outcome = Outcome::GateFailure {
                step: "read".into(),
                detail: format!("cannot read '{}': {e}", cli.script.display()),
            };
            return summary;
        }
    };
    let name = cli
        .script
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "script.rhai".into());
    let mut host = host(cli);
    let outcome = match host.run(&source, &name) {
        Ok(outcome) => outcome,
        Err(error) => {
            summary.outcome = classify(&error.kind, error.render(), *error.payload.clone());
            return summary;
        }
    };
    summary = summary.with_elapsed(outcome.elapsed);
    summary.printed = outcome.printed;
    let Some(mut model) = outcome.model else {
        summary.outcome = Outcome::GateFailure {
            step: "model".into(),
            detail: "the script created no model (call model() and keep it in a variable)".into(),
        };
        return summary;
    };
    for (key, value) in &cli.params {
        if let Err(e) = model
            .set_param(key, *value)
            .and_then(|_| model.replay().map(|_| ()))
        {
            let payload = e.payload();
            let kind = payload["kind"]
                .as_str()
                .unwrap_or("ArgumentError")
                .to_string();
            summary.outcome = classify(&kind, format!("--param {key}: {e}"), payload);
            return summary;
        }
    }
    summary.report = Some(report_value(&model));
    summary.features = model.features().len();
    if let Some(row) = model.report().failures().next() {
        summary.outcome = Outcome::KernelDecline {
            kind: "FeatureFailed".into(),
            rendered: format!("feature '{}' ({}) failed", row.label, row.operation),
            payload: serde_json::to_value(row).unwrap_or(serde_json::Value::Null),
        };
        return summary;
    }
    if let Some(path) = &cli.pkn {
        if let Err(e) = model.save(path) {
            summary.outcome = Outcome::GateFailure {
                step: "pkn".into(),
                detail: e.to_string(),
            };
            return summary;
        }
        summary.outputs.push(path.clone());
    }
    if cli.step.is_some() || cli.png.is_some() {
        let part = match model.to_part() {
            Ok(part) => part,
            Err(e) => {
                summary.outcome = Outcome::GateFailure {
                    step: "validate".into(),
                    detail: e.to_string(),
                };
                return summary;
            }
        };
        if let Some(path) = &cli.step {
            if let Err(e) = part.export_step(path) {
                summary.outcome = Outcome::GateFailure {
                    step: "step".into(),
                    detail: e.to_string(),
                };
                return summary;
            }
            summary.outputs.push(path.clone());
        }
        if let Some(path) = &cli.png {
            if let Err(e) =
                crate::render::render_png(part.brep(), path, cli.width, cli.height, cli.edges)
            {
                summary.outcome = Outcome::GateFailure {
                    step: "png".into(),
                    detail: e,
                };
                return summary;
            }
            summary.outputs.push(path.clone());
        }
    }
    summary
}
