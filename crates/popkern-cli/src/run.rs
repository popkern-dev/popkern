//! One run: script under budgets → parameter overrides → report →
//! exports. Every failure lands in the summary's outcome; nothing
//! panics, nothing is half-written (an export that fails leaves the
//! previous file untouched because it fails before writing).

use popkern_rhai::{HostOptions, ScriptHost};

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

/// The exit code for a summary.
pub fn exit_code(summary: &RunSummary) -> u8 {
    summary.outcome.exit_code()
}

fn host(cli: &Cli) -> ScriptHost {
    ScriptHost::with_options(HostOptions {
        max_operations: cli.max_ops,
        timeout: cli.timeout,
        ..HostOptions::default()
    })
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
            summary.outcome = popkern_rhai::outcome::from_run_error(&error);
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
            summary.outcome =
                popkern_rhai::outcome::from_script_error(&e, &format!("--param {key}"));
            return summary;
        }
    }
    summary.report = Some(popkern_rhai::report_value(&model));
    summary.features = model.features().len();
    if let Some(outcome) = popkern_rhai::outcome::first_failure(&model) {
        summary.outcome = outcome;
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
        // Whole-model view: every alive body (a multi-body document's
        // STEP and thumbnail carry all of it, not the active body).
        let part = match model.to_part_all() {
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
            let rendered = match cli.zebra {
                Some(normals) => crate::render::render_zebra_png(
                    part.brep(),
                    path,
                    cli.width,
                    cli.height,
                    cli.stripes,
                    normals,
                    cli.facet_angle_deg,
                )
                .map(|z| summary.zebra = Some(z)),
                None => {
                    crate::render::render_png(part.brep(), path, cli.width, cli.height, cli.edges)
                }
            };
            if let Err(e) = rendered {
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
