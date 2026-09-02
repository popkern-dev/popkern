//! How a hosted run ended — one class partition shared by every host
//! (the `popkern` command, the browser playground, the daemon): the
//! §10 error classes are sorted into *the script was wrong*, *the
//! kernel declined a feature*, and *a gate or export failed after the
//! script ran*, and each partition owns one process exit code.

use serde::Serialize;
use serde_json::Value;

use popkern_script::Model;

/// How a run ended.
#[derive(Debug, Clone, Serialize)]
#[serde(tag = "outcome", rename_all = "snake_case")]
pub enum Outcome {
    /// The script ran and every requested export succeeded.
    Served,
    /// The script itself failed (syntax, argument, budget).
    ScriptError {
        /// The §10 class (`SyntaxError`, `ArgumentError`, `BudgetExceeded`, …).
        kind: String,
        /// The rendered error (class, position, message, remediation).
        rendered: String,
        /// The §10 payload.
        payload: Value,
    },
    /// A feature declined in the kernel.
    KernelDecline {
        /// The §10 class (`FeatureFailed`, `BlendDeclined`, …).
        kind: String,
        /// The rendered error.
        rendered: String,
        /// The §10 payload.
        payload: Value,
    },
    /// A gate or export failed after the script ran.
    GateFailure {
        /// Which step (`assert_watertight`, `step`, `png`, `pkn`, `param`, `read`, `model`).
        step: String,
        /// The error text.
        detail: String,
    },
}

/// The §10 classes that mean the kernel declined a feature rather
/// than the script being wrong.
pub const KERNEL_CLASSES: &[&str] = &[
    "FeatureFailed",
    "BlendDeclined",
    "BooleanFailed",
    "SweepFailed",
    "ShellFailed",
    "DraftFailed",
    "ResolutionFailed",
    "CoreError",
];

/// The §10 classes that are gates or I/O.
pub const GATE_CLASSES: &[&str] = &["NotWatertight", "ValidationFailed", "IoFailed"];

/// Sort a §10 class into its partition.
pub fn classify(kind: &str, rendered: String, payload: Value) -> Outcome {
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

/// The outcome of a failed host run.
pub fn from_run_error(error: &crate::ScriptRunError) -> Outcome {
    classify(&error.kind, error.render(), *error.payload.clone())
}

/// The outcome of a recorder error raised outside a script (a
/// parameter override, an export): classified by its §10 payload.
pub fn from_script_error(error: &popkern_script::ScriptError, context: &str) -> Outcome {
    let payload = error.payload();
    let kind = payload["kind"]
        .as_str()
        .unwrap_or("ArgumentError")
        .to_string();
    classify(&kind, format!("{context}: {error}"), payload)
}

impl Outcome {
    /// The process exit code of this partition: `0` served, `2` the
    /// script itself, `3` a kernel decline, `4` a gate or export.
    pub fn exit_code(&self) -> u8 {
        match self {
            Outcome::Served => 0,
            Outcome::ScriptError { .. } => 2,
            Outcome::KernelDecline { .. } => 3,
            Outcome::GateFailure { .. } => 4,
        }
    }

    /// Whether the run served.
    pub fn is_served(&self) -> bool {
        matches!(self, Outcome::Served)
    }

    /// The §10 class, when the outcome carries one.
    pub fn kind(&self) -> Option<&str> {
        match self {
            Outcome::Served => None,
            Outcome::ScriptError { kind, .. } | Outcome::KernelDecline { kind, .. } => Some(kind),
            Outcome::GateFailure { step, .. } => Some(step),
        }
    }
}

/// The replay report (spec §11) as JSON with the derived counts every
/// host shows beside it: `is_clean`, `failures`, `downgraded_blends`.
pub fn report_value(model: &Model) -> Value {
    let report = model.report();
    let mut value = serde_json::to_value(report).unwrap_or(Value::Null);
    if let Value::Object(map) = &mut value {
        map.insert("is_clean".into(), report.is_clean().into());
        map.insert("failures".into(), report.failures().count().into());
        map.insert(
            "downgraded_blends".into(),
            report.downgraded_blends().count().into(),
        );
        // Per-row derived facts a panel shows directly: wall time in
        // milliseconds and the blend's downgrade verdict (measured
        // continuity below the request), which the record only
        // exposes as a method.
        if let Some(Value::Array(rows)) = map.get_mut("features") {
            for (row, feature) in rows.iter_mut().zip(&report.features) {
                let Value::Object(row) = row else { continue };
                row.insert(
                    "elapsed_ms".into(),
                    feature
                        .elapsed
                        .map(|d| Value::from(d.as_secs_f64() * 1000.0))
                        .unwrap_or(Value::Null),
                );
                if let (Some(Value::Object(blend)), Some(evidence)) =
                    (row.get_mut("blend"), feature.blend.as_ref())
                {
                    blend.insert("downgraded".into(), evidence.downgraded().into());
                }
            }
        }
    }
    value
}

/// The kernel-decline outcome for the first failed feature of the last
/// replay, if any (a feature that failed inside a replay after a
/// parameter override is a decline even though no script raised).
pub fn first_failure(model: &Model) -> Option<Outcome> {
    model
        .report()
        .failures()
        .next()
        .map(|row| Outcome::KernelDecline {
            kind: "FeatureFailed".into(),
            rendered: format!("feature '{}' ({}) failed", row.label, row.operation),
            payload: serde_json::to_value(row).unwrap_or(Value::Null),
        })
}
