//! The spec-§10 error taxonomy as script values. Every recorder error
//! crosses into Rhai as a runtime error whose value is a map with
//! `kind`, `message`, and the class's payload fields; the engine
//! stamps the script position of the call that raised it, and
//! `try`/`catch` receives the same map.

use popkern_history::{HistoryError, ReplayError};
use popkern_part::PartError;
use popkern_script::ScriptError;
use rhai::{Dynamic, EvalAltResult, Map, Position};
use serde_json::{Value, json};

/// Build a script runtime error from a §10 payload.
pub fn raise(payload: Value) -> Box<EvalAltResult> {
    let value =
        rhai::serde::to_dynamic(&payload).unwrap_or_else(|_| Dynamic::from(payload.to_string()));
    Box::new(EvalAltResult::ErrorRuntime(value, Position::NONE))
}

/// A skin-level `ArgumentError` (raised before the kernel runs).
pub fn argument_error(argument: &str, reason: impl Into<String>) -> Box<EvalAltResult> {
    let reason = reason.into();
    raise(json!({
        "kind": "ArgumentError",
        "message": format!("argument '{argument}': {reason}"),
        "argument": argument,
        "reason": reason,
    }))
}

/// The §10 payload of a recorder error.
pub fn payload(error: &ScriptError) -> Value {
    let message = error.to_string();
    let mut value = match error {
        ScriptError::FeatureFailed {
            label,
            operation,
            error,
        } => json!({
            "kind": "FeatureFailed",
            "feature": { "label": label, "operation": operation },
            "error": serde_json::to_value(error.as_ref()).unwrap_or(Value::Null),
        }),
        ScriptError::EmptySelection {
            selection,
            kind,
            available,
        } => json!({
            "kind": "EmptySelection",
            "selection": selection,
            "entity": kind.to_string(),
            "available": available,
        }),
        ScriptError::SelectionCardinality {
            selection,
            kind,
            expected,
            got,
        } => json!({
            "kind": "SelectionCardinality",
            "selection": selection,
            "entity": kind.to_string(),
            "expected": expected,
            "got": got,
        }),
        ScriptError::Argument { argument, reason } => json!({
            "kind": "ArgumentError",
            "argument": argument,
            "reason": reason,
        }),
        ScriptError::NoActiveBody => json!({
            "kind": "ArgumentError",
            "argument": "body",
            "reason": "no active body: record a sketch and a solid feature first",
        }),
        ScriptError::BodyNotAlive(id) => json!({
            "kind": "TimelineError",
            "reason": "body_not_found",
            "body": id,
        }),
        ScriptError::NotWatertight { open_edges } => json!({
            "kind": "NotWatertight",
            "open_edges": open_edges,
        }),
        ScriptError::Validation { detail } => json!({
            "kind": "ValidationFailed",
            "detail": detail,
            "issues": [],
        }),
        ScriptError::Io { path, detail } => json!({
            "kind": "IoFailed",
            "path": path,
            "detail": detail,
        }),
        ScriptError::History(e) => history_payload(e),
        ScriptError::Replay(e) => replay_payload(e),
        ScriptError::Part(e) => part_payload(e),
    };
    if let Value::Object(map) = &mut value {
        map.insert("message".into(), Value::String(message));
    }
    value
}

fn history_payload(error: &HistoryError) -> Value {
    match error {
        HistoryError::DuplicateLabel(label) => {
            json!({ "kind": "LabelError", "label": label, "reason": "duplicate" })
        }
        HistoryError::EmptyLabel => json!({ "kind": "LabelError", "label": "", "reason": "empty" }),
        HistoryError::LabelNotFound(label) => {
            json!({ "kind": "LabelError", "label": label, "reason": "not_found" })
        }
        HistoryError::ExpressionCycle(ids) => json!({
            "kind": "ParameterError", "reason": "cycle", "parameters": ids.len(),
        }),
        HistoryError::ExpressionUnresolved(_) => {
            json!({ "kind": "ParameterError", "reason": "unknown_reference" })
        }
        HistoryError::ParameterInUse(ids) => json!({
            "kind": "ParameterError", "reason": "referenced", "referenced_by": ids.len(),
        }),
        HistoryError::ReadOnlyParameter => {
            json!({ "kind": "ParameterError", "reason": "read_only" })
        }
        HistoryError::DuplicateParameterName(name) => {
            json!({ "kind": "ParameterError", "reason": "duplicate", "parameter": name })
        }
        HistoryError::DimensionMismatch {
            param,
            expected,
            got,
        } => json!({
            "kind": "ParameterError", "reason": "dimension_mismatch",
            "parameter": format!("{param:?}"),
            "expected": format!("{expected:?}"), "got": format!("{got:?}"),
        }),
        HistoryError::SuppressIfNotDimensionless(kind) => json!({
            "kind": "ParameterError", "reason": "dimension_mismatch",
            "expected": "Dimensionless", "got": format!("{kind:?}"),
        }),
        HistoryError::FeatureNotFound => {
            json!({ "kind": "TimelineError", "reason": "feature_not_found" })
        }
        HistoryError::InvalidPosition(index) => {
            json!({ "kind": "TimelineError", "reason": "invalid_position", "index": index })
        }
        HistoryError::DependencyCycle => {
            json!({ "kind": "TimelineError", "reason": "dependency_cycle" })
        }
        HistoryError::DependencyAfterInsert(what) => json!({
            "kind": "TimelineError", "reason": "dependency_after_insertion", "dependency": what,
        }),
        HistoryError::BodyNotFound => {
            json!({ "kind": "TimelineError", "reason": "body_not_found" })
        }
        HistoryError::EmptySeedFeatures => {
            json!({ "kind": "TimelineError", "reason": "seed_empty" })
        }
        HistoryError::SeedBodyMismatch => {
            json!({ "kind": "TimelineError", "reason": "seed_body_mismatch" })
        }
        HistoryError::SeedNotContiguous => {
            json!({ "kind": "TimelineError", "reason": "seed_not_contiguous" })
        }
        HistoryError::SeedMissingSketch => {
            json!({ "kind": "TimelineError", "reason": "seed_missing_sketch" })
        }
        HistoryError::SeedMissingMode => {
            json!({ "kind": "TimelineError", "reason": "seed_missing_mode" })
        }
        HistoryError::SeedAlreadyOwned => {
            json!({ "kind": "TimelineError", "reason": "seed_already_owned" })
        }
        HistoryError::ConfigNotFound => {
            json!({ "kind": "TimelineError", "reason": "config_not_found" })
        }
        HistoryError::InheritanceCycle => {
            json!({ "kind": "TimelineError", "reason": "config_inheritance_cycle" })
        }
        HistoryError::DuplicateConfigName(name) => json!({
            "kind": "LabelError", "label": name, "reason": "duplicate",
        }),
        HistoryError::SerializeError(detail) => json!({
            "kind": "DocumentError", "reason": "serialization", "detail": detail,
        }),
        HistoryError::UnsupportedVersion {
            file_version,
            min_supported,
            max_supported,
        } => json!({
            "kind": "DocumentError", "reason": "unsupported_version",
            "file_version": file_version, "min": min_supported, "max": max_supported,
        }),
        HistoryError::UnsupportedNamingScheme {
            file_scheme,
            supported,
        } => json!({
            "kind": "DocumentError", "reason": "unsupported_naming_scheme",
            "file_scheme": file_scheme, "supported": supported,
        }),
    }
}

fn replay_payload(error: &ReplayError) -> Value {
    match error {
        ReplayError::FeatureFailed { label, error } => json!({
            "kind": "FeatureFailed",
            "feature": { "label": label, "operation": "" },
            "error": serde_json::to_value(error.as_ref()).unwrap_or(Value::Null),
        }),
        ReplayError::ResolutionFailed(name) => json!({
            "kind": "ResolutionFailed",
            "name": format!("{name:?}"),
        }),
        ReplayError::CoreError(e) => json!({ "kind": "CoreError", "detail": e.to_string() }),
        ReplayError::HistoryError(e) => history_payload(e),
    }
}

fn part_payload(error: &PartError) -> Value {
    match error {
        PartError::Blend { selection, source } => {
            let blend = popkern_history::BlendDeclined::project(source, None);
            let mut value = serde_json::to_value(&blend).unwrap_or(Value::Null);
            if let Value::Object(map) = &mut value {
                map.insert("kind".into(), "BlendDeclined".into());
                map.insert("selection".into(), Value::String(selection.clone()));
            }
            value
        }
        PartError::Boolean { operation, source } => json!({
            "kind": "BooleanFailed", "operation": operation, "detail": source.to_string(),
        }),
        PartError::Sweep { feature, source } => json!({
            "kind": "SweepFailed", "feature": feature, "detail": source.to_string(),
        }),
        PartError::SmoothMesh { stage, detail } => json!({
            "kind": "NotApplicable", "detail": format!("smooth mesh failed at {stage}: {detail}"),
        }),
        PartError::EmptySelection {
            selection,
            entity,
            available,
        } => json!({
            "kind": "EmptySelection", "selection": selection, "entity": entity, "available": available,
        }),
        PartError::Validation { feature, source } => json!({
            "kind": "ValidationFailed", "feature": feature, "detail": source.to_string(), "issues": [],
        }),
        PartError::Core { feature, source } => json!({
            "kind": "CoreError", "feature": feature, "detail": source.to_string(),
        }),
        PartError::Mesh { source } => json!({ "kind": "CoreError", "detail": source.to_string() }),
        PartError::NotWatertight { open_edges } => {
            json!({ "kind": "NotWatertight", "open_edges": open_edges })
        }
        PartError::Io { path, source } => {
            json!({ "kind": "IoFailed", "path": path, "detail": source.to_string() })
        }
        PartError::NotApplicable(detail) => json!({ "kind": "NotApplicable", "detail": detail }),
        PartError::Cancelled => json!({ "kind": "Cancelled" }),
    }
}

/// Convert a recorder error into the script runtime error.
pub fn script_error(error: ScriptError) -> Box<EvalAltResult> {
    raise(payload(&error))
}

/// `Result` adaptor: `recorder_call().map_err(script_error)?`.
pub trait IntoScriptResult<T> {
    /// Map the recorder error into a script runtime error.
    fn script(self) -> Result<T, Box<EvalAltResult>>;
}

impl<T> IntoScriptResult<T> for Result<T, ScriptError> {
    fn script(self) -> Result<T, Box<EvalAltResult>> {
        self.map_err(script_error)
    }
}

/// Read the error map back out of a runtime error (host side).
pub fn error_map(error: &EvalAltResult) -> Option<Map> {
    match error {
        EvalAltResult::ErrorRuntime(value, _) => value.clone().try_cast::<Map>(),
        _ => None,
    }
}
