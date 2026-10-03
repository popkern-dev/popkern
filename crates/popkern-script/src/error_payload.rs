//! The spec-§10 payload of every recorder error as JSON: one
//! projection shared by every skin (Rhai error maps, Python exception
//! attributes, MCP `tool-error` bodies), so the classes and field
//! names never drift between front-ends. Kernel spec paths in
//! `InvalidParameter` payloads are translated to the surface's
//! parameter names, with the kernel path kept beside them.

use popkern_history::{HistoryError, ReplayError};
use popkern_part::PartError;
use serde_json::{Value, json};

use crate::error::ScriptError;

impl ScriptError {
    /// The §10 payload: `kind`, `message` (this error's `Display`),
    /// and the class's fields. Floats are carried as numbers, never
    /// formatted into text.
    pub fn payload(&self) -> Value {
        let message = self.to_string();
        let mut value = match self {
            ScriptError::FeatureFailed {
                label,
                operation,
                error,
                report,
            } => json!({
                "kind": "FeatureFailed",
                "feature": { "label": label, "operation": operation },
                "error": translate_parameters(
                    serde_json::to_value(error.as_ref()).unwrap_or(Value::Null),
                ),
                "report": serde_json::to_value(report.as_ref()).unwrap_or(Value::Null),
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
}

/// The surface parameter name for a kernel spec path
/// (`specs[0].radius` → `radius`, `options.profile.rho` → `profile`,
/// `dist_b` → `distance_2`). Unknown paths map to their last segment.
pub fn surface_parameter(kernel_path: &str) -> String {
    let path = kernel_path.trim();
    let leaf = path.rsplit('.').next().unwrap_or(path);
    let leaf = leaf.split('[').next().unwrap_or(leaf);
    match leaf {
        "edge" => "edges".to_string(),
        "start" if path.contains("radius") => "radius".to_string(),
        "end" if path.contains("radius") => "law".to_string(),
        "radius_end" => "law".to_string(),
        "rho" => "profile".to_string(),
        "dist_a" => "distance".to_string(),
        "dist_b" => "distance_2".to_string(),
        "angle_rad" => "angle".to_string(),
        other => other.to_string(),
    }
}

/// Rewrite every `InvalidParameter` record in a payload: `parameter`
/// becomes the surface name and `kernel_path` keeps the kernel's.
fn translate_parameters(mut value: Value) -> Value {
    fn walk(value: &mut Value) {
        match value {
            Value::Object(map) => {
                let invalid = ["kind", "blend_kind"]
                    .iter()
                    .any(|key| map.get(*key).and_then(Value::as_str) == Some("InvalidParameter"));
                if invalid && let Some(Value::String(path)) = map.get("parameter").cloned() {
                    map.insert("kernel_path".into(), Value::String(path.clone()));
                    map.insert("parameter".into(), Value::String(surface_parameter(&path)));
                }
                for child in map.values_mut() {
                    walk(child);
                }
            }
            Value::Array(items) => items.iter_mut().for_each(walk),
            _ => {}
        }
    }
    walk(&mut value);
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
        HistoryError::AssemblyDocument => {
            json!({ "kind": "DocumentError", "reason": "assembly_document" })
        }
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
            "error": translate_parameters(
                serde_json::to_value(error.as_ref()).unwrap_or(Value::Null),
            ),
            "report": Value::Null,
        }),
        ReplayError::ResolutionFailed(name) => json!({
            "kind": "ResolutionFailed",
            "name": format!("{name:?}"),
            "selection": Value::Null,
            "resolved": 0,
            "expected": 1,
        }),
        ReplayError::CoreError(e) => json!({ "kind": "CoreError", "detail": e.to_string() }),
        ReplayError::HistoryError(e) => history_payload(e),
    }
}

fn part_payload(error: &PartError) -> Value {
    match error {
        PartError::Blend { selection, source } => {
            let blend = popkern_history::BlendDeclined::project(source, None);
            let mut value =
                translate_parameters(serde_json::to_value(&blend).unwrap_or(Value::Null));
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
        PartError::Fastener { stage, source } => json!({
            "kind": "FastenerFailed", "stage": stage, "detail": source.to_string(),
        }),
        PartError::NotApplicable(detail) => json!({ "kind": "NotApplicable", "detail": detail }),
        PartError::Cancelled => json!({ "kind": "Cancelled" }),
    }
}
