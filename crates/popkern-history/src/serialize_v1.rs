// Structural lint marker: included from serialize.rs via `mod serialize_v1`.
//! V1 → V2 migration helpers for the `.pkn` document format.
//!
//! V1 stored feature dimensions as bare `f64` values.  V2 replaces them with
//! [`ParameterId`] handles backed by a full [`ParameterTable`].
//!
//! This module is private to [`crate::serialize`].

use serde_json::Value;

use crate::FeatureTree;
use crate::error::HistoryError;
use crate::feature::{BoolOpKind, Feature, ShellDirectionKind, SketchPlaneDesc, Vec3Desc};
use crate::node::BodyId;
use crate::parameter::ParameterId;
use crate::unit::Unit;
use popkern_core::naming::EntityName;

/// Migrate a v1 JSON tree (bare `f64` dimensions) into a v2 [`FeatureTree`].
///
/// V1 stored feature dimensions as bare `f64` values directly in the JSON.
/// V2 replaces them with [`ParameterId`] handles backed by a
/// [`crate::parameter::ParameterTable`].
///
/// The migration uses a purely JSON-value approach to avoid coupling to the
/// exact SlotMap wire format:
///
/// 1. Walk the v1 `features` slots in timeline order.
/// 2. For each feature, extract its dimension `f64` fields and create
///    corresponding parameter entries in a fresh tree.
/// 3. Return the reconstructed v2 [`FeatureTree`].
///
/// Because no real v1 `.pkn` files exist in the wild (the crate is new), this
/// is primarily a safety net.  If the v1 JSON is structurally invalid a
/// [`HistoryError::SerializeError`] is returned with a descriptive message.
pub(super) fn migrate_v1(tree_value: Value) -> Result<FeatureTree, HistoryError> {
    // Build a fresh tree to act as the accumulator.
    let mut new_tree = FeatureTree::new();

    // Copy simple scalar fields that are version-agnostic.
    if let Some(v) = tree_value.get("next_op_id").and_then(Value::as_u64) {
        new_tree.next_op_id = v;
    }
    if let Some(v) = tree_value.get("rollback_index").and_then(Value::as_u64) {
        new_tree.rollback_index = Some(v as usize);
    }
    if let Some(v) = tree_value.get("frozen_index").and_then(Value::as_u64) {
        new_tree.frozen_index = Some(v as usize);
    }

    // The SlotMap serializes as an array of slots.  Each occupied slot has a
    // non-null `value` field containing the `FeatureNode` JSON.
    let Some(slots) = tree_value.get("features").and_then(Value::as_array) else {
        return Err(HistoryError::SerializeError(
            "v1 migration: missing 'features' array".to_string(),
        ));
    };

    // Collect occupied nodes keyed by their slot index.
    let occupied: Vec<(usize, &Value)> = slots
        .iter()
        .enumerate()
        .filter_map(|(i, slot)| {
            slot.get("value")
                .filter(|v| !v.is_null())
                .map(|node| (i, node))
        })
        .collect();

    // The timeline is a Vec of SlotMap keys serialized as `{idx, version}`.
    // Follow the timeline ordering so features are re-inserted in the correct
    // evaluation order.
    let timeline_indices: Vec<u32> =
        if let Some(tl) = tree_value.get("timeline").and_then(Value::as_array) {
            tl.iter()
                .filter_map(|k| k.get("idx").and_then(Value::as_u64).map(|v| v as u32))
                .collect()
        } else {
            return Err(HistoryError::SerializeError(
                "v1 migration: missing 'timeline' array".to_string(),
            ));
        };

    // Build a map from slot idx → node JSON for O(1) lookup.
    let slot_map: std::collections::HashMap<u32, &Value> = occupied
        .iter()
        .map(|(i, node)| (*i as u32, *node))
        .collect();

    // Process features in timeline order.
    for slot_idx in &timeline_indices {
        let node = slot_map.get(slot_idx).ok_or_else(|| {
            HistoryError::SerializeError(format!(
                "v1 migration: timeline references missing slot {slot_idx}"
            ))
        })?;

        let label = node
            .get("label")
            .and_then(Value::as_str)
            .unwrap_or("Unnamed")
            .to_string();

        let feature_json = node.get("feature").ok_or_else(|| {
            HistoryError::SerializeError("v1 migration: node missing 'feature' field".to_string())
        })?;

        let v2_feature = convert_v1_feature_json(feature_json, &label, &mut new_tree)?;
        new_tree.add_feature(v2_feature, &label);
    }

    new_tree.rebuild_indices();
    Ok(new_tree)
}

/// Allocate a fresh parameter in `tree` with a unique name.
///
/// Tries the bare `base_name` first, then appends a numeric suffix until an
/// unused name is found.
fn alloc_param(
    tree: &mut FeatureTree,
    base_name: &str,
    value: f64,
    unit: Unit,
) -> Result<ParameterId, HistoryError> {
    if let Ok(id) = tree.parameters_mut().add_global(base_name, value, unit) {
        return Ok(id);
    }
    for n in 1u32..=10_000 {
        let candidate = format!("{base_name}_{n}");
        if let Ok(id) = tree.parameters_mut().add_global(&candidate, value, unit) {
            return Ok(id);
        }
    }
    Err(HistoryError::SerializeError(
        "v1 migration: could not allocate unique parameter name".to_string(),
    ))
}

/// Convert a v1 feature JSON value to a v2 [`Feature`].
///
/// Dimension fields that were `f64` in v1 are replaced by new parameter IDs
/// inserted into `tree`'s parameter table.
fn convert_v1_feature_json(
    feature_json: &Value,
    label: &str,
    tree: &mut FeatureTree,
) -> Result<Feature, HistoryError> {
    // Serde tagged enums produce `{"Variant": {...}}` for struct variants
    // and plain strings `"Variant"` for unit variants.
    let (variant, data) = if let Some(obj) = feature_json.as_object() {
        if let Some((k, v)) = obj.iter().next().filter(|_| obj.len() == 1) {
            (k.as_str(), v.clone())
        } else {
            return Err(HistoryError::SerializeError(
                "v1 migration: unexpected feature JSON shape".to_string(),
            ));
        }
    } else if let Some(s) = feature_json.as_str() {
        (s, Value::Object(Default::default()))
    } else {
        return Err(HistoryError::SerializeError(
            "v1 migration: feature JSON must be an object or string".to_string(),
        ));
    };

    let mut counter = 0u32;
    let mut alloc = |value: f64, unit: Unit| -> Result<ParameterId, HistoryError> {
        counter += 1;
        let name = format!("_migrated_{label}_{counter}");
        alloc_param(tree, &name, value, unit)
    };

    let err = |msg: &str| -> HistoryError {
        HistoryError::SerializeError(format!("v1 migration: {msg}"))
    };

    let get_f64 = |obj: &Value, field: &str| -> Result<f64, HistoryError> {
        obj.get(field)
            .and_then(Value::as_f64)
            .ok_or_else(|| err(&format!("missing f64 field '{field}'")))
    };

    let get_entity = |obj: &Value, field: &str| -> Result<EntityName, HistoryError> {
        let v = obj
            .get(field)
            .ok_or_else(|| err(&format!("missing field '{field}'")))?;
        serde_json::from_value(v.clone())
            .map_err(|e| err(&format!("bad EntityName '{field}': {e}")))
    };

    let get_entities = |obj: &Value, field: &str| -> Result<Vec<EntityName>, HistoryError> {
        let v = obj
            .get(field)
            .ok_or_else(|| err(&format!("missing field '{field}'")))?;
        serde_json::from_value(v.clone())
            .map_err(|e| err(&format!("bad EntityName array '{field}': {e}")))
    };

    let get_vec3 = |obj: &Value, field: &str| -> Result<Vec3Desc, HistoryError> {
        let v = obj
            .get(field)
            .ok_or_else(|| err(&format!("missing field '{field}'")))?;
        serde_json::from_value(v.clone()).map_err(|e| err(&format!("bad Vec3Desc '{field}': {e}")))
    };

    let get_arr3 = |obj: &Value, field: &str| -> Result<[f64; 3], HistoryError> {
        let v = obj
            .get(field)
            .ok_or_else(|| err(&format!("missing field '{field}'")))?;
        serde_json::from_value(v.clone()).map_err(|e| err(&format!("bad [f64;3] '{field}': {e}")))
    };

    match variant {
        "Sketch" => {
            let plane =
                convert_v1_plane_json(data.get("plane").unwrap_or(&Value::Null), label, tree)?;
            let sketch_data = data
                .get("sketch_data")
                .and_then(Value::as_str)
                .unwrap_or("")
                .to_string();
            Ok(Feature::Sketch { plane, sketch_data })
        }
        "Extrude" => {
            let face = get_entity(&data, "face")?;
            let direction = get_vec3(&data, "direction")?;
            let distance = alloc(get_f64(&data, "distance")?, Unit::Meters)?;
            Ok(Feature::Extrude {
                face,
                direction,
                distance,
            })
        }
        "Fillet" => {
            let edges = get_entities(&data, "edges")?;
            let radius = alloc(get_f64(&data, "radius")?, Unit::Meters)?;
            Ok(Feature::Fillet { edges, radius })
        }
        "Chamfer" => {
            let edges = get_entities(&data, "edges")?;
            let distance = alloc(get_f64(&data, "distance")?, Unit::Meters)?;
            Ok(Feature::Chamfer { edges, distance })
        }
        "Boolean" => {
            let operation: BoolOpKind = serde_json::from_value(
                data.get("operation")
                    .cloned()
                    .ok_or_else(|| err("missing 'operation'"))?,
            )
            .map_err(|e| err(&format!("bad BoolOpKind: {e}")))?;
            let tool_body: BodyId = serde_json::from_value(
                data.get("tool_body")
                    .cloned()
                    .ok_or_else(|| err("missing 'tool_body'"))?,
            )
            .map_err(|e| err(&format!("bad BodyId: {e}")))?;
            Ok(Feature::Boolean {
                operation,
                tool_body,
            })
        }
        "Shell" => {
            let faces_to_remove = get_entities(&data, "faces_to_remove")?;
            let thickness = alloc(get_f64(&data, "thickness")?, Unit::Meters)?;
            let direction: ShellDirectionKind = serde_json::from_value(
                data.get("direction")
                    .cloned()
                    .ok_or_else(|| err("missing 'direction'"))?,
            )
            .map_err(|e| err(&format!("bad ShellDirectionKind: {e}")))?;
            Ok(Feature::Shell {
                faces_to_remove,
                thickness,
                direction,
            })
        }
        "Draft" => {
            let faces = get_entities(&data, "faces")?;
            let angle = alloc(get_f64(&data, "angle")?, Unit::Radians)?;
            let neutral_plane = get_entity(&data, "neutral_plane")?;
            Ok(Feature::Draft {
                faces,
                angle,
                neutral_plane,
            })
        }
        "Sweep" => {
            let face = get_entity(&data, "face")?;
            let path_edge = get_entity(&data, "path_edge")?;
            let twist = alloc(get_f64(&data, "twist")?, Unit::Radians)?;
            Ok(Feature::Sweep {
                face,
                path_edge,
                twist,
            })
        }
        "Loft" => {
            let sections = get_entities(&data, "sections")?;
            let closed = data.get("closed").and_then(Value::as_bool).unwrap_or(false);
            Ok(Feature::Loft { sections, closed })
        }
        "Pipe" => {
            let path_edge = get_entity(&data, "path_edge")?;
            let radius = alloc(get_f64(&data, "radius")?, Unit::Meters)?;
            Ok(Feature::Pipe { path_edge, radius })
        }
        "Revolve" => {
            let face = get_entity(&data, "face")?;
            let axis_origin = get_arr3(&data, "axis_origin")?;
            let axis_direction = get_arr3(&data, "axis_direction")?;
            let angle = alloc(get_f64(&data, "angle")?, Unit::Radians)?;
            Ok(Feature::Revolve {
                face,
                axis_origin,
                axis_direction,
                angle,
            })
        }
        other => Err(err(&format!("unknown v1 feature variant '{other}'"))),
    }
}

/// Convert a v1 `SketchPlaneDesc` JSON value to a v2 [`SketchPlaneDesc`].
///
/// `OnFaceOffset` carries a bare `f64` offset in v1; migration creates a
/// parameter for it.
fn convert_v1_plane_json(
    plane_json: &Value,
    label: &str,
    tree: &mut FeatureTree,
) -> Result<SketchPlaneDesc, HistoryError> {
    let err = |msg: &str| HistoryError::SerializeError(format!("v1 migration: plane: {msg}"));

    // Unit variants are plain strings; struct variants are `{"Variant": {...}}`.
    if let Some(s) = plane_json.as_str() {
        return match s {
            "XY" => Ok(SketchPlaneDesc::XY),
            "XZ" => Ok(SketchPlaneDesc::XZ),
            "YZ" => Ok(SketchPlaneDesc::YZ),
            other => Err(err(&format!("unknown plane variant '{other}'"))),
        };
    }

    if let Some(obj) = plane_json.as_object() {
        if let Some(inner) = obj.get("OnFace") {
            let name: EntityName = serde_json::from_value(inner.clone())
                .map_err(|e| err(&format!("bad OnFace EntityName: {e}")))?;
            return Ok(SketchPlaneDesc::OnFace(name));
        }
        if let Some(inner) = obj.get("OnFaceOffset") {
            // V1: [EntityName, f64]
            let arr = inner
                .as_array()
                .ok_or_else(|| err("OnFaceOffset must be array"))?;
            if arr.len() != 2 {
                return Err(err("OnFaceOffset must have exactly 2 elements"));
            }
            let name: EntityName = serde_json::from_value(arr[0].clone())
                .map_err(|e| err(&format!("bad OnFaceOffset EntityName: {e}")))?;
            let offset_val = arr[1]
                .as_f64()
                .ok_or_else(|| err("OnFaceOffset offset must be f64"))?;
            let param_name = format!("_migrated_{label}_offset");
            let id = alloc_param(tree, &param_name, offset_val, Unit::Meters)?;
            return Ok(SketchPlaneDesc::OnFaceOffset(name, id));
        }
        if let Some(inner) = obj.get("Custom") {
            let origin: [f64; 3] =
                serde_json::from_value(inner.get("origin").cloned().unwrap_or(Value::Null))
                    .map_err(|e| err(&format!("bad Custom origin: {e}")))?;
            let normal: [f64; 3] =
                serde_json::from_value(inner.get("normal").cloned().unwrap_or(Value::Null))
                    .map_err(|e| err(&format!("bad Custom normal: {e}")))?;
            return Ok(SketchPlaneDesc::Custom { origin, normal });
        }
    }

    Err(err(&format!("unrecognized plane JSON: {plane_json}")))
}
