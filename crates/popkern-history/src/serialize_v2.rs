// Structural lint marker: included from serialize.rs via `mod serialize_v2`.
//! V2 → V3 migration helpers for the `.pkn` document format.
//!
//! V2 stored `Extrude` with `direction: Vec3Desc` and `distance: ParameterId`
//! (a flat `ParameterId`, not wrapped in `ExtrudeEnd`).  `Revolve` stored
//! `axis_origin: [f64; 3]`, `axis_direction: [f64; 3]`, and `angle:
//! ParameterId` as top-level fields instead of the new [`AxisDesc`] /
//! [`RevolveEnd`] wrappers.
//!
//! V3 replaces these with `end_condition: ExtrudeEnd` / `axis: AxisDesc` /
//! `end_condition: RevolveEnd` and adds optional `mode`, `draft_angle`, and
//! `second_direction` fields (all defaulting gracefully via `#[serde(default)]`).
//!
//! This module is private to [`crate::serialize`].

use serde_json::Value;

use crate::FeatureTree;
use crate::error::HistoryError;

/// Migrate a v2 JSON tree into a v3 [`FeatureTree`].
///
/// The tree JSON is structurally identical to v3 except that `Extrude` and
/// `Revolve` use the old flat field layout.  All other variants are
/// binary-compatible with v3 (new optional fields default via
/// `#[serde(default)]`).
///
/// Strategy: deserialize the outer tree skeleton normally, then walk the
/// `features` slot array and re-parse each feature, converting old
/// `Extrude`/`Revolve` shapes to their v3 equivalents in-place before
/// handing the whole value to the standard `serde_json::from_value` decoder.
///
/// # Errors
///
/// Returns [`HistoryError::SerializeError`] if the JSON is structurally
/// invalid.
pub(super) fn migrate_v2(tree_value: Value) -> Result<FeatureTree, HistoryError> {
    let mut tree_obj = tree_value
        .into_object()
        .map_err(|_| HistoryError::SerializeError("v2 migration: tree must be an object".into()))?;

    // Patch feature variants inside the SlotMap array.
    if let Some(Value::Array(slots)) = tree_obj.get_mut("features") {
        for slot in slots.iter_mut() {
            patch_slot(slot)?;
        }
    }

    let patched = Value::Object(tree_obj);
    let mut tree: FeatureTree = serde_json::from_value(patched)
        .map_err(|e| HistoryError::SerializeError(format!("v2 migration: {e}")))?;
    tree.rebuild_indices();
    Ok(tree)
}

/// Patch a single SlotMap slot in-place.
///
/// Occupied slots have a non-null `value` field containing a `FeatureNode`
/// whose `feature` key holds the enum JSON.  Vacant or null slots are left
/// untouched.
fn patch_slot(slot: &mut Value) -> Result<(), HistoryError> {
    let node = match slot
        .as_object_mut()
        .and_then(|o| o.get_mut("value"))
        .filter(|v| !v.is_null())
    {
        Some(n) => n,
        None => return Ok(()),
    };

    let feature_val = match node.as_object_mut().and_then(|o| o.get_mut("feature")) {
        Some(f) => f,
        None => return Ok(()),
    };

    patch_feature(feature_val)
}

/// Patch a single feature JSON value in-place.
///
/// Only `Extrude` and `Revolve` need structural changes; all other variants
/// are already v3-compatible.
///
/// Only visits a variant when the associated value is a non-null object, which
/// guards against corrupt documents where a key is present but set to `null`.
fn patch_feature(feature: &mut Value) -> Result<(), HistoryError> {
    let obj = match feature.as_object_mut() {
        Some(o) => o,
        None => return Ok(()),
    };

    // Only call patch_extrude / patch_revolve when the variant value is an
    // actual object.  `get_mut` returns `Some(&mut Value::Null)` for keys
    // whose value is null, so we must guard with `filter`.
    if let Some(extrude_val) = obj.get_mut("Extrude").filter(|v| v.is_object()) {
        patch_extrude(extrude_val)?;
    } else if let Some(revolve_val) = obj.get_mut("Revolve").filter(|v| v.is_object()) {
        patch_revolve(revolve_val)?;
    }
    Ok(())
}

/// Convert old v2 `Extrude` fields to v3 layout.
///
/// V2 had:
/// ```json
/// { "face": ..., "direction": ..., "distance": <ParameterId> }
/// ```
///
/// V3 needs:
/// ```json
/// { "face": ..., "direction": ..., "end_condition": {"Blind": {"distance": <ParameterId>}} }
/// ```
///
/// New optional fields (`mode`, `draft_angle`, `second_direction`) are absent
/// and will default via `#[serde(default)]`.
fn patch_extrude(extrude: &mut Value) -> Result<(), HistoryError> {
    let obj = extrude.as_object_mut().ok_or_else(|| {
        HistoryError::SerializeError("v2 migration: Extrude must be an object".into())
    })?;

    // Only migrate if old `distance` field is present (v2 shape).
    if let Some(distance_val) = obj.remove("distance") {
        let end_condition = serde_json::json!({ "Blind": { "distance": distance_val } });
        obj.insert("end_condition".to_string(), end_condition);

        // Wrap bare `direction` in `Some` if it's not already `null`.
        // In v2 `direction` was `Vec3Desc` (not Option); in v3 it is
        // `Option<Vec3Desc>`.  serde serializes `Some(v)` as just `v` and
        // `None` as `null` / absent, so we don't actually need to change the
        // JSON here — the outer `Option<Vec3Desc>` will deserialize from the
        // bare `Vec3Desc` value correctly when `direction` is present.
    }
    Ok(())
}

/// Convert old v2 `Revolve` fields to v3 layout.
///
/// V2 had:
/// ```json
/// { "face": ..., "axis_origin": [f64;3], "axis_direction": [f64;3], "angle": <ParameterId> }
/// ```
///
/// V3 needs:
/// ```json
/// {
///   "face": ...,
///   "axis": { "Literal": { "origin": [f64;3], "direction": [f64;3] } },
///   "end_condition": { "Blind": { "angle": <ParameterId> } }
/// }
/// ```
fn patch_revolve(revolve: &mut Value) -> Result<(), HistoryError> {
    let obj = revolve.as_object_mut().ok_or_else(|| {
        HistoryError::SerializeError("v2 migration: Revolve must be an object".into())
    })?;

    // Only migrate if the old flat fields are present.
    let has_old_fields = obj.contains_key("axis_origin")
        || obj.contains_key("axis_direction")
        || obj.contains_key("angle");

    if !has_old_fields {
        return Ok(());
    }

    let origin = obj.remove("axis_origin").ok_or_else(|| {
        HistoryError::SerializeError("v2 migration: Revolve missing 'axis_origin' field".into())
    })?;
    let direction = obj.remove("axis_direction").ok_or_else(|| {
        HistoryError::SerializeError("v2 migration: Revolve missing 'axis_direction' field".into())
    })?;
    let angle = obj.remove("angle").ok_or_else(|| {
        HistoryError::SerializeError("v2 migration: Revolve missing 'angle' field".into())
    })?;

    let axis = serde_json::json!({
        "Literal": { "origin": origin, "direction": direction }
    });
    let end_condition = serde_json::json!({ "Blind": { "angle": angle } });

    obj.insert("axis".to_string(), axis);
    obj.insert("end_condition".to_string(), end_condition);
    Ok(())
}

// ---------------------------------------------------------------------------
// Extension trait to make the Value API more ergonomic locally.
// ---------------------------------------------------------------------------

trait ValueExt {
    /// Consume `self` and return the inner `serde_json::Map`, or the original
    /// `Value` as an error if it's not an object.
    fn into_object(self) -> Result<serde_json::Map<String, Value>, Value>;
}

impl ValueExt for Value {
    fn into_object(self) -> Result<serde_json::Map<String, Value>, Value> {
        match self {
            Value::Object(m) => Ok(m),
            other => Err(other),
        }
    }
}
