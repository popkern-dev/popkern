//! Automatic boolean merge for Boss/Cut operations.
//!
//! When an operation's [`OperationMode`] is `Boss` or `Cut`, the tool body
//! produced by that feature should be automatically boolean-merged into the
//! host body.  This module provides the detection and execution logic.

use std::collections::HashMap;

use popkern_core::history::TrackingContext;
use popkern_core::naming_resolve::{Resolution, ResolvedEntity, resolve_name};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use slotmap::SlotMap;

use crate::error::ReplayError;
use crate::feature::{Feature, OperationMode, SketchPlaneDesc};
use crate::node::{BodyId, FeatureId, FeatureNode};
use popkern_core::naming::EntityName;

/// Extract the [`OperationMode`] from a [`Feature`] variant.
///
/// Returns [`OperationMode::NewBody`] for variants that do not carry a mode
/// field (Sketch, Fillet, Chamfer, Boolean, Shell, Draft, Reference, and
/// pattern variants).
pub(crate) fn feature_mode(feature: &Feature) -> OperationMode {
    match feature {
        Feature::Extrude { mode, .. }
        | Feature::Sweep { mode, .. }
        | Feature::Loft { mode, .. }
        | Feature::Pipe { mode, .. }
        | Feature::HelicalSweep { mode, .. }
        | Feature::Revolve { mode, .. } => *mode,
        Feature::Sketch { .. }
        | Feature::Fastener { .. }
        | Feature::TappedHole { .. }
        | Feature::ClearanceHole { .. }
        | Feature::RawBody { .. }
        | Feature::Fillet { .. }
        | Feature::Chamfer { .. }
        | Feature::Boolean { .. }
        | Feature::Transform { .. }
        | Feature::Shell { .. }
        | Feature::Thicken { .. }
        | Feature::OffsetSolid { .. }
        | Feature::Draft { .. }
        | Feature::Reference { .. }
        | Feature::LinearPattern { .. }
        | Feature::CircularPattern { .. }
        | Feature::CurvePattern { .. }
        | Feature::Mirror { .. } => OperationMode::NewBody,
    }
}

/// Detect the host body for an auto-merge operation.
///
/// Traces through the feature tree to find which body the given feature's
/// tool body should merge into:
///
/// 1. Get the feature's body ID.
/// 2. Find the Sketch feature that created this body (scan features for a
///    Sketch assigned to the same body).
/// 3. If the Sketch has `plane: SketchPlaneDesc::OnFace(face_name)`, resolve
///    which body owns that face.
/// 4. Return that body ID.
///
/// Returns `None` if no suitable host body can be found (e.g. the sketch is
/// on a reference plane rather than a model face, or no sketch exists for the
/// body).
pub(crate) fn detect_host_body(
    feature_id: FeatureId,
    features: &SlotMap<FeatureId, FeatureNode>,
    body_breps: &HashMap<BodyId, BRep>,
) -> Option<BodyId> {
    detect_host_body_among(
        feature_id,
        features,
        body_breps.iter().map(|(b, brep)| (*b, brep)),
    )
}

/// The body whose BRep resolves an on-face sketch plane's face — the
/// host a seed chain's on-face sketch (and the extrude after it) must
/// resolve against when replayed into a fresh body.
pub(crate) fn host_for_sketch_plane(
    plane: &SketchPlaneDesc,
    body_breps: &HashMap<BodyId, BRep>,
) -> Option<BodyId> {
    let face_name = match plane {
        SketchPlaneDesc::OnFace(name) | SketchPlaneDesc::OnFaceOffset(name, _) => name,
        _ => return None,
    };
    let mut bodies: Vec<&BodyId> = body_breps.keys().collect();
    bodies.sort_by_key(|b| b.0);
    bodies.into_iter().copied().find(|bid| {
        body_breps.get(bid).is_some_and(|brep| {
            matches!(
                resolve_name(face_name, brep),
                Resolution::Exact(ResolvedEntity::Face(_)) | Resolution::Split(_)
            )
        })
    })
}

/// The body whose BRep resolves an extrude's up-to / offset target
/// (surface or vertex) — the host such a feature is built against and,
/// for a Boss/Cut, merges into. `None` for blind ends.
pub(crate) fn host_for_extrude_target(
    feature: &Feature,
    own_body: BodyId,
    body_breps: &HashMap<BodyId, BRep>,
) -> Option<BodyId> {
    let Feature::Extrude {
        end_condition,
        second_direction,
        ..
    } = feature
    else {
        return None;
    };
    let target = |end: &crate::feature::ExtrudeEnd| -> Option<EntityName> {
        match end {
            crate::feature::ExtrudeEnd::UpToSurface { face }
            | crate::feature::ExtrudeEnd::OffsetFromSurface { face, .. } => Some(face.clone()),
            crate::feature::ExtrudeEnd::UpToVertex { vertex } => Some(vertex.clone()),
            _ => None,
        }
    };
    let name = target(end_condition).or_else(|| {
        second_direction
            .as_ref()
            .and_then(|s| target(&s.end_condition))
    })?;
    let mut bodies: Vec<&BodyId> = body_breps.keys().collect();
    bodies.sort_by_key(|b| b.0);
    bodies.into_iter().copied().find(|bid| {
        *bid != own_body
            && body_breps
                .get(bid)
                .is_some_and(|brep| !matches!(resolve_name(&name, brep), Resolution::NotFound))
    })
}

/// [`detect_host_body`] over any collection of live bodies.
pub(crate) fn detect_host_body_among<'a>(
    feature_id: FeatureId,
    features: &SlotMap<FeatureId, FeatureNode>,
    bodies: impl IntoIterator<Item = (BodyId, &'a BRep)>,
) -> Option<BodyId> {
    let node = features.get(feature_id)?;
    let body = node.body;

    // Find the Sketch feature that created this body.
    let sketch_face_name = features.values().find_map(|n| {
        if n.body != body {
            return None;
        }
        match &n.feature {
            Feature::Sketch {
                plane: SketchPlaneDesc::OnFace(face_name),
                ..
            } => Some(face_name.clone()),
            Feature::Sketch {
                plane: SketchPlaneDesc::OnFaceOffset(face_name, _),
                ..
            } => Some(face_name.clone()),
            _ => None,
        }
    })?;

    // Find which body owns the referenced face by checking all body BReps.
    for (bid, brep) in bodies {
        if bid == body {
            continue; // skip the feature's own body
        }
        match resolve_name(&sketch_face_name, brep) {
            Resolution::Exact(ResolvedEntity::Face(_)) | Resolution::Split(_) => return Some(bid),
            _ => {}
        }
    }

    None
}

/// Execute the boolean merge of a tool body into a host body.
///
/// Removes both bodies' BReps from `body_breps`, performs the boolean
/// operation (Union for Boss, Subtraction for Cut), and inserts the result
/// back under the host body key.  The tool body is consumed (not re-inserted).
///
/// # Errors
///
/// Returns [`ReplayError`] if either body's BRep is missing or the boolean
/// operation itself fails.
pub(crate) fn execute_auto_merge(
    mode: OperationMode,
    tool_body: BodyId,
    host_body: BodyId,
    body_breps: &mut HashMap<BodyId, BRep>,
    ctx: &mut TrackingContext,
    tolerance: &Tolerance,
) -> Result<(), ReplayError> {
    let tool_brep = body_breps
        .remove(&tool_body)
        .ok_or_else(|| ReplayError::FeatureFailed {
            label: "auto_merge".into(),
            error: format!("tool body {tool_body:?} has no BRep").into(),
        })?;

    let host_brep = body_breps
        .remove(&host_body)
        .ok_or_else(|| ReplayError::FeatureFailed {
            label: "auto_merge".into(),
            error: format!("host body {host_body:?} has no BRep").into(),
        })?;

    let bool_op = match mode {
        OperationMode::Boss => popkern_bool::BoolOperation::Union,
        OperationMode::Cut => popkern_bool::BoolOperation::Subtraction,
        OperationMode::NewBody => {
            // Should not be called with NewBody — put both bodies back.
            body_breps.insert(host_body, host_brep);
            body_breps.insert(tool_body, tool_brep);
            return Ok(());
        }
    };

    let result = popkern_bool::boolean(bool_op, host_brep, tool_brep, *tolerance, Some(ctx))
        .map_err(|e| ReplayError::FeatureFailed {
            label: "auto_merge".into(),
            error: e.to_string().into(),
        })?;

    body_breps.insert(host_body, result);
    // tool body is consumed — NOT re-inserted

    Ok(())
}
