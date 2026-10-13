//! Feature-replay engine for pattern instances.
//!
//! Provides helpers to clone and remap seed features for a specific pattern
//! instance, and to replay a chain of seed features against a fresh BRep to
//! produce a tool body.

use crate::working_bodies::WorkingBodies;
use std::collections::HashMap;

use popkern_core::history::TrackingContext;
use popkern_core::naming::OperationId;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::ReplayError;
use crate::feature::Feature;
use crate::parameter::ParameterTable;
use crate::replay;

/// Clone a seed feature for a specific pattern instance.
///
/// Walks every [`popkern_core::naming::EntityName`] reference inside the
/// feature and replaces `name.operation` values that appear in `op_remap`
/// with the mapped instance-specific [`OperationId`].
///
/// For `Sketch` features with `OnFace` planes, the face reference is kept
/// as-is: the chain replays at the seed's placement and the caller moves
/// the resulting tool by the instance transform.
pub(crate) fn clone_seed_feature(
    feature: &Feature,
    op_remap: &HashMap<OperationId, OperationId>,
) -> Feature {
    let mut cloned = feature.clone();
    remap_feature_entity_names(&mut cloned, op_remap);
    cloned
}

/// Replay a chain of seed features to produce a tool body BRep.
///
/// Creates a fresh [`BRep`] and executes each (feature, operation_id) pair in
/// order, returning the resulting BRep.
///
/// # Errors
///
/// Returns [`ReplayError`] if any feature in the chain fails.
pub(crate) fn replay_seed_chain(
    seed_features: &[(Feature, OperationId)],
    tolerance: &Tolerance,
    param_table: &ParameterTable,
    body_breps: &mut WorkingBodies,
    naming_scheme: u32,
) -> Result<BRep, ReplayError> {
    let mut brep = BRep::new(*tolerance);
    // The chain's on-face sketch names its host; the chain replays into
    // a fresh tool body, so the host is where planes, exact profiles and
    // up-to targets resolve.
    let host = seed_features.iter().find_map(|(f, _)| match f {
        crate::feature::Feature::Sketch { plane, .. } => {
            crate::auto_merge::host_for_sketch_plane(plane, body_breps)
        }
        _ => None,
    });
    for (feature, op_id) in seed_features {
        let mut ctx = TrackingContext::new(*op_id).with_naming_scheme(naming_scheme);
        let exact_profiles = crate::replay_exact_extrude::exact_profile_set(
            feature,
            seed_features.iter().map(|(f, op)| (f, *op)),
            host.and_then(|h| body_breps.get(&h)).unwrap_or(&brep),
            param_table,
        );
        replay::execute_feature(
            &mut brep,
            feature,
            &mut ctx,
            tolerance,
            body_breps,
            param_table,
            &exact_profiles,
            host,
            // The seed chain's evidence is discarded (only the geometry
            // is patterned), so certifying its seams is unread work.
            popkern_blend::ContinuityCertification::Deferred,
        )?;
    }
    Ok(brep)
}

// ---------------------------------------------------------------------------
// EntityName remapping
// ---------------------------------------------------------------------------

/// Remap all [`popkern_core::naming::EntityName`] references inside a feature.
///
/// For each `EntityName` whose `operation` field appears in `op_remap`, the
/// operation is replaced with the mapped value.
fn remap_feature_entity_names(feature: &mut Feature, op_remap: &HashMap<OperationId, OperationId>) {
    crate::feature_names::for_each_entity_name_mut(feature, &mut |name| {
        if let Some(new_op) = op_remap.get(&name.operation) {
            name.operation = *new_op;
        }
    });
}
