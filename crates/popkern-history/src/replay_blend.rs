//! Fillet and chamfer execution for [`super::execute_feature`]: the
//! blend operation plus its continuity evidence, certified inline or
//! deferred per the tree's policy (see [`crate::deferred_certification`]).

use popkern_blend::ContinuityCertification;
use popkern_core::history::TrackingContext;
use popkern_core::naming::EntityName;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::replay_fns::{chamfer_specs, fillet_specs};
use crate::error::ReplayError;
use crate::feature_blend::{ChamferKindDesc, FilletRecipeDesc};
use crate::parameter::{ParameterId, ParameterTable};
use crate::replay_report::{FeatureEvidence, face_census};
use crate::replay_side_effects::ExecuteSideEffects;
use crate::resolve::resolve_edges;

/// Execute a fillet feature in place on `brep`.
#[allow(clippy::too_many_arguments)]
pub(super) fn execute_fillet(
    brep: &mut BRep,
    edges: &[EntityName],
    radius: ParameterId,
    recipe: &FilletRecipeDesc,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
    param_table: &ParameterTable,
    certification: ContinuityCertification,
) -> Result<ExecuteSideEffects, ReplayError> {
    let edge_handles = resolve_edges(edges, brep)?;
    let edges_requested = edge_handles.len();
    let specs = fillet_specs(edge_handles, radius, recipe, param_table)?;
    let options = recipe.to_kernel_options();
    let pre_faces = face_census(brep);
    let result = popkern_blend::fillet_edges(brep, &specs, &options, *tolerance, Some(tracking))
        .map_err(|e| ReplayError::FeatureFailed {
            label: "fillet".into(),
            error: crate::failure::FeatureFailure::blend(&e, Some(brep)).into(),
        })?;
    let evidence = FeatureEvidence::blend(
        &result,
        &pre_faces,
        edges_requested,
        Some(options.continuity),
        certification,
    );
    *brep = result.brep;
    Ok(ExecuteSideEffects::reported(evidence))
}

/// Execute a chamfer feature in place on `brep`.
#[allow(clippy::too_many_arguments)]
pub(super) fn execute_chamfer(
    brep: &mut BRep,
    edges: &[EntityName],
    distance: ParameterId,
    kind: &ChamferKindDesc,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
    param_table: &ParameterTable,
    certification: ContinuityCertification,
) -> Result<ExecuteSideEffects, ReplayError> {
    let edge_handles = resolve_edges(edges, brep)?;
    let edges_requested = edge_handles.len();
    let specs = chamfer_specs(brep, edge_handles, distance, kind, param_table)?;
    let pre_faces = face_census(brep);
    let result = popkern_blend::chamfer_edges(
        brep,
        &specs,
        &popkern_blend::ChamferOptions::default(),
        *tolerance,
        Some(tracking),
    )
    .map_err(|e| ReplayError::FeatureFailed {
        label: "chamfer".into(),
        error: crate::failure::FeatureFailure::blend(&e, Some(brep)).into(),
    })?;
    let evidence =
        FeatureEvidence::blend(&result, &pre_faces, edges_requested, None, certification);
    *brep = result.brep;
    Ok(ExecuteSideEffects::reported(evidence))
}
