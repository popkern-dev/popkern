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
    let (specs, edges_requested) = if recipe.edge_radii.is_empty() {
        let edge_handles = resolve_edges(edges, brep)?;
        let count = edge_handles.len();
        (
            fillet_specs(edge_handles, radius, recipe, param_table)?,
            count,
        )
    } else {
        grouped_fillet_specs(brep, edges, radius, recipe, param_table)?
    };
    let mut options = recipe.to_kernel_options();
    options.corner_setback.ends = end_setbacks(brep, edges, recipe)?;
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

/// Kernel fillet specs for a fillet with radius groups (PK-602): every
/// grouped name at its group's radius, the rest at the feature radius.
/// Refuses a tree whose groups name an edge the feature omits, share an
/// edge, or combine with a radius law.
fn grouped_fillet_specs(
    brep: &BRep,
    edges: &[EntityName],
    radius: ParameterId,
    recipe: &FilletRecipeDesc,
    param_table: &ParameterTable,
) -> Result<(Vec<popkern_blend::FilletSpec>, usize), ReplayError> {
    let failed = |reason: String| ReplayError::FeatureFailed {
        label: "fillet".into(),
        error: reason.into(),
    };
    if recipe.radius_law != crate::feature_blend::FilletRadiusLaw::Constant {
        return Err(failed(
            "radius laws are per fillet; radius groups take constant radii".into(),
        ));
    }
    let mut grouped: std::collections::HashMap<&EntityName, usize> =
        std::collections::HashMap::new();
    for (g, group) in recipe.edge_radii.iter().enumerate() {
        for name in &group.edges {
            if !edges.contains(name) {
                return Err(failed(format!(
                    "radius group {} names edge {name:?} outside the fillet",
                    g + 1
                )));
            }
            if let Some(prev) = grouped.insert(name, g + 1) {
                return Err(failed(format!(
                    "edge {name:?} is in radius groups {prev} and {}",
                    g + 1
                )));
            }
        }
    }
    let base: Vec<EntityName> = edges
        .iter()
        .filter(|n| !grouped.contains_key(n))
        .cloned()
        .collect();
    let mut specs = Vec::new();
    let mut push = |names: &[EntityName], r: ParameterId| -> Result<(), ReplayError> {
        let value = super::replay_fns::resolve_param(param_table, r)?;
        for edge in resolve_edges(names, brep)? {
            specs.push(popkern_blend::FilletSpec {
                edge,
                radius: popkern_blend::FilletRadius::Constant(value),
            });
        }
        Ok(())
    };
    push(&base, radius)?;
    for group in &recipe.edge_radii {
        push(&group.edges, group.radius)?;
    }
    let count = specs.len();
    Ok((specs, count))
}

/// The kernel's per-end setbacks for the recipe's groups (PK-624): each
/// group edge's end at the corner it shares with one of the group's
/// `at` edges. Refuses a group naming an edge outside the fillet, an
/// edge that meets no `at` edge, or one that meets them at both ends.
fn end_setbacks(
    brep: &BRep,
    edges: &[EntityName],
    recipe: &FilletRecipeDesc,
) -> Result<Vec<popkern_blend::EdgeEndSetback>, ReplayError> {
    let failed = |reason: String| ReplayError::FeatureFailed {
        label: "fillet".into(),
        error: reason.into(),
    };
    let ends_of = |e: popkern_core::topology::EdgeHandle| {
        let mut ends: Vec<_> = brep
            .half_edges
            .values()
            .filter(|h| h.edge == e)
            .map(|h| h.origin)
            .collect();
        ends.sort();
        ends.dedup();
        ends
    };
    let mut out = Vec::new();
    for (g, group) in recipe.end_setbacks.iter().enumerate() {
        if let Some(name) = group
            .edges
            .iter()
            .chain(&group.at)
            .find(|n| !edges.contains(n))
        {
            return Err(failed(format!(
                "setback end group {g} names edge {name:?} outside the fillet"
            )));
        }
        let at = resolve_edges(&group.at, brep)?;
        for edge in resolve_edges(&group.edges, brep)? {
            let corners: Vec<_> = ends_of(edge)
                .into_iter()
                .filter(|v| at.iter().any(|&f| f != edge && ends_of(f).contains(v)))
                .collect();
            let [vertex] = corners[..] else {
                return Err(failed(format!(
                    "setback end group {g}: edge {edge:?} meets its `at` edges at {} ends, not one",
                    corners.len()
                )));
            };
            out.push(popkern_blend::EdgeEndSetback {
                vertex,
                edge,
                distance: group.distance,
            });
        }
    }
    Ok(out)
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
