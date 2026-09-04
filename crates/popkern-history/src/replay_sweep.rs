//! Sweep-family feature execution: sweep, loft, pipe, and helical
//! sweep arms of [`super::execute_feature`] (split per file-size rule).

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use std::collections::HashMap;

use super::replay_fns::{extract_edge_curve, extract_face_outer_curve, resolve_param};
use crate::error::ReplayError;
use crate::feature::AxisDesc;
use crate::node::BodyId;
use crate::parameter::{ParameterId, ParameterTable};
use crate::replay_report::FeatureEvidence;
use crate::replay_side_effects::ExecuteSideEffects;
use crate::resolve::{resolve_axis, resolve_edge, resolve_face};
use popkern_core::naming::EntityName;

/// The `Feature::Sweep` arm: profile face along a path edge.
#[allow(clippy::too_many_arguments)]
pub(super) fn execute_sweep(
    brep: &mut BRep,
    face: &EntityName,
    path_edge: &EntityName,
    twist: ParameterId,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
    body_breps: &HashMap<BodyId, BRep>,
    param_table: &ParameterTable,
    exact_profile: Option<&crate::replay_exact_extrude::ExactProfileWire>,
) -> Result<ExecuteSideEffects, ReplayError> {
    let tw = resolve_param(param_table, twist)?;
    let face_handle = resolve_face(face, brep)?;
    // The path edge usually lives on another body (the profile
    // sketch body has no path edges): resolve here first, then
    // across the other bodies — referenced, not consumed.
    let path_curve = if let Ok(e) = resolve_edge(path_edge, brep) {
        extract_edge_curve(brep, e)?
    } else if let Some(curve) = body_breps.values().find_map(|b| {
        resolve_edge(path_edge, b)
            .ok()
            .map(|e| extract_edge_curve(b, e))
    }) {
        curve?
    } else {
        return Err(ReplayError::FeatureFailed {
            label: "sweep".into(),
            error: format!("sweep path edge {path_edge:?} not found in any body").into(),
        });
    };
    // An analytic sketch sweeps its EXACT flat wire (the sweep
    // kernel's wire-contract form; arcs stay arcs) — the compiled
    // face's loop is only a vertex polyline (PK-408).
    let profile = match exact_profile {
        Some(exact) => popkern_sweep::SweepProfile::Wire(exact.flat.clone()),
        None => popkern_sweep::SweepProfile::Face {
            brep: Box::new(brep.clone()),
            face: face_handle,
        },
    };
    let options = popkern_sweep::SweepOptions {
        twist: tw,
        ..popkern_sweep::SweepOptions::default()
    };
    let result = popkern_sweep::sweep(&profile, &path_curve, &options, *tolerance, Some(tracking))
        .map_err(|e| ReplayError::FeatureFailed {
            label: "sweep".into(),
            error: crate::failure::FeatureFailure::sweep("sweep", &e).into(),
        })?;
    let evidence = FeatureEvidence::sweep(&result.warnings);
    *brep = result.brep;
    Ok(ExecuteSideEffects::reported(evidence))
}

/// The `Feature::Loft` arm: sections resolved cross-body, exact wires
/// where the sketches are analytic, other-body sections consumed.
pub(super) fn execute_loft(
    brep: &mut BRep,
    sections: &[EntityName],
    closed: bool,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
    body_breps: &HashMap<BodyId, BRep>,
    exact_profiles: &crate::replay_exact_extrude::ExactProfileSet,
) -> Result<ExecuteSideEffects, ReplayError> {
    // Sections are sketch faces, usually one per sketch body
    // (every sketch is its own body): resolve each name in the
    // executing body first, then across the other bodies —
    // names are operation-scoped, so the search is unambiguous.
    // A section found in another body consumes that body (the
    // boolean tool-body precedent): its geometry is spent
    // building the loft.
    let mut loft_sections = Vec::new();
    let mut consumed: Vec<BodyId> = Vec::new();
    for (i, section_name) in sections.iter().enumerate() {
        let (section_brep, owner, face_handle) = if let Ok(f) = resolve_face(section_name, brep) {
            (&*brep, None, f)
        } else if let Some((id, b, f)) = body_breps
            .iter()
            .find_map(|(id, b)| resolve_face(section_name, b).ok().map(|f| (*id, b, f)))
        {
            (b, Some(id), f)
        } else {
            return Err(ReplayError::FeatureFailed {
                label: "loft".into(),
                error: format!("loft section {i} ({section_name:?}) not found in any body").into(),
            });
        };
        // An analytic sketch lofts its EXACT world-space wire
        // (arcs stay arcs — the PK-408 rule); the compiled
        // face's loop is only a vertex polyline.
        let exact = exact_profiles.sections.get(i).and_then(|e| e.as_ref());
        let curve = match exact {
            Some(exact) => {
                let placement =
                    exact.plane.u_axis * exact.center[0] + exact.plane.v_axis * exact.center[1];
                popkern_core::geometry::NurbsCurve {
                    control_points: exact
                        .wire
                        .control_points
                        .iter()
                        .map(|p| *p + placement)
                        .collect(),
                    ..exact.wire.clone()
                }
            }
            None => extract_face_outer_curve(section_brep, face_handle)?,
        };
        if let Some(id) = owner
            && !consumed.contains(&id)
        {
            consumed.push(id);
        }
        loft_sections.push(popkern_sweep::LoftSection {
            curve,
            position: None,
        });
    }
    let options = popkern_sweep::LoftOptions {
        closed,
        ..popkern_sweep::LoftOptions::default()
    };
    let result = popkern_sweep::loft(&loft_sections, &options, *tolerance, Some(tracking))
        .map_err(|e| ReplayError::FeatureFailed {
            label: "loft".into(),
            error: crate::failure::FeatureFailure::sweep("loft", &e).into(),
        })?;
    let evidence = FeatureEvidence::sweep(&result.warnings);
    *brep = result.brep;
    Ok(ExecuteSideEffects {
        consumed_bodies: consumed,
        evidence,
    })
}

/// The `Feature::Pipe` arm: constant-radius tube along an edge of the
/// executing body (the recorder records the pipe on the path body).
pub(super) fn execute_pipe(
    brep: &mut BRep,
    path_edge: &EntityName,
    radius: ParameterId,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
    param_table: &ParameterTable,
) -> Result<ExecuteSideEffects, ReplayError> {
    let r = resolve_param(param_table, radius)?;
    let edge_handle = resolve_edge(path_edge, brep)?;
    let path_curve = extract_edge_curve(brep, edge_handle)?;
    let result = popkern_sweep::pipe(
        &path_curve,
        &popkern_sweep::PipeRadius::Constant(r),
        &popkern_sweep::SweepOptions::default(),
        *tolerance,
        Some(tracking),
    )
    .map_err(|e| ReplayError::FeatureFailed {
        label: "pipe".into(),
        error: crate::failure::FeatureFailure::sweep("pipe", &e).into(),
    })?;
    let evidence = FeatureEvidence::sweep(&result.warnings);
    *brep = result.brep;
    Ok(ExecuteSideEffects::reported(evidence))
}

/// Parameters of the `Feature::HelicalSweep` arm, resolved by the
/// caller from the feature's fields.
pub(super) struct HelicalSweepInputs<'a> {
    /// The profile face (stable name).
    pub face: &'a EntityName,
    /// The helix axis.
    pub axis: &'a AxisDesc,
    /// Radius at the helix start.
    pub start_radius: ParameterId,
    /// Radius at the helix end.
    pub end_radius: ParameterId,
    /// Axial advance per turn.
    pub pitch: ParameterId,
    /// Number of turns.
    pub turns: ParameterId,
    /// Phase about the axis.
    pub start_angle: ParameterId,
    /// Winding sense.
    pub left_handed: bool,
}

/// The `Feature::HelicalSweep` arm: profile face along a helix.
pub(super) fn execute_helical_sweep(
    brep: &mut BRep,
    inputs: HelicalSweepInputs<'_>,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
    param_table: &ParameterTable,
    exact_profile: Option<&crate::replay_exact_extrude::ExactProfileWire>,
) -> Result<ExecuteSideEffects, ReplayError> {
    let face_handle = resolve_face(inputs.face, brep)?;
    let (origin, direction) = resolve_axis(inputs.axis, brep)?;
    let spec = popkern_sweep::HelixSpec {
        axis_origin: origin,
        axis_dir: direction,
        start_radius: resolve_param(param_table, inputs.start_radius)?,
        end_radius: resolve_param(param_table, inputs.end_radius)?,
        pitch: popkern_sweep::HelixPitch::Constant(resolve_param(param_table, inputs.pitch)?),
        turns: resolve_param(param_table, inputs.turns)?,
        start_angle: resolve_param(param_table, inputs.start_angle)?,
        left_handed: inputs.left_handed,
    };
    // The profile contributes its shape; placement comes from
    // the helix (the kernel centres the profile on the path).
    // An analytic sketch sweeps its EXACT flat wire (PK-408).
    let profile = match exact_profile {
        Some(exact) => popkern_sweep::SweepProfile::Wire(exact.flat.clone()),
        None => popkern_sweep::SweepProfile::Face {
            brep: Box::new(brep.clone()),
            face: face_handle,
        },
    };
    let result = popkern_sweep::helical_sweep(
        &profile,
        &spec,
        &popkern_sweep::SweepOptions::default(),
        *tolerance,
        Some(tracking),
    )
    .map_err(|e| ReplayError::FeatureFailed {
        label: "helical_sweep".into(),
        error: crate::failure::FeatureFailure::sweep("helical_sweep", &e).into(),
    })?;
    let evidence = FeatureEvidence::sweep(&result.warnings);
    *brep = result.brep;
    Ok(ExecuteSideEffects::reported(evidence))
}
