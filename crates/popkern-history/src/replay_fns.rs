//! Helper functions for the feature execution engine.
//!
//! Extracted from `replay.rs` to keep each source file under 500 lines.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::error::ReplayError;
use crate::feature::{BoolOpKind, ShellDirectionKind, SketchPlaneDesc};
use crate::parameter::{ParameterId, ParameterTable};
use crate::resolve::resolve_sketch_plane;

/// Resolve a [`ParameterId`] to its `f64` value via the parameter table.
pub(super) fn resolve_param(table: &ParameterTable, id: ParameterId) -> Result<f64, ReplayError> {
    table
        .get(id)
        .map(|p| p.value)
        .ok_or_else(|| ReplayError::FeatureFailed {
            label: "parameter".into(),
            error: format!("parameter {id:?} not found"),
        })
}

/// Execute a sketch feature: deserialize the sketch, place it on a plane,
/// assign stable names, solve, and convert to BRep.
pub(super) fn execute_sketch(
    brep: &mut BRep,
    plane_desc: &SketchPlaneDesc,
    sketch_data: &str,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
    param_table: &ParameterTable,
) -> Result<(), ReplayError> {
    let plane = resolve_sketch_plane(plane_desc, brep, param_table)?;

    let mut sketch: popkern_constraint::Sketch =
        serde_json::from_str(sketch_data).map_err(|e| ReplayError::FeatureFailed {
            label: "sketch".into(),
            error: format!("failed to deserialize sketch: {e}"),
        })?;

    // Re-solve the sketch so parametric changes take effect on replay.
    sketch.solve().map_err(|e| ReplayError::FeatureFailed {
        label: "sketch".into(),
        error: format!("constraint solve failed: {e}"),
    })?;

    let mut sketch_on_plane = popkern_sketch::SketchOnPlane {
        sketch,
        plane,
        element_names: std::collections::HashMap::new(),
    };
    sketch_on_plane.assign_element_names(tracking.operation);

    let result = sketch_on_plane
        .to_brep(*tolerance, Some(tracking))
        .map_err(|e| ReplayError::FeatureFailed {
            label: "sketch".into(),
            error: e.to_string(),
        })?;

    *brep = result.brep;
    Ok(())
}

/// Convert a serialisable [`ShellDirectionKind`] to the runtime
/// [`popkern_shell::ShellDirection`].
pub(super) fn to_shell_direction(kind: &ShellDirectionKind) -> popkern_shell::ShellDirection {
    match kind {
        ShellDirectionKind::Inward => popkern_shell::ShellDirection::Inward,
        ShellDirectionKind::Outward => popkern_shell::ShellDirection::Outward,
        ShellDirectionKind::Midplane => popkern_shell::ShellDirection::Midplane,
    }
}

/// Convert a serialisable [`BoolOpKind`] to the runtime
/// [`popkern_bool::BoolOperation`].
pub(super) fn to_bool_operation(kind: BoolOpKind) -> popkern_bool::BoolOperation {
    match kind {
        BoolOpKind::Union => popkern_bool::BoolOperation::Union,
        BoolOpKind::Intersection => popkern_bool::BoolOperation::Intersection,
        BoolOpKind::Subtraction => popkern_bool::BoolOperation::Subtraction,
    }
}

/// Extract a [`NurbsCurve`] from an edge's curve geometry.
pub(super) fn extract_edge_curve(brep: &BRep, edge: EdgeHandle) -> Result<NurbsCurve, ReplayError> {
    let curve_on_edge =
        brep.edges[edge]
            .curve
            .as_ref()
            .ok_or_else(|| ReplayError::FeatureFailed {
                label: "sweep".into(),
                error: "path edge has no curve geometry".into(),
            })?;
    Ok(curve_on_edge
        .geometry
        .to_nurbs_in_range((curve_on_edge.t_start, curve_on_edge.t_end)))
}

/// Extract the outer loop of a face as a single [`NurbsCurve`].
///
/// Uses the first edge's curve from the face's outer loop. For simple
/// cross-section faces (single-curve boundary) this is sufficient; complex
/// multi-edge boundaries would need curve concatenation (future enhancement).
pub(super) fn extract_face_outer_curve(
    brep: &BRep,
    face: FaceHandle,
) -> Result<NurbsCurve, ReplayError> {
    let outer_loop = brep.faces[face].outer_loop;
    let first_he = brep.loops[outer_loop].half_edge;
    let first_edge = brep.half_edges[first_he].edge;
    extract_edge_curve(brep, first_edge)
}

/// Extract the plane origin and normal from a face's planar surface.
///
/// Returns an error if the face has no surface or its surface is not a plane.
pub(super) fn extract_plane_from_face(
    brep: &BRep,
    face: FaceHandle,
) -> Result<(Point3, Vec3), ReplayError> {
    let sof = brep.faces[face]
        .surface
        .as_ref()
        .ok_or_else(|| ReplayError::FeatureFailed {
            label: "draft".into(),
            error: "neutral plane face has no surface geometry".into(),
        })?;
    match &sof.geometry {
        SurfaceGeometry::Plane { origin, normal, .. } => {
            let n = if sof.orientation { *normal } else { -*normal };
            Ok((*origin, n))
        }
        _ => Err(ReplayError::FeatureFailed {
            label: "draft".into(),
            error: "neutral plane face is not planar".into(),
        }),
    }
}
