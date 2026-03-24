//! Topology construction -- apply fillet/chamfer geometry to the BRep.
//!
//! Clones the input BRep, splits perpendicular edges at the fillet offset
//! distance, splits faces along trim curves, merges consumed face regions
//! by removing the original fillet edge, and assigns the fillet/chamfer
//! surface geometry to the resulting face.

pub(crate) mod helpers;

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::analysis::EdgeAnalysis;
use crate::error::BlendFailure;

/// Geometry computed for one fillet/chamfer edge.
#[derive(Debug)]
pub struct EdgeGeometry {
    /// The edge analysis.
    pub analysis: EdgeAnalysis,
    /// Fillet or chamfer surface.
    pub surface: SurfaceGeometry,
    /// Trim curve on face 1 (3D).
    pub trim_3d_1: NurbsCurve,
    /// Trim curve on face 2 (3D).
    pub trim_3d_2: NurbsCurve,
    /// Pcurve on face 1.
    pub pcurve_1: NurbsCurve,
    /// Pcurve on face 2.
    pub pcurve_2: NurbsCurve,
}

/// Apply computed fillet/chamfer geometry to a BRep.
///
/// Clones the BRep, then for each edge geometry:
/// 1. Splits perpendicular edges at the fillet radius offset.
/// 2. Splits each adjacent face along corrected trim curves.
/// 3. Identifies consumed faces (portions eaten by the fillet).
/// 4. Merges consumed faces by removing the original fillet edge (`kef`).
/// 5. Assigns the fillet/chamfer surface to the merged face.
///
/// Returns the modified BRep or a list of per-edge failures.
pub fn apply_blend_geometry(
    brep: &BRep,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Result<BRep, Vec<(EdgeHandle, BlendFailure)>> {
    let mut result = brep.clone();
    let mut failures: Vec<(EdgeHandle, BlendFailure)> = Vec::new();

    for eg in edge_geometries {
        if let Err(reason) = apply_single_blend(&mut result, eg, tolerance) {
            failures.push((eg.analysis.edge, reason));
        }
    }

    if !failures.is_empty() && failures.len() == edge_geometries.len() {
        return Err(failures);
    }

    Ok(result)
}

/// Apply a single fillet/chamfer blend to the BRep.
///
/// Performs the four-phase topology modification for one edge:
/// 1. Split perpendicular edges at the offset distance (creates anchor vertices).
/// 2. Split faces along trim curves connecting the anchor vertices.
/// 3. Identify consumed faces (the portions that contain the fillet edge).
/// 4. Merge consumed faces via `kef` and assign fillet surface.
fn apply_single_blend(
    brep: &mut BRep,
    eg: &EdgeGeometry,
    tolerance: &Tolerance,
) -> Result<(), BlendFailure> {
    let fillet_edge = eg.analysis.edge;

    // Estimate the fillet radius from the trim curve geometry.
    let radius = helpers::estimate_radius_from_geometry(brep, fillet_edge, &eg.trim_3d_1);
    if radius < tolerance.linear {
        return Err(BlendFailure::InvalidTopology {
            reason: "estimated fillet radius is too small".into(),
        });
    }

    // Phase 0: Split perpendicular edges at the fillet offset distance.
    let split_info = helpers::split_perpendicular_edges_by_radius(
        brep,
        eg.analysis.face_1,
        eg.analysis.face_2,
        fillet_edge,
        radius,
        tolerance,
    )
    .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;

    // Compute pcurves for the corrected trim curves on each face.
    let pcurve_1 = helpers::pcurve_for_face(
        brep,
        eg.analysis.face_1,
        &brep.vertices[split_info.face_1_anchors.0].position,
        &brep.vertices[split_info.face_1_anchors.1].position,
    )
    .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;

    let pcurve_2 = helpers::pcurve_for_face(
        brep,
        eg.analysis.face_2,
        &brep.vertices[split_info.face_2_anchors.0].position,
        &brep.vertices[split_info.face_2_anchors.1].position,
    )
    .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;

    // Phase 1: Split each adjacent face along the corrected trim curves.
    let split_1 = split_face_with_trim(
        brep,
        eg.analysis.face_1,
        &pcurve_1,
        &split_info.trim_3d_1,
        tolerance,
    )?;

    let split_2 = split_face_with_trim(
        brep,
        eg.analysis.face_2,
        &pcurve_2,
        &split_info.trim_3d_2,
        tolerance,
    )?;

    // Phase 2: Identify which split face is consumed (contains the fillet edge).
    let (_kept_1, consumed_1) =
        helpers::identify_consumed_face(brep, split_1.face_a, split_1.face_b, fillet_edge)
            .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;

    let (_kept_2, consumed_2) =
        helpers::identify_consumed_face(brep, split_2.face_a, split_2.face_b, fillet_edge)
            .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;

    // Phase 3: Merge consumed faces by removing the original fillet edge.
    helpers::merge_consumed_faces(brep, consumed_1, consumed_2, fillet_edge, &eg.surface)
        .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;

    Ok(())
}

/// Split a face along a trim curve, handling errors gracefully.
fn split_face_with_trim(
    brep: &mut BRep,
    face: FaceHandle,
    pcurve: &NurbsCurve,
    trim_3d: &NurbsCurve,
    tolerance: &Tolerance,
) -> Result<popkern_core::topology::SplitFaceResult, BlendFailure> {
    brep.split_face(face, pcurve, trim_3d, tolerance)
        .map_err(|e| BlendFailure::InvalidTopology {
            reason: format!("face split failed: {e}"),
        })
}
