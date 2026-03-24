//! Topology construction — apply fillet/chamfer geometry to the BRep.
//!
//! Clones the input BRep, splits edges at trim-curve endpoints, splits faces
//! along trim curves, inserts new fillet/chamfer faces with the computed surface
//! geometry, and validates the result.

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
/// 1. Splits the fillet edge at the trim-curve endpoints to create anchor
///    vertices.
/// 2. Splits each adjacent face along its trim curve.
/// 3. Inserts a new face between the two trim curves carrying the fillet
///    surface.
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
/// Performs the three-phase topology modification for one edge:
/// split edges, split faces, insert the fillet face.
fn apply_single_blend(
    brep: &mut BRep,
    eg: &EdgeGeometry,
    tolerance: &Tolerance,
) -> Result<(), BlendFailure> {
    // Phase A: split fillet edge at trim-curve endpoints to create anchor vertices.
    let anchors = helpers::split_edge_at_trim_endpoints(
        brep,
        eg.analysis.edge,
        &eg.trim_3d_1,
        &eg.trim_3d_2,
        tolerance,
    )
    .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;

    // Phase B: split each adjacent face along its trim curve.
    let split_1 = split_face_with_trim(
        brep,
        eg.analysis.face_1,
        &eg.pcurve_1,
        &eg.trim_3d_1,
        tolerance,
    );
    let split_2 = split_face_with_trim(
        brep,
        eg.analysis.face_2,
        &eg.pcurve_2,
        &eg.trim_3d_2,
        tolerance,
    );

    // Phase C: insert the fillet face between the two trim curves.
    if let (Ok(s1), Ok(s2)) = (&split_1, &split_2) {
        helpers::insert_fillet_face(
            brep,
            &eg.surface,
            &anchors,
            &eg.trim_3d_1,
            &eg.trim_3d_2,
            &eg.pcurve_1,
            &eg.pcurve_2,
            s1,
            s2,
            tolerance,
        )
        .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;
    }

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
