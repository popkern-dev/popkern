//! Topology construction -- apply fillet/chamfer geometry to the BRep.
//!
//! Clones the input BRep, splits perpendicular edges at the fillet offset
//! distance, splits faces along trim curves, merges consumed face regions
//! by removing the original fillet edge, and assigns the fillet/chamfer
//! surface geometry to the resulting face.

pub(crate) mod helpers;
pub(crate) mod merge;

use std::collections::HashMap;

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::analysis::EdgeAnalysis;
use crate::error::{BlendFailure, BlendWarning};
use crate::types::VertexMode;

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

/// Tracks handle remapping after topology modifications.
///
/// When splitting edges, old edge handles become invalid and map to new
/// handles. This struct records those changes so subsequent operations
/// can look up the correct (remapped) handles.
#[derive(Debug, Default)]
struct HandleRemap {
    /// Maps an old edge handle to the two edges produced by splitting.
    edge_splits: HashMap<EdgeHandle, (EdgeHandle, EdgeHandle)>,
    /// Maps an old face handle to its replacement after splitting.
    face_remap: HashMap<FaceHandle, FaceHandle>,
}

impl HandleRemap {
    /// Maximum chain depth for resolving remapped handles.
    const MAX_DEPTH: usize = 64;

    /// Resolve an edge handle through the remap table.
    ///
    /// If the edge was split, returns the first child (edge_a). Follows
    /// the chain iteratively until a non-remapped handle is found, with
    /// a depth limit to prevent infinite loops.
    fn resolve_edge(&self, mut edge: EdgeHandle) -> EdgeHandle {
        for _ in 0..Self::MAX_DEPTH {
            if let Some(&(a, _b)) = self.edge_splits.get(&edge) {
                edge = a;
            } else {
                return edge;
            }
        }
        edge
    }

    /// Resolve a face handle through the remap table.
    fn resolve_face(&self, mut face: FaceHandle) -> FaceHandle {
        for _ in 0..Self::MAX_DEPTH {
            if let Some(&new_face) = self.face_remap.get(&face) {
                face = new_face;
            } else {
                return face;
            }
        }
        face
    }
}

/// Result of applying blend geometry, including warnings.
pub struct ApplyBlendResult {
    /// The modified BRep.
    pub brep: BRep,
    /// Warnings generated during vertex mode handling.
    pub warnings: Vec<(EdgeHandle, BlendWarning)>,
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
/// After all edges are processed, handles vertex blends according to
/// `vertex_mode`. For [`VertexMode::Mitred`], fillet faces terminate at
/// sharp corners (no extra face needed). For [`VertexMode::Smooth`], falls
/// back to mitred with a warning (smooth blending is planned for a future
/// release).
///
/// Tracks handle remapping across edges so that multi-edge fillets
/// correctly resolve handles invalidated by previous operations.
///
/// Returns the modified BRep or a list of per-edge failures.
pub fn apply_blend_geometry(
    brep: &BRep,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
    vertex_mode: &VertexMode,
) -> Result<ApplyBlendResult, Vec<(EdgeHandle, BlendFailure)>> {
    let mut result = brep.clone();
    let mut failures: Vec<(EdgeHandle, BlendFailure)> = Vec::new();
    let mut warnings: Vec<(EdgeHandle, BlendWarning)> = Vec::new();
    let mut remap = HandleRemap::default();

    for eg in edge_geometries {
        if let Err(reason) = apply_single_blend(&mut result, eg, tolerance, &mut remap) {
            failures.push((eg.analysis.edge, reason));
        }
    }

    if !failures.is_empty() && failures.len() == edge_geometries.len() {
        return Err(failures);
    }

    // Vertex blend handling: when Smooth mode is requested but not yet
    // implemented, fall back to Mitred (sharp corners) with a warning.
    // Mitred mode requires no extra topology — fillet faces already
    // terminate at the original vertex as sharp corners.
    if *vertex_mode == VertexMode::Smooth && edge_geometries.len() > 1 {
        // Smooth vertex blending not yet implemented; degrade to mitred.
        // Record a warning on the first edge so callers know.
        if let Some(eg) = edge_geometries.first() {
            warnings.push((
                eg.analysis.edge,
                BlendWarning::VertexModeFallback {
                    requested: VertexMode::Smooth,
                    actual: VertexMode::Mitred,
                },
            ));
        }
    }

    Ok(ApplyBlendResult {
        brep: result,
        warnings,
    })
}

/// Apply a single fillet/chamfer blend to the BRep.
///
/// Performs the four-phase topology modification for one edge:
/// 1. Split perpendicular edges at the offset distance (creates anchor vertices).
/// 2. Split faces along trim curves connecting the anchor vertices.
/// 3. Identify consumed faces (the portions that contain the fillet edge).
/// 4. Merge consumed faces via `kef` and assign fillet surface.
///
/// Updates `remap` with any handle changes caused by edge/face splits so
/// that subsequent blend operations can resolve invalidated handles.
fn apply_single_blend(
    brep: &mut BRep,
    eg: &EdgeGeometry,
    tolerance: &Tolerance,
    remap: &mut HandleRemap,
) -> Result<(), BlendFailure> {
    // Resolve handles through the remap table in case a prior blend
    // invalidated our edge or face handles.
    let fillet_edge = remap.resolve_edge(eg.analysis.edge);
    let face_1 = remap.resolve_face(eg.analysis.face_1);
    let face_2 = remap.resolve_face(eg.analysis.face_2);

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
        face_1,
        face_2,
        fillet_edge,
        radius,
        tolerance,
    )
    .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;

    // Record edge splits from Phase 0 for future remapping.
    for &(old_edge, edge_a, edge_b) in &split_info.edge_splits {
        remap.edge_splits.insert(old_edge, (edge_a, edge_b));
    }

    // Compute pcurves for the corrected trim curves on each face.
    let pcurve_1 = helpers::pcurve_for_face(
        brep,
        face_1,
        &brep.vertices[split_info.face_1_anchors.0].position,
        &brep.vertices[split_info.face_1_anchors.1].position,
    )
    .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;

    let pcurve_2 = helpers::pcurve_for_face(
        brep,
        face_2,
        &brep.vertices[split_info.face_2_anchors.0].position,
        &brep.vertices[split_info.face_2_anchors.1].position,
    )
    .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;

    // Phase 1: Split each adjacent face along the corrected trim curves.
    let split_1 = split_face_with_trim(brep, face_1, &pcurve_1, &split_info.trim_3d_1, tolerance)?;

    let split_2 = split_face_with_trim(brep, face_2, &pcurve_2, &split_info.trim_3d_2, tolerance)?;

    // Phase 2: Identify which split face is consumed (contains the fillet edge).
    let (kept_1, consumed_1) =
        helpers::identify_consumed_face(brep, split_1.face_a, split_1.face_b, fillet_edge)
            .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;

    let (kept_2, consumed_2) =
        helpers::identify_consumed_face(brep, split_2.face_a, split_2.face_b, fillet_edge)
            .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;

    // Record face remaps so future blends can find the surviving face.
    // Only remap if the kept face differs from the original.
    if kept_1 != face_1 {
        remap.face_remap.insert(face_1, kept_1);
    }
    if kept_2 != face_2 {
        remap.face_remap.insert(face_2, kept_2);
    }

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
