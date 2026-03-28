//! Topology construction -- apply fillet/chamfer geometry to the BRep.
//!
//! Processes each edge sequentially: splits perpendicular edges, computes
//! pcurves from anchor positions on the current face, splits faces along
//! trim curves, merges consumed faces via kef, and assigns the fillet
//! surface. Edge and face remap tables track handle invalidation across
//! edges. A post-pass corrects any surface misassignments from remapping.

pub(crate) mod blend_insert;
pub(crate) mod blend_topo;
pub(crate) mod helpers;
pub(crate) mod merge;
pub(crate) mod pcurve_attach;
pub(crate) mod surface_fix;

use std::collections::HashMap;

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

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
    /// Effective fillet/chamfer radius for perpendicular edge splitting.
    pub effective_radius: f64,
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
/// For each edge: resolves handles via remap tables, splits perpendicular
/// edges, computes pcurves from anchor positions, splits faces, identifies
/// consumed faces, merges them, and attaches pcurves. After all edges,
/// handles vertex blends and corrects any misassigned surfaces.
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

    // ── Phase A: ALL perpendicular edge splits ─────────────────────
    // Perp splits only modify edges (inserting vertices), not faces.
    // All can be done on the original face handles.
    struct PerpSplitData {
        face_1_anchors: (VertexHandle, VertexHandle),
        face_2_anchors: (VertexHandle, VertexHandle),
    }
    let mut perp_data: Vec<Option<PerpSplitData>> = Vec::new();

    for eg in edge_geometries {
        match helpers::split_perpendicular_edges_by_radius(
            &mut result,
            eg.analysis.face_1,
            eg.analysis.face_2,
            eg.analysis.edge,
            eg.effective_radius,
            tolerance,
        ) {
            Ok(info) => perp_data.push(Some(PerpSplitData {
                face_1_anchors: info.face_1_anchors,
                face_2_anchors: info.face_2_anchors,
            })),
            Err(e) => {
                failures.push((
                    eg.analysis.edge,
                    BlendFailure::InvalidTopology { reason: e },
                ));
                perp_data.push(None);
            }
        }
    }

    // ── Phase B: ALL face splits + merges ───────────────────────────
    // Now faces can be split. Each split may invalidate handles for
    // subsequent edges, so we track remaps.
    let mut face_remap: HashMap<FaceHandle, FaceHandle> = HashMap::new();

    for (idx, eg) in edge_geometries.iter().enumerate() {
        let Some(ref pd) = perp_data[idx] else {
            continue;
        };
        let face_1 = surface_fix::resolve(&face_remap, eg.analysis.face_1);
        let face_2 = surface_fix::resolve(&face_remap, eg.analysis.face_2);

        // Clip anchors to the current face boundary (handles crossing trims).
        let a1 =
            surface_fix::clip_anchors_to_face(&mut result, face_1, pd.face_1_anchors, tolerance);
        let p1a = result.vertices[a1.0].position;
        let p1b = result.vertices[a1.1].position;
        let pc1 = match helpers::pcurve_for_face(&result, face_1, &p1a, &p1b) {
            Ok(pc) => pc,
            Err(e) => {
                failures.push((
                    eg.analysis.edge,
                    BlendFailure::InvalidTopology { reason: e },
                ));
                continue;
            }
        };
        let t3d1 = helpers::line_nurbs(&p1a, &p1b);

        let a2 =
            surface_fix::clip_anchors_to_face(&mut result, face_2, pd.face_2_anchors, tolerance);
        let p2a = result.vertices[a2.0].position;
        let p2b = result.vertices[a2.1].position;
        let pc2 = match helpers::pcurve_for_face(&result, face_2, &p2a, &p2b) {
            Ok(pc) => pc,
            Err(e) => {
                failures.push((
                    eg.analysis.edge,
                    BlendFailure::InvalidTopology { reason: e },
                ));
                continue;
            }
        };
        let t3d2 = helpers::line_nurbs(&p2a, &p2b);

        let s1 = match split_face_with_trim(&mut result, face_1, &pc1, &t3d1, tolerance) {
            Ok(s) => s,
            Err(e) => {
                failures.push((eg.analysis.edge, e));
                continue;
            }
        };
        let s2 = match split_face_with_trim(&mut result, face_2, &pc2, &t3d2, tolerance) {
            Ok(s) => s,
            Err(e) => {
                failures.push((eg.analysis.edge, e));
                continue;
            }
        };

        let (kept_1, consumed_1) = match helpers::identify_consumed_face(
            &result,
            s1.face_a,
            s1.face_b,
            eg.analysis.edge,
        ) {
            Ok(kc) => kc,
            Err(e) => {
                failures.push((
                    eg.analysis.edge,
                    BlendFailure::InvalidTopology { reason: e },
                ));
                continue;
            }
        };
        let (kept_2, consumed_2) = match helpers::identify_consumed_face(
            &result,
            s2.face_a,
            s2.face_b,
            eg.analysis.edge,
        ) {
            Ok(kc) => kc,
            Err(e) => {
                failures.push((
                    eg.analysis.edge,
                    BlendFailure::InvalidTopology { reason: e },
                ));
                continue;
            }
        };

        if kept_1 != face_1 {
            face_remap.insert(face_1, kept_1);
        }
        if kept_2 != face_2 {
            face_remap.insert(face_2, kept_2);
        }

        match helpers::merge_consumed_faces(
            &mut result,
            consumed_1,
            consumed_2,
            eg.analysis.edge,
            &eg.surface,
        ) {
            Ok(fillet_face) => {
                let params = pcurve_attach::FilletPcurveParams {
                    fillet_face,
                    fillet_surface: &eg.surface,
                    trim_edges_1: &s1.new_edges,
                    trim_edges_2: &s2.new_edges,
                    eg_pcurve_1: &pc1,
                    eg_pcurve_2: &pc2,
                    trim_3d_1: &t3d1,
                    trim_3d_2: &t3d2,
                };
                pcurve_attach::attach_fillet_pcurves(&mut result, &params);
            }
            Err(e) => {
                failures.push((
                    eg.analysis.edge,
                    BlendFailure::InvalidTopology { reason: e },
                ));
            }
        }
    }

    if !failures.is_empty() && failures.len() == edge_geometries.len() {
        return Err(failures);
    }

    // Save non-planar surfaces keyed by face centroid position so we can
    // restore them after vertex blend (which destroys face handles via kef).
    use popkern_core::geometry::surface::SurfaceOnFace;
    let saved_fillet_surfaces: Vec<(popkern_core::geometry::Point3, SurfaceOnFace)> = result
        .faces
        .iter()
        .filter_map(|(fh, f)| {
            let sof = f.surface.as_ref()?;
            if matches!(sof.geometry, SurfaceGeometry::Plane { .. }) {
                return None;
            }
            let positions = result.face_outer_loop_positions(fh);
            if positions.is_empty() {
                return None;
            }
            let n = positions.len() as f64;
            let cx = positions.iter().map(|p| p.x).sum::<f64>() / n;
            let cy = positions.iter().map(|p| p.y).sum::<f64>() / n;
            let cz = positions.iter().map(|p| p.z).sum::<f64>() / n;
            Some((popkern_core::geometry::Point3::new(cx, cy, cz), sof.clone()))
        })
        .collect();

    // Vertex blend handling.
    if *vertex_mode == VertexMode::Smooth && edge_geometries.len() > 1 {
        let mut vertex_edges: HashMap<VertexHandle, Vec<usize>> = HashMap::new();
        for (idx, eg) in edge_geometries.iter().enumerate() {
            vertex_edges
                .entry(eg.analysis.vertex_start)
                .or_default()
                .push(idx);
            vertex_edges
                .entry(eg.analysis.vertex_end)
                .or_default()
                .push(idx);
        }
        for (vertex, edge_indices) in &vertex_edges {
            if edge_indices.len() >= 2 {
                match blend_insert::insert_vertex_blend(
                    &mut result,
                    *vertex,
                    edge_indices,
                    edge_geometries,
                    tolerance,
                ) {
                    Ok(()) => {}
                    Err(_) => {
                        let edge_idx = edge_indices[0];
                        if edge_idx < edge_geometries.len() {
                            warnings.push((
                                edge_geometries[edge_idx].analysis.edge,
                                BlendWarning::VertexModeFallback {
                                    requested: VertexMode::Smooth,
                                    actual: VertexMode::Mitred,
                                    vertex: Some(*vertex),
                                },
                            ));
                        }
                    }
                }
            }
        }
    }

    // Restore fillet surfaces destroyed by vertex blend's kef operations.
    // Match by face centroid proximity since handles change during kef.
    for (saved_centroid, saved_sof) in &saved_fillet_surfaces {
        for fh in result.faces.keys().collect::<Vec<_>>() {
            // Only restore to faces without a non-planar surface.
            if result.faces[fh]
                .surface
                .as_ref()
                .is_some_and(|s| !matches!(s.geometry, SurfaceGeometry::Plane { .. }))
            {
                continue;
            }
            let positions = result.face_outer_loop_positions(fh);
            if positions.is_empty() {
                continue;
            }
            let n = positions.len() as f64;
            let cx = positions.iter().map(|p| p.x).sum::<f64>() / n;
            let cy = positions.iter().map(|p| p.y).sum::<f64>() / n;
            let cz = positions.iter().map(|p| p.z).sum::<f64>() / n;
            let dist = ((cx - saved_centroid.x).powi(2)
                + (cy - saved_centroid.y).powi(2)
                + (cz - saved_centroid.z).powi(2))
            .sqrt();
            if dist < 0.2 {
                result.faces[fh].surface = Some(saved_sof.clone());
                break;
            }
        }
    }

    // Correct surface assignments corrupted by handle remapping.
    let fillet_surfaces: Vec<_> = edge_geometries
        .iter()
        .map(|eg| eg.surface.clone())
        .collect();
    surface_fix::correct_face_surfaces(&mut result, brep, &fillet_surfaces);

    Ok(ApplyBlendResult {
        brep: result,
        warnings,
    })
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
