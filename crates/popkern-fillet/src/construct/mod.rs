//! Topology construction -- apply fillet/chamfer geometry to the BRep.
//!
//! Clones the input BRep, splits perpendicular edges at the fillet offset
//! distance, splits faces along trim curves, merges consumed face regions
//! by removing the original fillet edge, and assigns the fillet/chamfer
//! surface geometry to the resulting face.

pub(crate) mod blend_insert;
pub(crate) mod blend_topo;
pub(crate) mod helpers;
pub(crate) mod merge;
pub(crate) mod pcurve_attach;

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
/// sharp corners (no extra face needed). For [`VertexMode::Smooth`],
/// inserts a G1-continuous blend face at each vertex where 2+ fillets
/// meet, falling back to mitred with a warning if insertion fails.
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

    // Vertex blend handling.
    //
    // Mitred mode requires no extra topology — fillet faces already
    // terminate at the original vertex as sharp corners.
    //
    // Smooth mode inserts a G1-continuous blend face at each vertex
    // where 2+ fillets meet, by carving corner triangles from surrounding
    // faces and merging them into one blend face.
    if *vertex_mode == VertexMode::Smooth && edge_geometries.len() > 1 {
        // Build vertex -> edge-index map to find multi-fillet vertices.
        let mut vertex_edges: HashMap<popkern_core::topology::VertexHandle, Vec<usize>> =
            HashMap::new();
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
                    Ok(()) => { /* smooth blend inserted successfully */ }
                    Err(_reason) => {
                        // Fall back to mitred with a warning.
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

    // Safety net: propagate surface geometry to any face that lost it
    // during topology modifications. A face that shares edges with a
    // face that has a surface should inherit that surface (unless it's a
    // fillet surface — in which case it was already assigned).
    propagate_missing_surfaces(&mut result);

    Ok(ApplyBlendResult {
        brep: result,
        warnings,
    })
}

/// Propagate surface geometry to faces that lost it during topology
/// modifications.
///
/// For each face without a surface, finds a neighboring face (sharing an
/// edge via half-edge twin) that has a surface with the same type as the
/// original (planar, cylindrical, etc.) and copies it. This handles cases
/// where face splits or `kef` operations leave a face without geometry.
fn propagate_missing_surfaces(brep: &mut BRep) {
    // Iterate until no more surfaces can be propagated (handles chains
    // of surfaceless faces that are only reachable transitively).
    let max_rounds = brep.faces.len();
    for _ in 0..max_rounds {
        let faces_without_surface: Vec<FaceHandle> = brep
            .faces
            .iter()
            .filter(|(_, f)| f.surface.is_none())
            .map(|(fh, _)| fh)
            .collect();

        if faces_without_surface.is_empty() {
            break;
        }

        let mut any_fixed = false;

        for fh in faces_without_surface {
            let outer_loop = brep.faces[fh].outer_loop;
            let start_he = match brep.loops.get(outer_loop) {
                Some(l) => l.half_edge,
                None => continue,
            };
            if !brep.half_edges.contains_key(start_he) {
                continue;
            }

            let mut cursor = start_he;
            let limit = brep.half_edges.len() + 1;
            let mut found_surface = None;

            for _ in 0..limit {
                let he = &brep.half_edges[cursor];
                let twin = he.twin;
                if brep.half_edges.contains_key(twin) {
                    let twin_loop = brep.half_edges[twin].loop_;
                    if brep.loops.contains_key(twin_loop) {
                        let neighbor_face = brep.loops[twin_loop].face;
                        if let Some(ref sof) = brep.faces[neighbor_face].surface {
                            found_surface = Some(sof.clone());
                            break;
                        }
                    }
                }
                cursor = he.next;
                if cursor == start_he {
                    break;
                }
            }

            // Also check inner loops for neighbor surfaces.
            if found_surface.is_none() {
                let inner_loops = brep.faces[fh].inner_loops.clone();
                for inner_lh in &inner_loops {
                    let il_start = match brep.loops.get(*inner_lh) {
                        Some(l) => l.half_edge,
                        None => continue,
                    };
                    if !brep.half_edges.contains_key(il_start) {
                        continue;
                    }
                    let mut ic = il_start;
                    for _ in 0..limit {
                        let he = &brep.half_edges[ic];
                        let twin = he.twin;
                        if brep.half_edges.contains_key(twin) {
                            let tl = brep.half_edges[twin].loop_;
                            if brep.loops.contains_key(tl) {
                                let nf = brep.loops[tl].face;
                                if let Some(ref sof) = brep.faces[nf].surface {
                                    found_surface = Some(sof.clone());
                                    break;
                                }
                            }
                        }
                        ic = he.next;
                        if ic == il_start {
                            break;
                        }
                    }
                    if found_surface.is_some() {
                        break;
                    }
                }
            }

            if let Some(surface) = found_surface {
                brep.faces[fh].surface = Some(surface);
                any_fixed = true;
            }
        }

        if !any_fixed {
            break;
        }
    }
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
    let fillet_face =
        helpers::merge_consumed_faces(brep, consumed_1, consumed_2, fillet_edge, &eg.surface)
            .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;

    // Phase 4: Attach pcurves to the fillet face and kept-face trim edges.
    let pcurve_params = pcurve_attach::FilletPcurveParams {
        fillet_face,
        fillet_surface: &eg.surface,
        trim_edges_1: &split_1.new_edges,
        trim_edges_2: &split_2.new_edges,
        eg_pcurve_1: &eg.pcurve_1,
        eg_pcurve_2: &eg.pcurve_2,
        trim_3d_1: &eg.trim_3d_1,
        trim_3d_2: &eg.trim_3d_2,
    };
    pcurve_attach::attach_fillet_pcurves(brep, &pcurve_params);

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
