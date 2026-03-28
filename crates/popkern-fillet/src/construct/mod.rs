//! Topology construction -- apply fillet/chamfer geometry to the BRep.
//!
//! Three-phase single-pass approach:
//!   Phase A: split ALL perpendicular edges (no face changes)
//!   Phase B: split ALL faces by trim curves (tracking sub-faces)
//!   Phase C: merge ALL consumed pairs via kef (assign fillet surfaces)
//!
//! A post-pass corrects any residual surface misassignments.

pub(crate) mod blend_insert;
pub(crate) mod blend_topo;
pub(crate) mod helpers;
pub(crate) mod merge;
pub(crate) mod pcurve_attach;
mod phase_split;
pub(crate) mod surface_fix;

use std::collections::HashMap;

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::analysis::EdgeAnalysis;
use crate::error::{BlendFailure, BlendWarning};
use crate::types::VertexMode;
use phase_split::{FaceTrim, TrimSplitResult, find_face_containing_vertices, split_face_with_trim};

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
/// Three-phase topology construction:
///   Phase A: split ALL perpendicular edges (no face changes)
///   Phase B: split ALL faces by trim curves (tracking sub-faces)
///   Phase C: merge ALL consumed pairs via kef (assign fillet surfaces)
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

    // ── Collect trim data from Phase A results ────────────────────
    // Use BTreeMap for deterministic face iteration order.
    let mut face_trims: std::collections::BTreeMap<FaceHandle, Vec<FaceTrim>> =
        std::collections::BTreeMap::new();

    for (idx, eg) in edge_geometries.iter().enumerate() {
        let Some(ref pd) = perp_data[idx] else {
            continue;
        };

        face_trims
            .entry(eg.analysis.face_1)
            .or_default()
            .push(FaceTrim {
                edge_idx: idx,
                side: 1,
                anchors: pd.face_1_anchors,
            });
        face_trims
            .entry(eg.analysis.face_2)
            .or_default()
            .push(FaceTrim {
                edge_idx: idx,
                side: 2,
                anchors: pd.face_2_anchors,
            });
    }

    // Sort trims within each face by edge_idx for deterministic ordering.
    for trims in face_trims.values_mut() {
        trims.sort_by_key(|t| (t.edge_idx, t.side));
    }

    // ── Phase B: Split ALL faces by their trim curves ─────────────
    let mut trim_results: Vec<TrimSplitResult> = Vec::new();

    for (original_face, trims) in &face_trims {
        // Track ALL sub-faces of this original face (both kept and
        // consumed). Consumed faces must remain searchable because
        // later trims may have anchor vertices on consumed sub-faces
        // (e.g., at cube corners where multiple fillet edges share
        // perpendicular edges). Phase C knows which faces are consumed
        // from TrimSplitResult.
        let mut sub_faces = vec![*original_face];

        for trim in trims {
            let target =
                find_face_containing_vertices(&result, &sub_faces, trim.anchors.0, trim.anchors.1);

            let Some(target) = target else {
                failures.push((
                    edge_geometries[trim.edge_idx].analysis.edge,
                    BlendFailure::InvalidTopology {
                        reason: format!(
                            "anchor vertices not found on any sub-face of {original_face:?}"
                        ),
                    },
                ));
                continue;
            };

            // Compute pcurve and 3D trim from actual vertex positions.
            // Done here (before kef) so the face surface identity is unambiguous.
            let p_a = result.vertices[trim.anchors.0].position;
            let p_b = result.vertices[trim.anchors.1].position;
            let pc = match helpers::pcurve_for_face(&result, target, &p_a, &p_b) {
                Ok(pc) => pc,
                Err(e) => {
                    failures.push((
                        edge_geometries[trim.edge_idx].analysis.edge,
                        BlendFailure::InvalidTopology { reason: e },
                    ));
                    continue;
                }
            };
            let t3d = helpers::line_nurbs(&p_a, &p_b);

            let split = match split_face_with_trim(&mut result, target, &pc, &t3d, tolerance) {
                Ok(s) => s,
                Err(e) => {
                    failures.push((edge_geometries[trim.edge_idx].analysis.edge, e));
                    continue;
                }
            };

            let (kept, consumed) = match helpers::identify_consumed_face(
                &result,
                split.face_a,
                split.face_b,
                edge_geometries[trim.edge_idx].analysis.edge,
            ) {
                Ok(kc) => kc,
                Err(e) => {
                    failures.push((
                        edge_geometries[trim.edge_idx].analysis.edge,
                        BlendFailure::InvalidTopology { reason: e },
                    ));
                    // Cannot determine which is consumed. Replace target
                    // with both halves so later trims can still find anchors.
                    if let Some(pos) = sub_faces.iter().position(|&f| f == target) {
                        sub_faces[pos] = split.face_a;
                        sub_faces.push(split.face_b);
                    }
                    continue;
                }
            };

            trim_results.push(TrimSplitResult {
                edge_idx: trim.edge_idx,
                side: trim.side,
                consumed,
                trim_edges: split.new_edges.clone(),
                pcurve: pc,
                trim_3d: t3d,
            });

            // Replace target with both sub-faces so later trims can
            // find anchor vertices on either half.
            if let Some(pos) = sub_faces.iter().position(|&f| f == target) {
                sub_faces[pos] = kept;
                sub_faces.push(consumed);
            }
        }
    }

    // ── Re-identify consumed faces ──────────────────────────────────
    // A consumed face from an earlier trim may have been split by a
    // later trim on the same original face. The stored `consumed` handle
    // may now refer to only part of the original strip. Re-identify by
    // finding which current face contains the fillet edge's half-edge.
    for tr in &mut trim_results {
        let fillet_edge = edge_geometries[tr.edge_idx].analysis.edge;
        if !result.faces.contains_key(tr.consumed) {
            continue; // Handle was recycled — skip, Phase C will fail gracefully.
        }
        // Check if the stored consumed face still contains the fillet edge.
        if helpers::identify_consumed_face_contains(&result, tr.consumed, fillet_edge) {
            continue; // Still correct.
        }
        // Search all faces for the one containing this fillet edge.
        for fh in result.faces.keys() {
            if helpers::identify_consumed_face_contains(&result, fh, fillet_edge) {
                tr.consumed = fh;
                break;
            }
        }
    }

    // ── Phase C: Merge consumed face pairs via kef ────────────────
    // Group by edge index: need side 1 + side 2 to merge.
    let mut edge_trims: std::collections::BTreeMap<usize, Vec<&TrimSplitResult>> =
        std::collections::BTreeMap::new();
    for tr in &trim_results {
        edge_trims.entry(tr.edge_idx).or_default().push(tr);
    }

    // Track which consumed faces were successfully merged for cleanup.
    let mut merged_consumed: std::collections::HashSet<FaceHandle> =
        std::collections::HashSet::new();

    for (idx, trims) in &edge_trims {
        let eg = &edge_geometries[*idx];
        let fillet_edge = eg.analysis.edge;

        let side_1 = trims.iter().find(|t| t.side == 1);
        let side_2 = trims.iter().find(|t| t.side == 2);

        let (Some(s1), Some(s2)) = (side_1, side_2) else {
            // Missing one side — consumed face(s) cleaned up below.
            continue;
        };

        match helpers::merge_consumed_faces(
            &mut result,
            s1.consumed,
            s2.consumed,
            fillet_edge,
            &eg.surface,
        ) {
            Ok(fillet_face) => {
                merged_consumed.insert(s1.consumed);
                merged_consumed.insert(s2.consumed);

                // Use pcurves computed during Phase B (before kef).
                let params = pcurve_attach::FilletPcurveParams {
                    fillet_face,
                    fillet_surface: &eg.surface,
                    trim_edges_1: &s1.trim_edges,
                    trim_edges_2: &s2.trim_edges,
                    eg_pcurve_1: &s1.pcurve,
                    eg_pcurve_2: &s2.pcurve,
                    trim_3d_1: &s1.trim_3d,
                    trim_3d_2: &s2.trim_3d,
                };
                pcurve_attach::attach_fillet_pcurves(&mut result, &params);
            }
            Err(e) => {
                failures.push((fillet_edge, BlendFailure::InvalidTopology { reason: e }));
            }
        }
    }

    // ── Phase C cleanup: undo orphaned consumed faces ─────────────
    // If a fillet edge only had one side split, its consumed face is
    // dangling. Undo the split by kef-ing the trim edge.
    for tr in &trim_results {
        if merged_consumed.contains(&tr.consumed) {
            continue;
        }
        if !result.faces.contains_key(tr.consumed) {
            continue;
        }
        for &trim_edge in &tr.trim_edges {
            if !result.edges.contains_key(trim_edge) {
                continue;
            }
            let he = result.edges[trim_edge].half_edge;
            let _ = result.kef(he); // Best-effort cleanup.
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
    // Only restore if the saved surface geometrically fits the face's vertices
    // (all vertices within 0.05 of the surface). This avoids the cross-matching
    // bug where centroid proximity caused wrong torus surfaces to be assigned
    // to planar faces.
    for (_saved_centroid, saved_sof) in &saved_fillet_surfaces {
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
            // Verify the saved surface actually fits this face's geometry.
            let max_dist: f64 = positions
                .iter()
                .filter_map(|p| {
                    popkern_core::geometry::uv_projection::project_point_to_uv(
                        &saved_sof.geometry,
                        p,
                    )
                    .map(|(u, v)| {
                        use popkern_core::geometry::surface::SurfaceEval;
                        p.distance_to(&saved_sof.geometry.point_at(u, v))
                    })
                })
                .fold(0.0_f64, f64::max);
            if max_dist < 0.05 {
                result.faces[fh].surface = Some(saved_sof.clone());
                break;
            }
        }
    }

    // Correct surface assignments — safety net.
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
