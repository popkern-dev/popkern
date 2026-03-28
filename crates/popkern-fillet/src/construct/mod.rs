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
    let mut edge_remap: HashMap<EdgeHandle, EdgeHandle> = HashMap::new();
    let mut face_remap: HashMap<FaceHandle, FaceHandle> = HashMap::new();

    for eg in edge_geometries {
        if let Err(reason) =
            apply_single_blend(&mut result, eg, tolerance, &mut edge_remap, &mut face_remap)
        {
            failures.push((eg.analysis.edge, reason));
        }
    }

    if !failures.is_empty() && failures.len() == edge_geometries.len() {
        return Err(failures);
    }

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

    // Correct surface assignments corrupted by handle remapping.
    let fillet_surfaces: Vec<_> = edge_geometries
        .iter()
        .map(|eg| eg.surface.clone())
        .collect();
    correct_face_surfaces(&mut result, brep, &fillet_surfaces);

    Ok(ApplyBlendResult {
        brep: result,
        warnings,
    })
}

/// Apply a single fillet/chamfer blend to the BRep.
fn apply_single_blend(
    brep: &mut BRep,
    eg: &EdgeGeometry,
    tolerance: &Tolerance,
    edge_remap: &mut HashMap<EdgeHandle, EdgeHandle>,
    face_remap: &mut HashMap<FaceHandle, FaceHandle>,
) -> Result<(), BlendFailure> {
    let fillet_edge = resolve(edge_remap, eg.analysis.edge);
    let face_1 = resolve(face_remap, eg.analysis.face_1);
    let face_2 = resolve(face_remap, eg.analysis.face_2);

    let radius = eg.effective_radius;
    if radius < tolerance.linear {
        return Err(BlendFailure::InvalidTopology {
            reason: "fillet radius is too small".into(),
        });
    }

    let split_info = helpers::split_perpendicular_edges_by_radius(
        brep,
        face_1,
        face_2,
        fillet_edge,
        radius,
        tolerance,
    )
    .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;

    for &(old, edge_a, _) in &split_info.edge_splits {
        edge_remap.insert(old, edge_a);
    }

    // Compute pcurves from anchor positions on current faces.
    let p1a = brep.vertices[split_info.face_1_anchors.0].position;
    let p1b = brep.vertices[split_info.face_1_anchors.1].position;
    let pc1 = helpers::pcurve_for_face(brep, face_1, &p1a, &p1b)
        .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;
    let t3d1 = helpers::line_nurbs(&p1a, &p1b);

    let p2a = brep.vertices[split_info.face_2_anchors.0].position;
    let p2b = brep.vertices[split_info.face_2_anchors.1].position;
    let pc2 = helpers::pcurve_for_face(brep, face_2, &p2a, &p2b)
        .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;
    let t3d2 = helpers::line_nurbs(&p2a, &p2b);

    let s1 = split_face_with_trim(brep, face_1, &pc1, &t3d1, tolerance)?;
    let s2 = split_face_with_trim(brep, face_2, &pc2, &t3d2, tolerance)?;

    let (kept_1, consumed_1) =
        helpers::identify_consumed_face(brep, s1.face_a, s1.face_b, fillet_edge)
            .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;
    let (kept_2, consumed_2) =
        helpers::identify_consumed_face(brep, s2.face_a, s2.face_b, fillet_edge)
            .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;

    if kept_1 != face_1 {
        face_remap.insert(face_1, kept_1);
    }
    if kept_2 != face_2 {
        face_remap.insert(face_2, kept_2);
    }

    let fillet_face =
        helpers::merge_consumed_faces(brep, consumed_1, consumed_2, fillet_edge, &eg.surface)
            .map_err(|e| BlendFailure::InvalidTopology { reason: e })?;

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
    pcurve_attach::attach_fillet_pcurves(brep, &params);

    Ok(())
}

/// Resolve a handle through a remap table with depth limit.
fn resolve<H: Copy + Eq + std::hash::Hash>(remap: &HashMap<H, H>, mut h: H) -> H {
    for _ in 0..64 {
        match remap.get(&h) {
            Some(&new) => h = new,
            None => return h,
        }
    }
    h
}

/// Correct face surface assignments after topology modifications.
///
/// Checks every face. If a face has no surface, or if its surface does
/// not match its boundary vertices (max vertex-to-surface distance
/// exceeds threshold), replaces it with the best matching surface from
/// the pool of original planar surfaces and fillet surfaces.
fn correct_face_surfaces(brep: &mut BRep, original: &BRep, fillet_surfaces: &[SurfaceGeometry]) {
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceOnFace};
    use popkern_core::geometry::uv_projection::project_point_to_uv;

    let mut candidates: Vec<SurfaceOnFace> = original
        .faces
        .iter()
        .filter_map(|(_, f)| f.surface.clone())
        .collect();
    for surf in fillet_surfaces {
        candidates.push(SurfaceOnFace {
            geometry: surf.clone(),
            orientation: true,
        });
    }

    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();

    for fh in face_handles {
        let positions = brep.face_outer_loop_positions(fh);
        if positions.is_empty() {
            continue;
        }

        // Skip faces that already have a non-planar surface. These
        // were explicitly assigned by merge_consumed_faces and should
        // not be replaced by the heuristic.
        if let Some(sof) = &brep.faces[fh].surface
            && !matches!(sof.geometry, SurfaceGeometry::Plane { .. })
        {
            continue;
        }

        // Check if current planar surface is geometrically correct.
        let current_ok = brep.faces[fh].surface.as_ref().is_some_and(|sof| {
            positions.iter().all(|p| {
                project_point_to_uv(&sof.geometry, p)
                    .map(|(u, v)| p.distance_to(&sof.geometry.point_at(u, v)))
                    .unwrap_or(f64::INFINITY)
                    < 0.01
            })
        });

        if current_ok {
            continue;
        }

        // Find the best matching surface. Prefer planes (tight match)
        // over fillet surfaces (looser match).
        let mut best: Option<(SurfaceOnFace, f64)> = None;
        for sof in &candidates {
            let max_dist: f64 = positions
                .iter()
                .map(|p| {
                    project_point_to_uv(&sof.geometry, p)
                        .map(|(u, v)| p.distance_to(&sof.geometry.point_at(u, v)))
                        .unwrap_or(f64::INFINITY)
                })
                .fold(0.0_f64, f64::max);

            let thresh = if matches!(sof.geometry, SurfaceGeometry::Plane { .. }) {
                0.01
            } else {
                0.15
            };
            if max_dist < thresh {
                match &best {
                    Some((_, bd)) if max_dist >= *bd => {}
                    _ => best = Some((sof.clone(), max_dist)),
                }
            }
        }

        if let Some((sof, _dist)) = best {
            brep.faces[fh].surface = Some(sof);
            continue;
        }

        // Last resort: fit a plane through face vertices. This handles
        // topologically malformed faces from handle remapping where
        // vertices span multiple original surfaces.
        if positions.len() >= 3 {
            let p0 = positions[0];
            let p1 = positions[1];
            let p2 = positions[2];
            let u = popkern_core::geometry::Vec3::new(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
            let v = popkern_core::geometry::Vec3::new(p2.x - p0.x, p2.y - p0.y, p2.z - p0.z);
            let n = u.cross(&v);
            let len = n.length();
            if len > 1e-12 {
                let normal = popkern_core::geometry::Vec3::new(n.x / len, n.y / len, n.z / len);
                let ul = u.length();
                let u_axis = popkern_core::geometry::Vec3::new(u.x / ul, u.y / ul, u.z / ul);
                let v_axis = normal.cross(&u_axis);
                brep.faces[fh].surface = Some(SurfaceOnFace {
                    geometry: SurfaceGeometry::Plane {
                        origin: p0,
                        normal,
                        u_axis,
                        v_axis,
                    },
                    orientation: true,
                });
            }
        }
    }
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
