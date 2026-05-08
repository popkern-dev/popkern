//! Apex edge resolution for the legacy A2 mutating builder.

use popkern_core::geometry::{CurveGeometry, CurveOnEdge, NurbsCurve};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::construct::blend_topo::find_halfedge_from_vertex_in_face;
use crate::fillet_surface::FilletSurfaceResult;

use super::construct_vertices::vertex_loop_face_handles;
use super::error::ApexBuildError;
use super::plan::ApexPlan;
use super::target_topology::TargetEdge;

pub(super) fn resolve_apex_edges(
    brep: &mut BRep,
    plan: &ApexPlan,
    apex_vertex_handles: &[VertexHandle],
    _fillet_surfaces: &[FilletSurfaceResult],
) -> Result<Vec<Option<EdgeHandle>>, ApexBuildError> {
    use super::target_topology::TargetEdgeRole;
    let mut out: Vec<Option<EdgeHandle>> =
        Vec::with_capacity(plan.target_topology.apex_edges.len());
    let mut cap_side_faces: std::collections::HashSet<FaceHandle> =
        std::collections::HashSet::new();
    for (i, target_edge) in plan.target_topology.apex_edges.iter().enumerate() {
        let v_origin = apex_vertex_handles
            .get(target_edge.origin_vertex_idx)
            .copied()
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "apex_edge[{i}] references origin_vertex_idx={} but plan has only {} apex_vertices",
                    target_edge.origin_vertex_idx,
                    apex_vertex_handles.len(),
                ),
            })?;
        let v_dest = apex_vertex_handles
            .get(target_edge.dest_vertex_idx)
            .copied()
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "apex_edge[{i}] references dest_vertex_idx={} but plan has only {} apex_vertices",
                    target_edge.dest_vertex_idx,
                    apex_vertex_handles.len(),
                ),
            })?;
        if matches!(target_edge.role, TargetEdgeRole::CapVsCyl { .. }) {
            if crate::a2_trace::enabled() {
                eprintln!(
                    "[a2-trace] apex_edge[{i}] (vidx {} -> {}) -> SKIPPED (CapVsCyl, \
                     §28.D)",
                    target_edge.origin_vertex_idx, target_edge.dest_vertex_idx,
                );
            }
            out.push(None);
            continue;
        }
        if let Some(eh) = find_edge_between(brep, v_origin, v_dest) {
            if crate::a2_trace::enabled() {
                eprintln!(
                    "[a2-trace] apex_edge[{i}] (vidx {} -> {}) -> FOUND existing edge \
                     {eh:?} between {v_origin:?} and {v_dest:?}",
                    target_edge.origin_vertex_idx, target_edge.dest_vertex_idx,
                );
            }
            out.push(Some(eh));
            continue;
        }
        if crate::a2_trace::enabled() {
            let origin_faces = vertex_loop_face_handles(brep, v_origin);
            let dest_faces = vertex_loop_face_handles(brep, v_dest);
            let shared: Vec<_> = origin_faces
                .iter()
                .filter(|f| dest_faces.contains(f))
                .copied()
                .collect();
            eprintln!(
                "[a2-trace] apex_edge[{i}] (vidx {} -> {}) NEED to create; \
                 v_origin={v_origin:?} on faces {origin_faces:?}; \
                 v_dest={v_dest:?} on faces {dest_faces:?}; \
                 SHARED faces = {shared:?}",
                target_edge.origin_vertex_idx, target_edge.dest_vertex_idx,
            );
        }
        let v_incident_faces = if matches!(target_edge.role, TargetEdgeRole::CapVsPlane { .. }) {
            faces_with_vertex_on_outer_loop(brep, plan.apex_vertex)
        } else {
            std::collections::HashSet::new()
        };
        let prefer_set: std::collections::HashSet<FaceHandle> = if !v_incident_faces.is_empty() {
            v_incident_faces
        } else {
            cap_side_faces.clone()
        };
        let (edge, mef_face_pair) =
            create_apex_edge_preferring(brep, v_origin, v_dest, target_edge, Some(&prefer_set))
                .map_err(|reason| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "apex_edge[{i}] (vertex_idx {} -> {}) create via mef failed: {reason}",
                        target_edge.origin_vertex_idx, target_edge.dest_vertex_idx,
                    ),
                })?;
        cap_side_faces.insert(mef_face_pair.0);
        cap_side_faces.insert(mef_face_pair.1);
        out.push(Some(edge));
    }
    Ok(out)
}

#[allow(dead_code)] // §28.D superseded; retained for diagnostic use.
fn adopt_fillet_end_arc(
    brep: &BRep,
    fillet_idx: usize,
    v_origin: VertexHandle,
    v_dest: VertexHandle,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Option<EdgeHandle> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let fillet = fillet_surfaces.get(fillet_idx)?;
    let SurfaceGeometry::Cylinder {
        origin: f_origin,
        axis: f_axis,
        radius: f_radius,
    } = &fillet.surface
    else {
        return None;
    };
    let f_axis_unit = f_axis.normalized()?;
    for (eh, edge) in &brep.edges {
        let primary = edge.half_edge;
        let he = brep.half_edges.get(primary)?;
        let twin_he = brep.half_edges.get(he.twin)?;
        let endpoints = (he.origin, twin_he.origin);
        if endpoints != (v_origin, v_dest) && endpoints != (v_dest, v_origin) {
            continue;
        }
        let face_a = brep.loops.get(he.loop_).map(|l| l.face);
        let face_b = brep.loops.get(twin_he.loop_).map(|l| l.face);
        for fh_opt in [face_a, face_b] {
            let Some(fh) = fh_opt else { continue };
            let Some(face) = brep.faces.get(fh) else {
                continue;
            };
            let Some(sof) = face.surface.as_ref() else {
                continue;
            };
            let SurfaceGeometry::Cylinder {
                origin: face_origin,
                axis: face_axis,
                radius: face_radius,
            } = &sof.geometry
            else {
                continue;
            };
            if (face_radius - f_radius).abs() > f_radius.abs() * 1e-4 + 1e-9 {
                continue;
            }
            let Some(face_axis_unit) = face_axis.normalized() else {
                continue;
            };
            if face_axis_unit.dot(&f_axis_unit).abs() < 1.0 - 1e-3 {
                continue;
            }
            let d = popkern_core::geometry::Vec3::new(
                face_origin.x - f_origin.x,
                face_origin.y - f_origin.y,
                face_origin.z - f_origin.z,
            );
            let along = d.dot(&f_axis_unit);
            let perp = popkern_core::geometry::Vec3::new(
                d.x - along * f_axis_unit.x,
                d.y - along * f_axis_unit.y,
                d.z - along * f_axis_unit.z,
            );
            if perp.length() > 1e-3 {
                continue;
            }
            return Some(eh);
        }
    }
    None
}

fn faces_with_vertex_on_outer_loop(
    brep: &BRep,
    vertex: VertexHandle,
) -> std::collections::HashSet<FaceHandle> {
    let mut out: std::collections::HashSet<FaceHandle> = std::collections::HashSet::new();
    for (fh, face) in &brep.faces {
        let Some(loop_data) = brep.loops.get(face.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut iters = 0;
        let max_iters = brep.half_edges.len() + 1;
        while let Some(he) = brep.half_edges.get(cursor) {
            if he.origin == vertex {
                out.insert(fh);
                break;
            }
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
    }
    out
}

fn find_edge_between(
    brep: &BRep,
    v_origin: VertexHandle,
    v_dest: VertexHandle,
) -> Option<EdgeHandle> {
    for (eh, edge) in &brep.edges {
        let primary = edge.half_edge;
        let Some(primary_he) = brep.half_edges.get(primary) else {
            continue;
        };
        let twin = primary_he.twin;
        let Some(twin_he) = brep.half_edges.get(twin) else {
            continue;
        };
        let endpoints = (primary_he.origin, twin_he.origin);
        if endpoints == (v_origin, v_dest) || endpoints == (v_dest, v_origin) {
            return Some(eh);
        }
    }
    None
}

fn find_shared_face_for_mef_preferring(
    brep: &BRep,
    v_origin: VertexHandle,
    v_dest: VertexHandle,
    prefer: Option<&std::collections::HashSet<FaceHandle>>,
) -> Option<(FaceHandle, HalfEdgeHandle, HalfEdgeHandle)> {
    let mut preferred: Option<(FaceHandle, HalfEdgeHandle, HalfEdgeHandle)> = None;
    let mut fallback: Option<(FaceHandle, HalfEdgeHandle, HalfEdgeHandle)> = None;
    for (fh, face) in &brep.faces {
        let outer_loop = face.outer_loop;
        let Some(loop_data) = brep.loops.get(outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut he_origin: Option<HalfEdgeHandle> = None;
        let mut he_dest: Option<HalfEdgeHandle> = None;
        let mut cursor = start;
        let max_iters = brep.half_edges.len() + 1;
        let mut iters = 0;
        while let Some(he) = brep.half_edges.get(cursor) {
            if he.origin == v_origin && he_origin.is_none() {
                he_origin = Some(cursor);
            }
            if he.origin == v_dest && he_dest.is_none() {
                he_dest = Some(cursor);
            }
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
        if let (Some(ho), Some(hd)) = (he_origin, he_dest) {
            if prefer.is_some_and(|s| s.contains(&fh)) {
                preferred = Some((fh, ho, hd));
                break;
            }
            if fallback.is_none() {
                fallback = Some((fh, ho, hd));
            }
        }
    }
    preferred.or(fallback)
}

fn create_apex_edge_preferring(
    brep: &mut BRep,
    v_origin: VertexHandle,
    v_dest: VertexHandle,
    target_edge: &TargetEdge,
    prefer: Option<&std::collections::HashSet<FaceHandle>>,
) -> Result<(EdgeHandle, (FaceHandle, FaceHandle)), String> {
    let (face, he_origin, he_dest) =
        find_shared_face_for_mef_preferring(brep, v_origin, v_dest, prefer).ok_or_else(|| {
            format!(
                "no face has both v_origin {v_origin:?} and v_dest {v_dest:?} on its outer loop"
            )
        })?;
    let he_origin_chk = find_halfedge_from_vertex_in_face(brep, face, v_origin)
        .map_err(|e| format!("origin half-edge lookup: {e}"))?;
    let he_dest_chk = find_halfedge_from_vertex_in_face(brep, face, v_dest)
        .map_err(|e| format!("dest half-edge lookup: {e}"))?;
    debug_assert_eq!(he_origin_chk, he_origin);
    debug_assert_eq!(he_dest_chk, he_dest);
    let mef_result = brep
        .mef(he_origin, he_dest)
        .map_err(|e| format!("mef on face {face:?}: {e:?}"))?;
    attach_curve_to_edge(brep, mef_result.edge, &target_edge.curve);
    Ok((mef_result.edge, (face, mef_result.face)))
}

pub(crate) fn attach_curve_to_edge(brep: &mut BRep, edge: EdgeHandle, curve: &NurbsCurve) {
    let p = curve.degree;
    let n = curve.control_points.len() - 1;
    let t_start = curve.knots[p];
    let t_end = curve.knots[n + 1];
    let geometry = CurveGeometry::Nurbs(curve.clone());
    let geometry = geometry
        .canonical_line_in_range((t_start, t_end), 1.0e-9)
        .unwrap_or(geometry);
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry,
        t_start,
        t_end,
    });
}
