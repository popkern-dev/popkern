use popkern_core::geometry::Point3;
use popkern_core::geometry::curve_inversion::invert_curve_at_point;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use super::super::constructive_topology;
use super::super::error::ApexBuildError;
use super::partition_cylinder_geometry::{
    cylinder_radial_vector, half_edge_end_vertex, point_is_not_edge_endpoint,
    point_lies_on_cylinder_face,
};
use super::plane_plane_cylinder_joint_points;

pub(super) fn materialize_support_boundary_generator_anchor(
    brep: &mut BRep,
    target_cylinder_face: FaceHandle,
    endpoint: VertexHandle,
    anchor: &SupportTopologyAnchorCandidate,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
    fillet_index: usize,
) -> Result<VertexHandle, ApexBuildError> {
    if let Some(existing) = coincident_vertex_on_face_outer_loop(
        brep,
        anchor.support_face,
        anchor.point,
        brep.tolerance.linear(),
    ) {
        if !point_lies_on_cylinder_face(
            brep,
            target_cylinder_face,
            brep.vertices[existing].position,
        ) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot reuse \
                     support boundary anchor {existing:?} for endpoint {endpoint:?}: it is not \
                     on target cylinder face {target_cylinder_face:?}"
                ),
            });
        }
        return Ok(existing);
    }
    let mut support_edges = Vec::new();
    if let Some(edge) = containing_outer_loop_edge_by_curve(brep, anchor.support_face, anchor.point)
    {
        support_edges.push((anchor.support_face, edge));
    }
    for support_face in allowed_support_faces {
        if *support_face == anchor.support_face {
            continue;
        }
        let Some(face_data) = brep.faces.get(*support_face) else {
            continue;
        };
        let Some(SurfaceGeometry::Plane { origin, normal, .. }) =
            face_data.surface.as_ref().map(|surface| &surface.geometry)
        else {
            continue;
        };
        if (anchor.point - *origin).dot(normal).abs() > brep.tolerance.linear().max(1.0e-9) {
            continue;
        }
        if let Some(edge) = containing_outer_loop_edge_by_curve(brep, *support_face, anchor.point) {
            support_edges.push((*support_face, edge));
        }
    }
    support_edges.sort_by_key(|(support_face, edge)| (*support_face, *edge));
    support_edges.dedup_by_key(|(_, edge)| *edge);
    if support_edges.is_empty() {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot materialize \
                 support-aware generator anchor for endpoint {endpoint:?}: computed point \
                 ({:.6},{:.6},{:.6}) lies on support face {:?} but not on a splittable \
                 curve-owned support boundary edge",
                anchor.point.x, anchor.point.y, anchor.point.z, anchor.support_face
            ),
        });
    }
    if std::env::var_os("POPKERN_A2_TRACE").is_some() {
        eprintln!(
            "[a2-trace §14.9/cylinder-anchor] support boundary promotion for endpoint \
             {endpoint:?} on {:?}: {:?}",
            target_cylinder_face, support_edges
        );
    }
    let tolerance = brep.tolerance;
    let mut primary_vertex = None;
    for (support_face, edge) in support_edges {
        let split = brep
            .split_edge_at_point(edge, anchor.point, &tolerance)
            .map_err(|err| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] failed to split \
                     support boundary edge {edge:?} on {support_face:?} for endpoint \
                     {endpoint:?}: {err:?}"
                ),
            })?;
        if let Some(primary) = primary_vertex {
            brep.weld_vertex_pair(primary, split.vertex)
                .map_err(|err| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] failed to \
                         weld support-promoted anchor {:?} into {:?}: {err:?}",
                        split.vertex, primary
                    ),
                })?;
        } else {
            primary_vertex = Some(split.vertex);
        }
    }
    let split_vertex = primary_vertex.ok_or_else(|| ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition FilletEnd[{fillet_index}] did not produce a support \
             boundary anchor for endpoint {endpoint:?}"
        ),
    })?;
    if !point_lies_on_cylinder_face(
        brep,
        target_cylinder_face,
        brep.vertices[split_vertex].position,
    ) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] split support \
                 boundary at {split_vertex:?}, but the resulting anchor is not on \
                 target cylinder face {target_cylinder_face:?}",
                split_vertex = split_vertex
            ),
        });
    }
    Ok(split_vertex)
}

#[derive(Clone)]
pub(super) struct SupportTopologyAnchorCandidate {
    pub(super) point: Point3,
    support_face: FaceHandle,
    endpoint_distance: f64,
}

pub(super) fn support_topology_generator_anchor_pair(
    brep: &BRep,
    face: FaceHandle,
    start_endpoint: Point3,
    end_endpoint: Point3,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
) -> Option<(
    SupportTopologyAnchorCandidate,
    SupportTopologyAnchorCandidate,
)> {
    let start_candidates =
        support_topology_anchor_candidates(brep, face, start_endpoint, allowed_support_faces);
    let end_candidates =
        support_topology_anchor_candidates(brep, face, end_endpoint, allowed_support_faces);
    let mut pairs = Vec::new();
    for start in &start_candidates {
        for end in &end_candidates {
            let direct_support = start.support_face == end.support_face;
            let transition_support = !plane_plane_cylinder_joint_points(
                brep,
                face,
                start.support_face,
                end.support_face,
            )
            .is_empty();
            if !direct_support && !transition_support {
                continue;
            }
            let score = start.endpoint_distance
                + end.endpoint_distance
                + if direct_support { 0.0 } else { 1.0e-6 };
            pairs.push((score, start.clone(), end.clone()));
        }
    }
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    pairs.into_iter().next().map(|(_, start, end)| (start, end))
}

fn support_topology_anchor_candidates(
    brep: &BRep,
    face: FaceHandle,
    endpoint: Point3,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
) -> Vec<SupportTopologyAnchorCandidate> {
    let Some(surface) = brep.faces.get(face).and_then(|face| face.surface.as_ref()) else {
        return Vec::new();
    };
    let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = surface.geometry
    else {
        return Vec::new();
    };
    let Some(radial) = cylinder_radial_vector(endpoint, origin, axis) else {
        return Vec::new();
    };
    if (radial.length() - radius).abs() > brep.tolerance.linear().max(radius.abs() * 1.0e-6) {
        return Vec::new();
    }
    let Some(radial_unit) = radial.normalized() else {
        return Vec::new();
    };
    let axis_len2 = axis.dot(&axis);
    if axis_len2 <= 1.0e-24 {
        return Vec::new();
    }
    let tol = brep.tolerance.linear().max(1.0e-9);
    let line_origin = origin + radial_unit * radius;
    let mut candidates = Vec::new();
    for (support_face, face_data) in brep.faces.iter() {
        if support_face == face {
            continue;
        }
        if !allowed_support_faces.contains(&support_face) {
            continue;
        }
        let Some(SurfaceGeometry::Plane {
            origin: plane_origin,
            normal,
            ..
        }) = face_data.surface.as_ref().map(|surface| &surface.geometry)
        else {
            continue;
        };
        let denom = axis.dot(normal);
        if denom.abs() <= 1.0e-12 {
            continue;
        }
        let lambda = (*plane_origin - line_origin).dot(normal) / denom;
        let point = line_origin + axis * lambda;
        if point.distance_to(&endpoint) <= tol {
            continue;
        }
        if !point_lies_on_cylinder_face(brep, face, point) {
            continue;
        }
        if (point - *plane_origin).dot(normal).abs() > tol {
            continue;
        }
        if candidates
            .iter()
            .any(|candidate: &SupportTopologyAnchorCandidate| {
                candidate.point.distance_to(&point) <= tol && candidate.support_face == support_face
            })
        {
            continue;
        }
        candidates.push(SupportTopologyAnchorCandidate {
            point,
            support_face,
            endpoint_distance: point.distance_to(&endpoint),
        });
    }
    candidates.sort_by(|a, b| {
        a.endpoint_distance
            .partial_cmp(&b.endpoint_distance)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    candidates
}

pub(super) fn generator_anchor_blocker(brep: &BRep, face: FaceHandle, endpoint: Point3) -> String {
    let Some(surface) = brep.faces.get(face).and_then(|face| face.surface.as_ref()) else {
        return "target face has no surface".to_string();
    };
    let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = surface.geometry
    else {
        return "target face is not a cylinder".to_string();
    };
    let Some(endpoint_radial) = cylinder_radial_vector(endpoint, origin, axis) else {
        return "endpoint has degenerate radial vector to cylinder axis".to_string();
    };
    let endpoint_surface_error = (endpoint_radial.length() - radius).abs();
    let Some(endpoint_unit) = endpoint_radial.normalized() else {
        return "endpoint radial vector cannot be normalized".to_string();
    };
    let Some(vertices) = constructive_topology::collect_outer_loop_vertices(brep, face) else {
        return "target cylinder face has no walkable outer loop".to_string();
    };
    let mut reports = vertices
        .into_iter()
        .filter_map(|vertex| {
            let point = brep.vertices.get(vertex)?.position;
            let radial = cylinder_radial_vector(point, origin, axis)?;
            let radial_distance = radial.length();
            let radial_unit = radial.normalized()?;
            let surface_error = (radial_distance - radius).abs();
            let radial_alignment = radial_unit.dot(&endpoint_unit);
            Some((vertex, surface_error, radial_alignment, point))
        })
        .collect::<Vec<_>>();
    reports.sort_by(|a, b| {
        let a_score = a.1 + (1.0 - a.2).abs() * radius.abs().max(1.0e-9);
        let b_score = b.1 + (1.0 - b.2).abs() * radius.abs().max(1.0e-9);
        a_score
            .partial_cmp(&b_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let nearest = reports
        .iter()
        .take(4)
        .map(|(vertex, surface_error, radial_alignment, point)| {
            format!(
                "{vertex:?}@({:.6},{:.6},{:.6}) surface_error={surface_error:.3e}, \
                 radial_alignment={radial_alignment:.6}",
                point.x, point.y, point.z
            )
        })
        .collect::<Vec<_>>()
        .join("; ");
    format!(
        "endpoint cylinder surface_error={endpoint_surface_error:.3e}; nearest radial loop \
         candidates [{}]. A commercial rebuild must materialize anchor vertices on the stored \
         cylinder surface before wiring generator trims.",
        nearest
    )
}

pub(super) fn outer_loop_path_between_vertices(
    brep: &BRep,
    face: FaceHandle,
    start_vertex: VertexHandle,
    end_vertex: VertexHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    let face_data = brep.faces.get(face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let mut cursor = loop_data.half_edge;
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        if half_edge.origin == start_vertex {
            let mut path = Vec::new();
            let mut path_cursor = cursor;
            for _ in 0..=brep.half_edges.len() {
                let path_half_edge = brep.half_edges.get(path_cursor)?;
                path.push(path_cursor);
                let end = half_edge_end_vertex(brep, path_cursor)?;
                if end == end_vertex {
                    return Some(path);
                }
                path_cursor = path_half_edge.next;
                if path_cursor == cursor {
                    return None;
                }
            }
            return None;
        }
        cursor = half_edge.next;
        if cursor == loop_data.half_edge {
            break;
        }
    }
    None
}

pub(super) fn path_off_cylinder_vertices(
    brep: &BRep,
    face: FaceHandle,
    path: &[HalfEdgeHandle],
) -> Vec<VertexHandle> {
    path.iter()
        .filter_map(|half_edge| {
            let vertex = brep.half_edges.get(*half_edge)?.origin;
            let point = brep.vertices.get(vertex)?.position;
            (!point_lies_on_cylinder_face(brep, face, point)).then_some(vertex)
        })
        .collect()
}

pub(super) fn vertex_is_on_face_outer_loop(
    brep: &BRep,
    face: FaceHandle,
    vertex: VertexHandle,
) -> bool {
    constructive_topology::collect_outer_loop_vertices(brep, face)
        .is_some_and(|vertices| vertices.contains(&vertex))
}

pub(super) fn coincident_vertex_on_face_outer_loop(
    brep: &BRep,
    face: FaceHandle,
    point: Point3,
    tolerance: f64,
) -> Option<VertexHandle> {
    let vertices = constructive_topology::collect_outer_loop_vertices(brep, face)?;
    let tol2 = tolerance * tolerance;
    vertices.into_iter().find(|vertex| {
        brep.vertices.get(*vertex).is_some_and(|vertex_data| {
            let delta = vertex_data.position - point;
            delta.dot(&delta) <= tol2
        })
    })
}

pub(super) fn containing_outer_loop_edge_by_curve(
    brep: &BRep,
    face: FaceHandle,
    point: Point3,
) -> Option<EdgeHandle> {
    let face_data = brep.faces.get(face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        let edge = brep.edges.get(half_edge.edge)?;
        if let Some(curve) = edge.curve.as_ref()
            && invert_curve_at_point(curve, point, &brep.tolerance).is_ok()
            && point_is_not_edge_endpoint(brep, half_edge.edge, point)
        {
            return Some(half_edge.edge);
        }
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    None
}
