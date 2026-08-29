//! Semantic support-side recovery for exact-boundary sector extraction.

use popkern_core::geometry::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

pub(super) fn same_oriented_cylinder_support_side(
    brep: &BRep,
    cap_face: FaceHandle,
    cap_half_edge: HalfEdgeHandle,
) -> Option<(HalfEdgeHandle, FaceHandle)> {
    same_oriented_support_side_matching(brep, cap_face, cap_half_edge, |surface| {
        matches!(surface.geometry, SurfaceGeometry::Cylinder { .. })
    })
}

pub(super) fn same_oriented_support_side(
    brep: &BRep,
    cap_face: FaceHandle,
    cap_half_edge: HalfEdgeHandle,
) -> Option<(HalfEdgeHandle, FaceHandle)> {
    same_oriented_support_side_matching(brep, cap_face, cap_half_edge, |_| true)
}

pub(super) fn coincident_support_side(
    brep: &BRep,
    cap_face: FaceHandle,
    cap_half_edge: HalfEdgeHandle,
) -> Option<(HalfEdgeHandle, FaceHandle)> {
    for (face, face_data) in brep.faces.iter() {
        if face == cap_face {
            continue;
        }
        let Some(sides) = collect_loop_half_edges(brep, face_data.outer_loop) else {
            continue;
        };
        for side in sides {
            if coincident_side(brep, cap_half_edge, side) {
                return Some((side, face));
            }
        }
    }
    None
}

fn same_oriented_support_side_matching(
    brep: &BRep,
    cap_face: FaceHandle,
    cap_half_edge: HalfEdgeHandle,
    surface_matches: impl Fn(&popkern_core::geometry::surface::SurfaceOnFace) -> bool,
) -> Option<(HalfEdgeHandle, FaceHandle)> {
    support_side_matching(
        brep,
        cap_face,
        cap_half_edge,
        surface_matches,
        same_oriented_side,
    )
}

fn support_side_matching(
    brep: &BRep,
    cap_face: FaceHandle,
    cap_half_edge: HalfEdgeHandle,
    surface_matches: impl Fn(&popkern_core::geometry::surface::SurfaceOnFace) -> bool,
    side_matches: impl Fn(&BRep, HalfEdgeHandle, HalfEdgeHandle) -> bool,
) -> Option<(HalfEdgeHandle, FaceHandle)> {
    for (face, face_data) in brep.faces.iter() {
        if face == cap_face {
            continue;
        }
        let Some(surface) = face_data.surface.as_ref() else {
            continue;
        };
        if !surface_matches(surface) {
            continue;
        }
        let Some(sides) = collect_loop_half_edges(brep, face_data.outer_loop) else {
            continue;
        };
        for side in sides {
            if side_matches(brep, cap_half_edge, side) {
                return Some((side, face));
            }
        }
    }
    None
}

fn same_oriented_side(brep: &BRep, cap_side: HalfEdgeHandle, support_side: HalfEdgeHandle) -> bool {
    let Some(cap_start) = brep.half_edges.get(cap_side).map(|data| data.origin) else {
        return false;
    };
    let Some(support_start) = brep.half_edges.get(support_side).map(|data| data.origin) else {
        return false;
    };
    if cap_start != support_start {
        return false;
    }
    let Some(cap_end) = half_edge_destination_vertex(brep, cap_side) else {
        return false;
    };
    let Some(support_end) = half_edge_destination_vertex(brep, support_side) else {
        return false;
    };
    cap_end == support_end
}

fn coincident_side(brep: &BRep, cap_side: HalfEdgeHandle, support_side: HalfEdgeHandle) -> bool {
    let Some(cap_start) = brep.half_edges.get(cap_side).map(|data| data.origin) else {
        return false;
    };
    let Some(support_start) = brep.half_edges.get(support_side).map(|data| data.origin) else {
        return false;
    };
    let Some(cap_end) = half_edge_destination_vertex(brep, cap_side) else {
        return false;
    };
    let Some(support_end) = half_edge_destination_vertex(brep, support_side) else {
        return false;
    };
    (cap_start == support_start && cap_end == support_end)
        || (cap_start == support_end && cap_end == support_start)
}

fn collect_loop_half_edges(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    let loop_data = brep.loops.get(loop_handle)?;
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    let mut cycle = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        cycle.push(cursor);
        cursor = brep.half_edges.get(cursor)?.next;
        if cursor == start {
            return Some(cycle);
        }
    }
    None
}

fn half_edge_destination_vertex(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<VertexHandle> {
    let data = brep.half_edges.get(half_edge)?;
    if let Some(twin) = brep.half_edges.get(data.twin) {
        return Some(twin.origin);
    }
    if !brep.boundary_edges.contains(&data.edge) {
        return brep
            .half_edges
            .get(data.next)
            .map(|next_half_edge| next_half_edge.origin);
    }
    let edge = brep.edges.get(data.edge)?;
    let Some(curve) = edge.curve.as_ref() else {
        return brep
            .half_edges
            .get(data.next)
            .map(|next_half_edge| next_half_edge.origin);
    };
    let start_position = brep.vertices.get(data.origin)?.position;
    let p0 = curve.geometry.point_at(curve.t_start);
    let p1 = curve.geometry.point_at(curve.t_end);
    let destination = if start_position.distance_to(&p0) <= start_position.distance_to(&p1) {
        p1
    } else {
        p0
    };
    let tolerance = brep.tolerance.linear().max(1.0e-9);
    brep.vertices
        .iter()
        .filter_map(|(vertex, data)| {
            let distance = data.position.distance_to(&destination);
            (distance <= tolerance).then_some((vertex, distance))
        })
        .min_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(vertex, _)| vertex)
}
