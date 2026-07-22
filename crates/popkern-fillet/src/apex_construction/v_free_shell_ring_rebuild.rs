//! Boundary-completion ring rebuild for winding reconciliation
//! (PK-153): a face contorted across both winding classes keeps its
//! trusted live sides and closes its ring by authoring twin halves
//! against coincident dead sides lying on its own surface.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle, HalfEdge, HalfEdgeHandle, VertexHandle};

use super::{side_face, side_span, trace};

/// Rebuild a contorted face's outer ring by boundary completion:
/// order the trusted sides into chains, then close the ring by
/// authoring twin halves against coincident dead sides lying on the
/// face's own surface. Returns seam vertices for coincident-record
/// merging, or None when the walk is ambiguous or fails to close.
pub(super) fn rebuild_face_ring(
    brep: &mut BRep,
    face: FaceHandle,
    trusted: &[HalfEdgeHandle],
) -> Option<Vec<VertexHandle>> {
    let tol = brep.tolerance.linear().max(1.0e-9) * 10.0;
    let loop_handle = brep.faces.get(face)?.outer_loop;
    let span_of = |brep: &BRep, side: HalfEdgeHandle| -> Option<(VertexHandle, VertexHandle)> {
        let data = brep.half_edges.get(side)?;
        let twin = brep.half_edges.get(data.twin)?;
        Some((data.origin, twin.origin))
    };
    let mut chain: Vec<HalfEdgeHandle> = Vec::new();
    let mut pool: Vec<HalfEdgeHandle> = trusted.to_vec();
    let first = pool.first().copied()?;
    pool.remove(0);
    chain.push(first);
    let (ring_start, mut open_end) = span_of(brep, first)?;
    let ring_start_pos = brep.vertices.get(ring_start)?.position;
    let mut seam_vertices = Vec::new();
    for _ in 0..brep.half_edges.len() + trusted.len() + 8 {
        let open_pos = brep.vertices.get(open_end)?.position;
        if open_pos.distance_to(&ring_start_pos) <= tol && chain.len() > 2 {
            break;
        }
        // Trusted continuation first.
        if let Some(index) = pool.iter().position(|side| {
            span_of(brep, *side).is_some_and(|(start, _)| {
                brep.vertices
                    .get(start)
                    .is_some_and(|v| v.position.distance_to(&open_pos) <= tol)
            })
        }) {
            let side = pool.remove(index);
            let (start, end) = span_of(brep, side)?;
            if start != open_end {
                seam_vertices.push(start);
                seam_vertices.push(open_end);
            }
            chain.push(side);
            open_end = end;
            continue;
        }
        // Dead-side completion: a dead side ending at the open
        // position whose span lies on this face's surface gets its
        // face-side twin authored here.
        let candidates: Vec<HalfEdgeHandle> = brep
            .half_edges
            .iter()
            .filter(|(_, data)| !brep.half_edges.contains_key(data.twin))
            .filter(|(handle, _)| side_face(brep, *handle) != Some(face))
            .filter(|(handle, _)| {
                let Some((start, end)) = side_span(brep, *handle) else {
                    return false;
                };
                end.distance_to(&open_pos) <= tol
                    && span_on_face_surface(brep, face, &start, &end, tol)
            })
            .map(|(handle, _)| handle)
            .collect();
        let [candidate] = candidates.as_slice() else {
            trace(&format!(
                "winding reconciliation: ring rebuild for {face:?} found {} dead \
                 continuation(s) at ({:.3},{:.3},{:.3})",
                candidates.len(),
                open_pos.x,
                open_pos.y,
                open_pos.z
            ));
            return None;
        };
        let dead = *candidate;
        let dead_data = brep.half_edges.get(dead)?;
        let dead_edge = dead_data.edge;
        let dead_origin = dead_data.origin;
        let dead_end = brep
            .half_edges
            .get(dead_data.next)
            .map(|next| next.origin)?;
        let new_side = brep.half_edges.insert(HalfEdge {
            origin: dead_end,
            twin: dead,
            next: Default::default(),
            prev: Default::default(),
            edge: dead_edge,
            loop_: loop_handle,
            pcurve: None,
        });
        brep.half_edges[dead].twin = new_side;
        brep.boundary_edges.remove(&dead_edge);
        if dead_end != open_end {
            seam_vertices.push(dead_end);
            seam_vertices.push(open_end);
        }
        chain.push(new_side);
        open_end = dead_origin;
    }
    let open_pos = brep.vertices.get(open_end)?.position;
    if open_pos.distance_to(&ring_start_pos) > tol || chain.len() < 3 {
        trace(&format!(
            "winding reconciliation: ring rebuild for {face:?} did not close ({} side(s))",
            chain.len()
        ));
        return None;
    }
    if open_end != ring_start {
        seam_vertices.push(open_end);
        seam_vertices.push(ring_start);
    }
    let n = chain.len();
    for k in 0..n {
        let side = chain[k];
        let next = chain[(k + 1) % n];
        brep.half_edges[side].next = next;
        brep.half_edges[next].prev = side;
        brep.half_edges[side].loop_ = loop_handle;
        brep.half_edges[side].pcurve = None;
    }
    brep.loops[loop_handle].half_edge = chain[0];
    Some(seam_vertices)
}

/// Both span endpoints lie on the face's surface (plane distance for
/// planes; other geometries decline so the walk stays conservative).
fn span_on_face_surface(
    brep: &BRep,
    face: FaceHandle,
    start: &Point3,
    end: &Point3,
    tol: f64,
) -> bool {
    let Some(sof) = brep.faces.get(face).and_then(|d| d.surface.as_ref()) else {
        return false;
    };
    let SurfaceGeometry::Plane { origin, normal, .. } = &sof.geometry else {
        return false;
    };
    let Some(n) = normal.normalized() else {
        return false;
    };
    let dist = |p: &Point3| ((*p - *origin).dot(&n)).abs();
    dist(start) <= tol && dist(end) <= tol
}
