//! Topology helpers for vertex-blend insertion.

use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::coherence::PendingCoherenceMarker;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::construct::EdgeGeometry;
use crate::construct::blend_topo::{
    find_face_trim_anchors, find_halfedge_from_vertex_in_face, identify_corner_face,
};

/// Phase 1: Split off corner regions from each face around V using
/// the per-face fillet trim endpoints as mef anchors.
pub(super) fn split_corner_triangles(
    brep: &mut BRep,
    vertex: VertexHandle,
    face_info: &[(FaceHandle, HalfEdgeHandle)],
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Result<Vec<FaceHandle>, String> {
    let mut corner_faces = Vec::new();

    for (face, _he_arriving_at_v) in face_info {
        let anchors = find_face_trim_anchors(brep, *face, vertex, edge_geometries, tolerance)?;
        if anchors.len() < 2 {
            return Err(format!(
                "face {face:?} has {} trim-endpoint anchors around vertex {vertex:?}, expected >=2",
                anchors.len()
            ));
        }
        let (anchor_a, anchor_b) = (anchors[0], anchors[1]);

        let he_at_a = find_halfedge_from_vertex_in_face(brep, *face, anchor_a)?;
        let he_at_b = find_halfedge_from_vertex_in_face(brep, *face, anchor_b)?;

        let loop_a = brep.half_edges[he_at_a].loop_;
        let loop_b = brep.half_edges[he_at_b].loop_;
        if loop_a != loop_b {
            return Err(format!(
                "anchor half-edges on face {face:?} are in different loops"
            ));
        }

        let mef_result = brep
            .mef(he_at_a, he_at_b)
            .map_err(|e| format!("mef failed on face {face:?}: {e}"))?;

        author_corner_split_edge_geometry(
            brep,
            *face,
            anchor_a,
            anchor_b,
            mef_result.edge,
            mef_result.half_edge,
            mef_result.half_edge_twin,
        );

        let corner = identify_corner_face(brep, *face, mef_result.face, vertex)?;
        corner_faces.push(corner);
    }

    Ok(corner_faces)
}

/// Stage C1: author the canonical edge curve and the matching
/// pcurves for the MEF-inserted corner-split edge. The new edge
/// runs across the parent face's surface from `anchor_a` to
/// `anchor_b`. For planar surfaces the canonical curve is a line
/// segment between the two endpoint positions and the two pcurves
/// are linear segments in the parent surface's `(u, v)` frame. For
/// non-planar surfaces (curved face support is left to a follow-up
/// substage) we record `PendingCoherenceMarker` entries so
/// downstream consumers must author the missing geometry before
/// steady state.
fn author_corner_split_edge_geometry(
    brep: &mut BRep,
    face: FaceHandle,
    anchor_a: VertexHandle,
    anchor_b: VertexHandle,
    edge: EdgeHandle,
    half_edge_ab: HalfEdgeHandle,
    half_edge_ba: HalfEdgeHandle,
) {
    let pos_a = brep.vertices[anchor_a].position;
    let pos_b = brep.vertices[anchor_b].position;
    let tol = brep.tolerance.linear();
    let surface_geometry = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone());

    match surface_geometry {
        Some(SurfaceGeometry::Plane {
            origin,
            u_axis,
            v_axis,
            ..
        }) => {
            if let Some(c) = CurveOnEdge::line_segment(pos_a, pos_b, tol.max(1.0e-12)) {
                brep.edges[edge].curve = Some(c);
                brep.clear_pending_authoring_for_edge(edge);
            } else {
                brep.record_pending_authoring(PendingCoherenceMarker::EdgeCurvePending {
                    edge,
                    owner_operator: "split_corner_triangles",
                });
            }
            let uv_a = plane_uv(pos_a, origin, u_axis, v_axis);
            let uv_b = plane_uv(pos_b, origin, u_axis, v_axis);
            brep.half_edges[half_edge_ab].pcurve = Some(Pcurve::linear(uv_a, uv_b, 0.0, 1.0));
            brep.half_edges[half_edge_ba].pcurve = Some(Pcurve::linear(uv_b, uv_a, 0.0, 1.0));
            brep.clear_pending_authoring_for_half_edge(half_edge_ab);
            brep.clear_pending_authoring_for_half_edge(half_edge_ba);
        }
        _ => {
            brep.record_pending_authoring(PendingCoherenceMarker::EdgeCurvePending {
                edge,
                owner_operator: "split_corner_triangles",
            });
            brep.record_pending_authoring(PendingCoherenceMarker::HalfEdgePcurvePending {
                half_edge: half_edge_ab,
                owner_operator: "split_corner_triangles",
            });
            brep.record_pending_authoring(PendingCoherenceMarker::HalfEdgePcurvePending {
                half_edge: half_edge_ba,
                owner_operator: "split_corner_triangles",
            });
        }
    }
}

fn plane_uv(p: Point3, origin: Point3, u_axis: Vec3, v_axis: Vec3) -> (f64, f64) {
    let d = p - origin;
    (d.dot(&u_axis), d.dot(&v_axis))
}

/// Normalize an under-resolved high-valence blend loop before the
/// N-sector Gregory/polar poke. When Euler merging collapses several
/// boundary intervals into fewer than `n_sides` topological edges,
/// split the longest live loop edge until the loop has `n_sides`
/// edges. The existing poke then performs its ordinary N-to-2N
/// midpoint split.
pub(super) fn normalize_underresolved_face_loop_to_n(
    brep: &mut BRep,
    face: FaceHandle,
    n_sides: usize,
) -> Result<usize, String> {
    if n_sides < 5 {
        return Ok(0);
    }
    let mut splits = 0usize;
    loop {
        let loop_hes = outer_loop_half_edges(brep, face)?;
        let count = loop_hes.len();
        if count >= n_sides {
            return Ok(splits);
        }
        let Some(edge) = longest_loop_edge(brep, &loop_hes) else {
            return Err(format!(
                "underresolved face-loop normalizer found no splittable edge at count {count}"
            ));
        };
        brep.split_edge(edge, 0.5)
            .map_err(|e| format!("underresolved face-loop split failed: {e}"))?;
        splits += 1;
        if splits > n_sides {
            return Err(format!(
                "underresolved face-loop normalizer exceeded split budget for N={n_sides}"
            ));
        }
    }
}

fn outer_loop_half_edges(brep: &BRep, face: FaceHandle) -> Result<Vec<HalfEdgeHandle>, String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| "face handle is invalid".to_string())?;
    let loop_data = brep
        .loops
        .get(face_data.outer_loop)
        .ok_or_else(|| "face outer loop is invalid".to_string())?;
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return Err("outer loop start half-edge is invalid".to_string());
    }
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..64 {
        out.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            return Ok(out);
        }
    }
    Err("outer loop walk exceeded safety cap".to_string())
}

fn longest_loop_edge(brep: &BRep, loop_hes: &[HalfEdgeHandle]) -> Option<EdgeHandle> {
    loop_hes
        .iter()
        .filter_map(|&he| {
            let he_data = brep.half_edges.get(he)?;
            let next = brep.half_edges.get(he_data.next)?;
            let a = brep.vertices.get(he_data.origin)?.position;
            let b = brep.vertices.get(next.origin)?.position;
            Some((he_data.edge, a.distance_to(&b)))
        })
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(edge, _)| edge)
}

/// Phase 2b: kev the remaining V-bridge edge (if any) to un-pinch
/// the merged face and remove V in one Euler step.
///
/// At entry, V may have zero or one incident edges. Zero means
/// Phase 2a completed cleanly (all V-edges were cross-face kefs);
/// one means there's a same-loop bridge whose both half-edges are
/// on the merged face. `kev` requires the edge's destination
/// vertex to be valence-1, so we pick the half-edge that ends at
/// V and call `kev(that)` to kill V plus the bridge together.
///
/// `pub(crate)` for §6c-DirectApex.2.e — same rationale as
/// [`kef_cross_face_v_edges`].
pub(crate) fn kev_v_bridge(brep: &mut BRep, vertex: VertexHandle) -> Result<(), String> {
    if !brep.vertices.contains_key(vertex) {
        return Ok(());
    }

    let outgoings: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, h)| h.origin == vertex)
        .map(|(h, _)| h)
        .collect();

    if outgoings.is_empty() {
        return Ok(());
    }
    if outgoings.len() != 1 {
        return Err(format!(
            "kev_v_bridge: V has valence {} (expected 0 or 1 for single-bridge teardown)",
            outgoings.len()
        ));
    }

    let he_out = outgoings[0];
    let he_in = brep.half_edges[he_out].twin;
    // kev kills the destination vertex of its argument; we want to
    // kill V, so pass the half-edge whose destination is V.
    debug_assert_eq!(brep.half_edges[brep.half_edges[he_in].twin].origin, vertex);
    let same_loop = brep.half_edges[he_in].loop_ == brep.half_edges[he_out].loop_;
    let is_only_edge = brep.half_edges[he_in].next == he_out
        && brep.half_edges[he_out].next == he_in
        && brep.half_edges[he_in].prev == he_out
        && brep.half_edges[he_out].prev == he_in;
    match brep.kev(he_in) {
        Ok(_) => Ok(()),
        Err(_) if same_loop && !is_only_edge => {
            brep.dissolve_same_loop_bridge(he_in)
                .map_err(|e| format!("same-loop bridge dissolve failed on V-bridge: {e}"))?;
            remove_orphaned_vertex(brep, vertex);
            Ok(())
        }
        Err(kev_error) => Err(format!("kev failed on V-bridge: {kev_error}")),
    }
}

/// Phase 3: Remove vertex V if it is no longer connected to any edges.
///
/// `pub(crate)` for §6c-DirectApex.2.e — same rationale as
/// [`kef_cross_face_v_edges`].
pub(crate) fn remove_orphaned_vertex(brep: &mut BRep, vertex: VertexHandle) {
    let v_has_edges = brep.half_edges.values().any(|he| he.origin == vertex);
    if !v_has_edges {
        brep.vertices.remove(vertex);
    }
}
