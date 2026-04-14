//! Closed fillet edge construction (§14 — plane-cylinder hole rims).
//!
//! When the fillet edge is a closed loop (vertex `start` == vertex
//! `end`, e.g. the rim of a cylindrical hole on a planar face),
//! the standard "split perpendicular edges + trim faces + KEF
//! merge" pipeline does not apply: there are no perpendicular
//! boundary edges to anchor on, and the fillet face is an annulus
//! bounded by two closed contact curves rather than a strip
//! between two anchor pairs.
//!
//! This module implements direct topology surgery for the closed
//! case:
//!
//! 1. Identify the two half-edges of the closed fillet edge,
//!    noting which loops they live on (one is typically a single-
//!    half-edge inner loop on a planar face, the other is part
//!    of a multi-half-edge loop on a curved face like a cylinder
//!    that has a seam).
//! 2. Compute the new contact-point positions on each support
//!    surface from the [`crate::construct::EdgeGeometry`]'s
//!    trim curves.
//! 3. Move the existing rim vertex (used by both the closed
//!    fillet edge AND any seam edges sharing that vertex) to the
//!    new position on the curved support's contact circle.
//! 4. Update the seam edge's curve geometry (if any) to the new
//!    shorter span.
//! 5. Create a new vertex on the planar support's contact circle.
//! 6. Create two new edges (one per contact circle) and four new
//!    half-edges (two per edge, one for each side).
//! 7. Splice the new half-edges into the support face loops in
//!    place of the old rim half-edges.
//! 8. Create the fillet face as an annulus with two boundary
//!    loops, one per new contact circle.
//! 9. Delete the old rim edge and its half-edges.
//!
//! Round 5c-3 scope: handles plane-cylinder closed fillets where
//! one side is a single-half-edge inner loop and the other side
//! is a multi-half-edge loop with at most one seam edge sharing
//! the rim vertex.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{
    BRep, Edge, EdgeHandle, Face, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, Vertex, VertexHandle,
};

use crate::construct::EdgeGeometry;
use crate::error::BlendFailure;

/// True when a fillet edge's start and end vertices are the same —
/// the topological signature of a closed fillet edge (e.g. a
/// circular hole rim).
pub(crate) fn is_closed_fillet_edge(brep: &BRep, edge: EdgeHandle) -> bool {
    let Some(edata) = brep.edges.get(edge) else {
        return false;
    };
    let Some(he) = brep.half_edges.get(edata.half_edge) else {
        return false;
    };
    let Some(twin) = brep.half_edges.get(he.twin) else {
        return false;
    };
    he.origin == twin.origin
}

/// Apply the closed-fillet topology surgery for a single
/// [`EdgeGeometry`] whose `analysis.edge` is a closed fillet edge.
///
/// Mutates `brep` in-place. Returns `Err` with the same
/// [`BlendFailure`] taxonomy the open path uses when the surgery
/// can't proceed (e.g. the precomputed geometry is inconsistent
/// with the BRep topology).
pub(crate) fn apply_closed_fillet(brep: &mut BRep, eg: &EdgeGeometry) -> Result<(), BlendFailure> {
    let fillet_edge = eg.analysis.edge;
    let edge_loc = crate::geometry_compute::edge_midpoint_by_handle(brep, fillet_edge);

    // 1. Locate the two half-edges of the closed edge and their loops.
    let edata = brep
        .edges
        .get(fillet_edge)
        .ok_or_else(|| BlendFailure::InvalidTopology {
            reason: format!("closed-fillet: edge {fillet_edge:?} not found"),
            location: edge_loc,
        })?;
    let he_a_h = edata.half_edge;
    let he_a =
        brep.half_edges
            .get(he_a_h)
            .cloned()
            .ok_or_else(|| BlendFailure::InvalidTopology {
                reason: format!("closed-fillet: half-edge {he_a_h:?} not found"),
                location: edge_loc,
            })?;
    let he_b_h = he_a.twin;
    let he_b =
        brep.half_edges
            .get(he_b_h)
            .cloned()
            .ok_or_else(|| BlendFailure::InvalidTopology {
                reason: format!("closed-fillet: twin half-edge {he_b_h:?} not found"),
                location: edge_loc,
            })?;

    let v_old = he_a.origin;
    if v_old != he_b.origin {
        return Err(BlendFailure::InvalidTopology {
            reason: format!("closed-fillet: edge {fillet_edge:?} is not closed (v_start != v_end)"),
            location: edge_loc,
        });
    }

    let loop_a = he_a.loop_;
    let loop_b = he_b.loop_;
    let face_a = brep
        .loops
        .get(loop_a)
        .ok_or_else(|| BlendFailure::InvalidTopology {
            reason: format!("closed-fillet: loop {loop_a:?} not found"),
            location: edge_loc,
        })?
        .face;
    let face_b = brep
        .loops
        .get(loop_b)
        .ok_or_else(|| BlendFailure::InvalidTopology {
            reason: format!("closed-fillet: loop {loop_b:?} not found"),
            location: edge_loc,
        })?
        .face;

    // 2. Pull the new contact-point positions from the precomputed
    //    trim curves. Both curves are closed, so start point ==
    //    end point — taking either is fine.
    let (t1_start, _) = eg.trim_3d_1.parameter_range();
    let (t2_start, _) = eg.trim_3d_2.parameter_range();
    let new_pos_a = eg.trim_3d_1.point_at(t1_start);
    let new_pos_b = eg.trim_3d_2.point_at(t2_start);

    // The geometry pipeline associates trim_3d_1 with face_1 and
    // trim_3d_2 with face_2 (per analysis.face_1/face_2). We need
    // to know which trim corresponds to which side of the rim
    // (loop_a vs loop_b). Use the analysis face IDs to map.
    let (pos_for_face_a, pos_for_face_b) = if face_a == eg.analysis.face_1 {
        (new_pos_a, new_pos_b)
    } else if face_a == eg.analysis.face_2 {
        (new_pos_b, new_pos_a)
    } else {
        return Err(BlendFailure::InvalidTopology {
            reason: format!(
                "closed-fillet: face {face_a:?} on rim is neither analysis.face_1 \
                 ({:?}) nor analysis.face_2 ({:?})",
                eg.analysis.face_1, eg.analysis.face_2
            ),
            location: edge_loc,
        });
    };
    let _ = face_b;

    // 3. Identify which side has a "single-half-edge loop" (the
    //    simple inner-loop case where the rim is an entire
    //    one-element loop) vs which side has a "multi-half-edge
    //    loop" (where the rim is one element among several, e.g.
    //    a cylinder with a seam edge).
    //
    //    The single-half-edge case: the loop contains exactly the
    //    rim half-edge with next == prev == self.
    //    The multi-half-edge case: rim's prev and next point at
    //    other half-edges in the loop chain.
    let a_is_single = he_a.next == he_a_h && he_a.prev == he_a_h;
    let b_is_single = he_b.next == he_b_h && he_b.prev == he_b_h;

    // For each side, identify any seam edge that shares the rim
    // vertex `v_old` and needs to be shortened so its endpoint
    // lands on the new contact circle. Walk the same loop looking
    // for half-edges other than the rim whose `origin` is `v_old`.
    let seam_a = if a_is_single {
        None
    } else {
        find_seam_edge_at_vertex(brep, loop_a, he_a_h, v_old)
    };
    let seam_b = if b_is_single {
        None
    } else {
        find_seam_edge_at_vertex(brep, loop_b, he_b_h, v_old)
    };

    // 4. Move the rim vertex to whichever side has the seam edge —
    //    that's the side where the existing vertex stays "owned"
    //    by another edge after the rim is removed. The other side
    //    gets a fresh new vertex.
    //
    //    For our plane-cylinder case: side A (cylinder) has a seam
    //    so v_old → new_pos_for_face_a. Side B (cube top) is a
    //    single-half-edge loop so a fresh vertex at
    //    new_pos_for_face_b.
    let (v_with_seam, pos_with_seam, v_fresh_pos) = if seam_a.is_some() {
        (v_old, pos_for_face_a, pos_for_face_b)
    } else if seam_b.is_some() {
        (v_old, pos_for_face_b, pos_for_face_a)
    } else {
        // Both sides single-half-edge loops — no seam to anchor
        // the existing vertex. Pick side A by convention; the
        // surgery still works because there's no seam to leave
        // dangling.
        (v_old, pos_for_face_a, pos_for_face_b)
    };

    brep.vertices[v_with_seam].position = pos_with_seam;

    // 5. If we moved a vertex that's an endpoint of a seam edge,
    //    update the seam edge's curve so its endpoint matches the
    //    new vertex position.
    if let Some(seam_edge) = seam_a.or(seam_b) {
        rebuild_seam_curve(brep, seam_edge, pos_with_seam);
    }

    // 6. Create a fresh vertex on the OTHER support's contact
    //    circle.
    let v_fresh = brep.vertices.insert(Vertex {
        position: v_fresh_pos,
        outgoing: HalfEdgeHandle::default(),
        name: None,
        tolerance: None,
    });

    // 7. Create the two new contact-circle edges. Each carries
    //    the corresponding trim_3d as its curve geometry.
    let edge_seam_side = brep.edges.insert(Edge {
        half_edge: HalfEdgeHandle::default(),
        curve: Some(popkern_core::geometry::curve::CurveOnEdge {
            geometry: popkern_core::geometry::curve::CurveGeometry::Nurbs(if seam_a.is_some() {
                eg.trim_3d_1.clone()
            } else {
                eg.trim_3d_2.clone()
            }),
            t_start: 0.0,
            t_end: 1.0,
        }),
        name: None,
        tolerance: None,
        canonical: Default::default(),
    });
    let edge_fresh_side = brep.edges.insert(Edge {
        half_edge: HalfEdgeHandle::default(),
        curve: Some(popkern_core::geometry::curve::CurveOnEdge {
            geometry: popkern_core::geometry::curve::CurveGeometry::Nurbs(if seam_a.is_some() {
                eg.trim_3d_2.clone()
            } else {
                eg.trim_3d_1.clone()
            }),
            t_start: 0.0,
            t_end: 1.0,
        }),
        name: None,
        tolerance: None,
        canonical: Default::default(),
    });

    // 8. Create the fillet face shell — we need this BEFORE the
    //    half-edges so each half-edge can reference its loop.
    //    Build with placeholder loop handles, then patch.
    let shell = brep.faces[face_a].shell;
    let fillet_face = brep.faces.insert(Face {
        outer_loop: LoopHandle::default(),
        inner_loops: vec![],
        shell,
        surface: Some(popkern_core::geometry::surface::SurfaceOnFace {
            geometry: eg.surface.clone(),
            orientation: true,
        }),
        name: None,
        tolerance: None,
    });

    let loop_fillet_outer = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: fillet_face,
    });
    let loop_fillet_inner = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: fillet_face,
    });

    brep.faces[fillet_face].outer_loop = loop_fillet_outer;
    brep.faces[fillet_face].inner_loops = vec![loop_fillet_inner];

    // 9. Create the four new half-edges. The "seam side" gets
    //    v_with_seam (the moved vertex); the "fresh side" gets
    //    v_fresh.
    //
    //    For a closed contact circle each half-edge has next ==
    //    prev == self when it lives on a single-half-edge loop
    //    (the fillet face's loops are always single-element).
    //    On the support side, the seam-side half-edge gets
    //    spliced into the existing multi-element loop and the
    //    fresh-side half-edge becomes the sole element of the
    //    single-element inner loop.
    let support_seam_loop = if seam_a.is_some() { loop_a } else { loop_b };
    let support_fresh_loop = if seam_a.is_some() { loop_b } else { loop_a };
    let support_seam_old_he = if seam_a.is_some() { he_a_h } else { he_b_h };
    let _support_fresh_old_he = if seam_a.is_some() { he_b_h } else { he_a_h };

    // Pre-fetch the prev/next of the seam-side rim half-edge so
    // we can splice the new one in.
    let (seam_old_prev, seam_old_next) = {
        let h = &brep.half_edges[support_seam_old_he];
        (h.prev, h.next)
    };

    let he_seam_support = brep.half_edges.insert(HalfEdge {
        origin: v_with_seam,
        twin: HalfEdgeHandle::default(),
        next: seam_old_next,
        prev: seam_old_prev,
        edge: edge_seam_side,
        loop_: support_seam_loop,
        pcurve: None,
    });
    let he_seam_fillet = brep.half_edges.insert(HalfEdge {
        origin: v_with_seam,
        twin: he_seam_support,
        next: HalfEdgeHandle::default(),
        prev: HalfEdgeHandle::default(),
        edge: edge_seam_side,
        loop_: loop_fillet_outer,
        pcurve: None,
    });
    brep.half_edges[he_seam_support].twin = he_seam_fillet;
    // The fillet-side half-edge is the sole element of its loop:
    // next/prev point back at itself.
    brep.half_edges[he_seam_fillet].next = he_seam_fillet;
    brep.half_edges[he_seam_fillet].prev = he_seam_fillet;

    let he_fresh_support = brep.half_edges.insert(HalfEdge {
        origin: v_fresh,
        twin: HalfEdgeHandle::default(),
        next: HalfEdgeHandle::default(),
        prev: HalfEdgeHandle::default(),
        edge: edge_fresh_side,
        loop_: support_fresh_loop,
        pcurve: None,
    });
    brep.half_edges[he_fresh_support].next = he_fresh_support;
    brep.half_edges[he_fresh_support].prev = he_fresh_support;

    let he_fresh_fillet = brep.half_edges.insert(HalfEdge {
        origin: v_fresh,
        twin: he_fresh_support,
        next: HalfEdgeHandle::default(),
        prev: HalfEdgeHandle::default(),
        edge: edge_fresh_side,
        loop_: loop_fillet_inner,
        pcurve: None,
    });
    brep.half_edges[he_fresh_support].twin = he_fresh_fillet;
    brep.half_edges[he_fresh_fillet].next = he_fresh_fillet;
    brep.half_edges[he_fresh_fillet].prev = he_fresh_fillet;

    // Update the new edges' half_edge back-references.
    brep.edges[edge_seam_side].half_edge = he_seam_support;
    brep.edges[edge_fresh_side].half_edge = he_fresh_support;

    // 10. Splice the seam-side new half-edge into the seam-side
    //     loop in place of the old rim half-edge.
    brep.half_edges[seam_old_prev].next = he_seam_support;
    brep.half_edges[seam_old_next].prev = he_seam_support;
    if brep.loops[support_seam_loop].half_edge == support_seam_old_he {
        brep.loops[support_seam_loop].half_edge = he_seam_support;
    }

    // 11. Replace the fresh-side single-half-edge inner loop's
    //     reference to point at the new half-edge.
    brep.loops[support_fresh_loop].half_edge = he_fresh_support;

    // 12. Set the fillet-face loops' half_edge pointers and update
    //     vertex outgoing references where needed.
    brep.loops[loop_fillet_outer].half_edge = he_seam_fillet;
    brep.loops[loop_fillet_inner].half_edge = he_fresh_fillet;
    brep.vertices[v_fresh].outgoing = he_fresh_support;
    if !brep
        .half_edges
        .contains_key(brep.vertices[v_with_seam].outgoing)
    {
        brep.vertices[v_with_seam].outgoing = he_seam_support;
    }

    // 13. Delete the old rim half-edges and edge.
    brep.half_edges.remove(he_a_h);
    brep.half_edges.remove(he_b_h);
    brep.edges.remove(fillet_edge);

    Ok(())
}

/// Find a half-edge in `loop_handle` other than `skip` whose
/// `origin` is `vertex` — i.e. an edge incident to the rim vertex
/// from a non-rim direction. For the cylinder case this returns
/// the seam edge half-edge.
fn find_seam_edge_at_vertex(
    brep: &BRep,
    loop_handle: LoopHandle,
    skip: HalfEdgeHandle,
    vertex: VertexHandle,
) -> Option<EdgeHandle> {
    let start = brep.loops.get(loop_handle)?.half_edge;
    let mut cursor = start;
    let mut iters = 0;
    loop {
        iters += 1;
        if iters > 10_000 {
            return None;
        }
        if cursor != skip {
            let h = brep.half_edges.get(cursor)?;
            if h.origin == vertex {
                return Some(h.edge);
            }
        }
        let next = brep.half_edges.get(cursor)?.next;
        cursor = next;
        if cursor == start {
            break;
        }
    }
    None
}

/// Replace a seam edge's curve geometry with a straight line that
/// preserves the canonical "forward" direction of the edge — the
/// curve starts at `edge.half_edge.origin`'s current position and
/// ends at the twin's `origin` current position. Used to update
/// the seam after one of its endpoint vertices moved to the new
/// contact-circle height.
///
/// The `_end_pos` parameter is kept for documentation but no
/// longer used — both vertex positions are read directly from the
/// BRep.
fn rebuild_seam_curve(brep: &mut BRep, edge: EdgeHandle, _end_pos: Point3) {
    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};

    let Some(edata) = brep.edges.get(edge) else {
        return;
    };
    let he = edata.half_edge;
    let Some(he_data) = brep.half_edges.get(he).cloned() else {
        return;
    };
    let Some(twin_data) = brep.half_edges.get(he_data.twin).cloned() else {
        return;
    };

    let v_start_pos = brep.vertices[he_data.origin].position;
    let v_end_pos = brep.vertices[twin_data.origin].position;
    let direction = Vec3::new(
        v_end_pos.x - v_start_pos.x,
        v_end_pos.y - v_start_pos.y,
        v_end_pos.z - v_start_pos.z,
    );
    if direction.length() < 1e-12 {
        return;
    }

    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: v_start_pos,
            direction,
        },
        t_start: 0.0,
        t_end: 1.0,
    });
}
