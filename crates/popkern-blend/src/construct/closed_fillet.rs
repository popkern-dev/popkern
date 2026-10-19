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
//! 5. Create a new vertex on the planar support's contact circle,
//!    two new edges (one per contact circle) with four half-edges,
//!    and splice them into the support face loops in place of the
//!    old rim half-edges.
//! 6. Create the fillet face as an annulus with two boundary
//!    loops, one per new contact circle; delete the old rim edge.
//!
//! scope: handles plane-cylinder closed fillets where
//! one side is a single-half-edge inner loop and the other side
//! is a multi-half-edge loop with at most one seam edge sharing
//! the rim vertex.

use popkern_core::geometry::Point3;
use popkern_core::topology::{
    BRep, Edge, EdgeHandle, Face, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, Vertex,
};
use popkern_math::PortableMath;

use super::closed_fillet_seam;
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
pub(crate) fn apply_closed_fillet(
    brep: &mut BRep,
    eg: &EdgeGeometry,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<popkern_core::topology::FaceHandle, BlendFailure> {
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

    // 2. Compute the new contact-point positions on each support.
    //    NOT the trim curve's own start point: for a closed rim the
    //    trim NURBS inherits its parameterization origin from the
    //    offset-surface intersection marcher's seed, which is an
    //    arbitrary, ULP-sensitive artifact (PK-285 measured: a 7e-15
    //    extrude-height difference moved the anchor from the seam
    //    angle to the diametrically opposite side, and the rebuilt
    //    seam line crossed the cylinder as a chord — check-9 caught
    //    it at deviation = radius). The surgery's anchor must be
    //    SEAM-CONSISTENT: the point on each contact circle nearest
    //    the old rim-seam vertex, so the shortened seam edge stays a
    //    generator of its surface.
    let v_old_pos = brep.vertices[v_old].position;
    let new_pos_a =
        closed_fillet_seam::seam_consistent_anchor(&eg.trim_circle_1, &eg.trim_3d_1, v_old_pos);
    let new_pos_b =
        closed_fillet_seam::seam_consistent_anchor(&eg.trim_circle_2, &eg.trim_3d_2, v_old_pos);

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
        closed_fillet_seam::find_seam_edge_at_vertex(brep, loop_a, he_a_h, v_old)
    };
    let seam_b = if b_is_single {
        None
    } else {
        closed_fillet_seam::find_seam_edge_at_vertex(brep, loop_b, he_b_h, v_old)
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

    let seam_face = brep.loops[if seam_a.is_some() { loop_a } else { loop_b }].face;
    let (pos_with_seam, v_fresh_pos) =
        super::closed_fillet_rail::closed_rim_anchors(brep, eg, seam_a.or(seam_b), seam_face)
            .unwrap_or((pos_with_seam, v_fresh_pos));
    brep.vertices[v_with_seam].position = pos_with_seam;

    // 5. If we moved a vertex that's an endpoint of a seam edge,
    //    update the seam edge's curve so its endpoint matches the
    //    new vertex position.
    if let Some(seam_edge) = seam_a.or(seam_b) {
        closed_fillet_seam::rebuild_seam_curve(brep, seam_edge);
    }

    // 6. Create a fresh vertex on the OTHER support's contact
    //    circle.
    let v_fresh = brep.vertices.insert(Vertex {
        position: v_fresh_pos,
        outgoing: HalfEdgeHandle::default(),
        name: None,
        tolerance: None,
    });

    // 7. Create the two new contact-circle edges. Closed-fillet
    //    contact curves are exact circles (rolling-ball contacts of
    //    a circular rim); carry them as analytic Arc geometry — the
    //    NURBS fit through the footprints has ~3e-3 interpolation
    //    error that fails every downstream surface-association gate
    //    (tessellation, export), and the angle parameterization
    //    pairs exactly with iso-linear pcurves on the supports.
    //    Fit failure falls back to the interpolated trim.
    if popkern_core::trace::trace_enabled("fillet") {
        eprintln!(
            "  [closed-fillet] seam_a={} circle_1={} circle_2={}",
            seam_a.is_some(),
            eg.trim_circle_1.is_some(),
            eg.trim_circle_2.is_some()
        );
    }
    let (trim_seam, trim_fresh, circle_seam, circle_fresh) = if seam_a.is_some() {
        (
            &eg.trim_3d_1,
            &eg.trim_3d_2,
            eg.trim_circle_1,
            eg.trim_circle_2,
        )
    } else {
        (
            &eg.trim_3d_2,
            &eg.trim_3d_1,
            eg.trim_circle_2,
            eg.trim_circle_1,
        )
    };
    // Each support's new canonical half-edge replaces that support's
    // rim half-edge: orient the circle to walk the way it did.
    let (old_seam_he, old_fresh_he) = if seam_a.is_some() {
        (he_a_h, he_b_h)
    } else {
        (he_b_h, he_a_h)
    };
    let circle_seam = &closed_fillet_seam::orient_circle_to_rim(brep, old_seam_he, circle_seam);
    let circle_fresh = &closed_fillet_seam::orient_circle_to_rim(brep, old_fresh_he, circle_fresh);
    let trim_curve_on_edge = |trim: &popkern_core::geometry::NurbsCurve,
                              circle: &Option<(Point3, popkern_core::geometry::Vec3, f64)>,
                              start_at: Point3,
                              old_he_axis: Option<popkern_core::geometry::Vec3>|
     -> popkern_core::geometry::CurveOnEdge {
        match circle {
            Some((center, normal, radius)) => {
                // Angle-parameterized analytic arc starting at the
                // edge's vertex. Angles measured in the Arc
                // evaluator's frame (circle_basis) so that
                // point_at(theta0) IS the vertex — see the exact_rim
                // anchor note in chamfer_curve.rs.
                let (e1, e2) = popkern_core::geometry::analytic_curves::circle_basis(*normal);
                let d = start_at - *center;
                let theta0 = d.dot(&e2).patan2(d.dot(&e1));
                if popkern_core::trace::trace_enabled("rim") {
                    eprintln!(
                        "  [rim-anchor] closed_fillet start_at=({:.3},{:.3},{:.3}) theta0={theta0:.4}",
                        start_at.x, start_at.y, start_at.z
                    );
                }
                popkern_core::geometry::CurveOnEdge {
                    geometry: popkern_core::geometry::curve::CurveGeometry::Arc {
                        center: *center,
                        normal: *normal,
                        radius: *radius,
                        start_angle: theta0,
                        end_angle: theta0 + std::f64::consts::TAU,
                    },
                    t_start: theta0,
                    t_end: theta0 + std::f64::consts::TAU,
                }
            }
            None => super::closed_fillet_rail::rail_on_edge(trim, start_at, old_he_axis),
        }
    };
    // ONE anchor authority per edge: each contact-circle curve is
    // anchored at ITS OWN vertex's position (pos_with_seam /
    // v_fresh_pos). The previous code passed pos_for_face_a/b —
    // selected by the analysis-face mapping — while trim_seam /
    // circle_seam are selected by the SEAM side; when the seam sits
    // on side B the two selectors disagree and the seam edge's arc
    // anchors at the OTHER circle's contact point (measured: a rim
    // trim carried the cap circle with the wall's start angle; the
    // principal-range arc→NURBS clamp collapsed both the vertex and
    // the crossed anchor to angle 0, hiding the crossing until the
    // exact conversion landed and it surfaced as a 0.25 mm
    // OpenLoop).
    let axis = |he| closed_fillet_seam::half_edge_winding_axis(brep, he);
    let seam_curve = trim_curve_on_edge(trim_seam, circle_seam, pos_with_seam, axis(old_seam_he));
    let fresh_curve = trim_curve_on_edge(trim_fresh, circle_fresh, v_fresh_pos, axis(old_fresh_he));
    let edge_seam_side = brep.insert_edge_tracked(Edge {
        half_edge: HalfEdgeHandle::default(),
        curve: Some(seam_curve),
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: Default::default(),
    });
    let edge_fresh_side = brep.insert_edge_tracked(Edge {
        half_edge: HalfEdgeHandle::default(),
        curve: Some(fresh_curve),
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: Default::default(),
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
        role: Default::default(),
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
    let support_fresh_old_he = if seam_a.is_some() { he_b_h } else { he_a_h };

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
    // Both supports seamed (a revolved groove's torus beside its wall
    // cylinder, PK-433): the fresh vertex takes over the second seam.
    if let Some(seam) = seam_a.and(seam_b) {
        let (rim, fresh) = (support_fresh_old_he, he_fresh_support);
        closed_fillet_seam::splice_fresh_seam_side(brep, rim, fresh, seam, v_old, v_fresh)
            .map_err(|reason| BlendFailure::InvalidTopology {
                reason,
                location: edge_loc,
            })?;
    }
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
    brep.vertices[v_with_seam].outgoing = he_seam_support;

    // 13. Delete the old rim half-edges and edge.
    brep.remove_half_edge_reanchored(he_a_h);
    brep.remove_half_edge_reanchored(he_b_h);
    brep.edges.remove(fillet_edge);

    // 13b. Pcurves for the four trim half-edges, angle-matched to
    //     the Arc edge curves: without them the STEP exporter's
    //     pcurve synthesis disagrees with the exact Arc and swaps it
    //     for a resampled B-spline (the wall face then imports
    //     broken in OpenCASCADE). Authored by the shared adaptive
    //     chart fit (`fit_chart_pcurve`), parameterized by the SAME
    //     angles as the Arc and mid-knot verified to the linear
    //     tolerance — the earlier fixed 17-sample polyline left the
    //     PLANE-side image with ~1.5e-2 chord sag, and a boolean
    //     trim endpoint solved against it landed 1.8e-5 off the rim
    //     arc, outside the split stage's vertex-insertion band
    //     (PK-292).
    for he in [
        he_seam_support,
        he_seam_fillet,
        he_fresh_support,
        he_fresh_fillet,
    ] {
        super::closed_fillet_rail_pcurve::author_trim_pcurve(brep, he, tolerance);
    }

    // 14. Seam the fillet face: its two boundary circles are
    //     disjoint loops, so the face charts as a full-wrap ring —
    //     unchartable for tessellation/export (the same reason a
    //     primitive cylinder's lateral carries a seam edge). Join
    //     the inner ring to the outer loop (mekr) along the strip's
    //     cross-section at the boundary vertices, which are aligned
    //     by construction: both trim curves start at the walker's
    //     spine start angle.
    {
        let he_outer = brep.loops[loop_fillet_outer].half_edge;
        let he_inner = brep.loops[loop_fillet_inner].half_edge;
        let p_out = brep.vertices[brep.half_edges[he_outer].origin].position;
        let p_in = brep.vertices[brep.half_edges[he_inner].origin].position;
        match brep.mekr(he_outer, he_inner) {
            Ok(mekr) => {
                // Seam curve, exact-first (PK-295): see
                // `closed_fillet_seam::author_strip_seam_curve` for
                // the exact torus tube arc and the UV-interpolated
                // polyline fallback (PK-285 branch alignment).
                if let Some(curve) =
                    closed_fillet_seam::author_strip_seam_curve(&eg.surface, p_out, p_in, tolerance)
                    && let Some(edge) = brep.edges.get_mut(mekr.edge)
                {
                    edge.curve = Some(curve);
                }
            }
            Err(e) => {
                // Structured decline: the strip stays valid (two
                // disjoint rings) but unchartable; tessellation
                // falls back to the vertex-blend path for this face.
                tracing_note_seam_decline(&format!("closed-fillet seam mekr declined: {e}"));
            }
        }
    }

    // PK-241: reconcile the provisional annulus flag (loops now wired).
    super::strip_orientation::reconcile_strip_orientation_boundary(brep, fillet_face);
    Ok(fillet_face)
}

/// Seam-decline observability without a tracing dependency: stderr
/// under the fillet trace env, silent otherwise.
fn tracing_note_seam_decline(msg: &str) {
    if popkern_core::trace::trace_enabled("fillet") {
        eprintln!("  [closed-fillet] {msg}");
    }
}
