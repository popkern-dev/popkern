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
pub(crate) fn apply_closed_fillet(
    brep: &mut BRep,
    eg: &EdgeGeometry,
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

    // 7. Create the two new contact-circle edges. Closed-fillet
    //    contact curves are exact circles (rolling-ball contacts of
    //    a circular rim); carry them as analytic Arc geometry — the
    //    NURBS fit through the footprints has ~3e-3 interpolation
    //    error that fails every downstream surface-association gate
    //    (tessellation, export), and the angle parameterization
    //    pairs exactly with iso-linear pcurves on the supports.
    //    Fit failure falls back to the interpolated trim.
    if std::env::var_os("PK_FILLET_TRACE").is_some() {
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
            &eg.trim_circle_1,
            &eg.trim_circle_2,
        )
    } else {
        (
            &eg.trim_3d_2,
            &eg.trim_3d_1,
            &eg.trim_circle_2,
            &eg.trim_circle_1,
        )
    };
    let trim_curve_on_edge = |trim: &popkern_core::geometry::NurbsCurve,
                              circle: &Option<(Point3, popkern_core::geometry::Vec3, f64)>,
                              start_at: Point3|
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
                let theta0 = d.dot(&e2).atan2(d.dot(&e1));
                if std::env::var_os("PK_RIM_TRACE").is_some() {
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
            None => popkern_core::geometry::CurveOnEdge {
                geometry: popkern_core::geometry::curve::CurveGeometry::Nurbs(trim.clone()),
                t_start: 0.0,
                t_end: 1.0,
            },
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
    let edge_seam_side = brep.insert_edge_tracked(Edge {
        half_edge: HalfEdgeHandle::default(),
        curve: Some(trim_curve_on_edge(trim_seam, circle_seam, pos_with_seam)),
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: Default::default(),
    });
    let edge_fresh_side = brep.insert_edge_tracked(Edge {
        half_edge: HalfEdgeHandle::default(),
        curve: Some(trim_curve_on_edge(trim_fresh, circle_fresh, v_fresh_pos)),
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

    // 13b. Pcurves for the four trim half-edges, angle-matched to
    //     the Arc edge curves: without them the STEP exporter's
    //     pcurve synthesis disagrees with the exact Arc and swaps it
    //     for a resampled B-spline (the wall face then imports
    //     broken in OpenCASCADE). Analytic inversion per sample; a
    //     degree-1 pcurve through 17 samples in each face's own
    //     chart, parameterized by the SAME angles as the Arc.
    {
        use popkern_core::geometry::curve::CurveEval;
        let author = |brep: &mut BRep, he: popkern_core::topology::HalfEdgeHandle| {
            let (edge_h, face_h) = {
                let d = &brep.half_edges[he];
                (d.edge, brep.loops[d.loop_].face)
            };
            let Some(curve) = brep.edges[edge_h].curve.clone() else {
                return;
            };
            let Some(sof) = brep.faces[face_h].surface.clone() else {
                return;
            };
            let ceval: &dyn CurveEval = &curve.geometry;
            let n = 17usize;
            let mut uv_pts = Vec::with_capacity(n);
            let mut prev_uv: Option<(f64, f64)> = None;
            for i in 0..n {
                let t = curve.t_start + (curve.t_end - curve.t_start) * i as f64 / (n - 1) as f64;
                let p = ceval.point_at(t);
                let uv = popkern_core::geometry::uv_projection::project_point_to_uv_analytic(
                    &sof.geometry,
                    &p,
                    prev_uv,
                );
                let Some(uv) = uv else { return };
                prev_uv = Some(uv);
                uv_pts.push(Point3::new(uv.0, uv.1, 0.0));
            }
            let mut knots = Vec::with_capacity(n + 2);
            knots.push(curve.t_start);
            for i in 0..n {
                knots.push(
                    curve.t_start + (curve.t_end - curve.t_start) * i as f64 / (n - 1) as f64,
                );
            }
            knots.push(curve.t_end);
            let pc = popkern_core::geometry::NurbsCurve {
                degree: 1,
                knots,
                control_points: uv_pts,
                weights: vec![1.0; n],
            };
            brep.half_edges[he].pcurve = Some(pc.into());
        };
        for he in [
            he_seam_support,
            he_seam_fillet,
            he_fresh_support,
            he_fresh_fillet,
        ] {
            author(brep, he);
        }
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
                // Seam curve: the strip's cross-section arc between
                // the two boundary vertices — a 17-point polyline
                // sampled from the FILLET SURFACE (chord sagitta
                // ~1e-7 at demo scales; exact-arc upgrade rides the
                // matched-pcurve follow-up recorded in trim.rs).
                let surface = &eg.surface;
                use popkern_core::geometry::surface::SurfaceEval;
                let uv_of = |p: &Point3| -> Option<(f64, f64)> {
                    popkern_core::geometry::uv_projection::project_point_to_uv(surface, p)
                };
                if let (Some(uv_a), Some(mut uv_b)) = (uv_of(&p_out), uv_of(&p_in)) {
                    // The two projections can land on DIFFERENT
                    // period branches (one endpoint at u = 2π, the
                    // other at u = 0): interpolating across that gap
                    // sweeps the "cross-section" a full lap around
                    // the chart, and the seam's pcurves become a
                    // diagonal through every u (measured: OCC flagged
                    // the strip wire self-intersecting). Align u_b to
                    // u_a's branch — the seam is a single station.
                    if let Some(period) = surface.parameter_periodicity().0 {
                        uv_b.0 -= period * ((uv_b.0 - uv_a.0) / period).round();
                    }
                    let n = 17usize;
                    let pts: Vec<Point3> = (0..n)
                        .map(|i| {
                            let f = i as f64 / (n - 1) as f64;
                            let u = uv_a.0 + (uv_b.0 - uv_a.0) * f;
                            let v = uv_a.1 + (uv_b.1 - uv_a.1) * f;
                            surface.point_at(u, v)
                        })
                        .collect();
                    // Endpoints exactly at the vertices.
                    let mut pts = pts;
                    pts[0] = p_out;
                    pts[n - 1] = p_in;
                    let mut knots = Vec::with_capacity(n + 2);
                    knots.push(0.0);
                    for i in 0..n {
                        knots.push(i as f64 / (n - 1) as f64);
                    }
                    knots.push(1.0);
                    if let Some(edge) = brep.edges.get_mut(mekr.edge) {
                        edge.curve = Some(popkern_core::geometry::CurveOnEdge {
                            geometry: popkern_core::geometry::CurveGeometry::Nurbs(
                                popkern_core::geometry::NurbsCurve {
                                    degree: 1,
                                    knots,
                                    control_points: pts,
                                    weights: vec![1.0; n],
                                },
                            ),
                            t_start: 0.0,
                            t_end: 1.0,
                        });
                    }
                }
            }
            Err(e) => {
                // Structured decline: the strip stays valid (two
                // disjoint rings) but unchartable; tessellation
                // falls back to the legacy path for this face.
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
    if std::env::var_os("PK_FILLET_TRACE").is_some() {
        eprintln!("  [closed-fillet] {msg}");
    }
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
