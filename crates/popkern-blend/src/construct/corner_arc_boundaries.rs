//! Corner-cap boundary arc reconstruction.
//!
//! The rolling-ball corner cap's true boundary curves are the
//! adjacent fillet strips' END CROSS-SECTIONS: circular arcs of the
//! fillet radius centered at the cap sphere's own center. They are
//! TANGENCY curves — no surface-intersection dispatcher can produce
//! them (the cap is tangent to its strips along them, and tangent
//! configurations are degenerate for intersection) — so they must be
//! reconstructed from the cap sphere itself. The carve phase leaves
//! straight chords in their place (2.5e-1 off-surface at demo
//! scales), which poisons tessellation and export association for
//! every face at the corner.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::BRep;
use popkern_math::PortableMath;

/// For each corner patch whose face is a rolling-ball Sphere,
/// replace chord boundary edges whose endpoints lie (within the
/// anchor-placement class, 5e-3) on the sphere with exact circular
/// arcs centered at the sphere center. Endpoint vertices are snapped
/// radially onto the sphere and adjacent chord-line edges rebuilt.
/// Returns the number of edges rewritten.
pub(super) fn arcify_corner_cap_boundaries(
    brep: &mut BRep,
    corner_patches: &[crate::types::CornerPatchInfo],
) -> usize {
    let mut rewritten = 0;
    for cp in corner_patches {
        let Some(face) = brep.faces.get(cp.face) else {
            continue;
        };
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        let SurfaceGeometry::Sphere { center, radius } = sof.geometry else {
            continue;
        };
        let start_he = brep.loops[face.outer_loop].half_edge;
        // Collect the loop's edges first (mutation below).
        let mut edges = Vec::new();
        let mut he = start_he;
        for _ in 0..brep.half_edges.len() + 1 {
            let d = &brep.half_edges[he];
            edges.push(d.edge);
            he = d.next;
            if he == start_he {
                break;
            }
        }
        for eh in edges {
            let Some(edge) = brep.edges.get(eh) else {
                continue;
            };
            let Some(curve) = edge.curve.as_ref() else {
                continue;
            };
            // Only chords: straight Lines or 2-point degree-1 NURBS.
            let is_chord = match &curve.geometry {
                CurveGeometry::Line { .. } => true,
                CurveGeometry::Nurbs(n) => n.degree == 1 && n.control_points.len() == 2,
                _ => false,
            };
            if !is_chord {
                continue;
            }
            let he0 = edge.half_edge;
            let (v1, v2) = {
                let d = &brep.half_edges[he0];
                (d.origin, brep.half_edges[d.twin].origin)
            };
            if v1 == v2 {
                continue;
            }
            let p1 = brep.vertices[v1].position;
            let p2 = brep.vertices[v2].position;
            let (d1, d2) = (p1.distance_to(&center), p2.distance_to(&center));
            let snap_gate = 5.0e-3;
            if (d1 - radius).abs() > snap_gate || (d2 - radius).abs() > snap_gate {
                continue;
            }
            if d1 < 1e-12 || d2 < 1e-12 {
                continue;
            }
            // Snap endpoints radially onto the sphere.
            let s1 = center + (p1 - center) * (radius / d1);
            let s2 = center + (p2 - center) * (radius / d2);
            crate::construct::edge_reconcile::snap_vertex_and_adjacent_chords(brep, v1, s1, eh);
            crate::construct::edge_reconcile::snap_vertex_and_adjacent_chords(brep, v2, s2, eh);
            // Exact minor arc centered at the sphere center,
            // verification-driven: the Arc evaluator derives its
            // frame from orthonormal_basis(normal), so rather than
            // reproducing that convention analytically (a reflected
            // frame silently traces the wrong sweep — measured as a
            // 2r loop gap on bottom-corner fragments), try both
            // normal directions and accept the candidate whose
            // endpoints round-trip and whose midpoint stays on the
            // sphere.
            let a = s1 - center;
            let b = s2 - center;
            let Some(plane_normal) = a.cross(&b).normalized() else {
                continue;
            };
            let mut accepted: Option<CurveOnEdge> = None;
            for arc_normal in [plane_normal, plane_normal * -1.0] {
                // Angles MUST be measured in the Arc evaluator's own
                // frame (circle_basis), not orthonormal_basis — the
                // two conventions differ and a reflected frame
                // traces the wrong sweep.
                let (f1, f2) = popkern_core::geometry::analytic_curves::circle_basis(arc_normal);
                let theta1 = a.dot(&f2).patan2(a.dot(&f1));
                let mut theta2 = b.dot(&f2).patan2(b.dot(&f1));
                while theta2 <= theta1 {
                    theta2 += std::f64::consts::TAU;
                }
                if theta2 - theta1 > std::f64::consts::PI + 1e-9 {
                    continue;
                }
                let cand = CurveOnEdge {
                    geometry: CurveGeometry::Arc {
                        center,
                        normal: arc_normal,
                        radius,
                        start_angle: theta1,
                        end_angle: theta2,
                    },
                    t_start: theta1,
                    t_end: theta2,
                };
                let ceval: &dyn CurveEval = &cand.geometry;
                let q0 = ceval.point_at(theta1);
                let q1 = ceval.point_at(theta2);
                let qm = ceval.point_at(0.5 * (theta1 + theta2));
                if q0.distance_to(&s1) < 1e-9
                    && q1.distance_to(&s2) < 1e-9
                    && (qm.distance_to(&center) - radius).abs() < 1e-9
                {
                    accepted = Some(cand);
                    break;
                }
            }
            let Some(new_curve) = accepted else {
                if popkern_core::trace::trace_enabled("corner-carve") {
                    eprintln!(
                        "  [arcify] REJECTED {eh:?}: |s1-C|-r {:.2e}, |s2-C|-r {:.2e}, chord {:.4}",
                        s1.distance_to(&center) - radius,
                        s2.distance_to(&center) - radius,
                        s1.distance_to(&s2)
                    );
                }
                continue;
            };
            if let Some(e) = brep.edges.get_mut(eh) {
                e.curve = Some(new_curve);
                rewritten += 1;
                // Stale pcurves would let later fix passes
                // resurrect the chord; recompute re-authors from
                // the arc.
                let he0 = e.half_edge;
                let tw = brep.half_edges[he0].twin;
                brep.half_edges[he0].pcurve = None;
                brep.half_edges[tw].pcurve = None;
            }
        }
    }
    rewritten
}

/// Route B corner surgery: merge each cap-boundary ANCHOR vertex
/// (a perp-split artifact at r·√2 from the cap center, off the
/// sphere) into its adjacent on-sphere CONTACT vertex via `kve`
/// edge contraction. The cap loop runs contact→anchor→contact…;
/// after contraction the surviving edges span contact→contact and
/// [`arcify_corner_cap_boundaries`] can emit exact sphere-centered
/// arcs. Incident straight edges (wall trim lines) are rebuilt to
/// the merged endpoint. Per-corner transactional: applied on a
/// candidate clone, accepted only when the structural invariants
/// hold; a declined corner keeps its chord topology and reports
/// through the trace.
pub(super) fn collapse_corner_anchor_chords(
    brep: &mut BRep,
    corner_patches: &[crate::types::CornerPatchInfo],
) -> usize {
    let mut contracted_total = 0;
    for cp in corner_patches {
        let Some(face) = brep.faces.get(cp.face) else {
            continue;
        };
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        let SurfaceGeometry::Sphere { center, radius } = sof.geometry else {
            continue;
        };
        let on_sphere = |p: Point3| (p.distance_to(&center) - radius).abs() < 1.0e-3;
        if popkern_core::trace::trace_enabled("corner-carve") {
            let start_he = brep.loops[face.outer_loop].half_edge;
            let mut he = start_he;
            let mut desc = Vec::new();
            for _ in 0..brep.half_edges.len() + 1 {
                let d = &brep.half_edges[he];
                let pp = brep.vertices[d.origin].position;
                desc.push(format!("{:.4}", pp.distance_to(&center)));
                he = d.next;
                if he == start_he {
                    break;
                }
            }
            eprintln!(
                "  [contact-collapse] cap {:?} loop |v-C|: {desc:?}",
                cp.face
            );
        }
        let mut candidate = brep.clone();
        let mut contracted = 0;
        // Iterate until no contact→anchor cap edge remains (each
        // contraction restructures the loop, so restart the walk).
        for _ in 0..64 {
            let Some(face) = candidate.faces.get(cp.face) else {
                break;
            };
            let start_he = candidate.loops[face.outer_loop].half_edge;
            let mut target: Option<popkern_core::topology::HalfEdgeHandle> = None;
            let mut he = start_he;
            for _ in 0..candidate.half_edges.len() + 1 {
                let d = &candidate.half_edges[he];
                let v_o = d.origin;
                let v_d = candidate.half_edges[d.twin].origin;
                let p_o = candidate.vertices[v_o].position;
                let p_d = candidate.vertices[v_d].position;
                let edge = d.edge;
                let is_chord = candidate.edges[edge]
                    .curve
                    .as_ref()
                    .map(|c| match &c.geometry {
                        CurveGeometry::Line { .. } => true,
                        CurveGeometry::Nurbs(n) => n.degree == 1 && n.control_points.len() == 2,
                        _ => false,
                    })
                    .unwrap_or(false);
                // Contract only chords from an on-sphere contact to
                // a near-corner anchor (bounded by 2r — fragment
                // chains reach farther and are left alone).
                if is_chord
                    && on_sphere(p_o)
                    && !on_sphere(p_d)
                    && p_d.distance_to(&center) < 2.0 * radius
                {
                    target = Some(he);
                    break;
                }
                he = d.next;
                if he == start_he {
                    break;
                }
            }
            let Some(he) = target else { break };
            match candidate.kve(he) {
                Ok(r) => {
                    rebuild_incident_straight_edges(&mut candidate, r.vertex);
                    contracted += 1;
                }
                Err(e) => {
                    if popkern_core::trace::trace_enabled("corner-carve") {
                        eprintln!("  [contact-collapse] kve declined at {:?}: {e:?}", cp.face);
                    }
                    break;
                }
            }
        }
        if contracted > 0 {
            if let Err(e) = candidate.validate()
                && popkern_core::trace::trace_enabled("corner-carve")
            {
                eprintln!(
                    "  [contact-collapse] corner {:?} candidate invalid: {e:?}",
                    cp.face
                );
                if let popkern_core::error::ValidationError::OpenLoop { loop_, .. } = e {
                    use popkern_core::geometry::curve::CurveEval;
                    if let Some(l) = candidate.loops.get(loop_) {
                        let start = l.half_edge;
                        let mut he = start;
                        for _ in 0..candidate.half_edges.len() + 1 {
                            let d = &candidate.half_edges[he];
                            if let Some(c) = candidate.edges[d.edge].curve.as_ref() {
                                let ceval: &dyn CurveEval = &c.geometry;
                                let p0 = ceval.point_at(c.t_start);
                                let p1 = ceval.point_at(c.t_end);
                                let kind = match &c.geometry {
                                    popkern_core::geometry::CurveGeometry::Line { .. } => "Line",
                                    popkern_core::geometry::CurveGeometry::Arc { .. } => "Arc",
                                    popkern_core::geometry::CurveGeometry::Nurbs(_) => "Nurbs",
                                    _ => "?",
                                };
                                eprintln!(
                                    "      edge {:?} {kind} ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                                    d.edge, p0.x, p0.y, p0.z, p1.x, p1.y, p1.z
                                );
                            }
                            he = d.next;
                            if he == start {
                                break;
                            }
                        }
                    }
                }
            }
            if candidate.validate().is_ok() {
                *brep = candidate;
                contracted_total += contracted;
            } else if popkern_core::trace::trace_enabled("corner-carve") {
                eprintln!(
                    "  [contact-collapse] corner {:?}: {contracted} contractions ROLLED BACK (validate failed)",
                    cp.face
                );
            }
        }
    }
    contracted_total
}

/// Rebuild every straight incident edge of `vertex` between its
/// current endpoint positions (post-contraction geometry repair —
/// the kve caller contract).
pub(super) fn rebuild_incident_straight_edges(
    brep: &mut BRep,
    vertex: popkern_core::topology::VertexHandle,
) {
    let incident: Vec<_> = brep
        .edges
        .iter()
        .filter(|(_, e)| {
            let he = e.half_edge;
            brep.half_edges.get(he).is_some_and(|h| h.origin == vertex)
                || brep
                    .half_edges
                    .get(he)
                    .and_then(|h| brep.half_edges.get(h.twin))
                    .is_some_and(|t| t.origin == vertex)
        })
        .map(|(h, _)| h)
        .collect();
    for eh in incident {
        let he = brep.edges[eh].half_edge;
        let a = brep.vertices[brep.half_edges[he].origin].position;
        let b = brep.vertices[brep.half_edges[brep.half_edges[he].twin].origin].position;
        let Some(curve) = brep.edges[eh].curve.clone() else {
            continue;
        };
        let new_curve = match &curve.geometry {
            CurveGeometry::Line { .. } => {
                let d = b - a;
                let len = d.length();
                let Some(dir) = (len > 1e-12).then(|| d.normalized()).flatten() else {
                    continue;
                };
                Some(CurveOnEdge {
                    geometry: CurveGeometry::Line {
                        origin: a,
                        direction: dir,
                    },
                    t_start: 0.0,
                    t_end: len,
                })
            }
            CurveGeometry::Nurbs(n) if n.degree == 1 && n.control_points.len() == 2 => {
                let d = b - a;
                let len = d.length();
                let Some(dir) = (len > 1e-12).then(|| d.normalized()).flatten() else {
                    continue;
                };
                Some(CurveOnEdge {
                    geometry: CurveGeometry::Line {
                        origin: a,
                        direction: dir,
                    },
                    t_start: 0.0,
                    t_end: len,
                })
            }
            CurveGeometry::Nurbs(n) => {
                // Endpoint-CP snap: move the control point nearest
                // the merged vertex onto it.
                use popkern_core::geometry::curve::CurveEval;
                let ceval: &dyn CurveEval = &curve.geometry;
                let start = ceval.point_at(curve.t_start);
                let end = ceval.point_at(curve.t_end);
                let vp = brep.vertices[vertex].position;
                let snap_first = start.distance_to(&vp) <= end.distance_to(&vp);
                let mut n2 = n.clone();
                if snap_first {
                    if let Some(cp) = n2.control_points.first_mut() {
                        *cp = vp;
                    }
                } else if let Some(cp) = n2.control_points.last_mut() {
                    *cp = vp;
                }
                Some(CurveOnEdge {
                    geometry: CurveGeometry::Nurbs(n2),
                    t_start: curve.t_start,
                    t_end: curve.t_end,
                })
            }
            CurveGeometry::Arc {
                center,
                normal,
                radius,
                ..
            } => {
                // Re-trim the arc between the (possibly moved)
                // endpoint angles, keeping center/normal/radius.
                let (f1, f2) =
                    popkern_core::geometry::analytic_surfaces::orthonormal_basis(*normal);
                let ang = |p: Point3| {
                    let d = p - *center;
                    d.dot(&f2).patan2(d.dot(&f1))
                };
                let ta = ang(a);
                let mut tb = ang(b);
                while tb <= ta {
                    tb += std::f64::consts::TAU;
                }
                // Keep the shorter sweep (the corner arcs are minor).
                if tb - ta > std::f64::consts::PI {
                    continue;
                }
                Some(CurveOnEdge {
                    geometry: CurveGeometry::Arc {
                        center: *center,
                        normal: *normal,
                        radius: *radius,
                        start_angle: ta,
                        end_angle: tb,
                    },
                    t_start: ta,
                    t_end: tb,
                })
            }
            _ => None,
        };
        if let Some(nc) = new_curve {
            brep.edges[eh].curve = Some(nc);
        }
        // Pcurves are stale for EVERY incident edge after the
        // endpoint move (later pipeline passes regenerate 3D curves
        // from pcurves and would resurrect the pre-merge position);
        // drop them so the global recompute re-authors from the new
        // 3D curves.
        let he0 = brep.edges[eh].half_edge;
        let tw = brep.half_edges[he0].twin;
        brep.half_edges[he0].pcurve = None;
        brep.half_edges[tw].pcurve = None;
    }
}
