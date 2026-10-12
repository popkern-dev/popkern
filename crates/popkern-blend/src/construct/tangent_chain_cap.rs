//! Tangent-cap ends of open tangent chains (PK-735,
//! `docs/architecture/fillet-open-tangent-chains.md`): the rim runs on
//! past the chain's end vertex V, tangent and unfilleted, so the strip
//! stops on its exact section at V and a planar CAP face in the section
//! plane closes the cut. The cap is bounded by the section arc X (P → Q)
//! and, on each support side, the edge from V to that side's contact:
//! the upper piece of the crease ("third edge") through V on that side,
//! split at its contact, or a new straight line on the support.
//!
//! [`prepare_cap`] solves the section and creates the cap's vertices,
//! edges and face before the strips are wired; [`wire_cap`] patches the
//! two support loops at V once the end strip's support half-edges have
//! replaced the chain edge, and closes the cap loop.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, EdgeHandle, Face, FaceHandle, HalfEdgeHandle, Loop, LoopHandle, ShellHandle, Vertex,
    VertexHandle,
};

use super::g1_junction::tangent_at_vertex;
use super::tangent_chain_geometry::{
    Link, ball_center, blank, chain_tangent, chain_tangent_at_origin, decline, exact_section,
    insert_edge, nearest, section_arc, segment_distance, surface_foot,
};
use crate::error::BlendFailure;

/// How a free end of an open tangent chain ends.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ChainEnd {
    /// Both rails end on the end vertex (PK-764): the strip runs out
    /// to zero width there.
    Vanishing,
    /// The rim continues tangent and unfilleted past the end vertex:
    /// the strip stops on its section, closed by a planar cap.
    TangentCap,
    /// The end vertex is a sharp corner: the strip runs out on the
    /// planar end face holding the corner's other two edges.
    EndFace,
}

/// One side of a cap: the edge from V to that side's contact.
#[derive(Clone, Copy, Debug)]
pub(super) enum CapSide {
    /// The upper piece of the crease through V, split at the contact.
    Third(EdgeHandle),
    /// A new line on the support; its two half-edges `(V → c, c → V)`.
    Line(HalfEdgeHandle, HalfEdgeHandle),
}

/// A prepared tangent cap at one end of an open chain.
#[derive(Clone, Copy, Debug)]
pub(super) struct Cap {
    /// The chain's end vertex (kept).
    pub(super) v: VertexHandle,
    /// The section contacts on the P and Q supports.
    pub(super) vp: VertexHandle,
    pub(super) vq: VertexHandle,
    /// The section arc's half-edges `P → Q` and `Q → P`.
    pub(super) x_pq: HalfEdgeHandle,
    pub(super) x_qp: HalfEdgeHandle,
    /// The section ball's centre.
    pub(super) center: Point3,
    /// The cap face: new, or the existing end face.
    pub(super) face: FaceHandle,
    /// Whether the strip ends on an existing END FACE: the corner
    /// vertex and both corner pieces go, the section takes their place
    /// in the end face's loop.
    pub(super) end_face: bool,
    pub(super) at_start: bool,
    pub(super) side_p: CapSide,
    pub(super) side_q: CapSide,
}

impl Cap {
    /// The section arc's half-edge on the strip's side.
    pub(super) fn strip_x(&self) -> HalfEdgeHandle {
        if self.at_start { self.x_pq } else { self.x_qp }
    }

    /// The edges the cap created or split (their faces need pcurves).
    pub(super) fn edges(&self, brep: &BRep) -> Vec<EdgeHandle> {
        let mut out = vec![brep.half_edges[self.x_pq].edge];
        for side in [self.side_p, self.side_q] {
            out.push(match side {
                CapSide::Third(e) => e,
                CapSide::Line(h, _) => brep.half_edges[h].edge,
            });
        }
        out
    }
}

fn fresh_vertex(brep: &mut BRep, position: Point3) -> VertexHandle {
    brep.vertices.insert(Vertex {
        position,
        outgoing: HalfEdgeHandle::default(),
        name: None,
        tolerance: None,
    })
}

/// A new straight edge on a support from `v` to a new contact vertex
/// at `c`, as a cap side: half-edges `(V → c, c → V)`, unlinked.
pub(super) fn support_line(brep: &mut BRep, v: VertexHandle, c: Point3) -> (VertexHandle, CapSide) {
    let vpos = brep.vertices[v].position;
    let vc = fresh_vertex(brep, c);
    let edge = insert_edge(
        brep,
        CurveOnEdge {
            geometry: CurveGeometry::Line {
                origin: vpos,
                direction: c - vpos,
            },
            t_start: 0.0,
            t_end: 1.0,
        },
    );
    let out = brep.half_edges.insert(blank(v, edge));
    let back = brep.half_edges.insert(blank(vc, edge));
    brep.half_edges[out].twin = back;
    brep.half_edges[back].twin = out;
    brep.edges[edge].half_edge = out;
    (vc, CapSide::Line(out, back))
}

fn support(brep: &BRep, he: HalfEdgeHandle) -> Option<SurfaceGeometry> {
    brep.faces[brep.loops[brep.half_edges[he].loop_].face]
        .surface
        .as_ref()
        .map(|s| s.geometry.clone())
}

/// Solve the section at the chain end (`at_start`: the origin of
/// `link.h`, else its destination), split the crease through it at its
/// contact and create the cap's contacts, edges and face. The support
/// loops are patched later by [`wire_cap`].
///
/// # Errors
///
/// A typed decline when the end is not a tangent cap (valence above 3,
/// no tangent continuation), the section has no ball, the crease
/// misses its contact, or a cap line would leave a curved support.
pub(super) fn prepare_cap(
    brep: &mut BRep,
    link: &Link<'_>,
    at_start: bool,
    shell: ShellHandle,
    band: f64,
    tolerance: &Tolerance,
) -> Result<Cap, BlendFailure> {
    let h = link.h;
    let tw = brep.half_edges[h].twin;
    let v = if at_start {
        brep.half_edges[h].origin
    } else {
        brep.half_edges[tw].origin
    };
    let vpos = brep.vertices[v].position;
    let e = brep.half_edges[h].edge;
    // V's neighbours on the P side (h's face) and the Q side.
    let (np, nq) = if at_start {
        (brep.half_edges[h].prev, brep.half_edges[tw].next)
    } else {
        (brep.half_edges[h].next, brep.half_edges[tw].prev)
    };
    let (a, b) = (brep.half_edges[np].edge, brep.half_edges[nq].edge);
    let valence = brep.half_edges.values().filter(|he| he.origin == v).count();
    let third = if a == b && valence == 2 && tangent_at_vertex(brep, e, a, v) {
        None
    } else if a != b && valence == 3 && tangent_at_vertex(brep, e, a, v) {
        Some((b, false))
    } else if a != b && valence == 3 && tangent_at_vertex(brep, e, b, v) {
        Some((a, true))
    } else {
        return Err(decline(
            format!("chain end {v:?} (valence {valence}) does not continue tangent"),
            vpos,
        ));
    };
    let tangent = if at_start {
        chain_tangent_at_origin(brep, h)
    } else {
        chain_tangent(brep, h)
    }
    .ok_or_else(|| decline(format!("chain end {v:?} has no tangent"), vpos))?;
    let (Some(sp), Some(sq)) = (support(brep, h), support(brep, tw)) else {
        return Err(decline(format!("chain end {v:?} lacks supports"), vpos));
    };
    // The radius at V along the edge's own direction.
    let forward = brep.half_edges[brep.edges[e].half_edge].origin == v;
    let r = if forward {
        link.eg.radius_at_start
    } else {
        link.eg.radius_at_end
    };
    let cp0 = nearest(link.p_trim, vpos);
    let cq0 = nearest(link.q_trim, vpos);
    let seed = ball_center(&sp, cp0, cq0, tangent, r).unwrap_or(vpos);
    let Some((center, cp, cq)) = exact_section(&sp, &sq, vpos, tangent, r, seed) else {
        return Err(decline(
            format!("no ball of radius {r} sits on both supports at chain end {v:?}"),
            vpos,
        ));
    };
    // Each side's contact vertex and its edge from V.
    let side = |brep: &mut BRep, on_p: bool, c: Point3| -> Result<_, BlendFailure> {
        if let Some((t, t_on_p)) = third
            && t_on_p == on_p
        {
            let off = match brep.edges[t].curve.as_ref() {
                Some(curve) => curve
                    .geometry
                    .closest_point(&c)
                    .map(|s| curve.geometry.point_at(s).distance_to(&c))
                    .unwrap_or(f64::INFINITY),
                None => segment_distance(brep, t, c),
            };
            if off > band {
                return Err(decline(
                    format!("crease {t:?} misses the cap section contact by {off:.3e}"),
                    vpos,
                ));
            }
            let split = brep
                .split_edge_at_point(t, c, tolerance)
                .map_err(|err| decline(format!("crease {t:?} split at its contact: {err}"), c))?;
            let ends = |brep: &BRep, edge: EdgeHandle| {
                let he = &brep.half_edges[brep.edges[edge].half_edge];
                [he.origin, brep.half_edges[he.twin].origin]
            };
            let upper = if ends(brep, split.edge_a).contains(&v) {
                split.edge_a
            } else {
                split.edge_b
            };
            return Ok((split.vertex, CapSide::Third(upper)));
        }
        let surface = if on_p { &sp } else { &sq };
        let mid = vpos + (c - vpos) * 0.5;
        let off = surface_foot(surface, mid).map_or(f64::INFINITY, |f| f.distance_to(&mid));
        if off > band {
            return Err(decline(
                format!("cap line at chain end {v:?} leaves its curved support by {off:.3e}"),
                vpos,
            ));
        }
        Ok(support_line(brep, v, c))
    };
    let (vp, side_p) = side(brep, true, cp)?;
    let (vq, side_q) = side(brep, false, cq)?;
    let (pp, pq) = (brep.vertices[vp].position, brep.vertices[vq].position);
    let arc = section_arc(pp, pq, center, tangent, r).ok_or_else(|| {
        decline(
            format!("no cap section arc of radius {r} at chain end {v:?}"),
            vpos,
        )
    })?;
    let x = insert_edge(brep, arc);
    let x_pq = brep.half_edges.insert(blank(vp, x));
    let x_qp = brep.half_edges.insert(blank(vq, x));
    brep.half_edges[x_pq].twin = x_qp;
    brep.half_edges[x_qp].twin = x_pq;
    brep.edges[x].half_edge = x_pq;
    let n = tangent
        .normalized()
        .ok_or_else(|| decline(format!("chain end {v:?} has no tangent"), vpos))?;
    let (u_axis, v_axis) = popkern_core::geometry::analytic_curves::circle_basis(n);
    let face = brep.faces.insert(Face {
        outer_loop: LoopHandle::default(),
        inner_loops: vec![],
        shell,
        surface: Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: vpos,
                normal: n,
                u_axis,
                v_axis,
            },
            orientation: true,
        }),
        name: None,
        tolerance: None,
        role: Default::default(),
    });
    Ok(Cap {
        v,
        vp,
        vq,
        x_pq,
        x_qp,
        center,
        face,
        end_face: false,
        at_start,
        side_p,
        side_q,
    })
}

/// Patch the support loops at the cap's V and close the cap loop.
/// `s_p` and `s_q` are the end strip's support half-edges (on P and Q),
/// already substituted for the chain edge.
///
/// # Errors
///
/// A typed decline when a support loop does not run through the split
/// crease as the preparation found it.
pub(super) fn wire_cap(
    brep: &mut BRep,
    cap: &Cap,
    s_p: HalfEdgeHandle,
    s_q: HalfEdgeHandle,
) -> Result<(), BlendFailure> {
    let at = brep.vertices[cap.v].position;
    let cap_x = if cap.at_start { cap.x_qp } else { cap.x_pq };
    let mut ring = vec![cap_x];
    for (side, s, on_p) in [(cap.side_p, s_p, true), (cap.side_q, s_q, false)] {
        // The support half-edge ends at the contact when it runs into
        // the cap end: P at the chain's end, Q at its start.
        let ends_at_contact = on_p != cap.at_start;
        let cap_he = match side {
            CapSide::Third(upper) => {
                let removed = if ends_at_contact {
                    brep.half_edges[s].next
                } else {
                    brep.half_edges[s].prev
                };
                if brep.half_edges[removed].edge != upper {
                    return Err(decline(
                        "support loop does not run through the split crease".into(),
                        at,
                    ));
                }
                let (pr, nx) = (brep.half_edges[removed].prev, brep.half_edges[removed].next);
                brep.half_edges[pr].next = nx;
                brep.half_edges[nx].prev = pr;
                let lp = brep.half_edges[removed].loop_;
                if brep.loops[lp].half_edge == removed {
                    brep.loops[lp].half_edge = nx;
                }
                // It changes faces: its pcurve was authored in the
                // support's chart (kept, it read the cap's winding
                // backwards — measured: uv area 9.3 for a 0.21 cap).
                brep.half_edges[removed].pcurve = None;
                removed
            }
            CapSide::Line(out, back) => {
                let lp = brep.half_edges[s].loop_;
                let (support_he, cap_he) = if ends_at_contact {
                    (back, out)
                } else {
                    (out, back)
                };
                let (pr, nx) = if ends_at_contact {
                    (s, brep.half_edges[s].next)
                } else {
                    (brep.half_edges[s].prev, s)
                };
                brep.half_edges[pr].next = support_he;
                brep.half_edges[support_he].prev = pr;
                brep.half_edges[support_he].next = nx;
                brep.half_edges[nx].prev = support_he;
                brep.half_edges[support_he].loop_ = lp;
                cap_he
            }
        };
        ring.push(cap_he);
    }
    if cap.end_face {
        return super::tangent_chain_face_end::close_end_face(brep, cap, cap_x, [ring[1], ring[2]]);
    }
    // Order the three half-edges head to tail.
    let dest = |brep: &BRep, he: HalfEdgeHandle| brep.half_edges[brep.half_edges[he].twin].origin;
    let first = ring[0];
    let mut ordered = vec![first];
    while ordered.len() < ring.len() {
        let tail = dest(brep, ordered[ordered.len() - 1]);
        let Some(&next) = ring
            .iter()
            .find(|&&he| !ordered.contains(&he) && brep.half_edges[he].origin == tail)
        else {
            return Err(decline("cap loop does not close".into(), at));
        };
        ordered.push(next);
    }
    if dest(brep, ordered[ordered.len() - 1]) != brep.half_edges[first].origin {
        return Err(decline("cap loop does not close".into(), at));
    }
    let lp = brep.loops.insert(Loop {
        half_edge: first,
        face: cap.face,
    });
    brep.faces[cap.face].outer_loop = lp;
    for (j, &he) in ordered.iter().enumerate() {
        let nx = ordered[(j + 1) % ordered.len()];
        brep.half_edges[he].next = nx;
        brep.half_edges[nx].prev = he;
        brep.half_edges[he].loop_ = lp;
    }
    orient_cap(brep, cap, &ordered);
    Ok(())
}

/// Orient the cap plane with its loop: the loop winds counter-clockwise
/// about the outward normal (Newell normal of the loop's points, the
/// section arc sampled at its middle).
fn orient_cap(brep: &mut BRep, cap: &Cap, ordered: &[HalfEdgeHandle]) {
    let mut points: Vec<Point3> = Vec::new();
    for &he in ordered {
        points.push(brep.vertices[brep.half_edges[he].origin].position);
        let edge = brep.half_edges[he].edge;
        if edge == brep.half_edges[cap.x_pq].edge
            && let Some(curve) = brep.edges[edge].curve.as_ref()
        {
            points.push(curve.geometry.point_at(0.5 * (curve.t_start + curve.t_end)));
        }
    }
    let mut newell = Vec3::new(0.0, 0.0, 0.0);
    for (j, p) in points.iter().enumerate() {
        let q = points[(j + 1) % points.len()];
        newell = newell + Vec3::new(p.x, p.y, p.z).cross(&Vec3::new(q.x, q.y, q.z));
    }
    if let Some(sof) = brep.faces[cap.face].surface.as_mut()
        && let SurfaceGeometry::Plane { normal, .. } = &sof.geometry
    {
        sof.orientation = normal.dot(&newell) > 0.0;
    }
}

/// Postcondition on a new cap: its orientation flag agrees with its
/// loop's winding in the plane chart (needs the cap's pcurves). The
/// serving boundary would reconcile a dishonest flag silently; the cap
/// is constructed honest, and a contradiction is a wiring defect.
///
/// # Errors
///
/// A typed decline when the winding contradicts the flag or cannot be
/// measured.
pub(super) fn check_cap_winding(brep: &BRep, cap: &Cap) -> Result<(), BlendFailure> {
    if cap.end_face {
        return Ok(());
    }
    let at = brep.vertices[cap.vp].position;
    let face = &brep.faces[cap.face];
    let Some(sof) = face.surface.as_ref() else {
        return Err(decline("cap without a surface".into(), at));
    };
    let area = brep
        .loop_signed_uv_area(face.outer_loop, &sof.geometry)
        .map_err(|reason| decline(format!("cap winding unmeasurable: {reason}"), at))?;
    if (area.signed > 0.0) != sof.orientation {
        return Err(decline(
            "cap orientation contradicts its loop winding".into(),
            at,
        ));
    }
    Ok(())
}
