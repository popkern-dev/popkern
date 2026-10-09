//! End-face ends of open tangent chains (PK-735,
//! `docs/architecture/fillet-open-tangent-chains.md`): the chain's end
//! vertex V is a sharp corner, and the strip runs out on the planar END
//! FACE F that holds V's other two edges `a` (P ∩ F) and `b` (Q ∩ F).
//! The strip's end is its section by F (exact circle or ellipse for a
//! cylinder strip, a verified fitted section otherwise) between the
//! rails' crossings of F; `a` and `b` are split there and the corner —
//! V and the two corner pieces — is cut away, the section taking its
//! place in F's loop.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle};

use super::tangent_chain_cap::{Cap, CapSide};
use super::tangent_chain_geometry::{
    Circle, Link, ball_center, blank, chain_tangent, chain_tangent_at_origin, decline,
    exact_section, insert_edge,
};
use crate::error::BlendFailure;

/// The strip axis must make at least this cosine with F's normal for
/// a closed section (as the per-edge strip ends, `strip_end_arcs`).
const OBLIQUE_SECTION_MIN_COS: f64 = 0.05;

/// Where the rail `curve` crosses the plane `(origin, n)`, nearest
/// `near`: bisection on its sign changes, or an end lying on the plane
/// within `band` (a rail ending square on a perpendicular end face).
fn rail_crossing(
    curve: &NurbsCurve,
    origin: Point3,
    n: Vec3,
    near: Point3,
    band: f64,
) -> Option<Point3> {
    let side = |p: Point3| (p - origin).dot(&n);
    let (lo, hi) = curve.parameter_range();
    const SAMPLES: usize = 64;
    let at = |i: usize| lo + (hi - lo) * i as f64 / SAMPLES as f64;
    let mut best: Option<Point3> = [curve.point_at(lo), curve.point_at(hi)]
        .into_iter()
        .filter(|p| side(*p).abs() <= band)
        .min_by(|p, q| p.distance_to(&near).total_cmp(&q.distance_to(&near)));
    for i in 0..SAMPLES {
        let (mut a, mut b) = (at(i), at(i + 1));
        let (fa, fb) = (side(curve.point_at(a)), side(curve.point_at(b)));
        if fa * fb > 0.0 {
            continue;
        }
        for _ in 0..80 {
            let m = 0.5 * (a + b);
            if side(curve.point_at(m)) * fa > 0.0 {
                a = m;
            } else {
                b = m;
            }
        }
        let p = curve.point_at(0.5 * (a + b));
        if best.is_none_or(|q| p.distance_to(&near) < q.distance_to(&near)) {
            best = Some(p);
        }
    }
    best
}

/// The point of `edge`'s curve nearest `p`, with its distance.
fn on_edge(brep: &BRep, edge: EdgeHandle, p: Point3) -> Option<(Point3, f64)> {
    let curve = brep.edges[edge].curve.as_ref()?;
    let s = curve.geometry.closest_point(&p)?;
    let q = curve.geometry.point_at(s);
    Some((q, q.distance_to(&p)))
}

/// Solve the strip's end on F at the chain end (`at_start`: the origin
/// of `link.h`) and split F's two corner edges at the rail crossings.
/// The corner is removed later, by [`close_end_face`] through
/// `wire_cap`.
///
/// # Errors
///
/// A typed decline when V is not a corner of one planar face (valence
/// 3, both other edges on F), F is near-parallel to the strip, a rail
/// does not cross F on its edge, or no verified section exists.
pub(super) fn prepare_face_end(
    brep: &mut BRep,
    link: &Link<'_>,
    at_start: bool,
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
    let (np, nq) = if at_start {
        (brep.half_edges[h].prev, brep.half_edges[tw].next)
    } else {
        (brep.half_edges[h].next, brep.half_edges[tw].prev)
    };
    let (a, b) = (brep.half_edges[np].edge, brep.half_edges[nq].edge);
    let valence = brep.half_edges.values().filter(|he| he.origin == v).count();
    let (tp, tq) = (brep.half_edges[np].twin, brep.half_edges[nq].twin);
    let face_of = |he: HalfEdgeHandle| brep.loops[brep.half_edges[he].loop_].face;
    let f = face_of(tp);
    let adjacent = if at_start {
        brep.half_edges[tq].next == tp
    } else {
        brep.half_edges[tp].next == tq
    };
    if valence != 3 || a == b || face_of(tq) != f || !adjacent {
        return Err(decline(
            format!("chain end {v:?} (valence {valence}) is not a corner of one end face"),
            vpos,
        ));
    }
    let Some(SurfaceGeometry::Plane { origin, normal, .. }) =
        brep.faces[f].surface.as_ref().map(|s| s.geometry.clone())
    else {
        return Err(decline(
            format!("chain end {v:?} runs out on a curved end face"),
            vpos,
        ));
    };
    let n = normal
        .normalized()
        .ok_or_else(|| decline("end face without a normal".into(), vpos))?;
    let tangent = if at_start {
        chain_tangent_at_origin(brep, h)
    } else {
        chain_tangent(brep, h)
    }
    .ok_or_else(|| decline(format!("chain end {v:?} has no tangent"), vpos))?;
    // Each rail's crossing of F, landed on F's edge on that side.
    let strip = &link.eg.surface;
    let crossing = |trim: (&NurbsCurve, Circle),
                    edge: EdgeHandle|
     -> Result<Point3, BlendFailure> {
        let raw = rail_crossing(trim.0, origin, n, vpos, band)
            .ok_or_else(|| decline(format!("a rail does not cross the end face at {v:?}"), vpos))?;
        let (on, off) = on_edge(brep, edge, raw)
            .ok_or_else(|| decline(format!("end-face edge {edge:?} has no curve"), vpos))?;
        if off > band {
            return Err(decline(
                format!("rail crosses the end face {off:.3e} off its edge {edge:?}"),
                vpos,
            ));
        }
        Ok(on)
    };
    let cp = crossing(link.p_trim, a)?;
    let cq = crossing(link.q_trim, b)?;
    // The section by F from the P crossing to the Q crossing.
    let (curve, cp, cq) = match strip {
        SurfaceGeometry::Cylinder {
            origin: axis_origin,
            axis,
            radius,
        } => {
            let denom = axis.dot(&n);
            if denom.abs() < OBLIQUE_SECTION_MIN_COS {
                return Err(decline(
                    format!("end face at {v:?} runs near-parallel to the strip"),
                    vpos,
                ));
            }
            let c0 = *axis_origin + *axis * ((origin - *axis_origin).dot(&n) / denom);
            let section = if denom.abs() >= 1.0 - 1e-9 {
                super::strip_end_arcs::strip_end_section_curves::verified_circle_section(
                    axis_origin,
                    axis,
                    *radius,
                    c0,
                    cp,
                    cq,
                    &origin,
                    &n,
                )
            } else {
                super::strip_end_arcs::strip_end_section_curves::verified_ellipse_section(
                    axis_origin,
                    axis,
                    *radius,
                    c0,
                    cp,
                    cq,
                    &origin,
                    &n,
                    denom,
                )
            };
            section
                .ok_or_else(|| decline(format!("no verified end-face section at {v:?}"), vpos))?
        }
        _ => {
            let fitted = super::hole_rim_terminal_section::strip_plane_section(
                strip,
                origin,
                n,
                cp,
                cq,
                tolerance.linear(),
            )
            .map_err(|e| decline(format!("end-face section at {v:?}: {e}"), vpos))?;
            let (lo, hi) = fitted.parameter_range();
            let curve = CurveOnEdge {
                geometry: CurveGeometry::Nurbs(fitted),
                t_start: lo,
                t_end: hi,
            };
            (curve, cp, cq)
        }
    };
    // Split F's corner edges at the crossings: the pieces through V are
    // the corner, cut away when the strip is wired.
    let mut split = |edge: EdgeHandle, at: Point3| -> Result<_, BlendFailure> {
        let res = brep
            .split_edge_at_point(edge, at, tolerance)
            .map_err(|e| decline(format!("end-face edge {edge:?} split: {e}"), at))?;
        brep.vertices[res.vertex].position = at;
        let has_v = |brep: &BRep, e: EdgeHandle| {
            let he = &brep.half_edges[brep.edges[e].half_edge];
            he.origin == v || brep.half_edges[he.twin].origin == v
        };
        let corner = if has_v(brep, res.edge_a) {
            res.edge_a
        } else {
            res.edge_b
        };
        Ok((res.vertex, CapSide::Third(corner)))
    };
    let (vp, side_p) = split(a, cp)?;
    let (vq, side_q) = split(b, cq)?;
    // The section edge, its curve running P → Q.
    let starts_at_p = curve.geometry.point_at(curve.t_start).distance_to(&cp)
        <= curve.geometry.point_at(curve.t_end).distance_to(&cp);
    let x = insert_edge(brep, curve);
    let x_pq = brep.half_edges.insert(blank(vp, x));
    let x_qp = brep.half_edges.insert(blank(vq, x));
    brep.half_edges[x_pq].twin = x_qp;
    brep.half_edges[x_qp].twin = x_pq;
    brep.edges[x].half_edge = if starts_at_p { x_pq } else { x_qp };
    // A ball centre for the strip's orientation: the section ball at V.
    let support = |he: HalfEdgeHandle| {
        brep.faces[brep.loops[brep.half_edges[he].loop_].face]
            .surface
            .as_ref()
            .map(|s| s.geometry.clone())
    };
    let r = link.eg.effective_radius;
    let center = match (support(h), support(tw)) {
        (Some(sp), Some(sq)) => {
            let seed = ball_center(&sp, cp, cq, tangent, r).unwrap_or(vpos);
            exact_section(&sp, &sq, vpos, tangent, r, seed).map_or(seed, |(c, _, _)| c)
        }
        _ => vpos,
    };
    Ok(Cap {
        v,
        vp,
        vq,
        x_pq,
        x_qp,
        center,
        face: f,
        end_face: true,
        at_start,
        side_p,
        side_q,
    })
}

/// Cut the corner away: in F's loop the two corner pieces (`removed`'s
/// twins, consecutive at V) are replaced by the section half-edge
/// `cap_x`; the corner pieces and V go.
///
/// # Errors
///
/// A typed decline when the corner pieces are not consecutive in F.
pub(super) fn close_end_face(
    brep: &mut BRep,
    cap: &Cap,
    cap_x: HalfEdgeHandle,
    removed: [HalfEdgeHandle; 2],
) -> Result<(), BlendFailure> {
    let at = brep.vertices[cap.v].position;
    let [t0, t1] = removed.map(|he| brep.half_edges[he].twin);
    let (first, second) = if brep.half_edges[t0].next == t1 {
        (t0, t1)
    } else if brep.half_edges[t1].next == t0 {
        (t1, t0)
    } else {
        return Err(decline(
            "end-face corner pieces are not consecutive".into(),
            at,
        ));
    };
    let (pr, nx) = (brep.half_edges[first].prev, brep.half_edges[second].next);
    let lp = brep.half_edges[first].loop_;
    brep.half_edges[pr].next = cap_x;
    brep.half_edges[cap_x].prev = pr;
    brep.half_edges[cap_x].next = nx;
    brep.half_edges[nx].prev = cap_x;
    brep.half_edges[cap_x].loop_ = lp;
    if brep.loops[lp].half_edge == first || brep.loops[lp].half_edge == second {
        brep.loops[lp].half_edge = cap_x;
    }
    for he in removed.into_iter().chain([t0, t1]) {
        let edge = brep.half_edges[he].edge;
        brep.remove_half_edge_reanchored(he);
        brep.edges.remove(edge);
    }
    brep.vertices.remove(cap.v);
    Ok(())
}
