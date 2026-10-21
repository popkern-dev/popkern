//! End-face ends of open tangent chains (PK-735,
//! `docs/architecture/fillet-open-tangent-chains.md`): the chain's end
//! vertex V is a sharp corner, and the strip runs out on the planar END
//! FACE F that holds V's other two edges `a` (P ∩ F) and `b` (Q ∩ F).
//! The strip's end is its section by F (exact circle or ellipse for a
//! cylinder strip, a verified fitted section otherwise) between the
//! rails' crossings of F; `a` and `b` are split there and the corner —
//! V and the two corner pieces — is cut away, the section taking its
//! place in F's loop.
//!
//! At a REFLEX corner (PK-798: a fin's root arc split by the housing
//! seam, so the arc runs as a tangent chain) one rail crosses F's plane
//! on its support past V, on the extension of F's edge on that side,
//! and the strip's end overhangs F. The corner then stays: a new
//! terminal DISC in F's plane closes the strip, bounded by the section,
//! the split piece of the other edge through V, and a new line on the
//! overhung support from V to its crossing. That is the disc the
//! per-edge pipeline builds for a lone root arc (PK-740).

use popkern_core::geometry::curve::{
    CurveEval, CurveGeometry, CurveOnEdge, closest_point_on_extension,
};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, Face, HalfEdgeHandle, LoopHandle};
use popkern_math::PortableMath;

use super::tangent_chain_cap::{Cap, CapSide};
use super::tangent_chain_geometry::{
    Circle, Link, ball_center, blank, chain_tangent, chain_tangent_at_origin, decline,
    exact_section, insert_edge, surface_foot,
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

/// The crossing of the exact `circle` with the plane `(origin, n)`
/// nearest `near`: `A cos θ + B sin θ = C` in the circle's basis. `None`
/// when the circle misses the plane or lies in it.
fn circle_crossing(
    circle: (Point3, Vec3, f64),
    origin: Point3,
    n: Vec3,
    near: Point3,
) -> Option<Point3> {
    let (center, normal, radius) = circle;
    let (e1, e2) = popkern_core::geometry::analytic_curves::circle_basis(normal);
    let (a, b) = (radius * e1.dot(&n), radius * e2.dot(&n));
    let c = (origin - center).dot(&n);
    let amp = (a * a + b * b).sqrt();
    if amp <= f64::EPSILON * radius.max(1.0) || c.abs() > amp {
        return None;
    }
    let phase = b.patan2(a);
    let spread = (c / amp).clamp(-1.0, 1.0).pacos();
    [phase + spread, phase - spread]
        .into_iter()
        .map(|t| center + e1 * (radius * t.pcos()) + e2 * (radius * t.psin()))
        .min_by(|p, q| p.distance_to(&near).total_cmp(&q.distance_to(&near)))
}

/// Where a rail's crossing of F lands relative to F's edge from V.
#[derive(Clone, Copy)]
enum Landing {
    /// On the edge: its corner piece through V is cut away.
    OnEdge(Point3),
    /// On the edge's extension past V: the strip overhangs F there.
    PastCorner(Point3),
}

impl Landing {
    fn point(self) -> Point3 {
        match self {
            Landing::OnEdge(p) | Landing::PastCorner(p) => p,
        }
    }
}

/// Land the crossing `raw` on `edge` (an edge of F through V): on its
/// trimmed range, or on its extension beyond V. `None` when it lies
/// off the edge's curve or past its far end.
fn landing(brep: &BRep, edge: EdgeHandle, raw: Point3, vpos: Point3, band: f64) -> Option<Landing> {
    let curve = brep.edges[edge].curve.as_ref()?;
    let s = curve.geometry.closest_point(&raw)?;
    let (lo, hi) = (
        curve.t_start.min(curve.t_end),
        curve.t_start.max(curve.t_end),
    );
    // Past the corner the rail lands on the edge curve's extension,
    // which the in-range projection does not reach (PK-907).
    let on_curve = closest_point_on_extension(&curve.geometry, &raw)
        .map(|e| curve.geometry.point_at(e))
        .into_iter()
        .chain([curve.geometry.point_at(s)])
        .min_by(|a, b| a.distance_to(&raw).total_cmp(&b.distance_to(&raw)))?;
    let clamped = curve.geometry.point_at(s.clamp(lo, hi));
    if clamped.distance_to(&raw) <= band {
        Some(Landing::OnEdge(clamped))
    } else if on_curve.distance_to(&raw) <= band && clamped.distance_to(&vpos) <= band {
        Some(Landing::PastCorner(on_curve))
    } else {
        None
    }
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
     -> Result<Landing, BlendFailure> {
        let raw = rail_crossing(trim.0, origin, n, vpos, band)
            .ok_or_else(|| decline(format!("a rail does not cross the end face at {v:?}"), vpos))?;
        // An exact circular rail crosses F exactly: the bisection runs
        // on its fitted trim, which strays ~1e-5 from the circle (the
        // strip's rail edge is the circle, so the contact must be too).
        let raw = trim
            .1
            .and_then(|circle| circle_crossing(circle, origin, n, raw))
            .unwrap_or(raw);
        landing(brep, edge, raw, vpos, band).ok_or_else(|| {
            decline(
                format!("rail crosses the end face off its edge {edge:?}"),
                vpos,
            )
        })
    };
    let (lp, lq) = (crossing(link.p_trim, a)?, crossing(link.q_trim, b)?);
    let reflex = match (lp, lq) {
        (Landing::OnEdge(_), Landing::OnEdge(_)) => false,
        (Landing::PastCorner(_), Landing::PastCorner(_)) => {
            return Err(decline(
                format!("both rails overhang the end face at {v:?}"),
                vpos,
            ));
        }
        _ => true,
    };
    let (cp, cq) = (lp.point(), lq.point());
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
    let split = |brep: &mut BRep, edge: EdgeHandle, at: Point3| -> Result<_, BlendFailure> {
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
    // An overhung side gets a new line on its support from V.
    let support_of = |brep: &BRep, he: HalfEdgeHandle| {
        brep.faces[brep.loops[brep.half_edges[he].loop_].face]
            .surface
            .as_ref()
            .map(|s| s.geometry.clone())
    };
    let overhang = |brep: &mut BRep, on: HalfEdgeHandle, c: Point3| {
        let surface = support_of(brep, on)
            .ok_or_else(|| decline(format!("chain end {v:?} lacks a support"), vpos))?;
        let mid = vpos + (c - vpos) * 0.5;
        let off = surface_foot(&surface, mid).map_or(f64::INFINITY, |f| f.distance_to(&mid));
        if off > band {
            return Err(decline(
                format!("overhang line at {v:?} leaves its support by {off:.3e}"),
                vpos,
            ));
        }
        Ok(super::tangent_chain_cap::support_line(brep, v, c))
    };
    let (vp, side_p) = match lp {
        Landing::OnEdge(_) => split(brep, a, cp)?,
        Landing::PastCorner(_) => overhang(brep, h, cp)?,
    };
    let (vq, side_q) = match lq {
        Landing::OnEdge(_) => split(brep, b, cq)?,
        Landing::PastCorner(_) => overhang(brep, tw, cq)?,
    };
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
    // A reflex end keeps F and closes the strip on a new disc in F's
    // plane (oriented with its loop when wired).
    let face = if reflex {
        let (u_axis, v_axis) = popkern_core::geometry::analytic_curves::circle_basis(n);
        let shell = brep.faces[f].shell;
        brep.faces.insert(Face {
            outer_loop: LoopHandle::default(),
            inner_loops: vec![],
            shell,
            surface: Some(SurfaceOnFace {
                geometry: SurfaceGeometry::Plane {
                    origin,
                    normal: n,
                    u_axis,
                    v_axis,
                },
                orientation: true,
            }),
            name: None,
            tolerance: None,
            role: Default::default(),
        })
    } else {
        f
    };
    Ok(Cap {
        v,
        vp,
        vq,
        x_pq,
        x_qp,
        center,
        face,
        end_face: !reflex,
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
