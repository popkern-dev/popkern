//! Geometry and topology helpers of the closed tangent-chain
//! construction (split from `tangent_chain_surgery` for the line cap):
//! junction sections, section arcs, support rails, strip orientation,
//! half-edge splicing.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::{
    BRep, Edge, EdgeHandle, FaceHandle, HalfEdge, HalfEdgeHandle, LoopHandle, VertexHandle,
};
use popkern_math::PortableMath;

pub(super) type Circle = Option<(Point3, Vec3, f64)>;

/// The chain direction along half-edge `h` at its origin.
pub(super) fn chain_tangent_at_origin(brep: &BRep, h: HalfEdgeHandle) -> Option<Vec3> {
    let he = &brep.half_edges[h];
    let t = super::g1_junction::edge_tangent(brep, he.edge, he.origin)?;
    let e = &brep.edges[he.edge];
    let along = if e.half_edge == h { 1.0 } else { -1.0 };
    Some(t * along)
}

/// The chain direction along half-edge `h` at its destination.
pub(super) fn chain_tangent(brep: &BRep, h: HalfEdgeHandle) -> Option<Vec3> {
    let he = &brep.half_edges[h];
    let to = brep.half_edges[he.twin].origin;
    let t = super::g1_junction::edge_tangent(brep, he.edge, to)?;
    // edge_tangent runs along the stored curve; orient it along h.
    let e = &brep.edges[he.edge];
    let along = if e.half_edge == h { 1.0 } else { -1.0 };
    Some(t * along)
}

/// The contact point nearest `p`: on the exact circle when the trim
/// is one (the fitted trim of a circular contact can stop short of a
/// junction — measured 1.4° on a rounded-rectangle corner), else on
/// the fitted trim.
pub(super) fn nearest(trim: (&NurbsCurve, Circle), p: Point3) -> Point3 {
    if let (_, Some((center, normal, radius))) = trim
        && let Some(n) = normal.normalized()
    {
        let d = p - center;
        let in_plane = d - n * d.dot(&n);
        if let Some(dir) = in_plane.normalized() {
            return center + dir * radius;
        }
    }
    trim.0
        .closest_point(&p)
        .map(|s| trim.0.point_at(s))
        .unwrap_or(p)
}

/// Distance from `p` to a curveless edge's chord.
pub(super) fn segment_distance(brep: &BRep, edge: EdgeHandle, p: Point3) -> f64 {
    let he = &brep.half_edges[brep.edges[edge].half_edge];
    let a = brep.vertices[he.origin].position;
    let b = brep.vertices[brep.half_edges[he.twin].origin].position;
    let d = b - a;
    let l2 = d.dot(&d).max(1e-300);
    let s = ((p - a).dot(&d) / l2).clamp(0.0, 1.0);
    (a + d * s).distance_to(&p)
}

/// Re-window the third edge's curve to its (moved) endpoint `v`.
pub(super) fn shorten_third(brep: &mut BRep, edge: EdgeHandle, v: VertexHandle) {
    let p = brep.vertices[v].position;
    let at_start = brep.half_edges[brep.edges[edge].half_edge].origin == v;
    if let Some(curve) = brep.edges[edge].curve.as_mut()
        && let Some(s) = curve.geometry.closest_point(&p)
    {
        if at_start {
            curve.t_start = s;
        } else {
            curve.t_end = s;
        }
    }
}

/// The trim between `a` and `b`, running in the chain direction `dir`
/// (an exact arc for a circular trim, else the fitted trim's window).
pub(super) fn trim_piece(
    trim: (&NurbsCurve, Circle),
    a: Point3,
    b: Point3,
    dir: Vec3,
) -> CurveOnEdge {
    if let (_, Some((center, normal, radius))) = trim {
        let ra = a - center;
        // Orient the circle so it turns from `a` along the chain.
        let n = if normal.cross(&ra).dot(&dir) >= 0.0 {
            normal
        } else {
            normal * -1.0
        };
        let (e1, e2) = popkern_core::geometry::analytic_curves::circle_basis(n);
        let ang = |p: Point3| {
            let d = p - center;
            d.dot(&e2).patan2(d.dot(&e1))
        };
        let a0 = ang(a);
        let mut sweep = (ang(b) - a0).rem_euclid(std::f64::consts::TAU);
        if sweep < 1e-12 {
            sweep = std::f64::consts::TAU;
        }
        return CurveOnEdge {
            geometry: CurveGeometry::Arc {
                center,
                normal: n,
                radius,
                start_angle: a0,
                end_angle: a0 + sweep,
            },
            t_start: a0,
            t_end: a0 + sweep,
        };
    }
    let curve = trim.0;
    let sa = curve.closest_point(&a).unwrap_or(0.0);
    let sb = curve.closest_point(&b).unwrap_or(1.0);
    if sa <= sb {
        CurveOnEdge {
            geometry: CurveGeometry::Nurbs(curve.clone()),
            t_start: sa,
            t_end: sb,
        }
    } else {
        let r = reversed(curve);
        let (lo, hi) = curve.parameter_range();
        CurveOnEdge {
            geometry: CurveGeometry::Nurbs(r),
            t_start: lo + hi - sa,
            t_end: lo + hi - sb,
        }
    }
}

pub(super) fn reversed(c: &NurbsCurve) -> NurbsCurve {
    let (lo, hi) = (c.knots[0], c.knots[c.knots.len() - 1]);
    NurbsCurve {
        degree: c.degree,
        knots: c.knots.iter().rev().map(|k| lo + hi - k).collect(),
        control_points: c.control_points.iter().rev().copied().collect(),
        weights: c.weights.iter().rev().copied().collect(),
    }
}

/// The exact junction section: the ball of radius `r` whose centre lies
/// in the plane through `v` normal to `n` and touches both supports —
/// solved by Newton on the two distance residuals from `guess` (the
/// fitted trims stop short of a junction: measured 0.26 at r = 6 on a
/// cylinder tee's saddle, so their ends cannot be the section). Returns
/// `(centre, contact on P, contact on Q)`.
pub(super) fn exact_section(
    sp: &SurfaceGeometry,
    sq: &SurfaceGeometry,
    v: Point3,
    n: Vec3,
    r: f64,
    guess: Point3,
) -> Option<(Point3, Point3, Point3)> {
    let n = n.normalized()?;
    let (e1, e2) = popkern_core::geometry::analytic_curves::circle_basis(n);
    let foot = |s: &SurfaceGeometry, c: Point3| -> Option<Point3> { surface_foot(s, c) };
    let resid = |a: f64, b: f64| -> Option<[f64; 2]> {
        let c = v + e1 * a + e2 * b;
        Some([
            foot(sp, c)?.distance_to(&c) - r,
            foot(sq, c)?.distance_to(&c) - r,
        ])
    };
    let g = guess - v;
    let (mut a, mut b) = (g.dot(&e1), g.dot(&e2));
    let h = r * 1e-6;
    for _ in 0..40 {
        let f = resid(a, b)?;
        if f[0].abs().max(f[1].abs()) < r * 1e-13 {
            break;
        }
        let fa = resid(a + h, b)?;
        let fb = resid(a, b + h)?;
        let j = [
            [(fa[0] - f[0]) / h, (fb[0] - f[0]) / h],
            [(fa[1] - f[1]) / h, (fb[1] - f[1]) / h],
        ];
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        if det.abs() < 1e-300 {
            return None;
        }
        a -= (f[0] * j[1][1] - f[1] * j[0][1]) / det;
        b -= (j[0][0] * f[1] - j[1][0] * f[0]) / det;
    }
    let f = resid(a, b)?;
    if f[0].abs().max(f[1].abs()) > r * 1e-9 {
        return None;
    }
    let c = v + e1 * a + e2 * b;
    Some((c, foot(sp, c)?, foot(sq, c)?))
}

/// The ball centre at a junction: in the plane normal to the chain
/// `tangent`, `r` from both contacts — of the two such points, the one
/// whose offset from `cp` runs along the P support's surface normal
/// there (the other is its mirror across the chord).
pub(super) fn ball_center(
    support_p: &SurfaceGeometry,
    cp: Point3,
    cq: Point3,
    tangent: Vec3,
    r: f64,
) -> Option<Point3> {
    let n = tangent.normalized()?;
    let chord = cq - cp;
    let half = chord.length() * 0.5;
    if half > r * (1.0 + 1e-9) || half < 1e-12 {
        return None;
    }
    let mid = cp + chord * 0.5;
    let w = n.cross(&chord).normalized()?;
    let h = (r * r - half * half).max(0.0).sqrt();
    let (u, v) = support_p.closest_point(&cp)?;
    let normal = support_p.normal_at(u, v).normalized()?;
    [mid + w * h, mid + w * (-h)].into_iter().max_by(|a, b| {
        let da = ((*a - cp) * (1.0 / r)).dot(&normal).abs();
        let db = ((*b - cp) * (1.0 / r)).dot(&normal).abs();
        da.total_cmp(&db)
    })
}

/// The straight section from `cp` to `cq` when it lies on BOTH strips'
/// surfaces (within `band`, sampled inside the segment): a chamfer's
/// cross-section is a line, and its plane or cone strips are ruled
/// along it. `None` for a fillet's curved section. The junction used
/// to author every section as a rolling-ball arc: a chamfered rounded
/// plate's corner cone met its side plane along a 45° line served as a
/// quarter circle 9.2e-3 off the cone, and check 9 declined the chamfer
/// (PK-429).
pub(super) fn section_line_on(
    cp: Point3,
    cq: Point3,
    surfaces: [&popkern_core::geometry::surface::SurfaceGeometry; 2],
    band: f64,
) -> Option<CurveOnEdge> {
    use popkern_core::geometry::surface::SurfaceEval;
    if cp.distance_to(&cq) <= band {
        return None;
    }
    let on_both = [0.25, 0.5, 0.75].into_iter().all(|f| {
        let p = cp + (cq - cp) * f;
        surfaces.iter().all(|s| {
            s.closest_point(&p)
                .is_some_and(|(u, v)| s.point_at(u, v).distance_to(&p) <= band)
        })
    });
    on_both.then(|| CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: cp,
            direction: cq - cp,
        },
        t_start: 0.0,
        t_end: 1.0,
    })
}

/// The short arc of the circle about `center` (in the plane normal to
/// `tangent`) from `cp` to `cq`: a fillet section subtends less than π.
pub(super) fn section_arc(
    cp: Point3,
    cq: Point3,
    center: Point3,
    tangent: Vec3,
    r: f64,
) -> Option<CurveOnEdge> {
    // The arc's plane holds the centre and both contacts (they are the
    // feet of the centre, `r` from it): a plane normal to the chain
    // tangent left the arc ends ~1e-6 off contacts lying a hair out of
    // it. The tangent only orients the normal.
    let (da, db) = (cp - center, cq - center);
    let mut n0 = da
        .cross(&db)
        .normalized()
        .or_else(|| tangent.normalized())?;
    if let Some(t) = tangent.normalized()
        && n0.dot(&t) < 0.0
    {
        n0 = n0 * -1.0;
    }
    for n in [n0, n0 * -1.0] {
        let (e1, e2) = popkern_core::geometry::analytic_curves::circle_basis(n);
        let ang = |p: Point3| {
            let d = p - center;
            d.dot(&e2).patan2(d.dot(&e1))
        };
        let a0 = ang(cp);
        let sweep = (ang(cq) - a0).rem_euclid(std::f64::consts::TAU);
        if sweep > 1e-12 && sweep < std::f64::consts::PI {
            return Some(CurveOnEdge {
                geometry: CurveGeometry::Arc {
                    center,
                    normal: n,
                    radius: r,
                    start_angle: a0,
                    end_angle: a0 + sweep,
                },
                t_start: a0,
                t_end: a0 + sweep,
            });
        }
    }
    None
}

/// Orient the strip outward: away from the ball centre on a convex
/// fillet (the ball is in the material), toward it on a concave one.
pub(super) fn orient_strip(
    brep: &mut BRep,
    face: FaceHandle,
    convex: bool,
    center: Point3,
    on_strip: Point3,
) {
    let Some(sof) = brep.faces[face].surface.as_mut() else {
        return;
    };
    let Some((u, v)) = sof.geometry.closest_point(&on_strip) else {
        return;
    };
    let normal = sof.geometry.normal_at(u, v);
    let away = on_strip - center;
    let outward = if convex { away } else { away * -1.0 };
    sof.orientation = normal.dot(&outward) > 0.0;
}

pub(super) fn blank(origin: VertexHandle, edge: EdgeHandle) -> HalfEdge {
    HalfEdge {
        origin,
        twin: HalfEdgeHandle::default(),
        next: HalfEdgeHandle::default(),
        prev: HalfEdgeHandle::default(),
        edge,
        loop_: LoopHandle::default(),
        pcurve: None,
    }
}

pub(super) fn insert_edge(brep: &mut BRep, curve: CurveOnEdge) -> EdgeHandle {
    brep.insert_edge_tracked(Edge {
        half_edge: HalfEdgeHandle::default(),
        curve: Some(curve),
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: Default::default(),
    })
}

/// Put `new` where `old` sits in its loop (same loop, prev, next).
pub(super) fn substitute(brep: &mut BRep, old: HalfEdgeHandle, new: HalfEdgeHandle) {
    let (prev, next, lp) = {
        let o = &brep.half_edges[old];
        (o.prev, o.next, o.loop_)
    };
    brep.half_edges[new].prev = prev;
    brep.half_edges[new].next = next;
    brep.half_edges[new].loop_ = lp;
    brep.half_edges[prev].next = new;
    brep.half_edges[next].prev = new;
    if brep.loops[lp].half_edge == old {
        brep.loops[lp].half_edge = new;
    }
}

/// `circle` when its points lie on `support` within 10 × the linear
/// tolerance (16 samples) AND it carries the fitted `trim` (the trim's
/// samples lie on it, within 5% of its radius: fitted trims stray
/// ~1e-3 from the exact circle, a misplaced circle by its whole offset),
/// else `None`: the fitted trim stands in. On an unbounded support the
/// first test alone cannot place the circle. A pocket-floor rim's wall
/// circle carried the cylinder origin's height (z = 0 for contacts at
/// z = 10) and lies on the infinite wall; only a failing multi-start
/// projection rejected it, until PK-737 made the projection exact.
pub(super) fn verified_circle(
    trim: &NurbsCurve,
    circle: Circle,
    support: Option<&SurfaceGeometry>,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Circle {
    let ((center, normal, radius), support) = (circle?, support?);
    let n = normal.normalized()?;
    let (e1, e2) = popkern_core::geometry::analytic_curves::circle_basis(n);
    let band = tolerance.linear().max(1e-12) * 10.0;
    let on = (0..16).all(|k| {
        let a = std::f64::consts::TAU * f64::from(k) / 16.0;
        let p = center + (e1 * a.pcos() + e2 * a.psin()) * radius;
        support
            .closest_point(&p)
            .map(|(u, v)| support.point_at(u, v).distance_to(&p) <= band)
            .unwrap_or(false)
    });
    // The fitted trim's samples on the circle.
    let carry = (band * 10.0).max(0.05 * radius);
    let (t0, t1) = trim.parameter_range();
    let carries = (0..=8).all(|k| {
        let q = trim.point_at(t0 + (t1 - t0) * f64::from(k) / 8.0);
        let d = q - center;
        let axial = d.dot(&n);
        let radial = (d - n * axial).length();
        axial.abs().phypot(radial - radius) <= carry
    });
    (on && carries).then_some((center, normal, radius))
}

/// The foot of `p` on `surface`: its closest point, exact for analytic
/// surfaces (PK-737; this helper once carried the closed forms itself,
/// working around the multi-start's 0.18 miss near a cylinder's axis).
pub(super) fn surface_foot(surface: &SurfaceGeometry, p: Point3) -> Option<Point3> {
    let (u, w) = surface.closest_point(&p)?;
    Some(surface.point_at(u, w))
}

/// One chain edge in chain order.
pub(super) struct Link<'a> {
    pub(super) eg: &'a super::EdgeGeometry,
    /// The half-edge running along the chain (origin `v_{i-1}`).
    pub(super) h: HalfEdgeHandle,
    pub(super) p_trim: (&'a NurbsCurve, Circle),
    pub(super) q_trim: (&'a NurbsCurve, Circle),
}

/// The third edge at a junction and its side.
#[derive(Clone, Copy)]
pub(super) struct Third {
    /// `None` at a valence-2 junction: both supports continue as the
    /// same faces across it (an edge split mid-run, as a per-edge
    /// corner fillet leaves its rim arc), so nothing is shortened.
    pub(super) edge: Option<EdgeHandle>,
    pub(super) on_p: bool,
}

/// A typed decline of the chain construction at `at`.
pub(super) fn decline(reason: String, at: Point3) -> crate::error::BlendFailure {
    crate::error::BlendFailure::InvalidTopology {
        reason: format!("tangent chain: {reason}"),
        location: Some(at),
    }
}
