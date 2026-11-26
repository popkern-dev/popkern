//! An edge-face pierce the edge-face solver cannot place is placed
//! where the records cross the edge (PK-1292, tangency design §4.8).
//!
//! Where an operand edge runs inside the zone of a tangent crossing
//! (`√(2·reach / min|λ|)`, where the two faces stay within the reach
//! of each other), the edge lies within the reach of the other face
//! along a stretch many tolerances long. Every point of that stretch
//! is a root at the resolution: the solver reports one of them, its
//! own closest approach, where no record passes (a cone's seam 1.0e-3
//! from a cylinder's tangent crossing: the hit 1.36e-3 from the
//! crossing, 8e-7 off the cylinder, while the crossing's arms cross
//! the seam elsewhere unsplit and the face's walk overwrapped its
//! period by the miss). The records of the edge's faces against the
//! face cross the edge transversally on the edge's own face, which is
//! well posed where the pierce is not (only beside a tangent crossing of
//! the edge's face with the hit's face, `beside_tangent_crossing`), and
//! they are the curves the
//! pierce must split: such a hit is replaced by the edge's crossings
//! with those records inside its stretch. A hit stays as it was where a
//! record passes it within the reach (records ending at an edge's exact
//! touch), where its stretch holds no record crossing (the edge touches
//! the face there without a curve meeting it), and where it is well
//! conditioned (its stretch under `ILL_CONDITIONED_STRETCH` reaches, or
//! running to an end of the edge).
//!
//! Tolerance: the linear reach of the boolean's tolerance (two places
//! within it are one), plus each record's own fit deviation where it
//! crosses the edge; no fixed constant is a length.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, closest_point_in_window};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::intersect::IntersectionRecord;
use crate::intersect_edge_face::EdgeFaceHit;

/// A stretch of the edge within the reach of the face longer than this
/// many reaches on either side of the hit makes the hit ill-conditioned
/// (a transversal pierce shallower than about one degree).
const ILL_CONDITIONED_STRETCH: f64 = 64.0;

/// Edge samples across the stretch when searching for record crossings.
const STRETCH_SAMPLES: usize = 256;

/// Replace each ill-conditioned edge-face hit by the edge's crossings
/// with the records of the edge's faces against the hit's face, inside
/// the stretch where the edge lies within the reach of that face.
pub(crate) fn respell_ill_conditioned_hits(
    hits: Vec<EdgeFaceHit>,
    records: &[IntersectionRecord],
    (a, b): (&BRep, &BRep),
    tolerance: &Tolerance,
) -> Vec<EdgeFaceHit> {
    let reach = tolerance.linear_reach();
    let mut out: Vec<EdgeFaceHit> = Vec::with_capacity(hits.len());
    for h in hits {
        let (owner, other) = if h.edge_of_a { (a, b) } else { (b, a) };
        let replaced = respell_one(&h, records, (owner, other), (a, b), reach);
        match replaced {
            Some(new) if !new.is_empty() => {
                if popkern_core::trace::trace_enabled("boolean") {
                    eprintln!(
                        "  [zone-hit] {:?} x {:?} at ({:.6},{:.6},{:.6}) -> {} record crossing(s)",
                        h.edge,
                        h.face,
                        h.point.x,
                        h.point.y,
                        h.point.z,
                        new.len()
                    );
                }
                for n in new {
                    let dup = out.iter().any(|o| {
                        o.edge == n.edge
                            && o.face == n.face
                            && o.edge_of_a == n.edge_of_a
                            && o.point.distance_to(&n.point) <= reach
                    });
                    if !dup {
                        out.push(n);
                    }
                }
            }
            _ => out.push(h),
        }
    }
    out
}

/// The record crossings standing for one hit, or `None` when the hit
/// is well conditioned (or its geometry is missing).
fn respell_one(
    h: &EdgeFaceHit,
    records: &[IntersectionRecord],
    (owner, other): (&BRep, &BRep),
    (a, b): (&BRep, &BRep),
    reach: f64,
) -> Option<Vec<EdgeFaceHit>> {
    let curve = owner.edges.get(h.edge)?.curve.as_ref()?;
    let surface = &other.faces.get(h.face)?.surface.as_ref()?.geometry;
    let (lo, hi) = (
        curve.t_start.min(curve.t_end),
        curve.t_start.max(curve.t_end),
    );
    let c = &curve.geometry;
    let gap = |t: f64| surface_gap(surface, &c.point_at(t));
    let stretch = within_reach_stretch(c, (lo, hi), h.t_edge, &gap, reach)?;
    let (f1, f2) = faces_of_edge(owner, h.edge)?;
    if ![f1, f2]
        .into_iter()
        .any(|f| beside_tangent_crossing((owner, f), (other, h.face), &h.point, reach))
    {
        return None;
    }
    let mut found = Vec::new();
    for r in records {
        let (edge_face, hit_face) = if h.edge_of_a {
            (r.face_a, r.face_b)
        } else {
            (r.face_b, r.face_a)
        };
        if hit_face != h.face || (edge_face != f1 && edge_face != f2) {
            continue;
        }
        for t in record_crossings(c, stretch, r, (a, b), reach) {
            let point = c.point_at(t);
            let Some(uv) = surface.closest_point(&point) else {
                continue;
            };
            found.push(EdgeFaceHit {
                edge: h.edge,
                face: h.face,
                point,
                t_edge: t,
                uv,
                edge_of_a: h.edge_of_a,
            });
        }
    }
    // A record passing the hit itself (records ending at an edge's exact
    // touch, or a pierce the solver placed well) keeps it: only a hit no
    // record passes is the solver's guess in the stretch.
    if found.iter().any(|n| n.point.distance_to(&h.point) <= reach) {
        return None;
    }
    Some(found)
}

/// The parameter interval round `t0` where the edge lies within the
/// reach of the face, when it runs past `ILL_CONDITIONED_STRETCH`
/// reaches on either side and leaves the reach again before both of
/// the edge's ends; `None` for a well-conditioned pierce or an edge
/// lying on the face to an end.
fn within_reach_stretch(
    c: &dyn CurveEval,
    (lo, hi): (f64, f64),
    t0: f64,
    gap: &dyn Fn(f64) -> f64,
    reach: f64,
) -> Option<(f64, f64)> {
    let speed = c.derivative_at(t0, 1).length();
    if speed.is_nan() || speed <= 0.0 || !gap(t0).is_finite() {
        return None;
    }
    let probe = ILL_CONDITIONED_STRETCH * reach / speed;
    let long_side = [-1.0, 1.0].into_iter().any(|s: f64| {
        let t = t0 + s * probe;
        (lo..=hi).contains(&t) && gap(t) <= reach
    });
    if !long_side {
        return None;
    }
    // Bounded inside the edge on both sides: an edge lying on the face
    // to its end (coplanar faces, an edge along a touch) is the
    // coincidence paths', not a pierce in a zone.
    let bound = |s: f64| -> Option<f64> {
        // Double outward while within the reach, then bisect the exit.
        let limit = if s < 0.0 { lo } else { hi };
        let (mut inside, mut step) = (t0, probe);
        loop {
            let t = inside + s * step;
            if (s < 0.0 && t <= limit) || (s > 0.0 && t >= limit) {
                if gap(limit) <= reach {
                    return None;
                }
                return Some(bisect_exit(inside, limit, gap, reach));
            }
            if gap(t) > reach {
                return Some(bisect_exit(inside, t, gap, reach));
            }
            inside = t;
            step *= 2.0;
        }
    };
    Some((bound(-1.0)?, bound(1.0)?))
}

/// The point between `inside` (within the reach) and `outside` (past
/// it) where the gap crosses the reach.
fn bisect_exit(mut inside: f64, mut outside: f64, gap: &dyn Fn(f64) -> f64, reach: f64) -> f64 {
    for _ in 0..60 {
        let m = 0.5 * (inside + outside);
        if gap(m) <= reach {
            inside = m;
        } else {
            outside = m;
        }
    }
    inside
}

/// The parameters inside `stretch` where the edge crosses record `r`:
/// local minima of the distance from the edge to the record, refined,
/// within the reach plus the record's fit deviation there.
fn record_crossings(
    c: &dyn CurveEval,
    (t_lo, t_hi): (f64, f64),
    r: &IntersectionRecord,
    (a, b): (&BRep, &BRep),
    reach: f64,
) -> Vec<f64> {
    let rc = &r.curve.curve_3d;
    let (r_lo, r_hi) = rc.parameter_range();
    let dist = |t: f64| -> (f64, Point3) {
        let p = c.point_at(t);
        let s = closest_point_in_window(rc, &p, r_lo, r_hi);
        let q = rc.point_at(s);
        (p.distance_to(&q), q)
    };
    let n = STRETCH_SAMPLES;
    let ts: Vec<f64> = (0..=n)
        .map(|k| t_lo + (t_hi - t_lo) * k as f64 / n as f64)
        .collect();
    let ds: Vec<f64> = ts.iter().map(|&t| dist(t).0).collect();
    let mut out: Vec<f64> = Vec::new();
    for k in 0..=n {
        let left = if k == 0 { f64::INFINITY } else { ds[k - 1] };
        let right = if k == n { f64::INFINITY } else { ds[k + 1] };
        if !(ds[k] <= left && ds[k] < right) {
            continue;
        }
        let (mut x0, mut x1) = (ts[k.saturating_sub(1)], ts[(k + 1).min(n)]);
        for _ in 0..80 {
            let m0 = x0 + (x1 - x0) / 3.0;
            let m1 = x1 - (x1 - x0) / 3.0;
            if dist(m0).0 <= dist(m1).0 {
                x1 = m1;
            } else {
                x0 = m0;
            }
        }
        let t = 0.5 * (x0 + x1);
        let (d, q) = dist(t);
        let deviation = crate::tolerance_zone::candidates::record_fit_deviation(a, b, r, false, &q);
        if d <= reach + deviation && !out.iter().any(|&u| (u - t).abs() <= f64::EPSILON * 16.0) {
            out.push(t);
        }
    }
    out
}

/// Distance from `p` to the surface (its closest point), infinite when
/// the projection fails.
fn surface_gap(surface: &SurfaceGeometry, p: &Point3) -> f64 {
    surface
        .closest_point(p)
        .map(|(u, v)| surface.point_at(u, v).distance_to(p))
        .filter(|d| d.is_finite())
        .unwrap_or(f64::INFINITY)
}

/// The two faces an edge bounds (one face twice for a seam).
fn faces_of_edge(brep: &BRep, eh: EdgeHandle) -> Option<(FaceHandle, FaceHandle)> {
    let he = brep.half_edges.get(brep.edges.get(eh)?.half_edge)?;
    let f1 = brep.loops.get(he.loop_)?.face;
    let tw = he.mate().and_then(|m| brep.half_edges.get(m))?;
    let f2 = brep.loops.get(tw.loop_)?.face;
    Some((f1, f2))
}

/// Whether the edge's face `f` and the hit's face meet at a tangent
/// crossing found from the hit (two branches of their intersection
/// crossing there), the zone the stretch belongs to. Faces tangent along
/// a curve (a boss a hair under tangency to a web) are a line contact,
/// one tangent curve at the resolution: their hits are not respelled.
fn beside_tangent_crossing(
    (fb, f): (&BRep, FaceHandle),
    (hb, hf): (&BRep, FaceHandle),
    near: &Point3,
    reach: f64,
) -> bool {
    use popkern_core::geometry::intersection::contact_classify::{
        ContactScale, SurfaceContactKind,
    };
    use popkern_core::geometry::intersection::contact_discovery::contact_from;
    let (Some((sa, wa)), Some((sb, wb))) = (window(fb, f), window(hb, hf)) else {
        return false;
    };
    let Some(seed) = popkern_core::geometry::uv_projection::project_point_to_uv(sa, near) else {
        return false;
    };
    let extent = (wa.0.1 - wa.0.0)
        .abs()
        .max((wa.1.1 - wa.1.0).abs())
        .max(1e-9);
    let size = sa
        .point_at(wa.0.0, wa.1.0)
        .distance_to(&sa.point_at(wa.0.1, wa.1.1));
    let scale = ContactScale {
        linear: reach,
        extent: size.max(reach),
    };
    contact_from((sa, wa), (sb, wb), seed, 1e-6 * extent, scale).is_some_and(|c| {
        matches!(
            c.classification.kind,
            SurfaceContactKind::CrossingPoint { .. }
        )
    })
}

/// A face's surface and its trim's parameter window.
type SurfaceWindow<'a> = (&'a SurfaceGeometry, ((f64, f64), (f64, f64)));

/// `face`'s surface and trim window.
fn window(b: &BRep, face: FaceHandle) -> Option<SurfaceWindow<'_>> {
    Some((
        &b.faces.get(face)?.surface.as_ref()?.geometry,
        b.face_trim_uv_window(face)?,
    ))
}
