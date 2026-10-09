//! Exact face-side test near curved trim boundaries (PK-747).
//!
//! `BRep::face_uv_trim_containment` classifies against a sampled
//! boundary polygon: six samples per pcurve. On a cycloid disc's face,
//! bounded by multi-span NURBS flank edges, the polygon's chords cut
//! past the concave flank by about the 0.02 mm gap to the ring pins. A
//! pin's surface point then read as lying on the disc face, and the
//! pair as contact at distance 0.
//!
//! When the polygon puts a point inside, this module re-decides near
//! the boundary with the boundary's own curves. It finds the nearest
//! boundary point and the direction across the boundary there, within
//! the surface. It then takes the material side from the polygon itself,
//! at two probes stepped across the boundary by twice the polygon's
//! sagitta: far enough that the polygon cannot misjudge them. The
//! test needs no loop-orientation convention.
//!
//! Points farther from every boundary edge than twice the polygon's
//! sagitta keep the polygon verdict, as do faces too thin for the
//! probes, where both probes read the same.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

/// Samples per edge for the coarse nearest-point search before the
/// closest-point solve.
const EDGE_SAMPLES: usize = 16;

/// Coarse-nearest boundary curves given an exact closest-point solve.
const EXACT_CANDIDATES: usize = 3;

/// Samples per pcurve in the trim polygon
/// (`popkern_core::validation`'s loop polygon: six, more on short
/// loops — six is the coarsest, so its sagitta bounds the error).
const POLYGON_SAMPLES: usize = 6;

/// Whether `x` (on `face`'s surface) is on the face's material side of
/// its nearest boundary edge, when that edge is near enough for the
/// polygon verdict to be in doubt; `None` when no boundary edge is that
/// near (the polygon verdict stands) or the geometry is missing.
pub(crate) fn material_side_near_boundary(
    brep: &BRep,
    face: FaceHandle,
    surface: &SurfaceGeometry,
    x: &Point3,
) -> Option<bool> {
    let f = brep.faces.get(face)?;
    // Coarse pass: every boundary curve's nearest sample.
    let mut coarse: Vec<(f64, &popkern_core::geometry::curve::CurveOnEdge, f64)> = Vec::new();
    for lh in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let start = brep.loops.get(lh)?.half_edge;
        let mut h = start;
        for _ in 0..=brep.half_edges.len() {
            let he = brep.half_edges.get(h)?;
            if let Some(c) = brep.edges.get(he.edge).and_then(|e| e.curve.as_ref()) {
                let (lo, hi) = (c.t_start.min(c.t_end), c.t_start.max(c.t_end));
                let mut near = (f64::INFINITY, lo);
                for k in 0..=EDGE_SAMPLES {
                    let t = lo + (hi - lo) * k as f64 / EDGE_SAMPLES as f64;
                    let d = c.geometry.point_at(t).distance_to(x);
                    if d < near.0 {
                        near = (d, t);
                    }
                }
                coarse.push((near.0, c, near.1));
            }
            h = he.next;
            if h == start {
                break;
            }
        }
    }
    coarse.sort_by(|a, b| a.0.total_cmp(&b.0));
    // Exact nearest point over the few coarse-nearest curves.
    let mut best: Option<(
        f64,
        Point3,
        Vec3,
        &popkern_core::geometry::curve::CurveOnEdge,
    )> = None;
    for &(_, c, t_near) in coarse.iter().take(EXACT_CANDIDATES) {
        let (lo, hi) = (c.t_start.min(c.t_end), c.t_start.max(c.t_end));
        let t = c
            .geometry
            .closest_point(x)
            .filter(|t| *t >= lo && *t <= hi)
            .unwrap_or(t_near);
        let p = c.geometry.point_at(t);
        let d = p.distance_to(x);
        if best.as_ref().is_none_or(|b| d < b.0) {
            best = Some((d, p, c.geometry.tangent_at(t), c));
        }
    }
    let (distance, at, tangent, curve) = best?;
    // How far the trim polygon's chords (POLYGON_SAMPLES per curve) can
    // stray from this curve: the band where its verdict may err.
    let (lo, hi) = (
        curve.t_start.min(curve.t_end),
        curve.t_start.max(curve.t_end),
    );
    let mut sagitta: f64 = 0.0;
    let mut span: f64 = 0.0;
    for k in 0..POLYGON_SAMPLES {
        let t0 = lo + (hi - lo) * k as f64 / POLYGON_SAMPLES as f64;
        let t1 = lo + (hi - lo) * (k + 1) as f64 / POLYGON_SAMPLES as f64;
        let (a, b) = (curve.geometry.point_at(t0), curve.geometry.point_at(t1));
        span = span.max(a.distance_to(&b));
        for j in 1..4 {
            let m = curve.geometry.point_at(t0 + (t1 - t0) * j as f64 / 4.0);
            sagitta = sagitta.max(distance_to_segment(&m, &a, &b));
        }
    }
    // Beyond twice the polygon's sagitta from the boundary its verdict
    // cannot be wrong: it stands.
    if distance > 2.0 * sagitta {
        return None;
    }
    let (u, v) = project_point_to_uv(surface, &at)?;
    let across = surface.normal_at(u, v).cross(&tangent).normalized()?;
    let inside = |q: Point3| -> Option<bool> {
        let (qu, qv) = project_point_to_uv(surface, &q)?;
        brep.face_uv_trim_containment(face, surface, qu, qv)
            .map(|c| c == popkern_core::validation::FaceUvContainment::Interior)
    };
    // The polygon strays from the curve by at most its sagitta: probes
    // twice that far across are classified correctly. (A chord-length
    // step jumped clean across a 6-sampled bore and read the face on
    // both sides.)
    let step = 2.0 * sagitta + 1e-9 * span;
    let plus = inside(at + across * step)?;
    let minus = inside(at - across * step)?;
    let material = match (plus, minus) {
        (true, false) => across,
        (false, true) => -across,
        _ => return None,
    };
    // A point on the boundary itself (within the solver scale of the
    // span) is on the face.
    Some((*x - at).dot(&material) >= -1e-9 * span.max(1e-12))
}

/// Distance from `p` to segment `a`–`b`.
fn distance_to_segment(p: &Point3, a: &Point3, b: &Point3) -> f64 {
    let ab = *b - *a;
    let len2 = ab.dot(&ab);
    if len2 <= 0.0 {
        return p.distance_to(a);
    }
    let t = ((*p - *a).dot(&ab) / len2).clamp(0.0, 1.0);
    p.distance_to(&(*a + ab * t))
}
