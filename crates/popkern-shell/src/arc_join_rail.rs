//! Offset images by the ORIGINAL face's normal field (PK-883).
//!
//! An offset surface is by definition `S + s·d·n`, so the offset image
//! of a point `p` on `S` is `p + s·d·n(p)` — the point where `p`'s
//! normal line meets the offset surface. Closest-point projection onto
//! the offset surface finds it only while `p` lies short of the offset
//! surface's focal point: past half the vanish distance (the wavy
//! plate's crest at d = 0.8, vanish ≈ 1) the image is a local MAXIMUM
//! of the distance across the crest, the projection lands on another
//! point, and rails, corner vertices and corner spheres sat up to
//! 2.6e-4 off their canals. [`offset_image`] solves the normal line ×
//! offset surface crossing by Newton from the projection's estimate,
//! needing no offset distance.

use crate::curved_recompute::project_to_surface;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, Vec3};

/// Newton steps for the normal-line crossing.
const IMAGE_STEPS: usize = 30;

/// The offset image of `p` (a point on `original`) on `offset`: the
/// crossing of `p`'s normal line with the offset surface, nearest the
/// projection's estimate. Falls back to the closest-point projection
/// when the crossing does not converge (or `original` has no normal
/// at `p`).
pub(crate) fn offset_image(
    original: &SurfaceGeometry,
    offset: &SurfaceGeometry,
    p: &Point3,
) -> Option<Point3> {
    let projected = project_to_surface(offset, p)?;
    Some(normal_crossing(original, offset, p, projected).unwrap_or(projected))
}

fn normal_crossing(
    original: &SurfaceGeometry,
    offset: &SurfaceGeometry,
    p: &Point3,
    projected: Point3,
) -> Option<Point3> {
    let (u0, v0) = project_point_to_uv(original, p)?;
    let n = original.normal_at(u0, v0).normalized()?;
    let mut t = (projected - *p).dot(&n);
    let scale = 1.0 + t.abs() + (*p - Point3::new(0.0, 0.0, 0.0)).length();
    if t.abs() <= 1e-12 * scale {
        return None;
    }
    // Seed beside the image: its projection is in the right basin.
    let (mut u, mut v) = project_point_to_uv(offset, &(*p + n * t))?;
    let m = n * -1.0;
    for _ in 0..IMAGE_STEPS {
        let r = offset.point_at(u, v) - (*p + n * t);
        if r.length() <= 1e-13 * scale {
            break;
        }
        let (su, sv) = (offset.partial_u(u, v), offset.partial_v(u, v));
        let det = su.dot(&sv.cross(&m));
        if det.abs() <= f64::MIN_POSITIVE {
            return None;
        }
        let rhs = r * -1.0;
        u += rhs.dot(&sv.cross(&m)) / det;
        v += su.dot(&rhs.cross(&m)) / det;
        t += su.dot(&sv.cross(&rhs)) / det;
    }
    let s = offset.point_at(u, v);
    (s.distance_to(&(*p + n * t)) <= 1e-10 * scale).then_some(s)
}

fn segment_distance(q: &Point3, a: &Point3, b: &Point3) -> f64 {
    let ab: Vec3 = *b - *a;
    let len2 = ab.dot(&ab);
    let t = if len2 > 0.0 {
        ((*q - *a).dot(&ab) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    q.distance_to(&(*a + ab * t))
}

/// The signed normal step that carries points of `original` onto
/// `offset`, solved at `probe` ([`offset_image`]'s crossing): an offset
/// surface displaces uniformly, so one solve serves a whole rail and
/// each further station is one original-normal evaluation.
pub(crate) fn normal_step(
    original: &SurfaceGeometry,
    offset: &SurfaceGeometry,
    probe: &Point3,
) -> Option<f64> {
    let n = unit_normal(original, probe)?;
    Some((offset_image(original, offset, probe)? - *probe).dot(&n))
}

/// `p` stepped `step` along `original`'s unit normal.
pub(crate) fn normal_station(original: &SurfaceGeometry, p: &Point3, step: f64) -> Option<Point3> {
    Some(*p + unit_normal(original, p)? * step)
}

fn unit_normal(surface: &SurfaceGeometry, p: &Point3) -> Option<Vec3> {
    let (u, v) = project_point_to_uv(surface, p)?;
    surface.normal_at(u, v).normalized()
}

/// Smallest spine-fraction interval [`adaptive_stations`] still splits.
const MIN_SPLIT: f64 = 1.0 / 65536.0;

/// Stations along `0..=1` whose chords hold the station curve within
/// `sag`: an interval splits while its midpoint or quarter-point image
/// sits off the chord by more than `sag`. Uniform stations under-sample
/// the fold: past half the vanish distance the image curve turns
/// sharply over the crest (129 uniform stations left 1.8e-4 of sag on
/// the wavy plate at d = 0.9), and the split concentrates there.
pub(crate) fn adaptive_stations(
    station: impl Fn(f64) -> Option<Point3>,
    sag: f64,
) -> Option<Vec<Point3>> {
    let mut out = vec![station(0.0)?];
    let mut stack: Vec<(f64, f64, Point3, Point3)> = Vec::new();
    let seeds: Vec<(f64, Point3)> = (0..=8)
        .map(|i| {
            let f = i as f64 / 8.0;
            station(f).map(|q| (f, q))
        })
        .collect::<Option<Vec<_>>>()?;
    for w in seeds.windows(2).rev() {
        stack.push((w[0].0, w[1].0, w[0].1, w[1].1));
    }
    while let Some((f0, f1, q0, q1)) = stack.pop() {
        let fm = 0.5 * (f0 + f1);
        let qm = station(fm)?;
        let off = [0.25, 0.75]
            .into_iter()
            .map(|k| station(f0 + k * (f1 - f0)).map(|q| segment_distance(&q, &q0, &q1)))
            .collect::<Option<Vec<_>>>()?
            .into_iter()
            .fold(segment_distance(&qm, &q0, &q1), f64::max);
        if off > sag && f1 - f0 > MIN_SPLIT {
            stack.push((fm, f1, qm, q1));
            stack.push((f0, fm, q0, qm));
        } else {
            out.push(q1);
        }
    }
    Some(out)
}
