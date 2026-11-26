//! A closed curve round a touch is the touch (W3 E cells, tangency
//! design §4.8).
//!
//! Where two surfaces touch at a point and curve apart (relative
//! curvature definite), a marcher can return a small closed curve round
//! the touch: the band where the gap is under the tolerance, a loop
//! about `√(2ε/λ)` across. It is no crossing: the gap's critical point
//! inside it lies within the tolerance, so at that resolution the
//! surfaces touch there and nowhere else. Read as a curve it split each
//! face of a boolean off a disc a tolerance deep, which the union kept
//! on both sides and the general policy glued rim to rim (CC-E 9 at
//! scale 1e-3: a loop 1.6e-4 across round the touch), and it made the
//! self-intersection check read the glued body's tangent faces as
//! crossing. A real small lens, a tolerance deep or deeper, keeps its
//! loop: its critical gap is outside the tangency band
//! (`Tolerance::tangency_band`, the one rule for a touch decided by
//! distance).

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval as _;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::intersection::contact_classify::{ContactScale, SurfaceContactKind};
use crate::geometry::intersection::contact_discovery::{ParamWindow, contact_from};
use crate::geometry::surface::{SurfaceEval as _, SurfaceGeometry};
use crate::geometry::uv_projection::project_point_to_uv;

/// Samples along the curve.
const SAMPLES: usize = 16;

/// Whether `curve`, an intersection of `s1` over `w1` and `s2` over
/// `w2`, is the tolerance band round a touch (see the module doc): every
/// point of it within twice `√(2·reach/λ)` of a definite touch whose gap
/// lies inside the tangency band of `linear` (λ the weaker eigenvalue of the relative
/// curvature there, which sets how far the band reaches). Closed, or cut open by a face's
/// window where the touch lies beside its boundary (CC-E 9: 8e-6 from
/// a cone's rim): a definite touch has no crossing that near it.
pub fn is_touch_band(
    curve: &CurveOnSurface,
    (s1, w1): (&SurfaceGeometry, ParamWindow),
    (s2, w2): (&SurfaceGeometry, ParamWindow),
    (linear, reach, extent): (f64, f64, f64),
) -> bool {
    let c = &curve.curve_3d;
    let (t0, t1) = c.parameter_range();
    let pts: Vec<Point3> = (0..=SAMPLES)
        .map(|k| c.point_at(t0 + (t1 - t0) * k as f64 / SAMPLES as f64))
        .collect();
    points_in_touch_band(&pts, (s1, w1), (s2, w2), (linear, reach, extent))
}

/// Whether `pts`, points where `s1` and `s2` meet within the tolerance,
/// all lie in the band round one definite touch (see [`is_touch_band`]):
/// a sampled meeting set read as crossing points (the validator's
/// cylinder-pair sampler) is a touch when they do.
pub fn points_in_touch_band(
    pts: &[Point3],
    (s1, w1): (&SurfaceGeometry, ParamWindow),
    (s2, w2): (&SurfaceGeometry, ParamWindow),
    (linear, reach, extent): (f64, f64, f64),
) -> bool {
    if pts.is_empty() {
        return false;
    }
    let n = pts.len() as f64;
    let centre = pts.iter().fold(Point3::new(0.0, 0.0, 0.0), |s, p| {
        Point3::new(s.x + p.x / n, s.y + p.y / n, s.z + p.z / n)
    });
    // The band is convex round the touch (a window cuts it, never
    // hollows it), so its points' centroid lies inside, where the
    // surfaces stand within the touch's critical gap of each other: a
    // crossing curve's centroid sits off one surface by its size (a
    // cylinder's circle on a plane, by the radius). Two projections
    // decide what the Newton search from five seeds would, for every
    // curve of every face pair (the six-axis arm's interference ran 4x
    // its CPU without it).
    let near = |s: &SurfaceGeometry| project_point_to_uv(s, &centre).map(|(u, v)| s.point_at(u, v));
    match (near(s1), near(s2)) {
        (Some(a), Some(b)) if a.distance_to(&b) <= 4.0 * reach => {}
        _ => return false,
    }
    let step = 1e-6 * (w1.0.1 - w1.0.0).abs().max(w1.1.1 - w1.1.0).max(1e-9);
    let scale = ContactScale {
        linear,
        extent: extent.max(linear),
    };
    // Newton from the middle, else from points of the curve: an
    // elongated band's middle can sit where the gap's Hessian misleads.
    let seeds = std::iter::once(centre).chain(pts.iter().step_by((pts.len() / 4).max(1)).copied());
    let touch = seeds
        .filter_map(|p| project_point_to_uv(s1, &p))
        .map(|seed| into_window(s1, w1, seed))
        .find_map(|seed| contact_from((s1, w1), (s2, w2), seed, step, scale));
    let Some(t) = touch else {
        return false;
    };
    // The band reaches √(2·reach/λ) from the touch along each principal
    // direction of the relative curvature λ: the weakest bounds it.
    let weakest = t
        .classification
        .eigenvalues
        .map(|[a, b]| a.abs().min(b.abs()))
        .filter(|l| *l > 0.0);
    let Some(weakest) = weakest else {
        return false;
    };
    let band = 2.0 * (2.0 * reach / weakest).sqrt();
    // A touch by distance lies inside the open tangency band (shallower
    // than the tolerance less a ten-thousandth): a lens one whole
    // tolerance deep is a secant, its rim kept, and a common that thin
    // is a sub-tolerance feature (the scripted thin lens, two balls
    // overlapping 1e-6, served empty without its warning when its rim
    // was read as a touch's band).
    let magnitude = t.point.x.abs().max(t.point.y.abs()).max(t.point.z.abs());
    let tangency = crate::tolerance::Tolerance::default()
        .with_linear(linear)
        .tangency_band(magnitude);
    t.gap.abs() < tangency
        && matches!(
            t.classification.kind,
            SurfaceContactKind::TangentPoint { .. }
        )
        && pts.iter().all(|p| p.distance_to(&t.point) <= band)
}

/// `uv` shifted by whole periods of `s` into window `w` where it fits:
/// a projection answers in the chart's base period, and a face's window
/// may sit a period over (a cylinder's u from 5.4 to 11.7).
fn into_window(s: &SurfaceGeometry, w: ParamWindow, (u, v): (f64, f64)) -> (f64, f64) {
    let (pu, pv) = crate::geometry::chart::ChartDescriptor::evaluation_periods_of(s);
    let shift = |x: f64, (lo, hi): (f64, f64), p: Option<f64>| match p {
        Some(p) if p > 0.0 && (x < lo || x > hi) => {
            let k = ((lo - x) / p).ceil();
            let y = x + k * p;
            if y <= hi { y } else { x }
        }
        _ => x,
    };
    (shift(u, w.0, pu), shift(v, w.1, pv))
}

/// The diagonal in space of `s` over `w`.
pub fn window_extent(s: &SurfaceGeometry, ((u0, u1), (v0, v1)): ParamWindow) -> f64 {
    s.point_at(u0, v0).distance_to(&s.point_at(u1, v1))
}

/// [`is_touch_band`] for an intersection of `s1` and `s2` over their
/// windows, judged at `tolerance`: a curve that is a touch's band.
pub fn is_touch_ring_of(
    isect: &crate::geometry::intersection::SurfaceSurfaceIntersection,
    (s1, w1): (&SurfaceGeometry, ParamWindow),
    (s2, w2): (&SurfaceGeometry, ParamWindow),
    tolerance: &crate::tolerance::Tolerance,
) -> bool {
    let crate::geometry::intersection::SurfaceSurfaceIntersection::Curve(c) = isect else {
        return false;
    };
    let extent = window_extent(s1, w1).max(window_extent(s2, w2));
    is_touch_band(
        c,
        (s1, w1),
        (s2, w2),
        (tolerance.linear(), tolerance.linear_reach(), extent),
    )
}

#[cfg(test)]
#[path = "touch_ring_tests.rs"]
mod tests;
