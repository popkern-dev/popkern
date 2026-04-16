//! Ray-vs-edge-curve intersection helpers used by the
//! distance+angle chamfer path
//! ([`crate::chamfer_2d_by_angle::chamfer_corner_2d_by_angle`]).
//!
//! Given a planar chamfer ray `t1 + r · ray_dir` and a supporting
//! edge's curve, each helper returns the ray parameter `r` at
//! the first forward crossing. Three curve types are supported:
//!
//! - **Line** — 2D Cramer solve via in-plane scalar cross
//!   product. Closed form, O(1).
//! - **Arc** — `|t1 + r · ray_dir − center|² = radius²` scalar
//!   quadratic. Closed form, O(1). Picks the smallest forward
//!   root on the supporting circle; the caller verifies the hit
//!   point lies on the edge's parameter range via
//!   [`crate::corner_2d_helpers::setback_from_corner`].
//! - **NURBS** — adaptive uniform sampling (64 points) to
//!   collect **all** sign-change brackets of
//!   `((C(t) − t1) × ray_dir) · plane_normal`, bisect each to
//!   refine, then return the **smallest positive** ray parameter.
//!   The multi-bracket approach guarantees the "nearest forward
//!   hit" semantics for exotic S-curves and self-intersecting
//!   projections, at the cost of one additional pass through
//!   the sample array.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle};

/// In-plane scalar cross product: `(u × v) · n`. Positive when
/// the rotation from `u` to `v` is CCW around `n`.
pub(crate) fn perp_cross(u: Vec3, v: Vec3, n: Vec3) -> f64 {
    u.cross(&v).dot(&n)
}

/// 2D Cramer solve for the ray parameter `r` where `t1 + r ·
/// ray_dir` meets the Line curve stored on `edge`. All vectors
/// live in the plane with normal `plane_normal`.
pub(crate) fn ray_line_edge_intersection(
    brep: &BRep,
    edge: EdgeHandle,
    t1: Point3,
    ray_dir: Vec3,
    plane_normal: Vec3,
) -> Result<f64, String> {
    let curve = brep.edges[edge]
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {edge:?} has no curve geometry"))?;
    let (origin, direction) = match &curve.geometry {
        CurveGeometry::Line { origin, direction } => (*origin, *direction),
        other => {
            return Err(format!(
                "ray_line_edge_intersection on non-line curve: {other:?}"
            ));
        }
    };
    let a = origin - t1;
    let denom = perp_cross(ray_dir, direction, plane_normal);
    if denom.abs() < 1.0e-12 {
        return Err("chamfer ray is parallel to edge_b line".into());
    }
    let r = perp_cross(a, direction, plane_normal) / denom;
    if r < -1.0e-9 {
        return Err(format!("ray-line intersection parameter negative ({r:.6})"));
    }
    Ok(r.max(0.0))
}

/// Solve `|t1 + r · ray_dir − center|² = radius²` for `r`,
/// picking the smallest non-negative root. Assumes `|ray_dir| =
/// 1`.
pub(crate) fn ray_arc_edge_intersection(
    brep: &BRep,
    edge: EdgeHandle,
    t1: Point3,
    ray_dir: Vec3,
) -> Result<f64, String> {
    let curve = brep.edges[edge]
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {edge:?} has no curve geometry"))?;
    let (center, radius) = match &curve.geometry {
        CurveGeometry::Arc { center, radius, .. } => (*center, *radius),
        other => {
            return Err(format!(
                "ray_arc_edge_intersection on non-arc curve: {other:?}"
            ));
        }
    };
    let u = t1 - center;
    let b = u.dot(&ray_dir);
    let c = u.dot(&u) - radius * radius;
    let disc = b * b - c;
    if disc < 0.0 {
        return Err(format!(
            "ray-arc intersection has no real solution (discriminant = {disc:.6})"
        ));
    }
    let sqrt_d = disc.sqrt();
    let r_near = -b - sqrt_d;
    let r_far = -b + sqrt_d;
    if r_near >= -1.0e-9 {
        Ok(r_near.max(0.0))
    } else if r_far >= -1.0e-9 {
        Ok(r_far.max(0.0))
    } else {
        Err(format!(
            "both ray-arc roots lie behind the chamfer ray origin (near={r_near:.6}, far={r_far:.6})"
        ))
    }
}

/// Find the first **forward** crossing of the chamfer ray with
/// the NURBS curve on `edge`. Collects every sign-change bracket
/// of `f(t) = ((C(t) − t1) × ray_dir) · plane_normal`, bisects
/// each, and returns the smallest non-negative ray parameter
/// across all refined hits.
///
/// This is the multi-bracket variant: robust against exotic
/// S-curves, inflection points, and near-tangential crossings
/// where a single-bracket search could pick up a "later" root
/// simply because the sampler encountered it first in parameter
/// order. The cost is one sweep of `SAMPLES` curve evaluations
/// plus up to ~60 bisection steps per bracket (typically 1–2
/// brackets in practice).
pub(crate) fn ray_nurbs_edge_intersection(
    brep: &BRep,
    edge: EdgeHandle,
    t1: Point3,
    ray_dir: Vec3,
    plane_normal: Vec3,
) -> Result<f64, String> {
    let curve = brep.edges[edge]
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {edge:?} has no curve geometry"))?;
    let nurbs = match &curve.geometry {
        CurveGeometry::Nurbs(n) => n,
        other => {
            return Err(format!(
                "ray_nurbs_edge_intersection on non-nurbs curve: {other:?}"
            ));
        }
    };

    // Signed planar distance from the ray line to `C(t)`.
    let f = |t: f64| -> f64 {
        let pt = nurbs.point_at(t);
        (pt - t1).cross(&ray_dir).dot(&plane_normal)
    };

    const SAMPLES: usize = 64;
    let t_lo = curve.t_start.min(curve.t_end);
    let t_hi = curve.t_start.max(curve.t_end);
    if (t_hi - t_lo).abs() < 1.0e-14 {
        return Err("ray-NURBS: degenerate edge parameter range".into());
    }
    let dt = (t_hi - t_lo) / (SAMPLES as f64);

    // Collect every (t, f(t)) sample up front so we can scan for
    // ALL sign-change brackets in one pass.
    let mut samples: Vec<(f64, f64)> = Vec::with_capacity(SAMPLES + 1);
    for i in 0..=SAMPLES {
        let t = t_lo + (i as f64) * dt;
        samples.push((t, f(t)));
    }

    // Every sign-change bracket is a candidate hit. Refine each
    // and compute the ray parameter, then pick the smallest
    // non-negative one.
    let mut best_r: Option<f64> = None;
    for window in samples.windows(2) {
        let (t_a, f_a) = window[0];
        let (t_b, f_b) = window[1];
        // Degenerate / exact-zero cases: treat either endpoint
        // as a direct hit.
        let hit_t = if f_a == 0.0 {
            Some(t_a)
        } else if f_b == 0.0 {
            Some(t_b)
        } else if f_a.signum() != f_b.signum() {
            // Real sign change — bisect.
            Some(refine_bracket(&f, t_a, t_b, f_a, f_b))
        } else {
            None
        };
        if let Some(t_hit) = hit_t {
            let pt_hit = nurbs.point_at(t_hit);
            let r = (pt_hit - t1).dot(&ray_dir);
            if r >= -1.0e-9 {
                let r_clamped = r.max(0.0);
                best_r = Some(match best_r {
                    Some(prev) if prev <= r_clamped => prev,
                    _ => r_clamped,
                });
            }
        }
    }

    best_r.ok_or_else(|| {
        "ray-NURBS: no forward intersection found (ray misses curve or only has backward hits)"
            .to_string()
    })
}

/// Bisect a bracket `[t_a, t_b]` with `f(t_a), f(t_b)` of
/// opposite sign until the interval shrinks below tolerance or
/// we hit the iteration cap. Returns the midpoint of the final
/// bracket.
fn refine_bracket(f: &dyn Fn(f64) -> f64, t_a: f64, t_b: f64, f_a: f64, f_b: f64) -> f64 {
    let mut lo = t_a;
    let mut hi = t_b;
    let mut f_lo = f_a;
    let _f_hi = f_b;
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        let f_mid = f(mid);
        if f_mid.abs() < 1.0e-12 || (hi - lo).abs() < 1.0e-14 {
            return mid;
        }
        if f_lo.signum() != f_mid.signum() {
            hi = mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
    }
    0.5 * (lo + hi)
}
