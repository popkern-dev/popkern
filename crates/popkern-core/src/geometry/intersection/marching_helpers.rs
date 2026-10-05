//! Small helpers for the marching loop: parameter wrapping, normal
//! estimation near singularities, tangent-plane basis construction,
//! and tangent-point direction search.
//!
//! Split out of `marching.rs` to keep the parent under the 500-line cap.

use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;
use std::f64::consts::TAU;

/// Whether a query window covers at least one full period of its
/// axis. The relative slack absorbs a full-circle face whose window
/// endpoints were computed a few ULPs short of the period.
fn window_covers_period(u_range: (f64, f64), period: f64) -> bool {
    let span = u_range.1 - u_range.0;
    period.is_finite() && period > 0.0 && span >= period * (1.0 - 1e-9)
}

/// Wrap a periodic parameter to the representative, modulo the TRUE
/// period, that lies in (or nearest to) the query window
/// `[u_range.0, u_range.1]`.
///
/// A wrap is only honest when it is an isometry of the surface: every
/// representative `u + k·period` is the same surface point. The
/// wrap therefore never uses any other modulus. Windows covering at
/// least one period (the boolean solvers pad 50%) land in
/// `[u_range.0, u_range.0 + period)`: wrapping modulo the padded span
/// let a cone trace run 1.5 laps with zero seam jumps (PK-460,
/// chamfered rod). Sub-period windows take the representative
/// nearest the window, so a sample just outside one edge stays just
/// outside THAT edge for [`exits_window`] to stop the march: the
/// former span-modulus wrap teleported it to the opposite edge — a
/// different surface point 0.2 away on a fillet torus corner — and
/// the march oscillated between the two preimages for its whole
/// 200k-step budget (PK-572).
///
/// Non-periodic axes (`period == None`), degenerate windows, and
/// non-finite inputs return `u` unchanged.
pub(super) fn wrap_parameter(u: f64, u_range: (f64, f64), period: Option<f64>) -> f64 {
    let Some(p) = period else {
        return u;
    };
    let span = u_range.1 - u_range.0;
    if span <= 0.0 || !u.is_finite() || !span.is_finite() || !p.is_finite() || p <= 0.0 {
        return u;
    }
    // Modular arithmetic instead of a while loop: bounded work for
    // extreme values.
    let wrapped = u_range.0 + (u - u_range.0).rem_euclid(p);
    if wrapped <= u_range.1 || window_covers_period(u_range, p) {
        return wrapped;
    }
    // In the gap (u_range.1, u_range.0 + p): pick the nearer edge's
    // representative.
    let below = wrapped - p;
    if u_range.0 - below < wrapped - u_range.1 {
        below
    } else {
        wrapped
    }
}

/// Whether a wrapped parameter has left its query window along an
/// axis the window bounds: every non-periodic axis, and periodic axes
/// whose window covers less than one period (a partial-arc face's
/// azimuth, a fillet torus's minor-angle strip). Full-period windows
/// on periodic axes never exit — the wrap keeps them inside.
pub(super) fn exits_window(u: f64, u_range: (f64, f64), period: Option<f64>) -> bool {
    if let Some(p) = period
        && window_covers_period(u_range, p)
    {
        return false;
    }
    u < u_range.0 || u > u_range.1
}

/// Estimate a surface normal near a singularity by averaging normals at offset points.
///
/// Samples normals at four offset positions (`±delta` in each parameter direction).
/// Returns the normalised average if at least one offset yields a non-degenerate
/// normal, otherwise returns `None`.
pub(super) fn estimate_normal_near_singularity(
    surface: &SurfaceGeometry,
    u: f64,
    v: f64,
    delta: f64,
) -> Option<Vec3> {
    let offsets = [(delta, 0.0), (-delta, 0.0), (0.0, delta), (0.0, -delta)];
    let mut sum = Vec3::zero();
    let mut count = 0usize;
    for (du, dv) in &offsets {
        let n = surface.normal_at(u + du, v + dv);
        if n.length() > 1e-10 {
            sum = sum + n;
            count += 1;
        }
    }
    if count > 0 {
        (sum * (1.0 / count as f64)).normalized()
    } else {
        None
    }
}

/// Build an orthonormal tangent basis in the plane perpendicular to `normal`.
///
/// Returns `(tu, tv)` where `tu` and `tv` are unit vectors perpendicular to
/// `normal` and to each other.
pub(super) fn tangent_basis(normal: &Vec3) -> (Vec3, Vec3) {
    let n = normal.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let axis = if n.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let tu = n
        .cross(&axis)
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let tv = n.cross(&tu);
    (tu, tv)
}

/// Find a march direction at a tangent point using perturbation sampling.
///
/// When `N1 × N2` is near-zero (surfaces are locally tangent), the standard
/// marching direction is undefined.  This function samples 8 directions in the
/// tangent plane of surface 1's normal and picks the one that minimises the
/// distance to `surface2.point_at(u2, v2)` — i.e., the direction that stays
/// closest to the second surface.
///
/// Returns `None` if the normal at `(u1, v1)` is degenerate or all sampled
/// directions are equidistant.
#[allow(clippy::too_many_arguments)]
pub(super) fn find_tangent_march_direction(
    surface1: &SurfaceGeometry,
    u1: f64,
    v1: f64,
    surface2: &SurfaceGeometry,
    u2: f64,
    v2: f64,
    point: &Point3,
    tolerance: &Tolerance,
) -> Option<Vec3> {
    let n1 = surface1.normal_at(u1, v1);
    if n1.length() < tolerance.numerical_floor() {
        return None;
    }
    let delta = tolerance.linear() * 10.0;
    let (tu, tv) = tangent_basis(&n1);

    let mut best_dir = Vec3::zero();
    let mut min_dist = f64::MAX;

    for i in 0..8usize {
        let angle = TAU * i as f64 / 8.0;
        let dir = tu * angle.pcos() + tv * angle.psin();
        let test_point = *point + dir * delta;
        let dist = (test_point - surface2.point_at(u2, v2)).length();
        if dist < min_dist {
            min_dist = dist;
            best_dir = dir;
        }
    }

    best_dir.normalized()
}

/// Wrap a march SEED's periodic coordinates onto the principal
/// branch of each window, by the surface's REAL period (analytic or
/// geometric chart closure — never the window span, which relocates
/// the seed geometrically on a sub-period patch: a grazing
/// strip × cap seed moved by the 0.6-rad window span and the tangent
/// contact vanished in one operand order). A seed expressed on
/// another branch (u = −0.09 on a period-1 closed chart) differs
/// from its first marched neighbour by a full period; the seam
/// splitter reads that as a spurious crossing, the closed-trace
/// rejoin glues the wrong-branch point onto the tail, and the
/// connector piece strands (PK-454: the z = 0.55 band grind served a
/// 24.6% volume lie).
#[allow(clippy::too_many_arguments)]
pub(super) fn principal_branch_seed(
    start: &super::MarchSample,
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
) -> super::MarchSample {
    let real = |analytic: Option<f64>, chart: Option<f64>| analytic.or(chart);
    let (a1u, a1v) = surface1.periodicity();
    let (c1u, c1v) = crate::geometry::chart_periods::geometric_chart_periods(surface1);
    let (a2u, a2v) = surface2.periodicity();
    let (c2u, c2v) = crate::geometry::chart_periods::geometric_chart_periods(surface2);
    let wrap = |x: f64, range: (f64, f64), p: Option<f64>| -> f64 {
        match p {
            Some(p) if p.is_finite() && p > 0.0 && x.is_finite() => {
                range.0 + (x - range.0).rem_euclid(p)
            }
            _ => x,
        }
    };
    let mut s = start.clone();
    s.u1 = wrap(s.u1, u1_range, real(a1u, c1u));
    s.v1 = wrap(s.v1, v1_range, real(a1v, c1v));
    s.u2 = wrap(s.u2, u2_range, real(a2u, c2u));
    s.v2 = wrap(s.v2, v2_range, real(a2v, c2v));
    s
}

/// Whether a forward march trace already closed on itself (endpoints
/// within twice the average step). Traces under the march topic.
pub(super) fn forward_trace_closed(forward: &[super::MarchSample], tolerance: &Tolerance) -> bool {
    forward.len() > 20
        && {
            let total_arc: f64 = forward
                .windows(2)
                .map(|w| w[0].point.distance_to(&w[1].point))
                .sum();
            let avg_step = total_arc / (forward.len() - 1) as f64;
            let tol = (avg_step * 2.0).max(tolerance.linear());
            let gap = forward[0]
                .point
                .distance_to(&forward[forward.len() - 1].point);
            if crate::trace::trace_enabled("march") {
                eprintln!(
                    "  [march] fwd {} pts arc {total_arc:.4} avg_step {avg_step:.5} tol {tol:.5} gap {gap:.5} -> closed={}",
                    forward.len(),
                    gap < tol
                );
            }
            gap < tol
        }
}

/// Recover a usable surface normal at a singular evaluation during the
/// march (split from `march_one_direction` under the 500-line file
/// budget): try direct estimation near the point, then walk `u`
/// forward with a halving step up to three times. Returns the normal
/// and the possibly reduced step size, or `None` when unrecoverable.
pub(super) fn recover_singular_normal(
    surface: &crate::geometry::surface::SurfaceGeometry,
    u: f64,
    v: f64,
    step_size: f64,
    tolerance: &crate::tolerance::Tolerance,
) -> Option<(crate::geometry::Vec3, f64)> {
    let delta = tolerance.linear() * 10.0;
    if let Some(estimated) = estimate_normal_near_singularity(surface, u, v, delta) {
        return Some((estimated, step_size));
    }
    let mut reduced_step = step_size * 0.5;
    for _attempt in 0..3 {
        if let Some(en) = estimate_normal_near_singularity(surface, u + reduced_step, v, delta) {
            return Some((en, reduced_step));
        }
        reduced_step *= 0.5;
    }
    None
}
