//! Station queries on center curves (PK-139): closest-point lookup
//! used by the per-station variable-radius center march.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;

/// Closest point on a NURBS curve to a target: dense sampling with
/// two rounds of local ternary refinement (analytic offset
/// intersections here are lines/arcs/low-degree fits — 64 samples
/// bracket the minimum comfortably).
pub(crate) fn closest_point_on_nurbs(curve: &NurbsCurve, target: &Point3) -> Point3 {
    let t = closest_t_on_nurbs(curve, target);
    curve.point_at(t)
}

/// Parameter of the closest point on a NURBS curve to a target (same
/// dense-sample + ternary refinement as [`closest_point_on_nurbs`]).
pub(crate) fn closest_t_on_nurbs(curve: &NurbsCurve, target: &Point3) -> f64 {
    let (t0, t1) = curve.parameter_range();
    let mut best_t = 0.5 * (t0 + t1);
    let mut best_d = f64::MAX;
    let n = 64usize;
    for i in 0..=n {
        let t = t0 + (t1 - t0) * i as f64 / n as f64;
        let d = curve.point_at(t).distance_to(target);
        if d < best_d {
            best_d = d;
            best_t = t;
        }
    }
    let mut lo = (best_t - (t1 - t0) / n as f64).max(t0);
    let mut hi = (best_t + (t1 - t0) / n as f64).min(t1);
    for _ in 0..40 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if curve.point_at(m1).distance_to(target) < curve.point_at(m2).distance_to(target) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    0.5 * (lo + hi)
}

/// Select the intersection curve whose midpoint is nearest to the given point.
pub(crate) fn select_nearest_curve(
    intersections: &[popkern_core::geometry::intersection::surface_surface::SurfaceSurfaceIntersection],
    target: &Point3,
) -> Result<NurbsCurve, String> {
    use popkern_core::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;

    let mut best: Option<(f64, NurbsCurve)> = None;

    for isect in intersections {
        if let SurfaceSurfaceIntersection::Curve(cos) = isect {
            let curve = &cos.curve_3d;
            let (t0, t1) = curve.parameter_range();
            let mid_pt = curve.point_at((t0 + t1) * 0.5);
            let dist = mid_pt.distance_to(target);
            let is_closer = best.as_ref().is_none_or(|(d, _)| dist < *d);
            if is_closer {
                best = Some((dist, curve.clone()));
            }
        }
    }

    best.map(|(_, c)| c)
        .ok_or_else(|| "no curve intersections found in offset result".into())
}

/// Try the rolling-ball contact solver at one spine sample.
///
/// Translates the legacy `(target, orientation_*)` interface into a
/// `ContactRequest`: seeds the initial UV guesses by projecting
/// `target` onto each surface (`closest_point`), maps `orientation_i`
/// to `BallSide::Outside` (matching popkern's "face normals point into
/// the solid" convention — the rolling ball is on the side opposite
/// the face normal, which is the `+normal_at` side for an outward-
/// pointing face surface normal), and uses the segment from `target`
/// toward the surface mid-point as the spine plane normal.
///
/// Returns `Some(center)` on success, `None` on any failure (parallel
/// normals, degenerate, out-of-bounds, etc.). The caller falls back
/// to the offset-intersect path on `None`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn try_contact_solver(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    radius: f64,
    target: &Point3,
    tolerance: &Tolerance,
    orientation_1: bool,
    orientation_2: bool,
    convexity: crate::analysis::Convexity,
) -> Option<Point3> {
    use crate::contact::{ContactOptions, ContactRequest, solve_rolling_ball_contact};

    // Seed each surface's initial UV from target via closest_point.
    let seed_uv0 = (surface1 as &dyn SurfaceEval).closest_point(target)?;
    let seed_uv1 = (surface2 as &dyn SurfaceEval).closest_point(target)?;

    // Ball side from edge convexity (PK-091): inscribed for convex,
    // void-side for concave.
    let (side_0, side_1) =
        crate::contact::ball_sides_for_edge(convexity, orientation_1, orientation_2);

    // Spine plane: pin the ball center close to `target`. Walk the
    // segment from target to the surface midpoint as the plane normal.
    // For one-shot use this is approximate; phase 1b walker passes
    // the actual spine tangent.
    let p0 = surface1.point_at(seed_uv0.0, seed_uv0.1);
    let p1 = surface2.point_at(seed_uv1.0, seed_uv1.1);
    let mid = popkern_core::geometry::Point3::new(
        (p0.x + p1.x) * 0.5,
        (p0.y + p1.y) * 0.5,
        (p0.z + p1.z) * 0.5,
    );
    let mut spine_plane_normal = mid - *target;
    if spine_plane_normal.length() < tolerance.numerical_floor() {
        // Target is at the contact midpoint; pick any tangent direction.
        // Use the cross product of the two face partials — well-defined
        // when the surfaces are not parallel, which we know they are
        // not since we're in the cone branch.
        let n0 = surface1.normal_at(seed_uv0.0, seed_uv0.1);
        let n1 = surface2.normal_at(seed_uv1.0, seed_uv1.1);
        spine_plane_normal = n0.cross(&n1);
        if spine_plane_normal.length() < tolerance.numerical_floor() {
            return None;
        }
    }

    let request = ContactRequest {
        surface_0: surface1,
        surface_1: surface2,
        side_0,
        side_1,
        radius,
        initial_uv0: seed_uv0,
        uv0_range: (u1_range, v1_range),
        initial_uv1: seed_uv1,
        uv1_range: (u2_range, v2_range),
        spine_point: *target,
        spine_plane_normal,
    };
    let solution =
        solve_rolling_ball_contact(&request, &ContactOptions::default(), tolerance, None).ok()?;
    Some(solution.center)
}
