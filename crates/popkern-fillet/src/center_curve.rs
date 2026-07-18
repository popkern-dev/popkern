//! Rolling-ball center curve computation for fillet operations.
//!
//! The center curve is the locus of the center of a rolling ball (sphere)
//! that is tangent to both adjacent surfaces. It is computed by offsetting
//! both surfaces inward by the fillet radius and intersecting the offsets.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::intersection::surface_surface::intersect_surfaces;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;

use crate::center_curve_helpers::{
    clamp_finite_range, clip_footprint_to_bounds, project_onto_surface,
};
use crate::offset::offset_surface;

// Re-export commonly used items from helpers so external callers
// (e.g. blend::surface tests) keep working with crate::center_curve::.
pub use crate::center_curve_helpers::fit_nurbs_through_points;
pub use crate::center_curve_walker::walker_based_center_curve;

/// Result of a center curve computation.
#[derive(Debug, Clone)]
pub struct CenterCurveResult {
    /// The rolling-ball center curve in 3D space.
    pub center_curve: NurbsCurve,
    /// Footprint points on the first surface (closest points from center curve).
    pub footprint_1: Vec<Point3>,
    /// Footprint points on the second surface (closest points from center curve).
    pub footprint_2: Vec<Point3>,
    /// Per-sample contact solver outputs from the spine walker, when
    /// the producer used the walker path. `None` for the legacy
    /// offset-intersect paths (which derive footprints from a
    /// surface intersection rather than per-sample contact solves).
    ///
    /// When populated, downstream code can route through
    /// `g2_cross_section::g2_sweep_from_contact_solutions` to build
    /// a degree-5 Bezier (G2) sweep instead of a circular-arc (G1)
    /// sweep. See §2.4 / §7 G2.
    pub contact_solutions: Option<Vec<crate::contact::ContactSolution>>,
}

/// Compute the center curve for a constant-radius fillet between two surfaces.
///
/// Offsets both surfaces inward by `radius`, intersects the offset surfaces,
/// and extracts the intersection curve nearest to `edge_midpoint`. Footprint
/// points are computed by projecting sampled center-curve points back onto
/// the original surfaces.
#[allow(clippy::too_many_arguments)]
pub fn constant_radius_center_curve(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    radius: f64,
    edge_midpoint: &Point3,
    tolerance: &Tolerance,
    orientation_1: bool,
    orientation_2: bool,
    convexity: crate::analysis::Convexity,
) -> Result<CenterCurveResult, String> {
    // Offset both surfaces toward the BALL side — inscribed in the
    // material for convex edges, void-side for concave — using the
    // same convexity-aware mapping as the contact solver.
    // `orientation_*` here are OUTWARD-honest (stored flag folded
    // with the measured shell sign by the caller); the former
    // `!orientation` flip encoded the legacy "face normals point
    // INTO the solid" convention that only held on primitive
    // constructor output.
    let (side_1, side_2) =
        crate::contact::ball_sides_for_edge(convexity, orientation_1, orientation_2);
    let inward_1 = matches!(side_1, crate::contact::BallSide::Inside);
    let inward_2 = matches!(side_2, crate::contact::BallSide::Inside);
    if std::env::var_os("PK_CONVEXITY_TRACE").is_some() {
        eprintln!(
            "    [offset-const] {convexity:?} o=({orientation_1},{orientation_2}) inward=({inward_1},{inward_2})"
        );
    }
    let offset1 = offset_surface(surface1, radius, inward_1, u1_range, v1_range, tolerance)?;
    let offset2 = offset_surface(surface2, radius, inward_2, u2_range, v2_range, tolerance)?;

    // 2. Compute parameter ranges for offset surfaces.
    let ((ou1_min, ou1_max), (ov1_min, ov1_max)) = offset1.parameter_range();
    let off1_u = clamp_finite_range((ou1_min, ou1_max), u1_range);
    let off1_v = clamp_finite_range((ov1_min, ov1_max), v1_range);
    let ((ou2_min, ou2_max), (ov2_min, ov2_max)) = offset2.parameter_range();
    let off2_u = clamp_finite_range((ou2_min, ou2_max), u2_range);
    let off2_v = clamp_finite_range((ov2_min, ov2_max), v2_range);

    // 3. Intersect offset surfaces.
    let result = intersect_surfaces(
        &offset1, off1_u, off1_v, &offset2, off2_u, off2_v, tolerance,
    )
    .map_err(|e| format!("offset intersection failed: {e}"))?;

    // §12.3 blend-on-blend fallback: if the analytic offset
    // surfaces don't intersect (common when one surface is a
    // Torus whose analytic intersection isn't supported by the
    // solver), retry by converting both to NURBS first. The
    // NURBS-NURBS intersection path handles arbitrary pairs via
    // Bezier decomposition + clipping.
    let result = if result.intersections.is_empty() {
        let nurbs1 = offset1.to_nurbs_in_range(off1_u, off1_v);
        let nurbs2 = offset2.to_nurbs_in_range(off2_u, off2_v);
        let n1_geo = popkern_core::geometry::surface::SurfaceGeometry::Nurbs(nurbs1.clone());
        let n2_geo = popkern_core::geometry::surface::SurfaceGeometry::Nurbs(nurbs2.clone());
        let n1_range = nurbs1.parameter_range();
        let n2_range = nurbs2.parameter_range();
        intersect_surfaces(
            &n1_geo, n1_range.0, n1_range.1, &n2_geo, n2_range.0, n2_range.1, tolerance,
        )
        .map_err(|e| format!("NURBS fallback offset intersection failed: {e}"))?
    } else {
        result
    };

    if result.intersections.is_empty() {
        return Err(
            "no intersection found between offset surfaces (analytic + NURBS fallback both empty)"
                .into(),
        );
    }

    // 4. Extract the curve nearest to the edge midpoint.
    let center_curve = select_nearest_curve(&result.intersections, edge_midpoint)?;

    // 5. Compute footprint points by projecting center-curve samples onto originals.
    let n_samples = 20;
    let (t_min, t_max) = center_curve.parameter_range();
    let mut footprint_1 = Vec::with_capacity(n_samples);
    let mut footprint_2 = Vec::with_capacity(n_samples);

    for i in 0..n_samples {
        let t = t_min + (t_max - t_min) * (i as f64) / ((n_samples - 1).max(1) as f64);
        let center_pt = center_curve.point_at(t);
        footprint_1.push(project_onto_surface(
            surface1, u1_range, v1_range, &center_pt,
        )?);
        footprint_2.push(project_onto_surface(
            surface2, u2_range, v2_range, &center_pt,
        )?);
    }

    // 6. Clip footprint points to face UV bounds so that trim curves
    //    do not extend beyond the face domain.
    footprint_1 = clip_footprint_to_bounds(surface1, u1_range, v1_range, &footprint_1);
    footprint_2 = clip_footprint_to_bounds(surface2, u2_range, v2_range, &footprint_2);

    Ok(CenterCurveResult {
        center_curve,
        footprint_1,
        footprint_2,
        contact_solutions: None,
    })
}

/// Compute the center curve for a variable-radius fillet between two surfaces.
///
/// Samples the edge at `n_samples` points, computes the local radius at each,
/// finds center points via Newton iteration, and fits a NURBS curve through them.
#[allow(clippy::too_many_arguments)]
pub fn variable_radius_center_curve(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    radius_start: f64,
    radius_end: f64,
    edge_anchor_at: &dyn Fn(f64) -> Point3,
    tolerance: &Tolerance,
    orientation_1: bool,
    orientation_2: bool,
    convexity: crate::analysis::Convexity,
) -> Result<CenterCurveResult, String> {
    parametric_radius_center_curve(
        surface1,
        u1_range,
        v1_range,
        surface2,
        u2_range,
        v2_range,
        |t| radius_start + (radius_end - radius_start) * t,
        edge_anchor_at,
        tolerance,
        orientation_1,
        orientation_2,
        convexity,
    )
}

/// Build a center curve for an arbitrary local-radius profile
/// `radius_at(t)` where `t` ranges over `[0, 1]` along the sharp
/// edge. Round 5l: ByPoints uses this with a piecewise-linear
/// closure so the full profile (not just the endpoint radii)
/// drives the rolling-ball center curve and footprints.
#[allow(clippy::too_many_arguments)]
pub fn parametric_radius_center_curve(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    radius_at: impl Fn(f64) -> f64,
    edge_anchor_at: &dyn Fn(f64) -> Point3,
    tolerance: &Tolerance,
    orientation_1: bool,
    orientation_2: bool,
    convexity: crate::analysis::Convexity,
) -> Result<CenterCurveResult, String> {
    let n_samples = 20usize;
    let mut center_points = Vec::with_capacity(n_samples);
    let mut footprint_1 = Vec::with_capacity(n_samples);
    let mut footprint_2 = Vec::with_capacity(n_samples);

    for i in 0..n_samples {
        let t = i as f64 / (n_samples - 1).max(1) as f64;
        let local_radius = radius_at(t);

        // Anchor each sample at ITS OWN station along the sharp edge.
        // A single shared anchor (the historic edge midpoint) made
        // every sample resolve the same nearest center: constant-law
        // samples all coincided ("degenerate: all points coincide")
        // and varying laws swept a radial spur instead of a curve
        // along the edge — the Variable ratio matrix completed with
        // OPEN bodies (PK-139 baseline).
        let anchor = edge_anchor_at(t);
        let center_pt = find_center_point_at_radius(
            surface1,
            u1_range,
            v1_range,
            surface2,
            u2_range,
            v2_range,
            local_radius,
            &anchor,
            tolerance,
            orientation_1,
            orientation_2,
            convexity,
        )?;

        let fp1 = project_onto_surface(surface1, u1_range, v1_range, &center_pt)?;
        let fp2 = project_onto_surface(surface2, u2_range, v2_range, &center_pt)?;

        center_points.push(center_pt);
        footprint_1.push(fp1);
        footprint_2.push(fp2);
    }

    let center_curve = fit_nurbs_through_points(&center_points)?;

    let footprint_1 = clip_footprint_to_bounds(surface1, u1_range, v1_range, &footprint_1);
    let footprint_2 = clip_footprint_to_bounds(surface2, u2_range, v2_range, &footprint_2);

    Ok(CenterCurveResult {
        center_curve,
        footprint_1,
        footprint_2,
        contact_solutions: None,
    })
}

// `walker_based_center_curve` lives in `center_curve_walker.rs` to
// keep this file under the 500-line structural cap. Re-exported below.

/// Piecewise-linear interpolation of a `ByPoints` radius profile.
/// `pts` must be sorted by parameter ascending. Out-of-range `t`
/// is clamped to the endpoint radii.
pub fn bypoints_local_radius(pts: &[(f64, f64)], t: f64) -> f64 {
    if pts.is_empty() {
        return 0.0;
    }
    if t <= pts[0].0 {
        return pts[0].1;
    }
    if t >= pts[pts.len() - 1].0 {
        return pts[pts.len() - 1].1;
    }
    for window in pts.windows(2) {
        let (t0, r0) = window[0];
        let (t1, r1) = window[1];
        if t >= t0 && t <= t1 {
            let span = (t1 - t0).max(1e-12);
            let frac = (t - t0) / span;
            return r0 + (r1 - r0) * frac;
        }
    }
    pts[pts.len() - 1].1
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Select the intersection curve whose midpoint is nearest to the given point.
fn select_nearest_curve(
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

/// Find the center point at a given radius for one spine sample.
///
/// **Routing:** if either surface is a [`SurfaceGeometry::Cone`] the
/// existing offset-intersect short-cut errors out (no analytic offset
/// for cones), so we route through the rolling-ball contact solver in
/// `crate::contact`. For every other surface pair we keep the
/// offset-intersect path as the primary because it produces analytic
/// surfaces (torus / sphere / cylinder) that downstream code can
/// recognise. On *unrelated* contact-solver failures (parallel
/// normals, etc.) we still fall through to offset-intersect so the
/// caller gets the historical error message.
///
/// This is phase 1a of the contact-solver integration. Phase 1b will
/// route more aggressively via a spine walker; phase 1d will let the
/// caller force the contact solver path even on analytic pairs.
#[allow(clippy::too_many_arguments)]
fn find_center_point_at_radius(
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
) -> Result<Point3, String> {
    // Cones don't have an analytic offset surface, so any cone in the
    // pair forces the contact solver path. For other surface pairs the
    // offset-intersect short-cut is exact and fast — keep it as the
    // primary route.
    let either_is_cone = matches!(surface1, SurfaceGeometry::Cone { .. })
        || matches!(surface2, SurfaceGeometry::Cone { .. });
    if either_is_cone
        && let Some(center) = try_contact_solver(
            surface1,
            u1_range,
            v1_range,
            surface2,
            u2_range,
            v2_range,
            radius,
            target,
            tolerance,
            orientation_1,
            orientation_2,
            convexity,
        )
    {
        return Ok(center);
    }

    // Offset toward the BALL side (convexity-aware, same truth as
    // the contact solver): inscribed for convex edges, void-side
    // for concave. The old `!orientation` flip was a legacy
    // convention that only encoded the concave case on
    // inward-natural primitive shells.
    let (side_1, side_2) =
        crate::contact::ball_sides_for_edge(convexity, orientation_1, orientation_2);
    if std::env::var_os("PK_CONVEXITY_TRACE").is_some() {
        eprintln!(
            "    [offset] {convexity:?} o=({orientation_1},{orientation_2}) sides=({side_1:?},{side_2:?})"
        );
    }
    let offset1 = offset_surface(
        surface1,
        radius,
        matches!(side_1, crate::contact::BallSide::Inside),
        u1_range,
        v1_range,
        tolerance,
    )?;
    let offset2 = offset_surface(
        surface2,
        radius,
        matches!(side_2, crate::contact::BallSide::Inside),
        u2_range,
        v2_range,
        tolerance,
    )?;

    let ((ou1_min, ou1_max), (ov1_min, ov1_max)) = offset1.parameter_range();
    let off1_u = clamp_finite_range((ou1_min, ou1_max), u1_range);
    let off1_v = clamp_finite_range((ov1_min, ov1_max), v1_range);
    let ((ou2_min, ou2_max), (ov2_min, ov2_max)) = offset2.parameter_range();
    let off2_u = clamp_finite_range((ou2_min, ou2_max), u2_range);
    let off2_v = clamp_finite_range((ov2_min, ov2_max), v2_range);

    let result = intersect_surfaces(
        &offset1, off1_u, off1_v, &offset2, off2_u, off2_v, tolerance,
    )
    .map_err(|e| format!("offset intersection failed: {e}"))?;

    // Find the intersection point nearest to the target: select the
    // nearest curve, then the closest point ON it. The historic
    // curve-MIDPOINT return ignored the target's station entirely —
    // every per-station center sample of a variable-radius march
    // received the same point (equal radii: "all points coincide";
    // varying radii: a radial spur swept into open bodies — the
    // PK-139 baseline).
    let curve = select_nearest_curve(&result.intersections, target)?;
    Ok(closest_point_on_nurbs(&curve, target))
}

/// Closest point on a NURBS curve to a target: dense sampling with
/// two rounds of local ternary refinement (analytic offset
/// intersections here are lines/arcs/low-degree fits — 64 samples
/// bracket the minimum comfortably).
fn closest_point_on_nurbs(curve: &NurbsCurve, target: &Point3) -> Point3 {
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
    curve.point_at(0.5 * (lo + hi))
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
fn try_contact_solver(
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

// Helper functions (clip, clamp, fit, project) are in center_curve_helpers.rs.

#[cfg(test)]
#[path = "center_curve_tests.rs"]
mod tests;
