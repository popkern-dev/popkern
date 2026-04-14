//! Fillet surface construction.
//!
//! Given a center curve and footprint points, builds the rolling-ball fillet
//! surface. For two-plane cases an analytic torus/cylinder section is returned;
//! otherwise a NURBS swept surface is constructed by skinning cross-section arcs.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};

use crate::center_curve::CenterCurveResult;

/// Result of a fillet surface construction.
#[derive(Debug, Clone)]
pub struct FilletSurfaceResult {
    /// The constructed fillet surface.
    pub surface: SurfaceGeometry,
    /// Trim curve on the first surface (3D curve along footprint_1).
    pub trim_curve_1: NurbsCurve,
    /// Trim curve on the second surface (3D curve along footprint_2).
    pub trim_curve_2: NurbsCurve,
}

/// Build the fillet surface from center-curve data and the two adjacent surfaces.
///
/// Attempts an analytic construction for the plane-plane case, falling back
/// to a NURBS swept surface for general surface pairs. The `is_constant_radius`
/// flag controls whether the analytic torus path is eligible: a torus has a
/// single minor radius, so it can only represent a constant-radius fillet —
/// passing `false` for a variable-radius spec forces the NURBS sweep path,
/// which honours the per-sample radii baked into `center_result.footprint_*`.
///
/// When `g2_cross_section` is `true` AND `center_result.contact_solutions`
/// is `Some` (i.e., the spine walker produced per-sample contact solutions),
/// the sweep is built via [`crate::g2_cross_section::g2_sweep_from_contact_solutions`]
/// — a degree-5 Bezier cross-section sweep that matches the surface
/// normal curvature at the boundary (G2). Otherwise the legacy
/// degree-1 arc sweep is used. The G2 path needs `orientation_*`
/// to translate the popkern face-normal convention into the
/// `BallSide` convention the contact solver uses.
#[allow(clippy::too_many_arguments)]
pub fn build_fillet_surface(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    center_result: &CenterCurveResult,
    radius: f64,
    arc_samples: usize,
    is_constant_radius: bool,
    g2_cross_section: bool,
    orientation_1: bool,
    orientation_2: bool,
) -> Result<FilletSurfaceResult, String> {
    // G2 path: route through the degree-5 Bezier sweep when the
    // caller asked for it AND the center-curve producer (the spine
    // walker) gave us per-sample contact solutions.
    if g2_cross_section
        && let Some(solutions) = center_result.contact_solutions.as_deref()
        && solutions.len() >= 2
    {
        use crate::contact::BallSide;
        use crate::g2_cross_section::g2_sweep_from_contact_solutions;

        // Translate popkern face orientations into BallSide. popkern
        // face normals point INTO the solid; the rolling ball sits
        // on the side OPPOSITE the face normal. So:
        //   orientation = true  → ball on -surface_normal side → BallSide::Inside
        //   orientation = false → ball on +surface_normal side → BallSide::Outside
        let side_0 = if orientation_1 {
            BallSide::Inside
        } else {
            BallSide::Outside
        };
        let side_1 = if orientation_2 {
            BallSide::Inside
        } else {
            BallSide::Outside
        };
        let nurbs = g2_sweep_from_contact_solutions(solutions, surface1, surface2, side_0, side_1)
            .map_err(|e| format!("G2 sweep failed: {e}"))?;
        let trim_1 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_1)?;
        let trim_2 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_2)?;
        return Ok(FilletSurfaceResult {
            surface: SurfaceGeometry::Nurbs(nurbs),
            trim_curve_1: trim_1,
            trim_curve_2: trim_2,
        });
    }

    // Try analytic torus fast paths for constant-radius fillets.
    // Variable-radius specs cannot use an analytic torus — its
    // single minor radius would silently average the per-sample
    // radii — so we must build a true NURBS sweep that follows
    // the varying footprints.
    if is_constant_radius {
        // §11.1 plane-plane → torus (dihedral fillet).
        if let (
            SurfaceGeometry::Plane { normal: n1, .. },
            SurfaceGeometry::Plane { normal: n2, .. },
        ) = (surface1, surface2)
            && let Some(result) = try_plane_plane_torus(n1, n2, center_result, radius)
        {
            return Ok(result);
        }
        // §11.1 plane-sphere → torus. Works on either (Plane,
        // Sphere) or (Sphere, Plane) ordering.
        if let Some(result) = try_plane_sphere_torus(surface1, surface2, center_result, radius) {
            return Ok(result);
        }
    }

    // General case: NURBS swept surface (G1 arc).
    build_nurbs_sweep(center_result, radius, arc_samples)
}

/// Try to build an analytic torus section for a plane-plane fillet.
///
/// Returns `None` if the geometry is degenerate (e.g., parallel planes).
fn try_plane_plane_torus(
    n1: &Vec3,
    n2: &Vec3,
    center_result: &CenterCurveResult,
    radius: f64,
) -> Option<FilletSurfaceResult> {
    // Compute the dihedral angle between the two planes.
    let cos_angle = n1.dot(n2).clamp(-1.0, 1.0);
    let angle = cos_angle.acos(); // angle between normals

    // The half-angle of the dihedral measured from the inside.
    let half_angle = angle * 0.5;
    let sin_half = half_angle.sin();

    // If planes are nearly parallel, can't build a torus fillet.
    if sin_half.abs() < 1e-10 {
        return None;
    }

    // Torus axis = cross product of the two normals (edge direction).
    let axis = n1.cross(n2);
    let axis = axis.normalized()?;

    // Center of the torus: evaluate the center curve at its midpoint.
    // The center curve is the locus of the rolling-ball center, so its
    // midpoint lies on the torus axis at the correct offset from the edge.
    // Using footprint points instead would inherit projection errors from
    // the offset-surface intersection step.
    use popkern_core::geometry::curve::CurveEval;
    let (t0, t1) = center_result.center_curve.parameter_range();
    let center = center_result.center_curve.point_at((t0 + t1) * 0.5);

    // For a fillet between two planes at angle `angle`:
    // Major radius = r / sin(half_angle) for a general dihedral.
    // For 90-degree planes (angle = PI/2), major_radius = r / sin(PI/4) = r * sqrt(2).
    // However, for a simple rolling-ball fillet the surface is always a
    // cylinder section when angle = PI/2, or a torus section otherwise.
    let major_radius = radius / sin_half;

    let surface = SurfaceGeometry::Torus {
        center,
        axis,
        major_radius,
        minor_radius: radius,
    };

    let trim_1 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_1).ok()?;
    let trim_2 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_2).ok()?;

    Some(FilletSurfaceResult {
        surface,
        trim_curve_1: trim_1,
        trim_curve_2: trim_2,
    })
}

/// Try to build an analytic torus for a plane-sphere fillet (§11.1).
///
/// **Geometry.** The locus of rolling-ball centers tangent to both a
/// plane and a sphere is the intersection of two offset surfaces: an
/// offset plane (parallel to the original at distance `r`) and an
/// offset sphere (concentric with the original at radius `R_s ± r`).
/// Two parallel planes intersected with a sphere give a circle, so
/// the center curve is a 3D circle lying in the offset plane. The
/// rolling ball swept along that circle generates an exact torus:
///
/// - **Axis** = plane normal (the circle's normal)
/// - **Center** = projection of the sphere center onto the offset
///   plane = centroid of the circular center curve
/// - **Major radius** = radius of the circular center curve
/// - **Minor radius** = the fillet radius
///
/// Rather than derive the torus center + major radius from raw
/// sphere/plane math (which needs sign-juggling over which side of
/// the plane the ball sits on and whether it is inside or outside
/// the sphere), we read them directly from the already-computed
/// `center_result.center_curve`: the centroid of N sampled points
/// *is* the torus center, and the mean distance from the centroid
/// *is* the major radius. The approach is analogous to the plane-
/// plane path above that pulls the torus center from the center
/// curve's midpoint.
///
/// Returns `None` for any surface pair that is not `(Plane, Sphere)`
/// or `(Sphere, Plane)`, or when the center curve is too short to
/// fit a meaningful circle (< 3 samples).
fn try_plane_sphere_torus(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    center_result: &CenterCurveResult,
    radius: f64,
) -> Option<FilletSurfaceResult> {
    // Match either ordering — the torus itself is symmetric in the
    // two supports.
    let plane_normal = match (surface1, surface2) {
        (SurfaceGeometry::Plane { normal, .. }, SurfaceGeometry::Sphere { .. })
        | (SurfaceGeometry::Sphere { .. }, SurfaceGeometry::Plane { normal, .. }) => *normal,
        _ => return None,
    };

    // Sample the center curve at N points and compute the centroid
    // + mean radial distance. Use N = 16 — dense enough to suppress
    // per-sample jitter, cheap to evaluate.
    use popkern_core::geometry::curve::CurveEval;
    let n = 16usize;
    let (t0, t1) = center_result.center_curve.parameter_range();
    let mut samples: Vec<Point3> = Vec::with_capacity(n);
    for i in 0..n {
        let t = t0 + (t1 - t0) * (i as f64) / ((n - 1).max(1) as f64);
        samples.push(center_result.center_curve.point_at(t));
    }
    if samples.len() < 3 {
        return None;
    }
    let cx = samples.iter().map(|p| p.x).sum::<f64>() / samples.len() as f64;
    let cy = samples.iter().map(|p| p.y).sum::<f64>() / samples.len() as f64;
    let cz = samples.iter().map(|p| p.z).sum::<f64>() / samples.len() as f64;
    let center = Point3::new(cx, cy, cz);

    let radii: Vec<f64> = samples.iter().map(|p| (*p - center).length()).collect();
    let major_radius = radii.iter().sum::<f64>() / radii.len() as f64;
    if !major_radius.is_finite() || major_radius <= 0.0 {
        return None;
    }

    // The center curve must be a CIRCLE — not just any 3D curve
    // that happens to average to a point. Enforce two conditions:
    //
    //   1. The sampled points lie on a plane with the expected
    //      normal `plane_normal` (no out-of-plane wobble). The
    //      max residual must be < 0.1% of the major radius.
    //
    //   2. The sampled points are equidistant from the centroid
    //      (no radial variation). The radial spread (max − min)
    //      must be < 10% of the major radius. The threshold is
    //      loose on purpose — a genuine circle sampled through a
    //      cubic NURBS fit can pick up a few-percent spread from
    //      the interpolation, while a non-circle like a straight
    //      line has > 100% spread (linear samples have min = 0
    //      at the centroid), so the two cases are well-separated.
    let plane_normal_n = plane_normal.normalized()?;
    let max_planar_residual = samples
        .iter()
        .map(|p| (*p - center).dot(&plane_normal_n).abs())
        .fold(0.0_f64, f64::max);
    if max_planar_residual > major_radius * 1e-3 {
        return None;
    }
    let radial_min = radii.iter().cloned().fold(f64::INFINITY, f64::min);
    let radial_max = radii.iter().cloned().fold(0.0_f64, f64::max);
    if (radial_max - radial_min) > major_radius * 0.1 {
        return None;
    }

    let surface = SurfaceGeometry::Torus {
        center,
        axis: plane_normal_n,
        major_radius,
        minor_radius: radius,
    };

    let trim_1 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_1).ok()?;
    let trim_2 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_2).ok()?;

    Some(FilletSurfaceResult {
        surface,
        trim_curve_1: trim_1,
        trim_curve_2: trim_2,
    })
}

/// Build a NURBS swept surface by skinning cross-section arcs along the center curve.
///
/// Each cross-section is a polyline arc from footprint_1[i] to footprint_2[i],
/// centered at the corresponding center-curve sample. The arcs are sampled
/// at `arc_samples` points and skinned into a tensor-product surface.
fn build_nurbs_sweep(
    center_result: &CenterCurveResult,
    radius: f64,
    arc_samples: usize,
) -> Result<FilletSurfaceResult, String> {
    let n_sections = center_result.footprint_1.len();
    if n_sections < 2 {
        return Err("need at least 2 cross-sections for sweep".into());
    }
    let n_arc = arc_samples.max(3); // at least 3 points per arc

    // Sample center curve at n_sections points.
    let (t_min, t_max) =
        popkern_core::geometry::curve::CurveEval::parameter_range(&center_result.center_curve);

    let mut control_points = Vec::with_capacity(n_sections * n_arc);

    for i in 0..n_sections {
        let t = t_min + (t_max - t_min) * (i as f64) / ((n_sections - 1).max(1) as f64);
        let center =
            popkern_core::geometry::curve::CurveEval::point_at(&center_result.center_curve, t);
        let fp1 = &center_result.footprint_1[i];
        let fp2 = &center_result.footprint_2[i];

        // Build the arc from fp1 to fp2 centered at `center`.
        let arc_pts = sample_arc(&center, fp1, fp2, radius, n_arc);
        control_points.extend_from_slice(&arc_pts);
    }

    // Build the NurbsSurface: U = along edge (n_sections), V = around arc (n_arc).
    let degree_u = 1.min(n_sections - 1);
    let degree_v = 1.min(n_arc - 1);

    let knots_u = build_clamped_uniform_knots(n_sections, degree_u);
    let knots_v = build_clamped_uniform_knots(n_arc, degree_v);
    let weights = vec![1.0; n_sections * n_arc];

    let nurbs = NurbsSurface {
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        control_points,
        num_u: n_sections,
        num_v: n_arc,
        weights,
    };

    let trim_1 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_1)?;
    let trim_2 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_2)?;

    Ok(FilletSurfaceResult {
        surface: SurfaceGeometry::Nurbs(nurbs),
        trim_curve_1: trim_1,
        trim_curve_2: trim_2,
    })
}

/// Sample a circular arc from `start` to `end` centered at `center`.
///
/// The arc lies in the plane defined by the three points. If `start` and
/// `end` are equidistant from `center` at approximately `radius`, the arc
/// sweeps the shorter angle between them.
fn sample_arc(
    center: &Point3,
    start: &Point3,
    end: &Point3,
    _radius: f64,
    n: usize,
) -> Vec<Point3> {
    let v1 = *start - *center;
    let v2 = *end - *center;
    let r1 = v1.length();
    let r2 = v2.length();

    // Handle degenerate cases.
    if r1 < 1e-15 || r2 < 1e-15 || n < 2 {
        return vec![*start; n];
    }

    // Normalize to unit vectors.
    let u = Vec3::new(v1.x / r1, v1.y / r1, v1.z / r1);
    let w_raw = v2 - u * v2.dot(&u);
    let w_len = w_raw.length();

    if w_len < 1e-15 {
        // start and end are nearly coincident or antipodal — just linearly interpolate.
        let mut pts = Vec::with_capacity(n);
        for i in 0..n {
            let t = i as f64 / (n - 1).max(1) as f64;
            pts.push(Point3::new(
                start.x + (end.x - start.x) * t,
                start.y + (end.y - start.y) * t,
                start.z + (end.z - start.z) * t,
            ));
        }
        return pts;
    }

    let v = Vec3::new(w_raw.x / w_len, w_raw.y / w_len, w_raw.z / w_len);

    // Compute sweep angle.
    let cos_theta = (v1.dot(&v2) / (r1 * r2)).clamp(-1.0, 1.0);
    let theta = cos_theta.acos();

    // Use average radius for the arc.
    let r = (r1 + r2) * 0.5;

    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 / (n - 1).max(1) as f64;
        let angle = theta * t;
        let x = r * angle.cos();
        let y = r * angle.sin();
        pts.push(Point3::new(
            center.x + u.x * x + v.x * y,
            center.y + u.y * x + v.y * y,
            center.z + u.z * x + v.z * y,
        ));
    }

    pts
}

/// Build clamped uniform knot vector for `n` control points at given degree.
fn build_clamped_uniform_knots(n: usize, degree: usize) -> Vec<f64> {
    let n_knots = n + degree + 1;
    let n_interior = n.saturating_sub(degree);
    let mut knots = Vec::with_capacity(n_knots);

    knots.extend(std::iter::repeat_n(0.0, degree + 1));
    for i in 1..n_interior {
        knots.push(i as f64 / n_interior as f64);
    }
    knots.extend(std::iter::repeat_n(1.0, degree + 1));

    knots
}

#[cfg(test)]
#[path = "fillet_surface_tests.rs"]
mod tests;
