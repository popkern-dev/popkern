//! Analytic fast-path surface producers for `build_fillet_surface`.
//!
//! Each `try_*` function matches a specific surface pair and returns
//! `Some(FilletSurfaceResult)` when it can produce an **exact**
//! analytic surface (torus, cylinder, …), or `None` otherwise so the
//! caller falls through to the generic NURBS sweep.
//!
//! Split out of `fillet_surface.rs` to keep both files under the
//! 500-line structural cap and to give the §11.1 analytic-path
//! catalogue a dedicated home as more cases land.
//!
//! ## §11.1 catalogue
//!
//! | Pair | Surface | Function |
//! |---|---|---|
//! | plane-plane (straight edge) | cylinder | [`try_plane_plane_analytic`] |
//! | plane-plane (curved edge) | torus | [`try_plane_plane_analytic`] |
//! | plane-sphere | torus | [`try_circular_center_torus`] |
//! | plane-cylinder (⊥ axis) | torus | [`try_circular_center_torus`] |
//! | plane-cone (⊥ axis) | torus | [`try_circular_center_torus`] |
//! | cylinder-cylinder (∥ axes) | cylinder | [`try_parallel_cylinder_cylinder`] |
//!
//! **Plane-cone reasoning.** A plane perpendicular to a cone's axis
//! cuts the cone in a circle; by rotational symmetry the rolling-ball
//! center curve is a 3D circle too, so the fillet surface is an
//! exact torus. Oblique plane-cone intersections (ellipse, parabola,
//! hyperbola) produce canal surfaces that popkern has no primitive
//! for, so the detector gates on `|plane_normal · cone_axis| > 0.999`
//! up front and leaves non-perpendicular cases to the generic NURBS
//! sweep.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};

use crate::center_curve::CenterCurveResult;
use crate::fillet_surface::FilletSurfaceResult;

/// Try to build an analytic surface for a plane-plane fillet.
///
/// Returns a **cylinder** when the center curve is approximately
/// a straight line (straight dihedral edge), or a **torus** when
/// it is curved. Returns `None` for degenerate cases (parallel
/// planes).
pub(crate) fn try_plane_plane_analytic(
    n1: &Vec3,
    n2: &Vec3,
    center_result: &CenterCurveResult,
    radius: f64,
) -> Option<FilletSurfaceResult> {
    let cos_angle = n1.dot(n2).clamp(-1.0, 1.0);
    let angle = cos_angle.acos();
    let half_angle = angle * 0.5;
    let sin_half = half_angle.sin();
    if sin_half.abs() < 1e-10 {
        return None;
    }

    use popkern_core::geometry::curve::CurveEval;
    let (t0, t1) = center_result.center_curve.parameter_range();
    let p_start = center_result.center_curve.point_at(t0);
    let p_end = center_result.center_curve.point_at(t1);
    let spine_mid = center_result.center_curve.point_at((t0 + t1) * 0.5);

    let chord = p_end - p_start;
    let chord_len = chord.length();

    // Straight-line test: midpoint deviation from chord < 1%.
    let is_straight = if chord_len > 1e-12 {
        let chord_dir = chord * (1.0 / chord_len);
        let to_mid = spine_mid - p_start;
        let proj = chord_dir * to_mid.dot(&chord_dir);
        (to_mid - proj).length() < 0.01 * chord_len
    } else {
        false
    };

    let surface = if is_straight {
        let axis = chord.normalized()?;
        SurfaceGeometry::Cylinder {
            origin: p_start,
            axis,
            radius,
        }
    } else {
        let axis = n1.cross(n2).normalized()?;
        let major_radius = radius / sin_half;
        // Derive the torus center from the center curve itself (the
        // spine is already on the correct side for both convex and
        // concave edges) instead of assuming the anti-bisector
        // placement, which holds only for convex arrangements.
        let center = circumcenter(&p_start, &spine_mid, &p_end)?;
        let measured_major = (p_start - center).length();
        if (measured_major - major_radius).abs() > 1.0e-3 * major_radius.max(1.0) {
            // The spine circle disagrees with the dihedral-derived
            // major radius: not an analytic torus configuration.
            return None;
        }
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius: radius,
        }
    };

    let trim_1 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_1).ok()?;
    let trim_2 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_2).ok()?;

    Some(FilletSurfaceResult {
        surface,
        trim_curve_1: trim_1,
        trim_curve_2: trim_2,
        continuity_downgrade: None,
    })
}

/// Circumcenter of three points (center of the circle through them);
/// `None` when they are collinear within numerical noise.
fn circumcenter(
    a: &popkern_core::geometry::Point3,
    b: &popkern_core::geometry::Point3,
    c: &popkern_core::geometry::Point3,
) -> Option<popkern_core::geometry::Point3> {
    let ab = *b - *a;
    let ac = *c - *a;
    let cross = ab.cross(&ac);
    let denom = 2.0 * cross.length_squared();
    if denom < 1.0e-18 {
        return None;
    }
    let offset = (cross.cross(&ab) * ac.length_squared() + ac.cross(&cross) * ab.length_squared())
        * (1.0 / denom);
    Some(*a + offset)
}

/// Try to build an analytic torus for a fillet whose center curve is
/// a 3D circle (§11.1 `plane-sphere`, `plane-cylinder`/`plane-cone`
/// with the plane perpendicular to the axis, and any other surface
/// pair with a circular center curve).
///
/// **Geometry.** Whenever the rolling-ball center locus is a 3D
/// circle, the resulting fillet surface is an exact torus:
///
/// - **Axis** = normal to the plane containing the circle
/// - **Center** = centroid of the circular center curve
/// - **Major radius** = radius of the circular center curve
/// - **Minor radius** = the fillet radius
///
/// For a **plane-sphere** fillet the center curve is always a circle
/// (intersection of a parallel plane with a concentric sphere).
///
/// For a **plane-cylinder** or **plane-cone** fillet the center curve
/// is a circle iff the plane is perpendicular to the curved surface's
/// axis. Oblique planes cut a cylinder in an ellipse and a cone in
/// an ellipse / parabola / hyperbola — none of which produce a
/// circular center curve — so the detector gates on
/// `|plane_normal · surface_axis| > 0.999` up front.
///
/// Rather than derive the torus center + major radius from raw
/// sphere/cylinder/cone/plane math (which needs sign-juggling over
/// which side of the plane the ball sits on), we read them directly
/// from the already-computed `center_result.center_curve`. The
/// approach mirrors [`try_plane_plane_analytic`] above.
pub(crate) fn try_circular_center_torus(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    center_result: &CenterCurveResult,
    radius: f64,
) -> Option<FilletSurfaceResult> {
    let plane_normal = match (surface1, surface2) {
        (SurfaceGeometry::Plane { normal, .. }, SurfaceGeometry::Sphere { .. })
        | (SurfaceGeometry::Sphere { .. }, SurfaceGeometry::Plane { normal, .. }) => *normal,
        (SurfaceGeometry::Plane { normal, .. }, SurfaceGeometry::Cylinder { axis, .. })
        | (SurfaceGeometry::Cylinder { axis, .. }, SurfaceGeometry::Plane { normal, .. }) => {
            let n_unit = normal.normalized()?;
            let a_unit = axis.normalized()?;
            if n_unit.dot(&a_unit).abs() < 0.999 {
                return None;
            }
            *normal
        }
        (SurfaceGeometry::Plane { normal, .. }, SurfaceGeometry::Cone { axis, .. })
        | (SurfaceGeometry::Cone { axis, .. }, SurfaceGeometry::Plane { normal, .. }) => {
            let n_unit = normal.normalized()?;
            let a_unit = axis.normalized()?;
            // Oblique plane-cone intersections are conic sections
            // other than circles (ellipse / parabola / hyperbola).
            // Their fillet surfaces are canal surfaces with non-
            // circular spines — popkern has no canal primitive in
            // `SurfaceGeometry`, so we leave those cases to the
            // generic NURBS sweep. The deferral is documented in
            // `docs/decisions/2026-04-14-oblique-plane-cone-canal-
            // deferred.md` with the trigger criteria for revisiting.
            if n_unit.dot(&a_unit).abs() < 0.999 {
                return None;
            }
            *normal
        }
        _ => return None,
    };

    use popkern_core::geometry::curve::CurveEval;
    let n = 16usize;
    let (t0, t1) = center_result.center_curve.parameter_range();
    // Closed spines are fitted as OPEN curves; the fit wiggles
    // radially near its end conditions (measured 5.3e-2 on an exact
    // circle whose interior sits at 1e-8). Judge circularity on the
    // interior 90 %; the torus parameters below come from the exact
    // spine, not from these samples.
    let closed_spine = center_result
        .center_curve
        .point_at(t0)
        .distance_to(&center_result.center_curve.point_at(t1))
        < (t1 - t0).abs().max(1.0) * 1e-3;
    let (s0, s1) = if closed_spine {
        (t0 + (t1 - t0) * 0.05, t1 - (t1 - t0) * 0.05)
    } else {
        (t0, t1)
    };
    let mut samples: Vec<Point3> = Vec::with_capacity(n);
    for i in 0..n {
        let t = s0 + (s1 - s0) * (i as f64) / ((n - 1).max(1) as f64);
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
    //   1. planarity — max out-of-plane residual < 0.1 % of the
    //      major radius.
    //   2. uniform radius — max − min radial distance < 10 % of the
    //      major radius. The second check is what catches a linear
    //      center curve, which has min = 0 (centroid lies on the
    //      line) and max > 0.
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

    // Torus parameters from the EXACT spine when a planar-support
    // footprint is available: the fitted center curve's wobble puts
    // the sampled center/major-radius ~1e-2 off the true surface.
    // The true ball-center circle is the planar-support footprint
    // lifted by the minor radius along the plane normal (sign toward
    // the other footprint) — exact to the contact-solver residual.
    let (center, major_radius) = {
        // The PLANE-side footprint is the one lying ON the plane
        // surface. (Selecting by surface order picked the sphere's
        // contact ring on plane–sphere fixtures; selecting by
        // planarity alone cannot discriminate — a sphere's contact
        // ring is also planar by symmetry.)
        let plane_origin = match (surface1, surface2) {
            (SurfaceGeometry::Plane { origin, .. }, _) => *origin,
            (_, SurfaceGeometry::Plane { origin, .. }) => *origin,
            _ => return None,
        };
        let dist_to_plane = |fp: &[Point3]| -> f64 {
            fp.iter()
                .map(|p| (*p - plane_origin).dot(&plane_normal_n).abs())
                .fold(0.0_f64, f64::max)
        };
        let planar_first =
            dist_to_plane(&center_result.footprint_1) <= dist_to_plane(&center_result.footprint_2);
        let (plane_fp, other_fp) = if planar_first {
            (&center_result.footprint_1, &center_result.footprint_2)
        } else {
            (&center_result.footprint_2, &center_result.footprint_1)
        };
        let lifted = if plane_fp.len() >= 4 && !other_fp.is_empty() {
            // Lift sign: toward the sampled spine centroid (always
            // available and approximately right — the refit only
            // sharpens center/radius). The other footprint's side is
            // ambiguous when the contacts are near-coplanar.
            let lift = if (center - plane_fp[0]).dot(&plane_normal_n) >= 0.0 {
                radius
            } else {
                -radius
            };
            let _ = &other_fp;
            let spine: Vec<Point3> = plane_fp
                .iter()
                .map(|p| *p + plane_normal_n * lift)
                .collect();
            fit_exact_circle_center(&spine)
        } else {
            None
        };
        lifted.unwrap_or((center, major_radius))
    };
    let surface = SurfaceGeometry::Torus {
        center,
        axis: plane_normal_n,
        major_radius,
        minor_radius: radius,
    };

    // Trim curves stay interpolation fits: their pcurves (fitted
    // downstream in trim.rs) share this parameterization, and an
    // exact rational circle (angle-quarter parameterized) cannot be
    // substituted without rebuilding matched pcurves per support
    // type — recorded PK-099 follow-up.
    let trim_1 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_1).ok()?;
    let trim_2 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_2).ok()?;

    Some(FilletSurfaceResult {
        surface,
        trim_curve_1: trim_1,
        trim_curve_2: trim_2,
        continuity_downgrade: None,
    })
}

#[path = "fillet_surface_analytic_cylinder_pair.rs"]
mod analytic_cylinder_pair;
pub(crate) use analytic_cylinder_pair::{fit_exact_circle_center, try_parallel_cylinder_cylinder};
