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
//! | cylinder-cylinder (∥ axes) | cylinder | [`try_cylindrical_analytic`] |
//! | plane-cylinder (∥ axis, axial edge) | cylinder | [`try_cylindrical_analytic`] |
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
use popkern_math::PortableMath;

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
    let angle = cos_angle.pacos();
    let half_angle = angle * 0.5;
    let sin_half = half_angle.psin();
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
    // circle whose interior sits at 1e-8). Judge on the interior
    // 90 %; the torus parameters below are derived, not sampled.
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
    // Selection samples: the fitted spine (cheap, wobbles up to
    // ~1e-3-relative between interpolation points — good enough to
    // pick the sign branch, not to accept it).
    let selection_samples: Vec<Point3> = (0..n)
        .map(|i| {
            let t = s0 + (s1 - s0) * (i as f64) / ((n - 1).max(1) as f64);
            center_result.center_curve.point_at(t)
        })
        .collect();
    if selection_samples.len() < 3 {
        return None;
    }

    // The torus is fully determined by the support pair + fillet
    // radius up to two branch signs (which side of the plane the
    // ball sits on × which side of the curved surface). Derive all
    // sign candidates EXACTLY and let the sampled spine select the
    // branch it actually lies on. This replaces the former
    // centroid-based circularity judgment, which structurally
    // rejected OPEN arc spines (a quarter arc's radial spread about
    // its sample centroid is ~60 % of its radius even when the arc
    // is exact — the §12.3 blend-on-blend decline class, PK-238),
    // and the lifted-footprint refit, whose 1e-6 relative
    // acceptance rejected fitted-spine wobble.
    let plane_normal_n = plane_normal.normalized()?;
    let plane_origin = match (surface1, surface2) {
        (SurfaceGeometry::Plane { origin, .. }, _) => *origin,
        (_, SurfaceGeometry::Plane { origin, .. }) => *origin,
        _ => return None,
    };
    let curved = if matches!(surface1, SurfaceGeometry::Plane { .. }) {
        surface2
    } else {
        surface1
    };
    // Candidate center circles `(center, major_radius)`; the circle
    // plane offsets the support plane by ±radius (ball tangent to
    // the plane), the major radius offsets the curved support by
    // ±radius along its own normal field.
    let mut candidates: Vec<(Point3, f64, f64)> = Vec::with_capacity(4);
    for s_plane in [1.0_f64, -1.0] {
        let q = plane_origin + plane_normal_n * (s_plane * radius);
        let on_circle_plane =
            |p: Point3| -> Point3 { p - plane_normal_n * (p - q).dot(&plane_normal_n) };
        for s_curved in [1.0_f64, -1.0] {
            match curved {
                SurfaceGeometry::Cylinder {
                    origin,
                    radius: r_cyl,
                    ..
                } => {
                    let major = r_cyl + s_curved * radius;
                    if major > 0.0 {
                        candidates.push((on_circle_plane(*origin), major, s_plane));
                    }
                }
                SurfaceGeometry::Sphere {
                    center: sphere_center,
                    radius: r_sphere,
                } => {
                    let h = (*sphere_center - q).dot(&plane_normal_n);
                    let rr = r_sphere + s_curved * radius;
                    let major_sq = rr * rr - h * h;
                    if major_sq > 0.0 {
                        candidates.push((
                            on_circle_plane(*sphere_center),
                            major_sq.sqrt(),
                            s_plane,
                        ));
                    }
                }
                SurfaceGeometry::Cone {
                    apex,
                    axis,
                    half_angle,
                } => {
                    let a = axis.normalized()?;
                    let denom = a.dot(&plane_normal_n);
                    if denom.abs() < 1e-9 {
                        continue;
                    }
                    let t_axis = (q - *apex).dot(&plane_normal_n) / denom;
                    if t_axis <= 0.0 {
                        // Circle plane on the empty side of the apex.
                        continue;
                    }
                    let rho = t_axis * half_angle.ptan();
                    let major = rho + s_curved * radius / half_angle.pcos();
                    if major > 0.0 {
                        candidates.push((*apex + a * t_axis, major, s_plane));
                    }
                }
                _ => return None,
            }
        }
    }

    // Branch selection on the fitted spine (nearest candidate),
    // then ACCEPTANCE on exact data: rolling-ball centers from the
    // walker's contact solutions when recorded, else the support
    // plane's raw footprint points lifted ±radius along the plane
    // normal (both are ball centers, exact to the solver residual —
    // free of the spine fit's interpolation wobble). Anything that
    // fails the verified candidate — an oblique section arc, a
    // wobbling footprint, a straight spine — declines to the
    // generic NURBS sweep.
    let deviation = |points: &[Point3], center: &Point3, major: f64| -> f64 {
        points
            .iter()
            .map(|p| {
                let w = *p - *center;
                let axial = w.dot(&plane_normal_n);
                let radial = (w - plane_normal_n * axial).length();
                ((radial - major).powi(2) + axial * axial).sqrt()
            })
            .fold(0.0_f64, f64::max)
    };
    let (center, major_radius, s_plane) = candidates
        .into_iter()
        .map(|(c, major, s_plane)| {
            let dev = deviation(&selection_samples, &c, major);
            (c, major, s_plane, dev)
        })
        .min_by(|a, b| a.3.total_cmp(&b.3))
        .map(|(c, major, s_plane, _)| (c, major, s_plane))?;
    let verify_samples: Vec<Point3> = match center_result.contact_solutions.as_deref() {
        Some(solutions) if solutions.len() >= 3 => solutions.iter().map(|s| s.center).collect(),
        _ => {
            // Plane-side footprint: the one lying ON the plane.
            // (Planarity alone cannot discriminate — a sphere's
            // contact ring is planar by symmetry too.)
            let dist_to_plane = |fp: &[Point3]| -> f64 {
                fp.iter()
                    .map(|p| (*p - plane_origin).dot(&plane_normal_n).abs())
                    .fold(0.0_f64, f64::max)
            };
            let plane_fp = if dist_to_plane(&center_result.footprint_1)
                <= dist_to_plane(&center_result.footprint_2)
            {
                &center_result.footprint_1
            } else {
                &center_result.footprint_2
            };
            if plane_fp.len() < 3 {
                return None;
            }
            plane_fp
                .iter()
                .map(|p| *p + plane_normal_n * (s_plane * radius))
                .collect()
        }
    };
    let spine_deviation = deviation(&verify_samples, &center, major_radius);
    let spine_gate = (major_radius + radius) * 1.0e-4;
    if popkern_core::trace::trace_enabled("analytic-intersection") {
        eprintln!(
            "    [torus-arm] best candidate center ({:.4},{:.4},{:.4}) major {major_radius:.4} dev {spine_deviation:.3e} gate {spine_gate:.3e}",
            center.x, center.y, center.z
        );
    }
    if spine_deviation > spine_gate {
        return None;
    }

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
pub(crate) use analytic_cylinder_pair::{fit_exact_circle_center, try_cylindrical_analytic};
