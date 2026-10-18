//! Line-cylinder intersection algorithm.

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::intersection::{SurfaceCurveIntersection, coordinate_scale};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;

/// Intersect a parametric line `P + t*D` with an infinite cylinder.
///
/// The cylinder is defined by an origin `O`, unit axis `A`, and radius `R`.
/// Surface parameters are:
/// - `u`: azimuthal angle in radians, computed via `atan2` in the orthonormal
///   frame produced by [`orthonormal_basis`].
/// - `v`: signed height along the axis from `O`.
///
/// Returns only intersections whose `t`, `u`, and `v` parameters all fall
/// within the supplied ranges (inclusive within tolerance).
///
/// # Degenerate case
/// When the line is parallel to the axis (`|D_perp|` is below epsilon) the
/// line either lies on the cylinder surface or runs along the interior/exterior.
/// Both cases produce a continuum of intersections (overlap); for now this
/// function returns an empty `Vec` and defers to the general solver.
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) fn intersect_line_cylinder(
    line_origin: &Point3,
    line_dir: &Vec3,
    t_range: (f64, f64),
    cyl_origin: &Point3,
    cyl_axis: &Vec3,
    cyl_radius: f64,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<SurfaceCurveIntersection> {
    // --- Step 1: normalise axis ---
    let a_unit = match cyl_axis.normalized() {
        Some(u) => u,
        None => return Vec::new(),
    };

    // --- Step 2: project line onto the plane perpendicular to axis ---
    let dp = *line_origin - *cyl_origin; // vector from cyl origin to line origin
    let dp_dot_a = dp.dot(&a_unit);
    let dp_perp = dp - a_unit * dp_dot_a;

    let d_dot_a = line_dir.dot(&a_unit);
    let d_perp = *line_dir - a_unit * d_dot_a;

    // --- Step 3: solve quadratic |dp_perp + t*D_perp|^2 = R^2 ---
    let a_coef = d_perp.dot(&d_perp);
    let b_coef = 2.0 * dp_perp.dot(&d_perp);
    let c_coef = dp_perp.dot(&dp_perp) - cyl_radius * cyl_radius;

    // Degenerate: line is parallel to axis (D_perp ≈ 0).
    // The line either lies wholly on the surface (overlap) or never intersects.
    // Defer to the general solver in both cases.
    if a_coef.abs() < tolerance.linear() * tolerance.linear() {
        return Vec::new();
    }

    // Tangency by DISTANCE, as the plane x cylinder solver decides it:
    // the line's closest approach to the axis, h, against the radius
    // within the linear tolerance. The former discriminant band
    // (tol^2 * a) admitted only h within ~tol^2/(8r) of r — a web edge
    // 1e-7 inside a radius-33 boss read as two crossings 5e-3 apart
    // while the web face's plane read the same contact as tangent, and
    // the boolean stranded the record ends between them (PK-789).
    let t_foot = -b_coef / (2.0 * a_coef);
    let h = (dp_perp + d_perp * t_foot).length();
    // A line OUTSIDE the radius touches nothing a record passes
    // through (the plane x cylinder solver emits no curve there either),
    // so it misses; one strictly less than a tolerance INSIDE is a
    // tangent contact. The band is open at a whole tolerance: a contact
    // exactly that deep is the secant the stadium fixtures serve on (a
    // closed band read them tangent and left 18 curveless edges).
    let t_candidates: Vec<f64> = if h > cyl_radius {
        // No real intersection — line misses the cylinder.
        return Vec::new();
    } else if h > cyl_radius
        - tolerance.tangency_band(coordinate_scale(&[*line_origin, *cyl_origin]))
    {
        // Tangent within tolerance — one contact at the foot point.
        vec![t_foot]
    } else {
        // Secant — two distinct solutions.
        let sqrt_disc = (b_coef * b_coef - 4.0 * a_coef * c_coef).max(0.0).sqrt();
        vec![
            (-b_coef - sqrt_disc) / (2.0 * a_coef),
            (-b_coef + sqrt_disc) / (2.0 * a_coef),
        ]
    };

    // Build the orthonormal basis for u-angle computation.
    let (perp1, perp2) = orthonormal_basis(a_unit);

    let mut results = Vec::new();

    for t in t_candidates {
        // --- filter t ---
        if t < t_range.0 - tolerance.parametric_base()
            || t > t_range.1 + tolerance.parametric_base()
        {
            continue;
        }

        // --- Step 4: compute intersection point and surface parameters ---
        let q = *line_origin + *line_dir * t;
        let q_rel = q - *cyl_origin; // Q - O

        // Height along axis: v = (Q - O) · A_unit
        let v_param = q_rel.dot(&a_unit);

        // Azimuthal angle: u = atan2((Q-O)·perp2, (Q-O)·perp1)
        let local_perp1 = q_rel.dot(&perp1);
        let local_perp2 = q_rel.dot(&perp2);
        let u_param = local_perp2.patan2(local_perp1);

        // Normalise u into the u_range window (same pattern as line_circle).
        let u_norm = normalise_angle(u_param, u_range);

        // --- filter u ---
        if u_norm < u_range.0 - tolerance.angular() || u_norm > u_range.1 + tolerance.angular() {
            continue;
        }

        // --- filter v ---
        if v_param < v_range.0 - tolerance.linear() || v_param > v_range.1 + tolerance.linear() {
            continue;
        }

        results.push(SurfaceCurveIntersection::Point {
            point: q,
            t,
            u: u_norm,
            v: v_param,
        });
    }

    results
}

/// Shift `angle` (from `atan2`) into the window `[range.0, range.0 + TAU)` so
/// that a full-circle range `[0, TAU)` accepts every angle.
fn normalise_angle(angle: f64, range: (f64, f64)) -> f64 {
    let tau = std::f64::consts::TAU;
    let mut a = angle;
    while a < range.0 - f64::EPSILON {
        a += tau;
    }
    while a >= range.0 + tau - f64::EPSILON {
        a -= tau;
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    fn default_tol() -> Tolerance {
        Tolerance::default()
    }

    /// Cylinder centred at origin, axis along Z, radius 1, unbounded in v.
    fn unit_cylinder() -> (Point3, Vec3, f64) {
        (Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0), 1.0)
    }

    // -----------------------------------------------------------------------
    // Secant: line passes through two points on the cylinder
    // -----------------------------------------------------------------------
    #[test]
    fn secant_two_intersections() {
        let (o, ax, r) = unit_cylinder();
        // Line along X at height z=0: P=(-2,0,0) D=(1,0,0)
        let hits = intersect_line_cylinder(
            &Point3::new(-2.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (-10.0, 10.0),
            &o,
            &ax,
            r,
            (0.0, TAU),
            (-10.0, 10.0),
            &default_tol(),
        );
        assert_eq!(hits.len(), 2, "expected 2 secant hits, got {}", hits.len());

        // Both hit points should lie on the cylinder surface.
        for hit in &hits {
            if let SurfaceCurveIntersection::Point { point, .. } = hit {
                let radial = (point.x.powi(2) + point.y.powi(2)).sqrt();
                assert!(
                    (radial - 1.0).abs() < 1e-9,
                    "point not on cylinder: radial={radial}"
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // Tangent: line just grazes the cylinder (discriminant ≈ 0)
    // -----------------------------------------------------------------------
    /// Tangency is a distance band (PK-789): a line up to a tolerance
    /// inside the radius touches once, at its foot point; one outside
    /// misses. A radius-33 cylinder and a line 1e-7 inside used to give
    /// two crossings 5e-3 apart.
    #[test]
    fn a_line_within_tolerance_inside_touches_once_and_outside_misses() {
        let (o, ax) = (Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0));
        let at = |y: f64| {
            intersect_line_cylinder(
                &Point3::new(-50.0, y, 0.0),
                &Vec3::new(1.0, 0.0, 0.0),
                (0.0, 100.0),
                &o,
                &ax,
                33.0,
                (0.0, TAU),
                (-10.0, 10.0),
                &default_tol(),
            )
        };
        let inside = at(33.0 - 1e-7);
        assert_eq!(inside.len(), 1, "{inside:?}");
        let SurfaceCurveIntersection::Point { point, .. } = &inside[0] else {
            panic!("point contact expected: {inside:?}");
        };
        assert!(point.x.abs() < 1e-9, "the foot point: {point:?}");
        assert!(
            at(33.0 + 1e-7).is_empty(),
            "a line outside the radius misses"
        );
        assert_eq!(
            at(33.0 - 1e-3).len(),
            2,
            "a real secant still crosses twice"
        );
    }

    #[test]
    fn tangent_one_intersection() {
        let (o, ax, r) = unit_cylinder();
        // Line tangent to unit cylinder at (0,1,0): P=(-2,1,0) D=(1,0,0)
        let hits = intersect_line_cylinder(
            &Point3::new(-2.0, 1.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (-10.0, 10.0),
            &o,
            &ax,
            r,
            (0.0, TAU),
            (-10.0, 10.0),
            &default_tol(),
        );
        assert_eq!(hits.len(), 1, "expected 1 tangent hit, got {}", hits.len());
    }

    // -----------------------------------------------------------------------
    // Miss: line does not intersect
    // -----------------------------------------------------------------------
    #[test]
    fn miss_no_intersection() {
        let (o, ax, r) = unit_cylinder();
        // Line at y=2 (outside radius 1): P=(-2,2,0) D=(1,0,0)
        let hits = intersect_line_cylinder(
            &Point3::new(-2.0, 2.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (-10.0, 10.0),
            &o,
            &ax,
            r,
            (0.0, TAU),
            (-10.0, 10.0),
            &default_tol(),
        );
        assert!(hits.is_empty(), "expected no hits, got {}", hits.len());
    }

    // -----------------------------------------------------------------------
    // Line along axis: degenerate — returns empty
    // -----------------------------------------------------------------------
    #[test]
    fn line_along_axis_returns_empty() {
        let (o, ax, r) = unit_cylinder();
        // Line along the cylinder axis itself: P=(0,0,-5) D=(0,0,1)
        let hits = intersect_line_cylinder(
            &Point3::new(0.0, 0.0, -5.0),
            &Vec3::new(0.0, 0.0, 1.0),
            (-10.0, 10.0),
            &o,
            &ax,
            r,
            (0.0, TAU),
            (-10.0, 10.0),
            &default_tol(),
        );
        assert!(
            hits.is_empty(),
            "parallel-to-axis line should return empty (degenerate)"
        );
    }

    // -----------------------------------------------------------------------
    // Line parallel to axis at radius R: tangent overlap — returns empty
    // -----------------------------------------------------------------------
    #[test]
    fn line_parallel_at_radius_returns_empty() {
        let (o, ax, r) = unit_cylinder();
        // Line at (1,0,z), parallel to Z: lies exactly on the cylinder surface.
        let hits = intersect_line_cylinder(
            &Point3::new(1.0, 0.0, -5.0),
            &Vec3::new(0.0, 0.0, 1.0),
            (-10.0, 10.0),
            &o,
            &ax,
            r,
            (0.0, TAU),
            (-10.0, 10.0),
            &default_tol(),
        );
        assert!(
            hits.is_empty(),
            "line on cylinder surface (overlap) should return empty"
        );
    }

    // -----------------------------------------------------------------------
    // Line parallel to axis outside radius: returns empty
    // -----------------------------------------------------------------------
    #[test]
    fn line_parallel_outside_radius_returns_empty() {
        let (o, ax, r) = unit_cylinder();
        // Line at (2,0,z) — parallel to Z but outside the cylinder.
        let hits = intersect_line_cylinder(
            &Point3::new(2.0, 0.0, -5.0),
            &Vec3::new(0.0, 0.0, 1.0),
            (-10.0, 10.0),
            &o,
            &ax,
            r,
            (0.0, TAU),
            (-10.0, 10.0),
            &default_tol(),
        );
        assert!(
            hits.is_empty(),
            "line parallel outside cylinder should return empty"
        );
    }

    // -----------------------------------------------------------------------
    // t_range filtering: hits outside t_range are discarded
    // -----------------------------------------------------------------------
    #[test]
    fn t_range_filters_hits() {
        let (o, ax, r) = unit_cylinder();
        // Line P=(-2,0,0) D=(1,0,0) hits the unit cylinder where x=±1:
        //   x = -2 + t = -1  →  t = 1  (left side)
        //   x = -2 + t =  1  →  t = 3  (right side)
        // Restrict t_range to (2, 5) so only t=3 survives.
        let hits = intersect_line_cylinder(
            &Point3::new(-2.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (2.0, 5.0),
            &o,
            &ax,
            r,
            (0.0, TAU),
            (-10.0, 10.0),
            &default_tol(),
        );
        assert_eq!(hits.len(), 1, "only one hit should survive t_range filter");
        if let SurfaceCurveIntersection::Point { t, .. } = &hits[0] {
            assert!(*t > 2.0, "surviving t should be > 2, got {t}");
        }
    }

    // -----------------------------------------------------------------------
    // v_range filtering: hits outside height range are discarded
    // -----------------------------------------------------------------------
    #[test]
    fn v_range_filters_hits() {
        let (o, ax, r) = unit_cylinder();
        // Oblique line: hits cylinder at two heights; restrict v to (0, 10).
        // P=(-2,0,-3), D=(1,0,2) — at t=1 height=2*1-3=-1 (below 0), at t=3 height=3.
        let hits = intersect_line_cylinder(
            &Point3::new(-2.0, 0.0, -3.0),
            &Vec3::new(1.0, 0.0, 2.0),
            (-10.0, 10.0),
            &o,
            &ax,
            r,
            (0.0, TAU),
            (0.0, 10.0), // v must be >= 0
            &default_tol(),
        );
        // Only the hit at positive height survives.
        for hit in &hits {
            if let SurfaceCurveIntersection::Point { v, .. } = hit {
                assert!(*v >= -1e-9, "v={v} is below v_range.0=0");
            }
        }
    }

    // -----------------------------------------------------------------------
    // Surface parameters are geometrically consistent
    // -----------------------------------------------------------------------
    #[test]
    fn surface_params_consistent() {
        let (o, ax, r) = unit_cylinder();
        let hits = intersect_line_cylinder(
            &Point3::new(-2.0, 0.0, 3.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (-10.0, 10.0),
            &o,
            &ax,
            r,
            (0.0, TAU),
            (-10.0, 10.0),
            &default_tol(),
        );
        assert_eq!(hits.len(), 2);

        let (perp1, perp2) = orthonormal_basis(ax);
        for hit in &hits {
            if let SurfaceCurveIntersection::Point { point, u, v, .. } = hit {
                // Reconstruct from (u, v) using the cylinder parameterisation.
                let reconstructed = o + perp1 * (r * u.cos()) + perp2 * (r * u.sin()) + ax * *v;
                assert!(
                    point.distance_to(&reconstructed) < 1e-9,
                    "parameter reconstruction mismatch: {point:?} vs {reconstructed:?}"
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // Non-unit axis input (should be normalised internally)
    // -----------------------------------------------------------------------
    #[test]
    fn non_unit_axis_normalised() {
        // Same geometry as secant test but axis supplied un-normalised.
        let hits = intersect_line_cylinder(
            &Point3::new(-2.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (-10.0, 10.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 5.0), // length 5, not 1
            1.0,
            (0.0, TAU),
            (-10.0, 10.0),
            &default_tol(),
        );
        assert_eq!(hits.len(), 2);
    }

    #[test]
    fn zero_axis_returns_empty() {
        // Zero-length axis should return empty (line 39)
        let hits = intersect_line_cylinder(
            &Point3::new(-2.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (-10.0, 10.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 0.0),
            1.0,
            (0.0, TAU),
            (-10.0, 10.0),
            &default_tol(),
        );
        assert!(hits.is_empty());
    }

    #[test]
    fn u_range_filters_cylinder() {
        // Restrict u_range to exclude some azimuthal angles (line 110)
        let (o, ax, r) = unit_cylinder();
        // Line at y=0: hits at x=-1 (u=PI) and x=1 (u=0)
        let hits = intersect_line_cylinder(
            &Point3::new(-2.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (-10.0, 10.0),
            &o,
            &ax,
            r,
            (0.5, 2.5), // excludes u=0 and u=PI
            (-10.0, 10.0),
            &default_tol(),
        );
        assert!(hits.is_empty());
    }
}
