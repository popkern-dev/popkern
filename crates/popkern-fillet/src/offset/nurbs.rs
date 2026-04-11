//! NURBS-based surface offset with adaptive refinement.
//!
//! Used as a fallback when analytic offset is not available (e.g. cones,
//! free-form NURBS surfaces). Samples the surface on an adaptive grid,
//! computes exact offset points, and fits a degree-1 NURBS surface through
//! them.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsSurface, Point3, Vec3};
use popkern_core::tolerance::Tolerance;

/// Build a uniform clamped knot vector for `n_pts` control points at the
/// given `degree`, spanning `[start, end]`.
fn clamped_knots(n_pts: usize, degree: usize, start: f64, end: f64) -> Vec<f64> {
    let n_knots = n_pts + degree + 1;
    let n_interior = n_pts - degree; // >= 1
    let mut knots = Vec::with_capacity(n_knots);
    for _ in 0..=degree {
        knots.push(start);
    }
    for i in 1..n_interior {
        let t = start + (end - start) * (i as f64) / (n_interior as f64);
        knots.push(t);
    }
    for _ in 0..=degree {
        knots.push(end);
    }
    knots
}

/// Compute the offset point at a single (u, v) parameter.
///
/// Returns `point + sign * distance * normalized_normal`. Falls back to the
/// surface point itself when the normal degenerates to zero length.
fn offset_point(surface: &SurfaceGeometry, u: f64, v: f64, distance: f64, sign: f64) -> Point3 {
    let pt = surface.point_at(u, v);
    let raw_normal = surface.normal_at(u, v);
    let n = raw_normal.normalized().unwrap_or(Vec3::new(0.0, 0.0, 0.0));
    pt + n * (sign * distance)
}

/// Compute an offset of an arbitrary surface using adaptive NURBS sampling.
///
/// The algorithm samples the surface on a uniform grid in UV space, offsets
/// each sample by `distance` along the surface normal, then builds a
/// degree-1 NURBS surface through the offset points. The grid is refined
/// (doubled) until the maximum deviation at inter-sample midpoints is within
/// `tolerance.linear() / 4`, up to a 64x64 ceiling.
pub fn offset_nurbs(
    surface: &SurfaceGeometry,
    distance: f64,
    inward: bool,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<SurfaceGeometry, String> {
    let sign = if inward { -1.0 } else { 1.0 };
    let stage_tol = tolerance.linear() / 4.0;
    let degree: usize = 1;

    let mut n_u: usize = 8;
    let mut n_v: usize = 8;

    loop {
        // Sample offset points on a uniform grid.
        let u_params: Vec<f64> = (0..n_u)
            .map(|i| u_range.0 + (u_range.1 - u_range.0) * (i as f64) / ((n_u - 1) as f64))
            .collect();
        let v_params: Vec<f64> = (0..n_v)
            .map(|j| v_range.0 + (v_range.1 - v_range.0) * (j as f64) / ((n_v - 1) as f64))
            .collect();

        let control_points: Vec<Point3> = u_params
            .iter()
            .flat_map(|&u| {
                v_params
                    .iter()
                    .map(move |&v| offset_point(surface, u, v, distance, sign))
            })
            .collect();

        let weights = vec![1.0_f64; n_u * n_v];
        let knots_u = clamped_knots(n_u, degree, u_range.0, u_range.1);
        let knots_v = clamped_knots(n_v, degree, v_range.0, v_range.1);

        let nurbs = NurbsSurface {
            degree_u: degree,
            degree_v: degree,
            knots_u,
            knots_v,
            control_points,
            num_u: n_u,
            num_v: n_v,
            weights,
        };

        // Check deviation at midpoints between sample positions.
        let mut max_dev: f64 = 0.0;
        for i in 0..(n_u - 1) {
            for j in 0..(n_v - 1) {
                let u_mid = (u_params[i] + u_params[i + 1]) / 2.0;
                let v_mid = (v_params[j] + v_params[j + 1]) / 2.0;

                let exact = offset_point(surface, u_mid, v_mid, distance, sign);
                let approx = nurbs.point_at(u_mid, v_mid);
                let dev = exact.distance_to(&approx);
                if dev > max_dev {
                    max_dev = dev;
                }
            }
        }

        if max_dev <= stage_tol || (n_u >= 64 && n_v >= 64) {
            return Ok(SurfaceGeometry::Nurbs(nurbs));
        }

        // Refine: double the grid.
        n_u = (n_u * 2).min(64);
        n_v = (n_v * 2).min(64);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::{Point3, Vec3};
    use popkern_core::tolerance::Tolerance;

    #[test]
    fn nurbs_offset_plane_z() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let tol = Tolerance::default();
        let result = offset_nurbs(&plane, 3.0, false, (-1.0, 1.0), (-1.0, 1.0), &tol).unwrap();
        // Evaluate at centre — should be at z ≈ 3.0
        match &result {
            SurfaceGeometry::Nurbs(n) => {
                let p = n.point_at(0.0, 0.0);
                assert!((p.z - 3.0).abs() < 0.01, "expected z≈3.0, got z={:.6}", p.z);
            }
            _ => panic!("expected NURBS"),
        }
    }

    #[test]
    fn nurbs_offset_cone_succeeds() {
        let cone = SurfaceGeometry::Cone {
            apex: Point3::new(0.0, 0.0, 1.0),
            axis: Vec3::new(0.0, 0.0, -1.0),
            half_angle: std::f64::consts::FRAC_PI_4,
        };
        let tol = Tolerance::default();
        let result = offset_nurbs(&cone, 0.1, false, (0.1, 1.0), (0.1, 1.0), &tol);
        assert!(result.is_ok());
    }
}
