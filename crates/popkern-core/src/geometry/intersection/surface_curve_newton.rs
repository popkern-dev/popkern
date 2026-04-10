//! Newton-iteration refinement for surface-curve intersections.
//!
//! Solves the 3×3 system `F(t,u,v) = curve(t) - surface(u,v) = 0` via Cramer's rule.
//! The Jacobian columns are:
//!   - ∂F/∂t =  curve'(t)
//!   - ∂F/∂u = -surface_u(u,v)
//!   - ∂F/∂v = -surface_v(u,v)

use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

/// Refine an approximate surface-curve intersection using Newton-Raphson.
///
/// Solves `F(t,u,v) = curve(t) - surface(u,v) = 0` via Cramer's rule on the
/// 3×3 linear system formed from the Jacobian of `F` at each iteration.
///
/// Returns `Some((t, u, v, point))` if the iteration converges within the
/// linear tolerance, or `None` if the determinant becomes singular or the
/// maximum iteration count is reached.
pub(crate) fn newton_surface_curve(
    curve: &CurveGeometry,
    t_approx: f64,
    surface: &SurfaceGeometry,
    u_approx: f64,
    v_approx: f64,
    tolerance: &Tolerance,
) -> Option<(f64, f64, f64, Point3)> {
    let max_iter = 20;
    let mut t = t_approx;
    let mut u = u_approx;
    let mut v = v_approx;

    for _ in 0..max_iter {
        let pc = curve.point_at(t);
        let ps = surface.point_at(u, v);

        // Residual vector: F = curve(t) - surface(u,v)
        let f = pc - ps;

        if f.length() < tolerance.linear() {
            return Some((t, u, v, pc));
        }

        // Jacobian columns:
        //   dc  =  ∂F/∂t =  curve'(t)
        //   su  =  ∂F/∂u = -surface_u(u,v)   (negate because F = curve - surface)
        //   sv  =  ∂F/∂v = -surface_v(u,v)
        let dc = curve.derivative_at(t, 1);
        let su = surface.partial_u(u, v) * -1.0;
        let sv = surface.partial_v(u, v) * -1.0;

        // det = dc · (su × sv)  — the scalar triple product of the Jacobian columns
        let su_cross_sv = su.cross(&sv);
        let det = dc.dot(&su_cross_sv);

        if det.abs() < 1e-30 {
            return None;
        }

        let inv = 1.0 / det;

        // Cramer's rule: replace each column with the residual f to get each update.
        //   Δt = det([f, su, sv]) / det  =  f · (su × sv) / det
        //   Δu = det([dc, f, sv]) / det  =  dc · (f × sv) / det
        //   Δv = det([dc, su, f]) / det  =  dc · (su × f) / det
        t -= f.dot(&su_cross_sv) * inv;
        u -= dc.dot(&f.cross(&sv)) * inv;
        v -= dc.dot(&su.cross(&f)) * inv;
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::CurveGeometry;
    use crate::geometry::surface::SurfaceGeometry;
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;

    // -----------------------------------------------------------------------
    // Test 1: Line through a plane — should converge to the exact intersection.
    //
    // Plane z=1 (origin=(0,0,1), normal=(0,0,1), u_axis=(1,0,0), v_axis=(0,1,0)).
    // Line along z-axis: origin=(0,0,0), direction=(0,0,1).
    // Exact intersection: t=1, u=0, v=0, point=(0,0,1).
    // -----------------------------------------------------------------------
    #[test]
    fn line_plane_converges_to_exact() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 1.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let line = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(0.0, 0.0, 1.0),
        };
        let tol = Tolerance::default();

        let result = newton_surface_curve(&line, 0.8, &plane, 0.1, 0.1, &tol);
        assert!(result.is_some(), "Expected convergence for line-plane");
        let (t, u, v, point) = result.unwrap();
        assert!((t - 1.0).abs() < 1e-10, "t should be 1.0, got {t}");
        assert!(u.abs() < 1e-10, "u should be 0.0, got {u}");
        assert!(v.abs() < 1e-10, "v should be 0.0, got {v}");
        assert!(
            (point.z - 1.0).abs() < 1e-10,
            "z should be 1.0, got {}",
            point.z
        );
    }

    // -----------------------------------------------------------------------
    // Test 2: Line parallel to a plane — determinant is zero, must return None.
    //
    // Plane z=0 (normal=(0,0,1)).
    // Line along x-axis at z=5: origin=(0,0,5), direction=(1,0,0).
    // The line is parallel to the plane and never intersects it; det=0.
    // -----------------------------------------------------------------------
    #[test]
    fn line_parallel_to_plane_returns_none() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        // Line at z=5, direction along x — parallel to the plane.
        let line = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 5.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        };
        let tol = Tolerance::default();

        let result = newton_surface_curve(&line, 0.0, &plane, 0.0, 0.0, &tol);
        assert!(
            result.is_none(),
            "Expected None for parallel line and plane"
        );
    }

    // -----------------------------------------------------------------------
    // Test 3: Line intersecting a sphere.
    //
    // Unit sphere at origin.
    // Line: origin=(2,0,0), direction=(-1,0,0) — passes through the equator.
    // Sphere parametrization: u = longitude, v = colatitude.
    //   point_at(u,v) = (sin(v)cos(u), sin(v)sin(u), cos(v))
    // First intersection at (1,0,0): t=1, u=0, v=π/2.
    // Start from an approximation close to that intersection.
    // -----------------------------------------------------------------------
    #[test]
    fn line_sphere_converges() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let line = CurveGeometry::Line {
            origin: Point3::new(2.0, 0.0, 0.0),
            direction: Vec3::new(-1.0, 0.0, 0.0),
        };
        let tol = Tolerance::default();

        // Approximate starting point near (1,0,0) → t≈1, u≈0, v≈π/2
        let v_approx = std::f64::consts::FRAC_PI_2;
        let result = newton_surface_curve(&line, 1.1, &sphere, 0.05, v_approx + 0.05, &tol);
        assert!(result.is_some(), "Expected convergence for line-sphere");
        let (t, _u, _v, point) = result.unwrap();

        // The intersection point should lie on the sphere surface.
        let dist_from_origin = (point.x * point.x + point.y * point.y + point.z * point.z).sqrt();
        // Convergence is within the default linear tolerance (1e-6); the actual
        // residual is well below that, so a check at 1e-5 is a safe bound.
        assert!(
            (dist_from_origin - 1.0).abs() < 1e-5,
            "Point not on sphere: dist={dist_from_origin}"
        );
        // And the curve parameter should place the point on the line.
        let line_pt = CurveGeometry::Line {
            origin: Point3::new(2.0, 0.0, 0.0),
            direction: Vec3::new(-1.0, 0.0, 0.0),
        };
        use crate::geometry::curve::CurveEval;
        let curve_pt = line_pt.point_at(t);
        assert!(
            (curve_pt.x - point.x).abs() < 1e-5,
            "Curve point x mismatch"
        );
    }
}
