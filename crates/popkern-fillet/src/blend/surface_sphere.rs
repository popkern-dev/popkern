//! Rolling-ball sphere corner blend: circumcenter fit over the
//! boundary samples with a radius-scaled acceptance window.
//!
//! Split from `blend/surface.rs` to honour the structural line cap.

use popkern_core::geometry::Point3;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;

use super::boundary::BlendBoundary;
use super::surface::BlendSurfaceResult;
use super::surface_math::solve_circumcenter_3x3;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ContinuityLevel, CornerMethod};

/// Try to build a spherical cap blend.
///
/// Computes the circumcenter of the boundary points — the point
/// equidistant from all of them. If the points lie on a sphere (max
/// deviation < tolerance), returns that sphere. This handles cube
/// corners with equal-radius perpendicular fillets exactly.
///
/// Returns `None` if the geometry is degenerate or the points don't
/// lie on a sphere.
pub(super) fn try_spherical_blend(
    _boundary: &BlendBoundary,
    _fillet_surfaces: &[FilletSurfaceResult],
    points: &[Point3],
    boundary_curves: &[NurbsCurve],
    tolerance: &Tolerance,
) -> Option<BlendSurfaceResult> {
    if points.len() < 4 {
        return None;
    }

    // Compute the circumcenter by solving the linear system:
    // 2*(p_j - p_0) · c = |p_j|² - |p_0|²  for j = 1..N
    //
    // This finds the point equidistant from all boundary points.
    // Use the first 3 non-degenerate point differences to solve for c,
    // then verify against all points.
    let p0 = points[0];
    let p0_sq = p0.x * p0.x + p0.y * p0.y + p0.z * p0.z;

    // Build a 3x4 augmented matrix [A | b] from the first 3 valid equations.
    let mut rows: Vec<[f64; 4]> = Vec::new();
    for p in points.iter().skip(1) {
        let dx = p.x - p0.x;
        let dy = p.y - p0.y;
        let dz = p.z - p0.z;
        if dx.abs() + dy.abs() + dz.abs() < tolerance.solver() {
            continue; // Skip near-duplicate points.
        }
        let rhs = p.x * p.x + p.y * p.y + p.z * p.z - p0_sq;
        rows.push([2.0 * dx, 2.0 * dy, 2.0 * dz, rhs]);
        if rows.len() == 3 {
            break;
        }
    }

    if rows.len() < 3 {
        return None; // Not enough independent equations.
    }

    // Solve 3x3 system by Gaussian elimination with partial pivoting.
    let (cx, cy, cz) = solve_circumcenter_3x3(&rows)?;
    let center = Point3::new(cx, cy, cz);

    // Compute radius as average distance from center to all points.
    let avg_radius: f64 =
        points.iter().map(|p| p.distance_to(&center)).sum::<f64>() / points.len() as f64;

    if avg_radius < tolerance.linear() {
        return None;
    }

    // Verify ALL boundary points are equidistant from center.
    let max_deviation: f64 = points
        .iter()
        .map(|p| (p.distance_to(&center) - avg_radius).abs())
        .fold(0.0_f64, f64::max);

    // Loosened from 10× to 100× after §5.3: boundary sample
    // points come off NURBS fits with ~1e-5 drift on unit-scale
    // geometry, and the 10× window was rejecting legitimate cube-
    // corner spheres. The sphere is an exact analytic surface —
    // if its max_deviation fits within 100× tol_linear, it's the
    // right answer and beats the polar-NURBS fallback on both
    // representation clarity and allocation cost.
    // SCALE-RELATIVE gate: the boundary samples come off NURBS fits
    // whose drift scales with the corner radius, so an absolute
    // window rejects exact cube-corner spheres at radii above ~1
    // (measured: uniform r=3 corner — sphere declined, the swept
    // fallback shipped a boundary lying 0.0712·r off the strips, and
    // check 9 caught it at serving; at r ≤ 1 the same corner is
    // machine-exact). The sphere is still verified as the analytic
    // answer, just against a fit-noise window that grows with r.
    let boundary_diag = {
        let (mut lo, mut hi) = (
            Point3::new(f64::MAX, f64::MAX, f64::MAX),
            Point3::new(f64::MIN, f64::MIN, f64::MIN),
        );
        for p in points {
            lo = Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        lo.distance_to(&hi)
    };
    let sphere_gate = (tolerance.linear() * 100.0).max(boundary_diag * 1.0e-3);
    if std::env::var_os("PK_SPHEREFIT_TRACE").is_some() {
        eprintln!(
            "[spherefit] rolling-ball sphere candidate r={avg_radius:.4} max_dev={max_deviation:.4e} gate={sphere_gate:.4e}"
        );
    }
    if max_deviation > sphere_gate {
        return None;
    }

    let _ = CornerMethod::RollingBallSphere; // classification tag
    let surface = SurfaceGeometry::Sphere {
        center,
        radius: avg_radius,
    };

    Some(BlendSurfaceResult {
        surface,
        boundary_curves: boundary_curves.to_vec(),
        method: {
            if std::env::var_os("PK_SPHEREFIT_TRACE").is_some() {
                eprintln!("[corner-method] rolling-ball");
            }
            CornerMethod::RollingBallSphere
        },
        continuity: ContinuityLevel::G1,
    })
}
