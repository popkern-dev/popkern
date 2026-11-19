//! Line-sphere intersection algorithm.

use crate::geometry::intersection::SurfaceCurveIntersection;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;

/// Intersect a line `P + t*D` with a sphere of given center `C` and radius `R`.
///
/// Solves `|P + t*D - C|² = R²`, a quadratic in `t`.
///
/// Surface parameters follow the Z-axis spherical convention used in `analytic_surfaces`:
/// - `u` = longitude = `atan2(Q.y, Q.x)`, normalised to `[0, 2π)`
/// - `v` = colatitude = `acos(Q.z / R)` in `[0, π]`
///
/// where `Q = intersection_point - center`.
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) fn intersect_line_sphere(
    line_origin: &Point3,
    line_dir: &Vec3,
    t_range: (f64, f64),
    sphere_center: &Point3,
    sphere_radius: f64,
    sphere_axis: &Vec3,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<SurfaceCurveIntersection> {
    // W = line_origin - sphere_center
    let w: Vec3 = *line_origin - *sphere_center;
    let d = *line_dir;
    let r = sphere_radius;

    // |W + t*D|² = R²
    // t²*(D·D) + 2t*(W·D) + (W·W - R²) = 0
    let a = d.dot(&d);
    let b = 2.0 * w.dot(&d);
    let c = w.dot(&w) - r * r;

    let disc = b * b - 4.0 * a * c;
    let eps2 = tolerance.linear() * tolerance.linear();
    if a < eps2 {
        // Degenerate: zero-length direction — no parameterised intersection
        return Vec::new();
    }
    // The line's distance h from the centre decides, as a length:
    // h² = R² − disc / 4a. Judged on `disc` against ε² (an area times
    // |D|⁴), rounding in the coordinates read a tangent line as a
    // crossing: a thousand units from the origin `disc` carried 2e-11
    // of noise, and the touch became two hits 2e-6 apart along the
    // line (an ice-cream cone's seam ruling, tangent to its ball at the
    // rim; the union's loops then failed to close, PK-1133). A tangent
    // amplifies a gap error δ to √(2Rδ) along the line; the tolerance
    // is the gap's, its band open at a whole tolerance as every
    // distance-decided touch's (`Tolerance::tangency_band`).
    let h = (r * r - disc / (4.0 * a)).max(0.0).sqrt();
    let scale = super::coordinate_scale(&[*line_origin, *sphere_center]) + r;
    let candidates: Vec<f64> = if h > r + tolerance.linear() {
        return Vec::new();
    } else if r - h < tolerance.tangency_band(scale) {
        // Tangent: the foot of the perpendicular from the centre.
        vec![-b / (2.0 * a)]
    } else {
        let sqrt_disc = disc.sqrt();
        vec![(-b - sqrt_disc) / (2.0 * a), (-b + sqrt_disc) / (2.0 * a)]
    };

    let tau = std::f64::consts::TAU;

    let mut results = Vec::new();

    for t in candidates {
        if t < t_range.0 - tolerance.linear() || t > t_range.1 + tolerance.linear() {
            continue;
        }

        let point = *line_origin + d * t;
        // In the sphere's chart frame (identity for a +Z pole).
        let q: Vec3 = crate::geometry::analytic_surfaces::sphere_to_local(
            sphere_axis,
            point - *sphere_center,
        );

        // Colatitude v = acos(q.z / R), clamped to handle floating-point noise
        let cos_v = (q.z / r).clamp(-1.0, 1.0);
        let v = cos_v.pacos(); // in [0, π]

        // Longitude u = atan2(q.y, q.x), normalised to [0, 2π)
        let u_raw = q.y.patan2(q.x);
        let u = ((u_raw % tau) + tau) % tau;

        // Filter by surface v_range
        if v < v_range.0 - tolerance.angular() || v > v_range.1 + tolerance.angular() {
            continue;
        }

        // Filter by surface u_range (handle full-circle)
        if !u_in_range(u, u_range, tolerance.angular()) {
            continue;
        }

        results.push(SurfaceCurveIntersection::Point { point, t, u, v });
    }

    // Sort by t parameter
    results.sort_by(|a, b| {
        a.t_val()
            .partial_cmp(&b.t_val())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    results
}

/// Returns true if `u` is within the given u_range (handles full-circle wrap-around).
fn u_in_range(u: f64, u_range: (f64, f64), angular_tol: f64) -> bool {
    let tau = std::f64::consts::TAU;
    let (lo, hi) = u_range;
    if hi - lo >= tau - angular_tol {
        return true;
    }
    u >= lo - angular_tol && u <= hi + angular_tol
}

/// Helper to extract `t` from a `SurfaceCurveIntersection::Point`.
trait TVal {
    fn t_val(&self) -> f64;
}

impl TVal for SurfaceCurveIntersection {
    fn t_val(&self) -> f64 {
        match self {
            SurfaceCurveIntersection::Point { t, .. } => *t,
            SurfaceCurveIntersection::Overlap { t_range, .. } => t_range.0,
        }
    }
}

#[cfg(test)]
#[path = "line_sphere_tests.rs"]
mod tests;
