//! Line-cone intersection algorithm.

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::intersection::{SurfaceCurveIntersection, coordinate_scale};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;

/// Intersect a line `P + t*D` with a cone defined by apex, unit axis `A`, and half-angle `α`.
///
/// The cone surface satisfies `(Q·A)² - cos²(α)*|Q|² = 0` where `Q = P + t*D - apex`.
/// Only the forward nappe (`Q·A > 0`) is considered.
///
/// Returns surface parameters:
/// - `v = Q·A` (height along axis from apex)
/// - `u = atan2(Q·perp2, Q·perp1)` (azimuth angle, normalised to `[0, 2π)`)
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) fn intersect_line_cone(
    line_origin: &Point3,
    line_dir: &Vec3,
    t_range: (f64, f64),
    cone_apex: &Point3,
    cone_axis: &Vec3,
    cone_half_angle: f64,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<SurfaceCurveIntersection> {
    // W = line_origin - cone_apex
    let w = *line_origin - *cone_apex;
    let d = *line_dir;
    let a = *cone_axis;

    let cos_alpha = cone_half_angle.pcos();
    let cos2 = cos_alpha * cos_alpha;

    // Quadratic coefficients from expanding (Q·A)² - cos²(α)*|Q|² = 0
    // Q = W + t*D
    // (W·A + t*(D·A))² - cos²(α) * |W + t*D|² = 0
    let da = d.dot(&a);
    let wa = w.dot(&a);
    let dd = d.dot(&d);
    let wd = w.dot(&d);
    let ww = w.dot(&w);

    let qa = da * da - cos2 * dd;
    let qb = 2.0 * (da * wa - cos2 * wd);
    let qc = wa * wa - cos2 * ww;

    let disc = qb * qb - 4.0 * qa * qc;

    // Threshold for treating as degenerate quadratic
    let eps2 = tolerance.linear() * tolerance.linear();

    let candidates: Vec<f64> = if qa < 0.0 && qa.abs() >= eps2 {
        crossing_or_touch(&w, &d, cone_half_angle, (qa, qb, qc), cone_apex, tolerance)
    } else if qa.abs() < eps2 {
        // Degenerate: linear in t (or no solution)
        if qb.abs() < eps2 {
            // qc must also be ~0 for a solution; no isolated t
            return Vec::new();
        }
        vec![-qc / qb]
    } else if disc < -eps2 {
        // No real roots
        return Vec::new();
    } else if disc <= eps2 {
        // Tangent — one root
        vec![-qb / (2.0 * qa)]
    } else {
        let sqrt_disc = disc.sqrt();
        vec![
            (-qb - sqrt_disc) / (2.0 * qa),
            (-qb + sqrt_disc) / (2.0 * qa),
        ]
    };

    let (perp1, perp2) = orthonormal_basis(a);
    let tau = std::f64::consts::TAU;

    let mut results = Vec::new();

    for t in candidates {
        // Filter by line t_range
        if t < t_range.0 - tolerance.linear() || t > t_range.1 + tolerance.linear() {
            continue;
        }

        let q = w + d * t; // Q = P + t*D - apex

        // Forward nappe filter: v = Q·A must be > 0 (within tolerance)
        let v = q.dot(&a);
        if v < -tolerance.linear() {
            continue;
        }

        // Filter by surface v_range
        if v < v_range.0 - tolerance.linear() || v > v_range.1 + tolerance.linear() {
            continue;
        }

        // Compute azimuth u from projection onto the orthonormal basis
        let q_perp1 = q.dot(&perp1);
        let q_perp2 = q.dot(&perp2);
        let u_raw = q_perp2.patan2(q_perp1);
        // Normalise u to [0, 2π)
        let u = ((u_raw % tau) + tau) % tau;

        // Filter by surface u_range
        if !u_in_range(u, u_range, tolerance.angular()) {
            continue;
        }

        let point = *cone_apex + a * v + perp1 * q_perp1 + perp2 * q_perp2;

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

/// The line parameters where a line meeting the cone across its
/// outside (`qa < 0`: the cone's form `F = (Q·A)² − cos²α|Q|²` peaks
/// along the line, positive inside) meets it: none when the peak stays
/// outside past rounding, the peak alone when it reaches inside by less
/// than the tangency band (a touch), else the two crossings. Decided by
/// distance, as the cylinder decides it (PK-789): the peak's depth is
/// `F / |∇F|`, `|∇F| = 2|Q| sin α cos α` on the cone. The discriminant
/// band it replaces (`linear²` in the form's units, which grow as the
/// square of the scale) read a plate edge touching a kilometre dimple's
/// rim as a miss and the cut declined (PK-1143).
fn crossing_or_touch(
    w: &Vec3,
    d: &Vec3,
    half_angle: f64,
    (qa, qb, qc): (f64, f64, f64),
    apex: &Point3,
    tolerance: &Tolerance,
) -> Vec<f64> {
    let t_peak = -qb / (2.0 * qa);
    let q = *w + *d * t_peak;
    let peak = qa * t_peak * t_peak + qb * t_peak + qc;
    let grad = 2.0 * q.length() * half_angle.psin() * half_angle.pcos();
    let scale = coordinate_scale(&[*apex + *w, *apex]);
    if grad <= f64::EPSILON * scale {
        // The peak at the apex: no tangent plane to measure from.
        let disc = qb * qb - 4.0 * qa * qc;
        return if disc < 0.0 {
            Vec::new()
        } else {
            let root = disc.sqrt();
            vec![(-qb + root) / (2.0 * qa), (-qb - root) / (2.0 * qa)]
        };
    }
    let depth = peak / grad;
    if depth < -64.0 * f64::EPSILON * scale {
        Vec::new()
    } else if depth < tolerance.tangency_band(scale) {
        vec![t_peak]
    } else {
        let root = (qb * qb - 4.0 * qa * qc).max(0.0).sqrt();
        vec![(-qb + root) / (2.0 * qa), (-qb - root) / (2.0 * qa)]
    }
}

/// Returns true if `u` is within the given u_range (handles full-circle wrap-around).
fn u_in_range(u: f64, u_range: (f64, f64), angular_tol: f64) -> bool {
    let tau = std::f64::consts::TAU;
    let (lo, hi) = u_range;
    // Full circle
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
#[path = "line_cone_tests.rs"]
mod tests;
