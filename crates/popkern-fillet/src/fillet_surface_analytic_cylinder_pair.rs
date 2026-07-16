//! Analytic fillet surfaces for parallel cylinder-cylinder edges and
//! exact circle fitting.

use crate::center_curve::CenterCurveResult;
use crate::fillet_surface::FilletSurfaceResult;
use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;

/// Try to build an analytic cylinder for a fillet between two cylinders
/// with **parallel axes** (§11.1).
///
/// **Geometry.** Two parallel infinite cylinders intersect in 0, 1, or
/// 2 straight lines, parallel to the common axis direction. Along one
/// such intersection line the rolling ball of radius `r` sits at a
/// fixed 2D position (x_c, y_c) in the cross-section plane, determined
/// by the circle-circle tangency problem `|center − axis_A| = R_A ± r`
/// and similarly for cylinder B. That position is constant in the
/// axial direction, so the center curve is a straight line parallel to
/// both cylinder axes and the fillet surface swept by the ball is an
/// exact cylinder:
///
/// - **Axis** = cylinder_axis_direction (shared by both supports)
/// - **Origin** = any point on the center curve
/// - **Radius** = the fillet radius
///
/// The detector reads the center-curve direction from sampling the
/// already-computed `center_result.center_curve` and bails if:
///
/// - the cylinders are not parallel (dot < 0.999);
/// - the center curve is not a straight line (max lateral residual
///   > 0.1 % of the curve span); or
/// - the center-curve direction is not parallel to the cylinder axes
///   (dot < 0.999) — which also catches cylinder-cylinder cases where
///   the axes happen to be parallel but the center curve isn't
///   (degenerate setups).
pub(crate) fn try_parallel_cylinder_cylinder(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    center_result: &CenterCurveResult,
    radius: f64,
) -> Option<FilletSurfaceResult> {
    let (axis1, axis2) = match (surface1, surface2) {
        (
            SurfaceGeometry::Cylinder { axis: a1, .. },
            SurfaceGeometry::Cylinder { axis: a2, .. },
        ) => (*a1, *a2),
        _ => return None,
    };
    let a1_n = axis1.normalized()?;
    let a2_n = axis2.normalized()?;
    if a1_n.dot(&a2_n).abs() < 0.999 {
        return None;
    }

    use popkern_core::geometry::curve::CurveEval;
    let n = 16usize;
    let (t0, t1) = center_result.center_curve.parameter_range();
    let samples: Vec<Point3> = (0..n)
        .map(|i| {
            let t = t0 + (t1 - t0) * (i as f64) / ((n - 1).max(1) as f64);
            center_result.center_curve.point_at(t)
        })
        .collect();
    if samples.len() < 3 {
        return None;
    }

    let p0 = samples[0];
    let p_last = samples[samples.len() - 1];
    let span = p_last - p0;
    let span_len = span.length();
    if span_len < 1e-12 {
        return None;
    }
    let dir = span.normalized()?;

    // The center curve direction must be parallel to the cylinder
    // axes. This is the load-bearing check — without it, a non-
    // straight center curve could still have a matching span
    // direction if the endpoints happen to line up.
    if dir.dot(&a1_n).abs() < 0.999 {
        return None;
    }

    // All intermediate samples must lie on the line through p0 in
    // direction `dir`. Max perpendicular residual < 0.1 % of the
    // span length. A straight line has residual ≈ 0 (floating-point
    // noise only); a curved center curve (e.g. from two cylinders
    // with non-parallel axes that accidentally passed the axis
    // check) will have larger residual and bail.
    let max_residual = samples
        .iter()
        .map(|p| {
            let v = *p - p0;
            let along = v.dot(&dir);
            let perp = v - dir * along;
            perp.length()
        })
        .fold(0.0_f64, f64::max);
    if max_residual > span_len * 1e-3 {
        return None;
    }

    let surface = SurfaceGeometry::Cylinder {
        origin: p0,
        axis: a1_n,
        radius,
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

/// Kåsa least-squares circle through near-circular planar points:
/// center, radius. Returns `None` unless the points are planar and
/// on the fitted circle to 1e-6 relative — intended for rolling-ball
/// contact footprints and spines, which are exact to the solver
/// residual. (The centroid is NOT the center when samples are
/// unevenly spaced; adaptive walker steps put 2e-2 of apparent
/// radial spread on an exact circle.)
pub(crate) fn fit_exact_circle_center(points: &[Point3]) -> Option<(Point3, f64)> {
    if points.len() < 4 {
        return None;
    }
    let n = points.len() as f64;
    let c0 = Point3::new(
        points.iter().map(|p| p.x).sum::<f64>() / n,
        points.iter().map(|p| p.y).sum::<f64>() / n,
        points.iter().map(|p| p.z).sum::<f64>() / n,
    );
    let mut nx = 0.0;
    let mut ny = 0.0;
    let mut nz = 0.0;
    for i in 0..points.len() {
        let a = points[i];
        let b = points[(i + 1) % points.len()];
        nx += (a.y - b.y) * (a.z + b.z);
        ny += (a.z - b.z) * (a.x + b.x);
        nz += (a.x - b.x) * (a.y + b.y);
    }
    let normal = popkern_core::geometry::Vec3::new(nx, ny, nz).normalized()?;
    let (e1, e2) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(normal);
    let xy: Vec<(f64, f64)> = points
        .iter()
        .map(|p| {
            let d = *p - c0;
            (d.dot(&e1), d.dot(&e2))
        })
        .collect();
    let (mut sx, mut sy, mut sxx, mut syy, mut sxy, mut sxz, mut syz, mut sz) =
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for &(x, y) in &xy {
        let z = x * x + y * y;
        sx += x;
        sy += y;
        sxx += x * x;
        syy += y * y;
        sxy += x * y;
        sxz += x * z;
        syz += y * z;
        sz += z;
    }
    let m = [[sxx, sxy, sx], [sxy, syy, sy], [sx, sy, n]];
    let rhs = [sxz, syz, sz];
    let det3 = |m: &[[f64; 3]; 3]| -> f64 {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d0 = det3(&m);
    if d0.abs() < 1e-30 {
        return None;
    }
    let solve_col = |col: usize| -> f64 {
        let mut mm = m;
        for (row, mrow) in mm.iter_mut().enumerate() {
            mrow[col] = rhs[row];
        }
        det3(&mm) / d0
    };
    let (a, b, c2) = (solve_col(0), solve_col(1), solve_col(2));
    let (cx2, cy2) = (0.5 * a, 0.5 * b);
    let r = (c2 + cx2 * cx2 + cy2 * cy2).max(0.0).sqrt();
    if r < 1e-12 {
        return None;
    }
    let center = c0 + e1 * cx2 + e2 * cy2;
    let gate = r * 1e-6;
    for p in points {
        let d = *p - center;
        if (d.length() - r).abs() > gate || d.dot(&normal).abs() > gate {
            return None;
        }
    }
    Some((center, r))
}
