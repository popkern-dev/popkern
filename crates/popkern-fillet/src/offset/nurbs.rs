//! NURBS-based surface offset with certified-by-measurement fitting.
//!
//! Used as a fallback when no closed-form offset exists (free-form
//! NURBS, Gregory, polyhedral spline). Exact offset points are sampled
//! on a uniform grid and approximated by a degree-3 least-squares
//! B-spline ([`super::nurbs_fit`]); the achieved accuracy is then
//! *measured* on an independent dense grid (finer than, and offset
//! from, the fit samples) and reported to the caller, who gates on it.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;

use super::error::OffsetError;
use super::nurbs_fit::fit_degree3_surface;

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

fn uniform_params(range: (f64, f64), m: usize) -> Vec<f64> {
    (0..m)
        .map(|i| range.0 + (range.1 - range.0) * (i as f64) / ((m - 1) as f64))
        .collect()
}

/// Maximum deviation between the exact offset and the fitted surface on
/// an independent dense grid: `2m+1` per direction, which contains both
/// the fit-sample parameters and every midpoint between them.
fn measure_deviation(
    surface: &SurfaceGeometry,
    distance: f64,
    sign: f64,
    u_range: (f64, f64),
    v_range: (f64, f64),
    fitted: &popkern_core::geometry::NurbsSurface,
    m: usize,
) -> f64 {
    let dense_u = uniform_params(u_range, 2 * m + 1);
    let dense_v = uniform_params(v_range, 2 * m + 1);
    let mut max_dev: f64 = 0.0;
    for &u in &dense_u {
        for &v in &dense_v {
            let exact = offset_point(surface, u, v, distance, sign);
            let approx = fitted.point_at(u, v);
            let dev = exact.distance_to(&approx);
            if dev > max_dev {
                max_dev = dev;
            }
        }
    }
    max_dev
}

/// Compute an offset of an arbitrary surface by degree-3 least-squares
/// fitting, reporting the measured deviation of the final fit.
///
/// The exact offset is sampled on a uniform grid (twice as dense as the
/// control net), fit by a degree-3 B-spline, and the fit's maximum
/// deviation from the exact offset is measured on an independent dense
/// grid. The control net is refined (doubled — bounded knot refinement
/// on uniform clamped knots) until the measured deviation is within
/// `tolerance.linear() / 4`, up to a 64×64 control-net ceiling.
///
/// Returns the surface together with the measured maximum deviation on
/// the final grid. At the ceiling the deviation may exceed the stage
/// tolerance — the caller decides whether that is acceptable (the
/// surface never silently claims to be in tolerance).
///
/// # Errors
///
/// [`OffsetError::FitFailed`] when the least-squares system degenerates.
pub fn offset_nurbs_measured(
    surface: &SurfaceGeometry,
    distance: f64,
    inward: bool,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<(SurfaceGeometry, f64), OffsetError> {
    let sign = if inward { -1.0 } else { 1.0 };
    let stage_tol = tolerance.linear() / 4.0;

    let mut n: usize = 8;
    loop {
        let m = 2 * n;
        let u_params = uniform_params(u_range, m);
        let v_params = uniform_params(v_range, m);
        let samples: Vec<Point3> = u_params
            .iter()
            .flat_map(|&u| {
                v_params
                    .iter()
                    .map(move |&v| offset_point(surface, u, v, distance, sign))
            })
            .collect();

        let fitted = fit_degree3_surface(&u_params, &v_params, &samples, n, n, u_range, v_range)?;
        let max_dev = measure_deviation(surface, distance, sign, u_range, v_range, &fitted, m);

        if max_dev <= stage_tol || n >= 64 {
            return Ok((SurfaceGeometry::Nurbs(fitted), max_dev));
        }
        n *= 2;
    }
}

/// Compute an offset of an arbitrary surface by degree-3 least-squares
/// fitting.
///
/// Convenience wrapper over [`offset_nurbs_measured`] that discards the
/// measured deviation. Quality-gated callers use the measured form.
pub fn offset_nurbs(
    surface: &SurfaceGeometry,
    distance: f64,
    inward: bool,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<SurfaceGeometry, OffsetError> {
    offset_nurbs_measured(surface, distance, inward, u_range, v_range, tolerance)
        .map(|(surface, _deviation)| surface)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::{NurbsSurface, Point3, Vec3};
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
        // A plane offset is itself a plane — the degree-3 fit must
        // reproduce it to numerical precision, not the old degree-1
        // "close enough at the samples" level.
        match &result {
            SurfaceGeometry::Nurbs(n) => {
                let p = n.point_at(0.3, -0.7);
                assert!(
                    (p.z - 3.0).abs() < 1e-9,
                    "expected z≈3.0, got z={:.12}",
                    p.z
                );
            }
            _ => panic!("expected NURBS"),
        }
    }

    #[test]
    fn curved_nurbs_offset_measured_within_tolerance() {
        // A genuinely curved base surface: the fit must land within the
        // stage tolerance and REPORT its measured deviation.
        let bump = NurbsSurface {
            degree_u: 2,
            degree_v: 2,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 0.0, 0.4),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 0.5, 0.4),
                Point3::new(0.5, 0.5, 0.8),
                Point3::new(1.0, 0.5, 0.4),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.5, 1.0, 0.4),
                Point3::new(1.0, 1.0, 0.0),
            ],
            num_u: 3,
            num_v: 3,
            weights: vec![1.0; 9],
        };
        let surface = SurfaceGeometry::Nurbs(bump);
        let tol = Tolerance::default();
        let (fitted, dev) =
            offset_nurbs_measured(&surface, 0.1, false, (0.0, 1.0), (0.0, 1.0), &tol).unwrap();
        assert!(dev.is_finite() && dev >= 0.0);
        assert!(
            dev <= tol.linear() / 4.0,
            "measured deviation {dev:.3e} exceeds stage tolerance"
        );
        // Independent spot check away from grid alignment.
        match &fitted {
            SurfaceGeometry::Nurbs(n) => {
                let exact = offset_point(&surface, 0.377, 0.613, 0.1, 1.0);
                let approx = n.point_at(0.377, 0.613);
                assert!(exact.distance_to(&approx) <= tol.linear());
            }
            _ => panic!("expected NURBS"),
        }
    }
}
