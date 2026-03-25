//! Loft surface generation through cross-sections.
//!
//! Interpolates a NURBS surface through a sequence of section curves using
//! chord-length parameterization and column-wise curve fitting.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::fitting::fit_curve;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3};
use popkern_core::tolerance::Tolerance;

use crate::error::SweepError;
use crate::surface::loft_continuity::{ContinuityEnd, apply_continuity};
use crate::types::{LoftOptions, LoftSection};

/// Number of uniform samples for resampling incompatible sections.
const RESAMPLE_COUNT: usize = 32;

/// Result of loft surface generation.
#[derive(Debug)]
pub struct LoftSurfaceResult {
    /// The lofted NURBS surface.
    pub surface: NurbsSurface,
    /// Start cap curve (first section).
    pub start_cap: NurbsCurve,
    /// End cap curve (last section).
    pub end_cap: NurbsCurve,
}

/// Generate a loft surface through a sequence of cross-sections.
///
/// Sections are made compatible (same number of control points), then
/// chord-length parameterization determines v-parameters for each section.
/// Column-wise curve fitting through corresponding control points produces
/// the v-direction of the resulting surface.
///
/// # Errors
///
/// Returns `SweepError::IncompatibleSections` if fewer than 2 sections
/// are provided or if section compatibility cannot be achieved.
pub fn generate_loft_surface(
    sections: &[LoftSection],
    options: &LoftOptions,
    tolerance: &Tolerance,
) -> Result<LoftSurfaceResult, SweepError> {
    // Validate section count.
    if sections.len() < 2 {
        return Err(SweepError::IncompatibleSections(format!(
            "need at least 2 sections, got {}",
            sections.len()
        )));
    }

    // Step 1: Make sections compatible (same control point count).
    let compatible = make_sections_compatible(sections, tolerance)?;

    // Step 2: Compute v-parameterization via chord-length.
    let v_params = compute_v_params(&compatible, options);

    // Step 3: Column interpolation to build the surface.
    let mut surface = interpolate_columns(&compatible, &v_params)?;

    // Step 4: Apply G1/G2 continuity constraints.
    apply_continuity(
        &mut surface,
        &options.start_continuity,
        ContinuityEnd::Start,
    )?;
    apply_continuity(&mut surface, &options.end_continuity, ContinuityEnd::End)?;

    let start_cap = compatible[0].clone();
    let end_cap = compatible[compatible.len() - 1].clone();

    Ok(LoftSurfaceResult {
        surface,
        start_cap,
        end_cap,
    })
}

/// Make all sections compatible by ensuring they have the same number of
/// control points.
///
/// If all sections already share the same control point count, they are
/// returned as-is. Otherwise, all sections are resampled at uniform
/// parameter intervals to produce polyline sections with matching point
/// counts.
fn make_sections_compatible(
    sections: &[LoftSection],
    tolerance: &Tolerance,
) -> Result<Vec<NurbsCurve>, SweepError> {
    let curves: Vec<&NurbsCurve> = sections.iter().map(|s| &s.curve).collect();

    // Check if all sections already have the same control point count.
    let first_count = curves[0].control_points.len();
    let all_same = curves.iter().all(|c| c.control_points.len() == first_count);

    if all_same {
        // Also check if all have the same degree — if so, use as-is.
        let first_degree = curves[0].degree;
        let same_degree = curves.iter().all(|c| c.degree == first_degree);

        if same_degree {
            // Check knot vectors are compatible (same length, near-equal values).
            let first_knots = &curves[0].knots;
            let knots_match = curves.iter().all(|c| {
                c.knots.len() == first_knots.len()
                    && c.knots
                        .iter()
                        .zip(first_knots.iter())
                        .all(|(a, b)| (a - b).abs() < tolerance.linear)
            });

            if knots_match {
                return Ok(curves.iter().map(|c| (*c).clone()).collect());
            }
        }
    }

    // Sections are incompatible — resample all to uniform polylines.
    resample_sections(&curves)
}

/// Resample each section curve at N uniform parameter values, producing
/// degree-1 polyline curves with identical point counts.
fn resample_sections(curves: &[&NurbsCurve]) -> Result<Vec<NurbsCurve>, SweepError> {
    let n = RESAMPLE_COUNT;
    let mut result = Vec::with_capacity(curves.len());

    for curve in curves {
        let (t_min, t_max) = curve.parameter_range();
        let mut points = Vec::with_capacity(n);

        for i in 0..n {
            let frac = i as f64 / (n - 1).max(1) as f64;
            let t = t_min + frac * (t_max - t_min);
            points.push(curve.point_at(t));
        }

        // Build a degree-1 polyline through these points.
        let mut knots = Vec::with_capacity(n + 1 + 1);
        knots.push(0.0); // clamped start
        for i in 0..n {
            knots.push(i as f64 / (n - 1).max(1) as f64);
        }
        knots.push(1.0); // clamped end
        // Actually for degree 1 with n CPs: knot vector has n + 2 entries, clamped.
        let knots = build_clamped_knots(n, 1);
        let weights = vec![1.0; n];

        result.push(NurbsCurve {
            degree: 1,
            knots,
            control_points: points,
            weights,
        });
    }

    Ok(result)
}

/// Compute chord-length v-parameterization from section centroids.
fn compute_v_params(sections: &[NurbsCurve], options: &LoftOptions) -> Vec<f64> {
    let n = sections.len();
    if n == 0 {
        return vec![];
    }
    if n == 1 {
        return vec![0.0];
    }

    let centroids: Vec<Point3> = sections.iter().map(compute_centroid).collect();

    let mut cumulative = vec![0.0_f64; n];
    for i in 1..n {
        let dist = centroids[i].distance_to(&centroids[i - 1]);
        cumulative[i] = cumulative[i - 1] + dist;
    }

    let total = cumulative[n - 1];
    if total < 1e-30 {
        // All centroids coincide — use uniform spacing.
        return (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
    }

    let mut params: Vec<f64> = cumulative.iter().map(|c| c / total).collect();

    if options.closed {
        // For closed loft, the last param wraps to 1.0 (period).
        if let Some(last) = params.last_mut() {
            *last = 1.0;
        }
    }

    params
}

/// Interpolate columns of control points in the v-direction to build
/// the loft surface.
///
/// For each control point index j across all sections, fits a degree-3
/// curve through those points at the given v-parameters. The resulting
/// control points form the surface grid.
fn interpolate_columns(
    sections: &[NurbsCurve],
    v_params: &[f64],
) -> Result<NurbsSurface, SweepError> {
    let num_sections = sections.len();
    let num_u = sections[0].control_points.len();

    // Determine v-direction degree: min(3, num_sections - 1).
    let v_degree = 3.min(num_sections - 1);

    // For exact interpolation: num_control_points == num_data_points.
    // But fit_curve requires num_control_points <= num_data_points
    // and num_control_points >= degree + 1.
    let num_v_cp = num_sections;
    if num_v_cp < v_degree + 1 {
        return Err(SweepError::IncompatibleSections(format!(
            "need at least {} sections for degree-{} interpolation, got {}",
            v_degree + 1,
            v_degree,
            num_sections
        )));
    }

    // Fit a curve through each column of control points.
    let mut first_knots_v: Option<Vec<f64>> = None;
    let mut surface_cps = Vec::with_capacity(num_u * num_v_cp);
    let mut surface_weights = Vec::with_capacity(num_u * num_v_cp);

    for j in 0..num_u {
        // Collect j-th control point from each section.
        let column_points: Vec<Point3> = sections.iter().map(|s| s.control_points[j]).collect();

        let fitted = fit_curve(&column_points, v_params, v_degree, num_v_cp)
            .map_err(|e| SweepError::IncompatibleSections(format!("column fit failed: {e}")))?;

        // All columns must produce the same knot vector.
        match &first_knots_v {
            Some(knots) => {
                if knots.len() != fitted.knots.len() {
                    return Err(SweepError::IncompatibleSections(
                        "inconsistent v-knot vectors across columns".into(),
                    ));
                }
            }
            None => {
                first_knots_v = Some(fitted.knots.clone());
            }
        }

        // The fitted curve's control points become row j of the surface grid.
        // Surface layout is row-major: cp(i_u, i_v) = cps[i_u * num_v + i_v].
        surface_cps.extend_from_slice(&fitted.control_points);
        surface_weights.extend_from_slice(&fitted.weights);
    }

    let knots_v = first_knots_v
        .ok_or_else(|| SweepError::IncompatibleSections("no columns to interpolate".into()))?;

    // u-direction uses the section's knot vector and degree.
    let degree_u = sections[0].degree;
    let knots_u = sections[0].knots.clone();

    Ok(NurbsSurface {
        degree_u,
        degree_v: v_degree,
        knots_u,
        knots_v,
        control_points: surface_cps,
        num_u,
        num_v: num_v_cp,
        weights: surface_weights,
    })
}

/// Compute centroid of a curve's control points.
fn compute_centroid(curve: &NurbsCurve) -> Point3 {
    let n = curve.control_points.len();
    if n == 0 {
        return Point3::new(0.0, 0.0, 0.0);
    }
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sz = 0.0;
    for p in &curve.control_points {
        sx += p.x;
        sy += p.y;
        sz += p.z;
    }
    let inv = 1.0 / n as f64;
    Point3::new(sx * inv, sy * inv, sz * inv)
}

/// Build a clamped uniform knot vector for `n` control points at given degree.
fn build_clamped_knots(n: usize, degree: usize) -> Vec<f64> {
    let m = n + degree + 1;
    let mut knots = vec![0.0; degree + 1];

    let n_interior = m.saturating_sub(2 * (degree + 1));
    for i in 1..=n_interior {
        knots.push(i as f64 / (n_interior + 1) as f64);
    }

    knots.extend(std::iter::repeat_n(1.0, degree + 1));
    knots
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::surface::SurfaceEval;

    /// Create a unit circle as a rational NURBS (degree 2, 9 control points).
    fn unit_circle_at_z(z: f64, radius: f64) -> NurbsCurve {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let r = radius;
        NurbsCurve {
            degree: 2,
            knots: vec![
                0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
            ],
            control_points: vec![
                Point3::new(r, 0.0, z),
                Point3::new(r, r, z),
                Point3::new(0.0, r, z),
                Point3::new(-r, r, z),
                Point3::new(-r, 0.0, z),
                Point3::new(-r, -r, z),
                Point3::new(0.0, -r, z),
                Point3::new(r, -r, z),
                Point3::new(r, 0.0, z),
            ],
            weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
        }
    }

    fn make_section(curve: NurbsCurve) -> LoftSection {
        LoftSection {
            curve,
            position: None,
        }
    }

    #[test]
    fn loft_two_identical_circles() {
        let sections = vec![
            make_section(unit_circle_at_z(0.0, 1.0)),
            make_section(unit_circle_at_z(1.0, 1.0)),
        ];
        let options = LoftOptions::default();
        let tol = Tolerance::default();
        let result = generate_loft_surface(&sections, &options, &tol);
        assert!(result.is_ok(), "loft failed: {:?}", result.err());

        let result = result.unwrap(); // ok in tests
        let surf = &result.surface;

        // Should have 9 control points in u (from circle) and 2 in v.
        assert_eq!(surf.num_u, 9);
        assert_eq!(surf.num_v, 2);

        // Evaluate at the midpoint: should be approximately on a unit circle at z=0.5.
        let mid = surf.point_at(0.0, 0.5);
        assert!((mid.z - 0.5).abs() < 0.3, "expected z~0.5, got {}", mid.z);
        let r = (mid.x * mid.x + mid.y * mid.y).sqrt();
        assert!((r - 1.0).abs() < 0.3, "expected radius~1.0, got {}", r);
    }

    #[test]
    fn loft_circle_to_larger_circle() {
        let sections = vec![
            make_section(unit_circle_at_z(0.0, 1.0)),
            make_section(unit_circle_at_z(1.0, 2.0)),
        ];
        let options = LoftOptions::default();
        let tol = Tolerance::default();
        let result = generate_loft_surface(&sections, &options, &tol);
        assert!(result.is_ok(), "loft failed: {:?}", result.err());

        let result = result.unwrap(); // ok in tests
        let surf = &result.surface;

        // At v=0, should be on the smaller circle.
        let p0 = surf.point_at(0.0, 0.0);
        let r0 = (p0.x * p0.x + p0.y * p0.y).sqrt();
        assert!((r0 - 1.0).abs() < 0.3, "expected r~1.0 at v=0, got {}", r0);

        // At v=1, should be on the larger circle.
        let p1 = surf.point_at(0.0, 1.0);
        let r1 = (p1.x * p1.x + p1.y * p1.y).sqrt();
        assert!((r1 - 2.0).abs() < 0.3, "expected r~2.0 at v=1, got {}", r1);
    }

    #[test]
    fn loft_three_sections() {
        let sections = vec![
            make_section(unit_circle_at_z(0.0, 1.0)),
            make_section(unit_circle_at_z(1.0, 1.5)),
            make_section(unit_circle_at_z(2.0, 1.0)),
        ];
        let result =
            generate_loft_surface(&sections, &LoftOptions::default(), &Tolerance::default());
        let surf = &result.unwrap().surface; // ok in tests
        assert_eq!(surf.num_u, 9);
        assert_eq!(surf.num_v, 3);
        assert!(surf.point_at(0.0, 0.0).z.abs() < 0.3);
        assert!((surf.point_at(0.0, 1.0).z - 2.0).abs() < 0.3);
    }

    #[test]
    fn loft_too_few_sections_errors() {
        let opts = LoftOptions::default();
        let tol = Tolerance::default();
        // Zero sections.
        let r0 = generate_loft_surface(&[], &opts, &tol);
        assert!(r0.unwrap_err().to_string().contains("at least 2"));
        // One section.
        let r1 = generate_loft_surface(&[make_section(unit_circle_at_z(0.0, 1.0))], &opts, &tol);
        assert!(r1.unwrap_err().to_string().contains("at least 2"));
    }

    #[test]
    fn loft_incompatible_sections_resampled() {
        // Circle (9 CPs, degree 2) vs line segment (2 CPs, degree 1).
        let circle = unit_circle_at_z(0.0, 1.0);
        let line = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(-1.0, 0.0, 1.0), Point3::new(1.0, 0.0, 1.0)],
            weights: vec![1.0, 1.0],
        };
        let sections = vec![make_section(circle), make_section(line)];
        let options = LoftOptions::default();
        let tol = Tolerance::default();
        let result = generate_loft_surface(&sections, &options, &tol);
        // Should succeed after resampling.
        assert!(result.is_ok(), "loft failed: {:?}", result.err());

        let result = result.unwrap(); // ok in tests
        // Both should have been resampled to RESAMPLE_COUNT points.
        assert_eq!(result.surface.num_u, RESAMPLE_COUNT);
    }

    #[test]
    fn loft_caps_match_sections() {
        let sections = vec![
            make_section(unit_circle_at_z(0.0, 1.0)),
            make_section(unit_circle_at_z(3.0, 1.0)),
        ];
        let result =
            generate_loft_surface(&sections, &LoftOptions::default(), &Tolerance::default());
        let result = result.unwrap(); // ok in tests
        assert!(result.start_cap.point_at(0.0).z.abs() < 0.01);
        assert!((result.end_cap.point_at(0.0).z - 3.0).abs() < 0.01);
    }

    #[test]
    fn v_params_uniform_for_coincident_centroids() {
        // Two circles at the same location.
        let sections = vec![
            unit_circle_at_z(0.0, 1.0),
            unit_circle_at_z(0.0, 1.0),
            unit_circle_at_z(0.0, 1.0),
        ];
        let options = LoftOptions::default();
        let params = compute_v_params(&sections, &options);
        assert_eq!(params.len(), 3);
        assert!((params[0]).abs() < 1e-10);
        assert!((params[1] - 0.5).abs() < 1e-10);
        assert!((params[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn build_clamped_knots_basic() {
        let knots = build_clamped_knots(5, 2);
        // Length = 5 + 2 + 1 = 8
        assert_eq!(knots.len(), 8);
        assert_eq!(knots[0], 0.0);
        assert_eq!(knots[1], 0.0);
        assert_eq!(knots[2], 0.0);
        assert_eq!(knots[7], 1.0);
        assert_eq!(knots[6], 1.0);
        assert_eq!(knots[5], 1.0);
    }
}
