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
    generate_loft_surface_at_stations(sections, options, tolerance, None)
}

/// Loft with the v-parameterization supplied by the caller.
///
/// Span-wise lofting MUST use one shared station parameterization for
/// every span family: per-span chord-length params (each span's own
/// centroid chain) interpolate the shared boundary columns at
/// DIFFERENT parameters, so adjacent span surfaces disagreed along
/// their common boundary between sections — measured 9.1e-5 at a
/// propeller blade's girdle junctions, where the two boundary-crossing
/// refinements each landed exactly on their own span's boundary iso
/// and the merge graph minted disjoint junction nodes (PK-443). With
/// shared stations the boundary columns pose identical interpolation
/// problems on both sides, so adjacent spans share their boundary
/// curve exactly.
///
/// `stations` must be strictly increasing with one value per section;
/// `None` computes chord-length parameters from these sections'
/// centroids (the single-surface path).
///
/// # Errors
///
/// Returns `SweepError::IncompatibleSections` as
/// [`generate_loft_surface`] does, or when `stations` has the wrong
/// length.
pub fn generate_loft_surface_at_stations(
    sections: &[LoftSection],
    options: &LoftOptions,
    tolerance: &Tolerance,
    stations: Option<&[f64]>,
) -> Result<LoftSurfaceResult, SweepError> {
    // Validate section count.
    if sections.len() < 2 {
        return Err(SweepError::IncompatibleSections(format!(
            "need at least 2 sections, got {}",
            sections.len()
        )));
    }
    if let Some(stations) = stations
        && stations.len() != sections.len()
    {
        return Err(SweepError::IncompatibleSections(format!(
            "{} stations for {} sections",
            stations.len(),
            sections.len()
        )));
    }

    // Step 1: Make sections compatible (same control point count).
    let compatible = make_sections_compatible(sections, tolerance)?;

    // Step 2: v-parameterization — the caller's shared stations, else
    // chord-length from these sections' centroids.
    let v_params = match stations {
        Some(stations) => stations.to_vec(),
        None => compute_v_params(&compatible, options),
    };

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
                        .all(|(a, b)| (a - b).abs() < tolerance.linear())
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

/// Chord-length station parameters from the FULL section curves — the
/// shared v-parameterization for span-wise lofting (one station per
/// section, owned by the section, not by any span family). See
/// [`generate_loft_surface_at_stations`].
pub fn loft_station_params(sections: &[LoftSection], options: &LoftOptions) -> Vec<f64> {
    let curves: Vec<NurbsCurve> = sections.iter().map(|s| s.curve.clone()).collect();
    compute_v_params(&curves, options)
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
        // Interpolate the j-th column in HOMOGENEOUS space: rational
        // sections (arcs, circles) carry their shape in the weights,
        // and fitting the Euclidean control points alone stamped
        // weight 1.0 on every column — every arc span lofted to a
        // PARABOLIC sheet 6.1% of the arc radius off the true
        // surface, while the cap rim edges carried the exact rational
        // arcs (measured: 1.52e-3 at a blade root's 0.025 arc,
        // 3.66e-4 at its 0.006 tip, both 6.1%·r — PK-443). With the
        // square interpolation system both (w·P) and (w) reproduce
        // the section rows exactly, so every station curve of the
        // surface IS its section, and constant-weight columns stay
        // exactly rational in between.
        let column_w: Vec<f64> = sections.iter().map(|s| s.weights[j]).collect();
        let weighted_points: Vec<Point3> = sections
            .iter()
            .zip(&column_w)
            .map(|(s, &w)| {
                let p = s.control_points[j];
                Point3::new(p.x * w, p.y * w, p.z * w)
            })
            .collect();
        let weight_points: Vec<Point3> =
            column_w.iter().map(|&w| Point3::new(w, 0.0, 0.0)).collect();

        let fitted_wp = fit_curve(&weighted_points, v_params, v_degree, num_v_cp)
            .map_err(|e| SweepError::IncompatibleSections(format!("column fit failed: {e}")))?;
        let fitted_w = fit_curve(&weight_points, v_params, v_degree, num_v_cp)
            .map_err(|e| SweepError::IncompatibleSections(format!("weight fit failed: {e}")))?;

        // All columns must produce the same knot vector.
        match &first_knots_v {
            Some(knots) => {
                if knots.len() != fitted_wp.knots.len() {
                    return Err(SweepError::IncompatibleSections(
                        "inconsistent v-knot vectors across columns".into(),
                    ));
                }
            }
            None => {
                first_knots_v = Some(fitted_wp.knots.clone());
            }
        }

        // The fitted homogeneous points become row j of the surface
        // grid (row-major: cp(i_u, i_v) = cps[i_u * num_v + i_v]).
        for (wp, w_cp) in fitted_wp
            .control_points
            .iter()
            .zip(&fitted_w.control_points)
        {
            let w = w_cp.x;
            if w <= 1.0e-12 || !w.is_finite() {
                return Err(SweepError::IncompatibleSections(format!(
                    "column {j} weight interpolation left the rational cone (w = {w:.3e})"
                )));
            }
            surface_cps.push(Point3::new(wp.x / w, wp.y / w, wp.z / w));
            surface_weights.push(w);
        }
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
#[path = "loft_surface_tests.rs"]
mod tests;
