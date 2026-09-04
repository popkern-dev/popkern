//! Sweep surface generation via cross-section skinning.
//!
//! For each smooth profile span, transforms the profile control points
//! by each frame in the sequence (position + normal/binormal orientation),
//! applies scaling, and skins the result into a NURBS surface.

use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3};
use popkern_core::tolerance::Tolerance;

use crate::error::SweepError;
use crate::frame::{FrameSequence, PathFrame};
use crate::profile::AnalyzedProfile;
use crate::types::ScaleMethod;

/// Result of sweep surface generation.
#[derive(Debug)]
pub struct SweepSurfaceResult {
    /// One NURBS surface per smooth profile span.
    pub surfaces: Vec<NurbsSurface>,
    /// Start cap profile (transformed to path start).
    pub start_cap: Vec<NurbsCurve>,
    /// End cap profile (transformed to path end).
    pub end_cap: Vec<NurbsCurve>,
}

/// Generate sweep surface(s) from a profile and frame sequence.
///
/// Each profile span is transformed by every frame in the sequence,
/// producing a grid of control points that becomes a NURBS surface.
/// The u-direction follows the profile and the v-direction follows the path.
pub fn generate_sweep_surface(
    profile: &AnalyzedProfile,
    frames: &FrameSequence,
    scale: &ScaleMethod,
    _tolerance: &Tolerance,
) -> Result<SweepSurfaceResult, SweepError> {
    let n_frames = frames.frames.len();
    if n_frames < 2 {
        return Err(SweepError::InvalidPath(
            "need at least 2 frames for sweep".into(),
        ));
    }

    if profile.spans.is_empty() {
        return Err(SweepError::InvalidProfile("no profile spans".into()));
    }

    let total_arc = frames.s_values.last().copied().unwrap_or(1.0).max(1e-30);

    // The wire contract (`SweepProfile::Wire`): profile `x` → frame
    // normal, `y` → binormal, `z` dropped. A span with a `z` extent is
    // a world-space profile the caller forgot to express in the path
    // frame — it would skin to a collapsed, zero-normal face, so it is
    // refused here with the measured extent (PK-416).
    let z_extent = profile
        .spans
        .iter()
        .flat_map(|s| s.control_points.iter().map(|p| p.z))
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), z| {
            (lo.min(z), hi.max(z))
        });
    if z_extent.1 - z_extent.0 > _tolerance.linear() {
        return Err(SweepError::InvalidProfile(format!(
            "wire profile must lie in the path frame's plane (x → normal, y → binormal, z \
             constant); its z extent is {:.3e} — express the profile in the sketch plane's \
             frame and place the swept body afterwards",
            z_extent.1 - z_extent.0
        )));
    }

    let mut surfaces = Vec::with_capacity(profile.spans.len());

    for span in &profile.spans {
        let surface = skin_span(span, frames, scale, total_arc)?;
        surfaces.push(surface);
    }

    let start_cap = transform_profile_at_frame(&profile.spans, &frames.frames[0], scale, 0.0);
    let end_cap =
        transform_profile_at_frame(&profile.spans, &frames.frames[n_frames - 1], scale, 1.0);

    Ok(SweepSurfaceResult {
        surfaces,
        start_cap,
        end_cap,
    })
}

/// Skin a single profile span into a NURBS surface using the frame sequence.
fn skin_span(
    span: &NurbsCurve,
    frames: &FrameSequence,
    scale: &ScaleMethod,
    total_arc: f64,
) -> Result<NurbsSurface, SweepError> {
    let num_u = span.control_points.len();
    let num_v = frames.frames.len();

    if num_u == 0 {
        return Err(SweepError::InvalidProfile(
            "span has no control points".into(),
        ));
    }

    // Build the control grid by SKINNING: each profile control point's
    // row of frame images is interpolated as a v-curve at the frames'
    // arc-length parameters (shared across rows, so the rows stay a
    // valid tensor product; each row's weight is constant along v, so
    // Euclidean row interpolation is exact for the rational profile).
    // Using the frame images directly as control points only
    // approximated the sweep — a clamped uniform B-spline traverses
    // its last control legs as a near-chord, and a tapered helical
    // sweep visibly cut the corner of its final quarter-turn (PK-440:
    // the spring gallery's "blade"; measured 95° of swing in the last
    // 2.5% of v, with the tube centre 0.01 off the helix).
    let degree_v = compute_v_degree(num_v);
    let v_params: Vec<f64> = frames.s_values.iter().map(|s| s / total_arc).collect();

    let mut control_points = Vec::with_capacity(num_u * num_v);
    let mut weights = Vec::with_capacity(num_u * num_v);
    let mut knots_v: Option<Vec<f64>> = None;

    for (cp_idx, cp) in span.control_points.iter().enumerate() {
        let row: Vec<Point3> = frames
            .frames
            .iter()
            .enumerate()
            .map(|(frame_idx, frame)| {
                let s_frac = frames.s_values[frame_idx] / total_arc;
                let s_factor = evaluate_scale(scale, s_frac);
                transform_point(cp, frame, s_factor)
            })
            .collect();
        let fitted = popkern_core::geometry::fitting::fit_curve(&row, &v_params, degree_v, num_v)
            .map_err(|e| {
            SweepError::InvalidPath(format!("sweep skin interpolation failed: {e}"))
        })?;
        debug_assert_eq!(fitted.control_points.len(), num_v);
        if let Some(knots) = &knots_v {
            debug_assert_eq!(knots, &fitted.knots);
        } else {
            knots_v = Some(fitted.knots.clone());
        }
        control_points.extend(fitted.control_points);
        weights.extend(std::iter::repeat_n(span.weights[cp_idx], num_v));
    }

    let degree_u = span.degree;
    let knots_u = span.knots.clone();
    let knots_v = knots_v.unwrap_or_else(|| clamped_uniform_knots(num_v, degree_v));

    Ok(NurbsSurface {
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        control_points,
        num_u,
        num_v,
        weights,
    })
}

/// Transform a profile control point by a frame with scaling.
///
/// Profile points are in local coordinates where x maps to normal
/// and y maps to binormal. The z component is ignored (profile is 2D).
fn transform_point(cp: &Point3, frame: &PathFrame, scale: f64) -> Point3 {
    let x = cp.x;
    let y = cp.y;
    let offset_n = frame.normal * (scale * x);
    let offset_b = frame.binormal * (scale * y);
    frame.position + offset_n + offset_b
}

/// Evaluate scale factor at a given arc-length fraction.
fn evaluate_scale(method: &ScaleMethod, s_frac: f64) -> f64 {
    match method {
        ScaleMethod::None => 1.0,
        ScaleMethod::Uniform(pairs) => interpolate_pairs(pairs, s_frac),
    }
}

/// Linearly interpolate a value from sorted (parameter, value) pairs.
fn interpolate_pairs(pairs: &[(f64, f64)], t: f64) -> f64 {
    if pairs.is_empty() {
        return 1.0;
    }
    if pairs.len() == 1 || t <= pairs[0].0 {
        return pairs[0].1;
    }
    if t >= pairs[pairs.len() - 1].0 {
        return pairs[pairs.len() - 1].1;
    }

    let idx = pairs.partition_point(|&(s, _)| s < t);
    if idx == 0 {
        return pairs[0].1;
    }

    let (s0, v0) = pairs[idx - 1];
    let (s1, v1) = pairs[idx];
    let denom = s1 - s0;
    if denom.abs() < 1e-30 {
        return v0;
    }
    let frac = (t - s0) / denom;
    v0 + frac * (v1 - v0)
}

/// Choose v-direction degree based on the number of frames.
fn compute_v_degree(num_v: usize) -> usize {
    if num_v <= 2 {
        1
    } else if num_v <= 4 {
        num_v - 1
    } else {
        3
    }
}

/// Generate clamped uniform knot vector for `n` control points with given degree.
fn clamped_uniform_knots(n: usize, degree: usize) -> Vec<f64> {
    let m = n + degree + 1; // total number of knots
    let mut knots = vec![0.0; degree + 1];

    // Interior knots
    let n_interior = m.saturating_sub(2 * (degree + 1));
    for i in 1..=n_interior {
        knots.push(i as f64 / (n_interior + 1) as f64);
    }

    // Last (degree + 1) knots = 1.0
    knots.extend(std::iter::repeat_n(1.0, degree + 1));

    knots
}

/// Transform profile spans at a single frame (for cap generation).
fn transform_profile_at_frame(
    spans: &[NurbsCurve],
    frame: &PathFrame,
    scale: &ScaleMethod,
    s_frac: f64,
) -> Vec<NurbsCurve> {
    let s_factor = evaluate_scale(scale, s_frac);
    spans
        .iter()
        .map(|span| {
            let control_points = span
                .control_points
                .iter()
                .map(|cp| transform_point(cp, frame, s_factor))
                .collect();
            NurbsCurve {
                degree: span.degree,
                knots: span.knots.clone(),
                control_points,
                weights: span.weights.clone(),
            }
        })
        .collect()
}

#[cfg(test)]
#[path = "sweep_surface_tests.rs"]
mod tests;
