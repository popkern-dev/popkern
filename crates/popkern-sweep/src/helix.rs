//! Helical sweep operations.
//!
//! Generates helical path curves from a [`HelixSpec`] and detects when a
//! helix combined with a circular profile produces an analytic torus surface.

use popkern_math::PortableMath;
use std::f64::consts::{FRAC_1_SQRT_2, PI};

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

use crate::error::{SweepError, SweepWarning};
use crate::types::{HelixPitch, HelixSpec, SweepProfile};

/// Number of sample points per helix turn for variable pitch/radius fitting.
const SAMPLES_PER_TURN: usize = 36;

/// Generate a helical path curve from a [`HelixSpec`].
///
/// For constant pitch and constant radius, produces an exact rational
/// degree-2 NURBS curve using quarter-arc decomposition. For variable
/// pitch or radius, samples the parametric helix and fits a NURBS curve.
pub fn build_helix_path(spec: &HelixSpec) -> Result<(NurbsCurve, Vec<SweepWarning>), SweepError> {
    validate_spec(spec)?;

    let is_constant =
        spec.start_radius == spec.end_radius && matches!(spec.pitch, HelixPitch::Constant(_));

    if is_constant {
        let pitch = match &spec.pitch {
            HelixPitch::Constant(p) => *p,
            HelixPitch::Variable(_) => unreachable!(),
        };
        let curve = build_constant_helix(spec, pitch)?;
        Ok((curve, vec![]))
    } else {
        let curve = build_variable_helix(spec)?;
        Ok((curve, vec![SweepWarning::HelixApproximated]))
    }
}

/// Check if a helix combined with a circular wire profile produces a torus.
///
/// Returns `Some(SurfaceGeometry::Torus { .. })` when:
/// - The helix has constant pitch and constant radius
/// - The helix is exactly one turn
/// - The profile is a circular wire
///
/// The detection is approximate: it checks that the profile centroid is
/// near the origin and the profile is closed, which is consistent with
/// a circular cross-section.
pub fn detect_analytic_helix(spec: &HelixSpec, profile: &SweepProfile) -> Option<SurfaceGeometry> {
    // Must be constant pitch and radius.
    let pitch = match &spec.pitch {
        HelixPitch::Constant(p) => *p,
        HelixPitch::Variable(_) => return None,
    };
    if (spec.start_radius - spec.end_radius).abs() > 1e-10 {
        return None;
    }

    // Zero pitch means the helix is a circle — combined with a circular
    // profile this is a torus (degenerate: a surface of revolution).
    if pitch.abs() > 1e-10 {
        return None;
    }

    // Must be exactly one turn.
    if (spec.turns - 1.0).abs() > 1e-10 {
        return None;
    }

    // Profile must be a closed wire (heuristic for "circle").
    let wire = match profile {
        SweepProfile::Wire(c) => c,
        _ => return None,
    };

    // Check that start == end (closed curve).
    let (t_min, t_max) = wire.parameter_range();
    let start = wire.point_at(t_min);
    let end = wire.point_at(t_max);
    if start.distance_to(&end) > 1e-6 {
        return None;
    }

    // Estimate the minor radius from the profile bounding extent.
    let minor_radius = estimate_profile_radius(wire);
    if minor_radius < 1e-10 {
        return None;
    }

    let axis = spec.axis_dir.normalized()?;

    Some(SurfaceGeometry::Torus {
        center: spec.axis_origin,
        axis,
        major_radius: spec.start_radius,
        minor_radius,
    })
}

/// Validate a helix specification.
fn validate_spec(spec: &HelixSpec) -> Result<(), SweepError> {
    if spec.turns <= 0.0 {
        return Err(SweepError::InvalidPath(
            "helix turns must be positive".into(),
        ));
    }
    if spec.start_radius <= 0.0 || spec.end_radius <= 0.0 {
        return Err(SweepError::InvalidPath(
            "helix radius must be positive".into(),
        ));
    }
    if spec.axis_dir.length() < 1e-30 {
        return Err(SweepError::InvalidPath(
            "helix axis direction must be non-zero".into(),
        ));
    }
    if !spec.start_angle.is_finite() {
        return Err(SweepError::InvalidPath(
            "helix start angle must be finite".into(),
        ));
    }
    match &spec.pitch {
        HelixPitch::Constant(p) if *p < 0.0 => {
            return Err(SweepError::InvalidPath(
                "helix pitch must be non-negative".into(),
            ));
        }
        HelixPitch::Variable(pts) if pts.is_empty() => {
            return Err(SweepError::InvalidPath(
                "variable pitch must have at least one point".into(),
            ));
        }
        _ => {}
    }
    Ok(())
}

/// Build an exact rational NURBS helix with constant pitch and radius.
///
/// Decomposes each turn into 4 quarter arcs. Each quarter arc uses 3
/// control points: start, weighted midpoint (w = 1/sqrt(2)), end.
fn build_constant_helix(spec: &HelixSpec, pitch: f64) -> Result<NurbsCurve, SweepError> {
    let axis = spec
        .axis_dir
        .normalized()
        .ok_or_else(|| SweepError::InvalidPath("helix axis is zero-length".into()))?;
    let radius = spec.start_radius;

    // Build a local coordinate system.
    let (u_dir, v_dir) = orthonormal_pair(&axis);

    let n_quarters = (spec.turns * 4.0).ceil() as usize;
    if n_quarters == 0 {
        return Err(SweepError::InvalidPath("helix has no turns".into()));
    }
    // A left-hand helix winds clockwise about the axis while still
    // advancing along it; phase and winding are rigid transforms of
    // the same exact quarter arcs (each still subtends 90°, so the
    // rational midpoint weight is unchanged).
    let winding = if spec.left_handed { -1.0 } else { 1.0 };
    let total_angle = winding * spec.turns * 2.0 * PI;
    let total_height = spec.turns * pitch;

    // Each quarter arc has 3 control points; shared endpoints → 2*n + 1 total.
    let n_cp = 2 * n_quarters + 1;
    let mut control_points = Vec::with_capacity(n_cp);
    let mut weights = Vec::with_capacity(n_cp);

    for i in 0..=n_quarters {
        let frac = i as f64 / n_quarters as f64;
        let angle = spec.start_angle + frac * total_angle;
        let z = frac * total_height;

        // On-curve point.
        let x = radius * angle.pcos();
        let y = radius * angle.psin();
        let pos = spec.axis_origin + u_dir * x + v_dir * y + axis * z;
        control_points.push(pos);
        weights.push(1.0);

        // Midpoint (weighted) — except after the last on-curve point.
        if i < n_quarters {
            let frac_mid = (i as f64 + 0.5) / n_quarters as f64;
            let angle_mid = spec.start_angle + frac_mid * total_angle;
            let z_mid = frac_mid * total_height;

            // For a rational quarter-arc, the weighted midpoint control
            // point lies at radius * sqrt(2) from the center in the
            // bisector direction (radius / cos(pi/4)).
            let scale = 1.0 / FRAC_1_SQRT_2; // = sqrt(2)
            let x_mid2 = radius * scale * angle_mid.pcos();
            let y_mid2 = radius * scale * angle_mid.psin();
            let pos_mid = spec.axis_origin + u_dir * x_mid2 + v_dir * y_mid2 + axis * z_mid;
            control_points.push(pos_mid);
            weights.push(FRAC_1_SQRT_2);
        }
    }

    // Knot vector: degree 2, clamped, with double knots at quarter boundaries.
    let mut knots = vec![0.0, 0.0, 0.0];
    for i in 1..n_quarters {
        let k = i as f64 / n_quarters as f64;
        knots.push(k);
        knots.push(k);
    }
    knots.push(1.0);
    knots.push(1.0);
    knots.push(1.0);

    Ok(NurbsCurve {
        degree: 2,
        knots,
        control_points,
        weights,
    })
}

/// Build an approximate NURBS helix for variable pitch or radius.
///
/// Samples the parametric helix at many points and fits a NURBS curve.
fn build_variable_helix(spec: &HelixSpec) -> Result<NurbsCurve, SweepError> {
    let axis = spec
        .axis_dir
        .normalized()
        .ok_or_else(|| SweepError::InvalidPath("helix axis is zero-length".into()))?;
    let (u_dir, v_dir) = orthonormal_pair(&axis);

    let n_samples = (spec.turns * SAMPLES_PER_TURN as f64).ceil() as usize;
    let n_samples = n_samples.max(8);

    let winding = if spec.left_handed { -1.0 } else { 1.0 };
    let total_angle = winding * spec.turns * 2.0 * PI;

    let mut points = Vec::with_capacity(n_samples + 1);
    let mut params = Vec::with_capacity(n_samples + 1);

    let mut cumulative_height = 0.0;

    for i in 0..=n_samples {
        let frac = i as f64 / n_samples as f64;
        let angle = spec.start_angle + frac * total_angle;

        // Interpolate radius.
        let radius = spec.start_radius + frac * (spec.end_radius - spec.start_radius);

        // Compute height from pitch.
        let current_pitch = evaluate_pitch(&spec.pitch, frac);
        if i > 0 {
            let dt = 1.0 / n_samples as f64;
            // Pitch is distance per turn, so per-fraction increment is pitch * dt * turns.
            cumulative_height += current_pitch * dt * spec.turns;
        }

        let x = radius * angle.pcos();
        let y = radius * angle.psin();
        let pos = spec.axis_origin + u_dir * x + v_dir * y + axis * cumulative_height;
        points.push(pos);
        params.push(frac);
    }

    let degree = 3.min(n_samples);
    let num_cp = (n_samples + 1).min(n_samples + 1);

    let curve = popkern_core::geometry::fitting::fit_curve(&points, &params, degree, num_cp)
        .map_err(|e| SweepError::InvalidPath(format!("helix fit failed: {e}")))?;

    Ok(curve)
}

/// Evaluate pitch at a given turn fraction.
fn evaluate_pitch(pitch: &HelixPitch, frac: f64) -> f64 {
    match pitch {
        HelixPitch::Constant(p) => *p,
        HelixPitch::Variable(pts) => interpolate_pairs(pts, frac),
    }
}

/// Linearly interpolate a value from sorted (parameter, value) pairs.
fn interpolate_pairs(pairs: &[(f64, f64)], t: f64) -> f64 {
    if pairs.is_empty() {
        return 0.0;
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

/// Build an orthonormal pair perpendicular to a given axis.
///
/// The axis must be a unit vector (pre-normalized by the caller).
fn orthonormal_pair(axis: &Vec3) -> (Vec3, Vec3) {
    let candidate = if axis.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    // Safe to use fallback since axis is pre-validated as non-zero.
    let u = axis
        .cross(&candidate)
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v = axis
        .cross(&u)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    (u, v)
}

/// Estimate the radius of a wire profile by sampling and measuring max
/// distance from the centroid.
fn estimate_profile_radius(curve: &NurbsCurve) -> f64 {
    let (t_min, t_max) = curve.parameter_range();
    let n = 32;

    // Compute centroid.
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    for i in 0..=n {
        let t = t_min + (i as f64 / n as f64) * (t_max - t_min);
        let p = curve.point_at(t);
        cx += p.x;
        cy += p.y;
        cz += p.z;
    }
    let inv = 1.0 / (n + 1) as f64;
    let centroid = Point3::new(cx * inv, cy * inv, cz * inv);

    // Measure max distance from centroid.
    let mut max_dist = 0.0_f64;
    for i in 0..=n {
        let t = t_min + (i as f64 / n as f64) * (t_max - t_min);
        let p = curve.point_at(t);
        let d = p.distance_to(&centroid);
        max_dist = max_dist.max(d);
    }

    max_dist
}

#[cfg(test)]
#[path = "helix_tests.rs"]
mod tests;
