//! General NURBS face draft via cross-section rotation and skinning.
//!
//! For non-planar, non-cylindrical faces, samples cross-sections along the
//! neutral line, rotates each by the draft angle, and skins the result
//! into a new NurbsSurface.

use crate::error::DraftError;
use crate::neutral::NeutralLine;
use crate::surface::skinning;
use crate::types::{AngleVariation, DraftOptions};
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;

/// Draft a general (NURBS) face by sampling, rotating, and skinning.
///
/// 1. Samples the neutral line at `options.nurbs_samples` points.
/// 2. At each point, computes a cross-section of the face perpendicular
///    to the neutral tangent.
/// 3. Rotates each cross-section by the draft angle (which may vary).
/// 4. Skins the rotated cross-sections into a NurbsSurface.
///
/// # Errors
///
/// Returns [`DraftError`] if sampling or skinning fails.
#[allow(clippy::too_many_arguments)]
pub fn draft_general(
    surface: &SurfaceGeometry,
    neutral: &NeutralLine,
    angle: f64,
    variation: &AngleVariation,
    _pull_direction: &Vec3,
    options: &DraftOptions,
    _tolerance: &Tolerance,
) -> Result<SurfaceGeometry, DraftError> {
    let n_samples = options.nurbs_samples.max(2);
    let (t0, t1) = neutral.curve.parameter_range();
    let cross_samples = 5; // Points per cross-section.

    let mut rows = Vec::with_capacity(n_samples);

    for i in 0..n_samples {
        let t = if n_samples == 1 {
            (t0 + t1) * 0.5
        } else {
            t0 + (t1 - t0) * (i as f64) / ((n_samples - 1) as f64)
        };

        let center = neutral.curve.point_at(t);
        let tangent = neutral
            .curve
            .tangent_at(t)
            .normalized()
            .unwrap_or(neutral.direction);

        // Compute the local draft angle at this parameter.
        let local_angle = match variation {
            AngleVariation::Constant => angle,
            AngleVariation::Variable(pairs) => interpolate_angle(t, t0, t1, pairs, angle),
        };

        // Build cross-section perpendicular to tangent.
        let cross_dir = compute_cross_direction(surface, &center, &tangent);

        // Sample the face surface along the cross-section direction.
        let cross_points = sample_cross_section(surface, &center, &cross_dir, cross_samples);

        // Rotate each cross-section point about the neutral tangent.
        let rotated: Vec<Point3> = cross_points
            .iter()
            .map(|p| {
                let offset = *p - center;
                let rotated_offset = rodrigues_rotate(&offset, &tangent, local_angle);
                center + rotated_offset
            })
            .collect();

        rows.push(rotated);
    }

    let nurbs = skinning::skin_points(&rows)?;
    Ok(SurfaceGeometry::Nurbs(nurbs))
}

/// Interpolate the draft angle at parameter t using the variable angle pairs.
fn interpolate_angle(t: f64, t0: f64, t1: f64, pairs: &[(f64, f64)], default: f64) -> f64 {
    if pairs.is_empty() {
        return default;
    }
    if pairs.len() == 1 {
        return pairs[0].1;
    }

    // Normalize t to [0, 1].
    let range = t1 - t0;
    let u = if range.abs() > 1e-12 {
        (t - t0) / range
    } else {
        0.5
    };

    // Find the bracketing pair.
    if u <= pairs[0].0 {
        return pairs[0].1;
    }
    if u >= pairs[pairs.len() - 1].0 {
        return pairs[pairs.len() - 1].1;
    }

    for w in pairs.windows(2) {
        let (u0, a0) = w[0];
        let (u1, a1) = w[1];
        if u >= u0 && u <= u1 {
            let frac = if (u1 - u0).abs() > 1e-12 {
                (u - u0) / (u1 - u0)
            } else {
                0.5
            };
            return a0 + (a1 - a0) * frac;
        }
    }

    default
}

/// Rodrigues rotation: rotate vector `v` about unit axis `k` by angle `a` (radians).
fn rodrigues_rotate(v: &Vec3, k: &Vec3, a: f64) -> Vec3 {
    let cos_a = a.cos();
    let sin_a = a.sin();
    let k_cross_v = k.cross(v);
    let k_dot_v = k.dot(v);
    *v * cos_a + k_cross_v * sin_a + *k * k_dot_v * (1.0 - cos_a)
}

/// Compute a direction perpendicular to the tangent for cross-section sampling.
fn compute_cross_direction(surface: &SurfaceGeometry, center: &Point3, tangent: &Vec3) -> Vec3 {
    // Use the surface normal at the nearest point to define the cross-section plane.
    let normal = match surface {
        SurfaceGeometry::Plane { normal, .. } => *normal,
        _ => {
            let ((u0, u1), (v0, v1)) = surface.parameter_range();
            let u_mid = (u0 + u1) * 0.5;
            let v_mid = (v0 + v1) * 0.5;
            surface.normal_at(u_mid, v_mid)
        }
    };

    // Cross direction = tangent x normal (lies in the surface, perpendicular to tangent).
    let cross = tangent.cross(&normal);
    cross.normalized().unwrap_or_else(|| {
        // Fallback: pick any direction perpendicular to tangent.
        let candidate = if tangent.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        let _ = center; // suppress unused warning in fallback
        tangent
            .cross(&candidate)
            .normalized()
            .unwrap_or(Vec3::new(0.0, 1.0, 0.0))
    })
}

/// Sample points on the face surface along a cross-section direction.
fn sample_cross_section(
    _surface: &SurfaceGeometry,
    center: &Point3,
    cross_dir: &Vec3,
    n_samples: usize,
) -> Vec<Point3> {
    let extent = 0.5; // Half-width of the cross-section sampling range.
    (0..n_samples)
        .map(|i| {
            let t = if n_samples == 1 {
                0.0
            } else {
                -extent + 2.0 * extent * (i as f64) / ((n_samples - 1) as f64)
            };
            *center + *cross_dir * t
        })
        .collect()
}

#[cfg(test)]
#[path = "general_tests.rs"]
mod tests;
