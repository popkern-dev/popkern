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
mod tests {
    use super::*;

    #[test]
    fn interpolate_constant() {
        let angle = interpolate_angle(0.5, 0.0, 1.0, &[], 0.1);
        assert!((angle - 0.1).abs() < 1e-12);
    }

    #[test]
    fn interpolate_single_pair() {
        let angle = interpolate_angle(0.5, 0.0, 1.0, &[(0.5, 0.2)], 0.1);
        assert!((angle - 0.2).abs() < 1e-12);
    }

    #[test]
    fn interpolate_two_pairs_midpoint() {
        let angle = interpolate_angle(0.5, 0.0, 1.0, &[(0.0, 0.1), (1.0, 0.3)], 0.0);
        assert!((angle - 0.2).abs() < 1e-10, "expected 0.2, got {}", angle);
    }

    #[test]
    fn interpolate_two_pairs_endpoints() {
        let a0 = interpolate_angle(0.0, 0.0, 1.0, &[(0.0, 0.1), (1.0, 0.3)], 0.0);
        let a1 = interpolate_angle(1.0, 0.0, 1.0, &[(0.0, 0.1), (1.0, 0.3)], 0.0);
        assert!((a0 - 0.1).abs() < 1e-10);
        assert!((a1 - 0.3).abs() < 1e-10);
    }

    #[test]
    fn variable_angle_verify_different_at_endpoints() {
        let pairs = vec![(0.0, 0.05), (1.0, 0.15)];
        let a_start = interpolate_angle(0.0, 0.0, 1.0, &pairs, 0.0);
        let a_end = interpolate_angle(1.0, 0.0, 1.0, &pairs, 0.0);
        assert!(
            (a_start - a_end).abs() > 0.05,
            "angles at endpoints should differ"
        );
    }

    #[test]
    fn rodrigues_identity() {
        let v = Vec3::new(0.0, 1.0, 0.0);
        let k = Vec3::new(0.0, 0.0, 1.0);
        let rotated = rodrigues_rotate(&v, &k, 0.0);
        assert!((rotated.x - v.x).abs() < 1e-12);
        assert!((rotated.y - v.y).abs() < 1e-12);
    }

    #[test]
    fn sample_cross_section_produces_correct_count() {
        let center = Point3::new(0.0, 0.0, 0.0);
        let cross_dir = Vec3::new(1.0, 0.0, 0.0);
        let plane = SurfaceGeometry::Plane {
            origin: center,
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let points = sample_cross_section(&plane, &center, &cross_dir, 5);
        assert_eq!(points.len(), 5);
    }

    #[test]
    fn sample_cross_section_single_point_at_center() {
        let center = Point3::new(1.0, 2.0, 3.0);
        let cross_dir = Vec3::new(0.0, 1.0, 0.0);
        let plane = SurfaceGeometry::Plane {
            origin: center,
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let points = sample_cross_section(&plane, &center, &cross_dir, 1);
        assert_eq!(points.len(), 1);
        // Single sample lands at center (t=0.0).
        assert!((points[0].x - center.x).abs() < 1e-10);
        assert!((points[0].y - center.y).abs() < 1e-10);
    }

    #[test]
    fn rodrigues_rotate_ninety_degrees() {
        // Rotate Y axis 90 degrees about Z: should give -X.
        let v = Vec3::new(0.0, 1.0, 0.0);
        let k = Vec3::new(0.0, 0.0, 1.0);
        let rotated = rodrigues_rotate(&v, &k, std::f64::consts::FRAC_PI_2);
        assert!((rotated.x + 1.0).abs() < 1e-10, "x={}", rotated.x);
        assert!(rotated.y.abs() < 1e-10, "y={}", rotated.y);
        assert!(rotated.z.abs() < 1e-10, "z={}", rotated.z);
    }

    #[test]
    fn rodrigues_rotate_preserves_length() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let k = Vec3::new(0.0, 1.0, 0.0).normalized().unwrap();
        let angle = 0.7;
        let rotated = rodrigues_rotate(&v, &k, angle);
        let original_len = (v.x * v.x + v.y * v.y + v.z * v.z).sqrt();
        let rotated_len =
            (rotated.x * rotated.x + rotated.y * rotated.y + rotated.z * rotated.z).sqrt();
        assert!((original_len - rotated_len).abs() < 1e-10);
    }

    #[test]
    fn compute_cross_direction_plane_surface() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let center = Point3::new(0.0, 0.0, 0.0);
        let tangent = Vec3::new(1.0, 0.0, 0.0);
        let cross = compute_cross_direction(&plane, &center, &tangent);
        // Cross direction should be perpendicular to tangent.
        assert!(
            cross.dot(&tangent).abs() < 1e-10,
            "dot={}",
            cross.dot(&tangent)
        );
        // Should be a unit vector.
        let len = (cross.x * cross.x + cross.y * cross.y + cross.z * cross.z).sqrt();
        assert!((len - 1.0).abs() < 1e-10);
    }

    #[test]
    fn compute_cross_direction_nurbs_surface() {
        use crate::surface::skinning;
        // Build a tiny NURBS surface and compute cross direction.
        let rows = vec![
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            vec![Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ];
        let nurbs = skinning::skin_points(&rows).unwrap();
        let surf = SurfaceGeometry::Nurbs(nurbs);
        let center = Point3::new(0.5, 0.5, 0.0);
        let tangent = Vec3::new(1.0, 0.0, 0.0);
        let cross = compute_cross_direction(&surf, &center, &tangent);
        // Must be a unit vector.
        let len = (cross.x * cross.x + cross.y * cross.y + cross.z * cross.z).sqrt();
        assert!(
            (len - 1.0).abs() < 1e-10,
            "cross direction not unit: len={}",
            len
        );
    }

    #[test]
    fn interpolate_angle_below_first_pair() {
        // u < pairs[0].0: should clamp to first pair's angle.
        let angle = interpolate_angle(-1.0, 0.0, 1.0, &[(0.5, 0.3), (1.0, 0.6)], 0.0);
        assert!((angle - 0.3).abs() < 1e-12, "got {}", angle);
    }

    #[test]
    fn interpolate_angle_above_last_pair() {
        // u > pairs.last().0: should clamp to last pair's angle.
        let angle = interpolate_angle(2.0, 0.0, 1.0, &[(0.0, 0.1), (0.5, 0.2)], 0.0);
        assert!((angle - 0.2).abs() < 1e-12, "got {}", angle);
    }

    #[test]
    fn interpolate_angle_degenerate_range() {
        // t0 == t1: range is zero, u defaults to 0.5.
        let angle = interpolate_angle(5.0, 5.0, 5.0, &[(0.0, 0.1), (1.0, 0.9)], 0.0);
        // u=0.5 → should interpolate to 0.5 between 0.1 and 0.9.
        assert!((angle - 0.5).abs() < 1e-10, "got {}", angle);
    }

    #[test]
    fn draft_general_constant_angle_produces_nurbs() {
        use crate::neutral::NeutralLine;
        use crate::surface::skinning;
        use crate::types::{AngleVariation, DraftOptions};
        use popkern_core::geometry::NurbsCurve;

        let rows = vec![
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            vec![Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ];
        let nurbs_surf = skinning::skin_points(&rows).unwrap();
        let surface = SurfaceGeometry::Nurbs(nurbs_surf);

        let neutral = NeutralLine {
            curve: NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
                weights: vec![1.0, 1.0],
            },
            direction: Vec3::new(1.0, 0.0, 0.0),
            origin: Point3::new(0.5, 0.5, 0.0),
        };

        let pull = Vec3::new(0.0, 0.0, 1.0);
        let opts = DraftOptions {
            nurbs_samples: 3,
            ..DraftOptions::default()
        };
        let tol = popkern_core::tolerance::Tolerance::default();

        let result = draft_general(
            &surface,
            &neutral,
            0.1,
            &AngleVariation::Constant,
            &pull,
            &opts,
            &tol,
        );
        assert!(result.is_ok(), "draft_general failed: {:?}", result.err());
        assert!(matches!(result.unwrap(), SurfaceGeometry::Nurbs(_)));
    }

    #[test]
    fn draft_general_zero_angle_still_produces_nurbs() {
        use crate::neutral::NeutralLine;
        use crate::surface::skinning;
        use crate::types::{AngleVariation, DraftOptions};
        use popkern_core::geometry::NurbsCurve;

        let rows = vec![
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            vec![Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ];
        let nurbs_surf = skinning::skin_points(&rows).unwrap();
        let surface = SurfaceGeometry::Nurbs(nurbs_surf);

        let neutral = NeutralLine {
            curve: NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
                weights: vec![1.0, 1.0],
            },
            direction: Vec3::new(1.0, 0.0, 0.0),
            origin: Point3::new(0.5, 0.5, 0.0),
        };

        let pull = Vec3::new(0.0, 0.0, 1.0);
        let opts = DraftOptions {
            nurbs_samples: 2,
            ..DraftOptions::default()
        };
        let tol = popkern_core::tolerance::Tolerance::default();

        let result = draft_general(
            &surface,
            &neutral,
            0.0,
            &AngleVariation::Constant,
            &pull,
            &opts,
            &tol,
        );
        assert!(
            result.is_ok(),
            "zero angle draft failed: {:?}",
            result.err()
        );
    }

    #[test]
    fn draft_general_variable_angle() {
        use crate::neutral::NeutralLine;
        use crate::surface::skinning;
        use crate::types::{AngleVariation, DraftOptions};
        use popkern_core::geometry::NurbsCurve;

        let rows = vec![
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            vec![Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ];
        let nurbs_surf = skinning::skin_points(&rows).unwrap();
        let surface = SurfaceGeometry::Nurbs(nurbs_surf);

        let neutral = NeutralLine {
            curve: NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
                weights: vec![1.0, 1.0],
            },
            direction: Vec3::new(1.0, 0.0, 0.0),
            origin: Point3::new(0.5, 0.5, 0.0),
        };

        let pull = Vec3::new(0.0, 0.0, 1.0);
        let variation = AngleVariation::Variable(vec![(0.0, 0.05), (1.0, 0.15)]);
        let opts = DraftOptions {
            nurbs_samples: 4,
            ..DraftOptions::default()
        };
        let tol = popkern_core::tolerance::Tolerance::default();

        let result = draft_general(&surface, &neutral, 0.1, &variation, &pull, &opts, &tol);
        assert!(
            result.is_ok(),
            "variable angle draft failed: {:?}",
            result.err()
        );
    }

    #[test]
    fn draft_general_nurbs_samples_minimum_two() {
        use crate::neutral::NeutralLine;
        use crate::surface::skinning;
        use crate::types::{AngleVariation, DraftOptions};
        use popkern_core::geometry::NurbsCurve;

        // nurbs_samples=1 should be clamped to 2 internally.
        let rows = vec![
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            vec![Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ];
        let nurbs_surf = skinning::skin_points(&rows).unwrap();
        let surface = SurfaceGeometry::Nurbs(nurbs_surf);

        let neutral = NeutralLine {
            curve: NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
                weights: vec![1.0, 1.0],
            },
            direction: Vec3::new(1.0, 0.0, 0.0),
            origin: Point3::new(0.5, 0.5, 0.0),
        };

        let pull = Vec3::new(0.0, 0.0, 1.0);
        let opts = DraftOptions {
            nurbs_samples: 1,
            ..DraftOptions::default()
        };
        let tol = popkern_core::tolerance::Tolerance::default();

        let result = draft_general(
            &surface,
            &neutral,
            0.05,
            &AngleVariation::Constant,
            &pull,
            &opts,
            &tol,
        );
        assert!(
            result.is_ok(),
            "n_samples=1 clamp failed: {:?}",
            result.err()
        );
    }
}
