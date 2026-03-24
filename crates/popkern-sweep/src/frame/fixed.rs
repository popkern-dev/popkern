//! Fixed-direction frame computation.

use crate::arc_length::ArcLengthTable;
use crate::error::SweepError;
use crate::frame::{FrameSequence, PathFrame};
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Vec3};

/// Compute frames using a fixed up-direction projected onto each tangent plane.
///
/// At each sample point, the fixed direction is projected onto the plane
/// perpendicular to the tangent. If the fixed direction is parallel to the
/// tangent, a `FrameFailed` error is returned.
pub fn fixed_frames(
    curve: &NurbsCurve,
    table: &ArcLengthTable,
    fixed_dir: Vec3,
    n_samples: usize,
    twist: f64,
) -> Result<FrameSequence, SweepError> {
    if n_samples == 0 {
        return Err(SweepError::InvalidPath("n_samples must be > 0".to_string()));
    }

    let fixed_normalized = fixed_dir
        .normalized()
        .ok_or_else(|| SweepError::FrameFailed {
            t: 0.0,
            reason: "fixed direction is zero vector".to_string(),
        })?;

    let total_length = table.total_length;
    let mut frames = Vec::with_capacity(n_samples + 1);
    let mut s_values = Vec::with_capacity(n_samples + 1);

    for i in 0..=n_samples {
        let s = total_length * i as f64 / n_samples as f64;
        let t = table.s_to_t(s);
        let pos = curve.point_at(t);
        let tangent = curve.tangent_at(t);

        if tangent.length() < 1e-12 {
            return Err(SweepError::FrameFailed {
                t,
                reason: "zero tangent".to_string(),
            });
        }

        // Project fixed direction onto plane perpendicular to tangent.
        let projected = fixed_normalized - tangent * fixed_normalized.dot(&tangent);
        let normal = projected
            .normalized()
            .ok_or_else(|| SweepError::FrameFailed {
                t,
                reason: "fixed direction parallel to path".to_string(),
            })?;

        let binormal = tangent.cross(&normal);
        let binormal = binormal.normalized().unwrap_or(Vec3::zero());

        // Apply twist.
        let fraction = if total_length > 1e-15 {
            s / total_length
        } else {
            0.0
        };
        let angle = twist * fraction;
        let (n_twisted, b_twisted) = if angle.abs() > 1e-15 {
            let cos_a = angle.cos();
            let sin_a = angle.sin();
            (
                normal * cos_a + binormal * sin_a,
                binormal * cos_a - normal * sin_a,
            )
        } else {
            (normal, binormal)
        };

        frames.push(PathFrame {
            position: pos,
            tangent,
            normal: n_twisted,
            binormal: b_twisted,
        });
        s_values.push(s);
    }

    Ok(FrameSequence { frames, s_values })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arc_length::ArcLengthTable;
    use popkern_core::geometry::Point3;

    fn xy_plane_line() -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        }
    }

    #[test]
    fn fixed_z_up_on_xy_plane() {
        let curve = xy_plane_line();
        let table = ArcLengthTable::build(&curve, 100);
        let seq = fixed_frames(&curve, &table, Vec3::new(0.0, 0.0, 1.0), 10, 0.0).unwrap();

        assert_eq!(seq.frames.len(), 11);

        for frame in &seq.frames {
            // Normal should be Z (projected fixed dir onto plane perp to X-tangent).
            assert!(
                (frame.normal.z - 1.0).abs() < 1e-6,
                "normal Z should be ~1: {:?}",
                frame.normal
            );
            // Normal should have no X component (tangent direction).
            assert!(
                frame.normal.x.abs() < 1e-6,
                "normal X should be ~0: {:?}",
                frame.normal
            );
        }
    }

    #[test]
    fn fixed_parallel_to_path_errors() {
        let curve = xy_plane_line();
        let table = ArcLengthTable::build(&curve, 100);
        // Fixed direction is along X, same as the tangent of the X-axis line.
        let result = fixed_frames(&curve, &table, Vec3::new(1.0, 0.0, 0.0), 10, 0.0);
        assert!(result.is_err());
    }

    #[test]
    fn fixed_zero_dir_errors() {
        let curve = xy_plane_line();
        let table = ArcLengthTable::build(&curve, 100);
        let result = fixed_frames(&curve, &table, Vec3::zero(), 10, 0.0);
        assert!(result.is_err());
    }

    #[test]
    fn zero_samples_error() {
        let curve = xy_plane_line();
        let table = ArcLengthTable::build(&curve, 100);
        let result = fixed_frames(&curve, &table, Vec3::new(0.0, 0.0, 1.0), 0, 0.0);
        assert!(result.is_err());
    }
}
