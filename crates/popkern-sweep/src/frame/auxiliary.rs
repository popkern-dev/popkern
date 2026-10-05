//! Auxiliary-curve frame computation.

use crate::arc_length::ArcLengthTable;
use crate::error::SweepError;
use crate::frame::{FrameSequence, PathFrame};
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Vec3};
use popkern_math::PortableMath;

/// Compute frames using an auxiliary (spine) curve to define the up-direction.
///
/// At each sample, the vector from the path point to the corresponding auxiliary
/// curve point is projected onto the plane perpendicular to the tangent to form
/// the frame normal.
pub fn auxiliary_frames(
    curve: &NurbsCurve,
    aux_curve: &NurbsCurve,
    table: &ArcLengthTable,
    n_samples: usize,
    twist: f64,
) -> Result<FrameSequence, SweepError> {
    if n_samples == 0 {
        return Err(SweepError::InvalidPath("n_samples must be > 0".to_string()));
    }

    let total_length = table.total_length;
    let (aux_t_min, aux_t_max) = aux_curve.parameter_range();
    let aux_range = aux_t_max - aux_t_min;

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

        // Evaluate auxiliary curve at proportional arc-length position.
        let fraction = if total_length > 1e-15 {
            s / total_length
        } else {
            0.0
        };
        let aux_t = aux_t_min + fraction * aux_range;
        let aux_pos = aux_curve.point_at(aux_t);

        // Up direction = auxiliary point - path point.
        let up = aux_pos - pos;

        // Project onto plane perpendicular to tangent.
        let projected = up - tangent * up.dot(&tangent);
        let normal = projected
            .normalized()
            .ok_or_else(|| SweepError::FrameFailed {
                t,
                reason: "auxiliary curve coincides with path or is parallel to tangent".to_string(),
            })?;

        let binormal = tangent.cross(&normal);
        let binormal = binormal.normalized().unwrap_or(Vec3::zero());

        // Apply twist.
        let angle = twist * fraction;
        let (n_twisted, b_twisted) = if angle.abs() > 1e-15 {
            let cos_a = angle.pcos();
            let sin_a = angle.psin();
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

    fn x_axis_line() -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        }
    }

    fn offset_line_z() -> NurbsCurve {
        // Line offset 1 unit in Z: (0,0,1) to (1,0,1)
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 1.0), Point3::new(1.0, 0.0, 1.0)],
            weights: vec![1.0, 1.0],
        }
    }

    #[test]
    fn auxiliary_offset_curve_tracks_guide() {
        let curve = x_axis_line();
        let aux = offset_line_z();
        let table = ArcLengthTable::build(&curve, 100);
        let seq = auxiliary_frames(&curve, &aux, &table, 10, 0.0).unwrap();

        assert_eq!(seq.frames.len(), 11);

        for frame in &seq.frames {
            // The auxiliary is offset in Z, so normal should point in +Z.
            assert!(
                (frame.normal.z - 1.0).abs() < 1e-6,
                "normal should point +Z: {:?}",
                frame.normal
            );
        }
    }

    #[test]
    fn auxiliary_coincident_errors() {
        let curve = x_axis_line();
        // Auxiliary curve is the same as the path.
        let aux = x_axis_line();
        let table = ArcLengthTable::build(&curve, 100);
        let result = auxiliary_frames(&curve, &aux, &table, 10, 0.0);
        assert!(result.is_err());
    }

    #[test]
    fn zero_samples_error() {
        let curve = x_axis_line();
        let aux = offset_line_z();
        let table = ArcLengthTable::build(&curve, 100);
        let result = auxiliary_frames(&curve, &aux, &table, 0, 0.0);
        assert!(result.is_err());
    }
}
