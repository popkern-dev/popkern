//! Frenet-Serret frame computation.

use crate::arc_length::ArcLengthTable;
use crate::error::{SweepError, SweepWarning};
use crate::frame::{FrameSequence, PathFrame};
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Vec3};

/// Compute Frenet-Serret frames along a curve.
///
/// The Frenet frame is defined by the tangent (T), normal (N = curvature direction),
/// and binormal (B = T x N). At points of zero curvature (straight segments,
/// inflection points), the Frenet frame is undefined. In such cases, the Bishop
/// frame is used as a fallback and a `FrenetFallback` warning is emitted.
pub fn frenet_frames(
    curve: &NurbsCurve,
    table: &ArcLengthTable,
    n_samples: usize,
    twist: f64,
) -> Result<(FrameSequence, Vec<SweepWarning>), SweepError> {
    if n_samples == 0 {
        return Err(SweepError::InvalidPath("n_samples must be > 0".to_string()));
    }

    let total_length = table.total_length;
    let mut frames = Vec::with_capacity(n_samples + 1);
    let mut s_values = Vec::with_capacity(n_samples + 1);
    let mut warnings = Vec::new();

    // Keep track of the last valid normal for fallback propagation.
    let mut last_normal: Option<Vec3> = None;
    let mut last_binormal: Option<Vec3> = None;

    for i in 0..=n_samples {
        let s = total_length * i as f64 / n_samples as f64;
        let t = table.s_to_t(s);
        let pos = curve.point_at(t);

        let d1 = curve.derivative_at(t, 1);
        let tangent = match d1.normalized() {
            Some(v) => v,
            None => {
                return Err(SweepError::FrameFailed {
                    t,
                    reason: "zero tangent".to_string(),
                });
            }
        };

        let d2 = curve.derivative_at(t, 2);
        // Remove tangential component: N_raw = d2 - (d2 . T) * T
        let n_raw = d2 - tangent * d2.dot(&tangent);

        let (normal, binormal, is_fallback) = if let Some(n) = n_raw.normalized() {
            let b = tangent.cross(&n);
            let b = b.normalized().unwrap_or(Vec3::zero());
            (n, b, false)
        } else {
            // Zero curvature: use last valid frame or pick arbitrary perpendicular.
            let (n, b) = match (last_normal, last_binormal) {
                (Some(ln), Some(lb)) => (ln, lb),
                _ => {
                    let n = pick_perpendicular_fallback(tangent);
                    let b = tangent.cross(&n);
                    let b = b.normalized().unwrap_or(Vec3::zero());
                    (n, b)
                }
            };
            (n, b, true)
        };

        if is_fallback {
            warnings.push(SweepWarning::FrenetFallback { t });
        }

        last_normal = Some(normal);
        last_binormal = Some(binormal);

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

    Ok((FrameSequence { frames, s_values }, warnings))
}

/// Pick a vector perpendicular to `tangent` for fallback.
fn pick_perpendicular_fallback(tangent: Vec3) -> Vec3 {
    let ax = tangent.x.abs();
    let ay = tangent.y.abs();
    let az = tangent.z.abs();

    let seed = if ax <= ay && ax <= az {
        Vec3::new(1.0, 0.0, 0.0)
    } else if ay <= az {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };

    let candidate = tangent.cross(&seed);
    candidate.normalized().unwrap_or(Vec3::new(0.0, 1.0, 0.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arc_length::ArcLengthTable;
    use popkern_core::geometry::Point3;

    fn unit_circle() -> NurbsCurve {
        use std::f64::consts::{FRAC_PI_2, PI, TAU};
        let w = 1.0_f64 / 2.0_f64.sqrt();
        NurbsCurve {
            degree: 2,
            knots: vec![
                0.0,
                0.0,
                0.0,
                FRAC_PI_2,
                FRAC_PI_2,
                PI,
                PI,
                PI + FRAC_PI_2,
                PI + FRAC_PI_2,
                TAU,
                TAU,
                TAU,
            ],
            control_points: vec![
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(-1.0, 1.0, 0.0),
                Point3::new(-1.0, 0.0, 0.0),
                Point3::new(-1.0, -1.0, 0.0),
                Point3::new(0.0, -1.0, 0.0),
                Point3::new(1.0, -1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
        }
    }

    fn straight_line() -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        }
    }

    #[test]
    fn circle_normal_points_to_center() {
        let curve = unit_circle();
        let table = ArcLengthTable::build(&curve, 400);
        let (seq, warnings) = frenet_frames(&curve, &table, 50, 0.0).unwrap();

        assert!(warnings.is_empty(), "circle should have no fallbacks");

        // For a unit circle centered at origin, the normal should point toward the center.
        // The Frenet normal for a circle points inward (toward center).
        for frame in &seq.frames[..seq.frames.len() - 1] {
            let to_center = Point3::new(0.0, 0.0, 0.0) - frame.position;
            let to_center_norm = to_center.normalized().unwrap_or(Vec3::zero());
            let dot = frame.normal.dot(&to_center_norm);
            // Normal should be roughly aligned with direction to center.
            assert!(
                dot.abs() > 0.5,
                "normal should point toward center: dot={dot}, normal={:?}, pos={:?}",
                frame.normal,
                frame.position
            );
        }
    }

    #[test]
    fn straight_line_all_fallbacks() {
        let curve = straight_line();
        let table = ArcLengthTable::build(&curve, 100);
        let (seq, warnings) = frenet_frames(&curve, &table, 10, 0.0).unwrap();

        assert_eq!(seq.frames.len(), 11);
        // A straight line has zero curvature everywhere, so all frames should be fallbacks.
        assert!(
            !warnings.is_empty(),
            "straight line should produce FrenetFallback warnings"
        );

        for w in &warnings {
            assert!(matches!(w, SweepWarning::FrenetFallback { .. }));
        }
    }

    #[test]
    fn zero_samples_error() {
        let curve = straight_line();
        let table = ArcLengthTable::build(&curve, 100);
        let result = frenet_frames(&curve, &table, 0, 0.0);
        assert!(result.is_err());
    }
}
