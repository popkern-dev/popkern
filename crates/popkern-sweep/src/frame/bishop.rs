//! Bishop (rotation-minimizing) frame computation.

use crate::arc_length::ArcLengthTable;
use crate::error::SweepError;
use crate::frame::{FrameSequence, PathFrame};
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Vec3};

/// Compute Bishop rotation-minimizing frames along a curve.
///
/// The Bishop frame avoids the torsion-induced spin of the Frenet frame and is
/// well-defined even at inflection points (zero curvature). The initial normal
/// is chosen perpendicular to the tangent, then parallel-transported along the
/// curve using Rodrigues' rotation formula.
///
/// For closed paths, holonomy correction is applied so that the frame matches
/// at the start and end.
pub fn bishop_frames(
    curve: &NurbsCurve,
    table: &ArcLengthTable,
    n_samples: usize,
    twist: f64,
    is_closed: bool,
) -> Result<FrameSequence, SweepError> {
    if n_samples == 0 {
        return Err(SweepError::InvalidPath("n_samples must be > 0".to_string()));
    }

    let total_length = table.total_length;
    let mut s_values = Vec::with_capacity(n_samples + 1);
    let mut tangents = Vec::with_capacity(n_samples + 1);
    let mut positions = Vec::with_capacity(n_samples + 1);

    // Sample positions and tangents at uniform arc-length intervals.
    for i in 0..=n_samples {
        let s = total_length * i as f64 / n_samples as f64;
        let t = table.s_to_t(s);
        let pos = curve.point_at(t);
        let tan = curve.tangent_at(t);

        if tan.length() < 1e-12 {
            return Err(SweepError::FrameFailed {
                t,
                reason: "zero tangent".to_string(),
            });
        }

        s_values.push(s);
        positions.push(pos);
        tangents.push(tan);
    }

    // Choose initial normal perpendicular to T0.
    let t0 = tangents[0];
    let initial_normal = pick_perpendicular(t0)?;
    let initial_binormal = t0.cross(&initial_normal);

    // Parallel transport using Rodrigues' rotation.
    let mut normals = Vec::with_capacity(n_samples + 1);
    let mut binormals = Vec::with_capacity(n_samples + 1);
    normals.push(initial_normal);
    binormals.push(initial_binormal);

    for i in 1..=n_samples {
        let t_prev = tangents[i - 1];
        let t_curr = tangents[i];
        let n_prev = normals[i - 1];
        let b_prev = binormals[i - 1];

        let (n_new, b_new) = rotate_frame(t_prev, t_curr, n_prev, b_prev);
        normals.push(n_new);
        binormals.push(b_new);
    }

    // For closed paths, correct holonomy.
    let holonomy_correction = if is_closed {
        compute_holonomy_correction(&normals, &tangents)
    } else {
        0.0
    };

    // Apply twist and holonomy correction.
    let mut frames = Vec::with_capacity(n_samples + 1);
    for i in 0..=n_samples {
        let fraction = if total_length > 1e-15 {
            s_values[i] / total_length
        } else {
            0.0
        };
        let angle = twist * fraction + holonomy_correction * fraction;
        let (n_twisted, b_twisted) =
            rotate_around_axis(normals[i], binormals[i], tangents[i], angle);

        frames.push(PathFrame {
            position: positions[i],
            tangent: tangents[i],
            normal: n_twisted,
            binormal: b_twisted,
        });
    }

    Ok(FrameSequence { frames, s_values })
}

/// Pick a vector perpendicular to `tangent`.
///
/// Uses the coordinate axis with the smallest component of `tangent`
/// (numerical stability), ORTHOGONALIZED against the tangent — not
/// `tangent × seed`. The distinction is user-visible: the initial
/// normal is where the profile's local +x lands, and the cross-
/// product form rotated every profile 90° about the path (a
/// rect(2, 1) extruded along +z came out spanning ±1 in Y — masked
/// for years by circles and squares). With the orthogonalized seed,
/// an extrusion along +z maps profile x→world x, y→world y.
fn pick_perpendicular(tangent: Vec3) -> Result<Vec3, SweepError> {
    let mut ax = tangent.x.abs();
    let mut ay = tangent.y.abs();
    let mut az = tangent.z.abs();

    // Snap noise components to zero before choosing the seed axis.
    // A numerically-differentiated tangent of an axis-aligned path
    // carries ULP-level trash in the other components (a lerp of two
    // identical coordinates is not always bit-exact), and comparing
    // two ~1e-15 values picks a data-dependent axis: the whole
    // profile frame flipped 90° when a sweep path's y-coordinate
    // changed in the 5th decimal (measured: a 22-gon deck built
    // rotated about its axis purely from its centroid value).
    let m = ax.max(ay).max(az);
    if m > 0.0 {
        let floor = m * 1.0e-9;
        if ax < floor {
            ax = 0.0;
        }
        if ay < floor {
            ay = 0.0;
        }
        if az < floor {
            az = 0.0;
        }
    }

    let seed = if ax <= ay && ax <= az {
        Vec3::new(1.0, 0.0, 0.0)
    } else if ay <= az {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };

    let candidate = seed - tangent * seed.dot(&tangent);
    candidate
        .normalized()
        .ok_or_else(|| SweepError::FrameFailed {
            t: 0.0,
            reason: "could not find perpendicular to tangent".to_string(),
        })
}

/// Rotate (normal, binormal) from `t_prev` tangent direction to `t_curr` using
/// Rodrigues' rotation formula.
fn rotate_frame(t_prev: Vec3, t_curr: Vec3, normal: Vec3, binormal: Vec3) -> (Vec3, Vec3) {
    let axis = t_prev.cross(&t_curr);
    let axis_len = axis.length();

    if axis_len < 1e-14 {
        // Tangents are (nearly) parallel; no rotation needed.
        return (normal, binormal);
    }

    let axis_norm = axis * (1.0 / axis_len);
    let cos_angle = t_prev.dot(&t_curr).clamp(-1.0, 1.0);
    let sin_angle = axis_len; // |a x b| = |a||b|sin(theta), both unit

    let n_new = rodrigues(normal, axis_norm, cos_angle, sin_angle);
    let b_new = rodrigues(binormal, axis_norm, cos_angle, sin_angle);

    // Re-normalize to combat drift.
    let n_out = n_new.normalized().unwrap_or(normal);
    let b_out = b_new.normalized().unwrap_or(binormal);

    (n_out, b_out)
}

/// Rodrigues' rotation: rotate `v` around `axis` by angle with given cos/sin.
fn rodrigues(v: Vec3, axis: Vec3, cos_a: f64, sin_a: f64) -> Vec3 {
    let dot = axis.dot(&v);
    let cross = axis.cross(&v);
    v * cos_a + cross * sin_a + axis * (dot * (1.0 - cos_a))
}

/// Rotate normal and binormal around the tangent axis by the given angle.
fn rotate_around_axis(normal: Vec3, binormal: Vec3, tangent: Vec3, angle: f64) -> (Vec3, Vec3) {
    if angle.abs() < 1e-15 {
        return (normal, binormal);
    }
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    let n = normal * cos_a + binormal * sin_a;
    let b = binormal * cos_a - normal * sin_a;
    let n = n.normalized().unwrap_or(normal);
    let b = b.normalized().unwrap_or(binormal);
    // Suppress unused variable warning for tangent (kept for API clarity).
    let _ = tangent;
    (n, b)
}

/// Compute the angle needed to correct holonomy for closed paths.
///
/// This is the angle between the initial normal and the transported normal at
/// the end, measured in the tangent plane at the start.
fn compute_holonomy_correction(normals: &[Vec3], tangents: &[Vec3]) -> f64 {
    let n_start = normals[0];
    let n_end = normals[normals.len() - 1];
    let t_start = tangents[0];

    // Project n_end onto the start tangent plane.
    let n_end_proj = n_end - t_start * n_end.dot(&t_start);
    let n_end_proj = match n_end_proj.normalized() {
        Some(v) => v,
        None => return 0.0,
    };

    // Signed angle from n_end_proj to n_start in the tangent plane.
    let cos_h = n_start.dot(&n_end_proj).clamp(-1.0, 1.0);
    let sin_h = n_start.cross(&n_end_proj).dot(&t_start);
    let holonomy = sin_h.atan2(cos_h);

    // We want to distribute -holonomy so that the final frame matches the start.
    -holonomy
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arc_length::ArcLengthTable;
    use popkern_core::geometry::Point3;

    /// Create a straight line NURBS from (0,0,0) to (1,0,0).
    fn straight_line() -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        }
    }

    /// Create a quarter-circle NURBS in the XY plane.
    fn quarter_circle() -> NurbsCurve {
        use std::f64::consts::FRAC_1_SQRT_2;
        NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            weights: vec![1.0, FRAC_1_SQRT_2, 1.0],
        }
    }

    /// Create a full circle NURBS in the XY plane.
    fn full_circle() -> NurbsCurve {
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

    #[test]
    fn straight_line_constant_frame() {
        let curve = straight_line();
        let table = ArcLengthTable::build(&curve, 100);
        let seq = bishop_frames(&curve, &table, 10, 0.0, false).unwrap();

        assert_eq!(seq.frames.len(), 11);

        // All normals should be the same (no rotation on a straight line).
        let n0 = seq.frames[0].normal;
        for frame in &seq.frames {
            let diff = (frame.normal - n0).length();
            assert!(diff < 1e-6, "normal changed on straight line: diff={diff}");
        }
    }

    #[test]
    fn quarter_circle_rotates_90() {
        let curve = quarter_circle();
        let table = ArcLengthTable::build(&curve, 200);
        let seq = bishop_frames(&curve, &table, 50, 0.0, false).unwrap();

        // Tangent at start should be roughly +Y, at end roughly -X.
        let t_start = seq.frames[0].tangent;
        let t_end = seq.frames.last().unwrap().tangent;

        // Start tangent should point roughly in +Y direction.
        assert!(
            t_start.y > 0.5,
            "start tangent should point +Y: {t_start:?}"
        );
        // End tangent should point roughly in -X direction.
        assert!(t_end.x < -0.5, "end tangent should point -X: {t_end:?}");

        // Angle between start and end tangent should be ~90 degrees.
        let angle = t_start.angle_to(&t_end);
        assert!(
            (angle - std::f64::consts::FRAC_PI_2).abs() < 0.1,
            "angle should be ~90 degrees: {angle}"
        );
    }

    #[test]
    fn closed_circle_holonomy_corrected() {
        let curve = full_circle();
        let table = ArcLengthTable::build(&curve, 400);
        let seq = bishop_frames(&curve, &table, 100, 0.0, true).unwrap();

        // First and last frame should match (closed path with holonomy correction).
        let first = &seq.frames[0];
        let last = seq.frames.last().unwrap();

        let n_diff = (first.normal - last.normal).length();
        let b_diff = (first.binormal - last.binormal).length();

        assert!(
            n_diff < 0.15,
            "normals should match on closed path: diff={n_diff}"
        );
        assert!(
            b_diff < 0.15,
            "binormals should match on closed path: diff={b_diff}"
        );
    }

    #[test]
    fn zero_samples_error() {
        let curve = straight_line();
        let table = ArcLengthTable::build(&curve, 100);
        let result = bishop_frames(&curve, &table, 0, 0.0, false);
        assert!(result.is_err());
    }
}
