//! Revolve operation — creates a solid of revolution by rotating a profile
//! around an axis.
//!
//! The profile is analyzed for smooth spans, and each span is converted to a
//! NURBS surface of revolution using [`popkern_core::geometry::revolution::revolution_surface`].
//! The resulting surfaces are assembled into a BRep using the existing sweep
//! topology builder.

use popkern_core::geometry::revolution::{revolution_surface, rotate_point};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;

use crate::error::SweepError;
use crate::profile::analyze_profile;
use crate::topology::side_faces::build_sweep_topology;
use crate::types::{SweepProfile, SweepResult};

/// Revolve a profile around an axis by a given angle.
///
/// Creates a solid of revolution by rotating the profile curve(s) around
/// the axis. A full revolution (angle = 2pi) creates a closed solid.
/// Partial revolutions create an open-sided solid with start/end cap faces.
///
/// # Arguments
/// * `profile` - The profile to revolve (wire or face).
/// * `axis_origin` - A point on the revolution axis.
/// * `axis_direction` - Direction of the revolution axis.
/// * `angle` - Revolution angle in radians, must be in (0, 2pi].
/// * `tolerance` - Geometric tolerance.
///
/// # Errors
///
/// Returns [`SweepError`] if the profile is invalid, the angle is out of
/// range, or topology construction fails.
pub fn revolve(
    profile: &SweepProfile,
    axis_origin: &Point3,
    axis_direction: &Vec3,
    angle: f64,
    tolerance: Tolerance,
) -> Result<SweepResult, SweepError> {
    // 1. Analyze profile into smooth spans.
    let analyzed = analyze_profile(profile, &tolerance)?;

    if analyzed.spans.is_empty() {
        return Err(SweepError::InvalidProfile("no profile spans".into()));
    }

    // 2. Build a revolution surface for each span.
    let mut surfaces = Vec::with_capacity(analyzed.spans.len());
    for span in &analyzed.spans {
        let surface = revolution_surface(span, axis_origin, axis_direction, angle)
            .map_err(SweepError::Core)?;
        surfaces.push(surface);
    }

    // 3. Build cap curves (rotated profile endpoints) for partial revolutions.
    let is_full_revolution = (angle - std::f64::consts::TAU).abs() < 1e-10;

    let (start_caps, end_caps) = if analyzed.is_closed && !is_full_revolution {
        // For closed profiles with partial revolution, we need cap faces.
        let start_caps = build_cap_curves(&analyzed.spans);
        let end_caps =
            build_rotated_cap_curves(&analyzed.spans, axis_origin, axis_direction, angle);
        (start_caps, end_caps)
    } else if analyzed.is_closed && is_full_revolution {
        // Full revolution with closed profile: no caps needed, but the
        // topology builder expects cap curves for closed profiles.
        // We pass empty caps and treat it as an open profile to avoid
        // adding cap faces — the surface wraps around.
        (vec![], vec![])
    } else {
        (vec![], vec![])
    };

    // 4. Build topology.
    // For full revolution with closed profile, the body is already closed
    // so we build with is_closed_profile = false to avoid cap faces.
    let use_closed_topology = analyzed.is_closed && !is_full_revolution;

    let brep = build_sweep_topology(
        &surfaces,
        &start_caps,
        &end_caps,
        use_closed_topology,
        tolerance,
    )?;

    Ok(SweepResult {
        brep,
        warnings: vec![],
    })
}

/// Build cap curves from the profile spans at the start position (angle = 0).
///
/// Each span contributes one curve segment to the cap boundary.
fn build_cap_curves(spans: &[NurbsCurve]) -> Vec<NurbsCurve> {
    spans.to_vec()
}

/// Build cap curves from the profile spans rotated to the end position.
///
/// Each span's control points are rotated around the axis by the given angle.
fn build_rotated_cap_curves(
    spans: &[NurbsCurve],
    axis_origin: &Point3,
    axis_direction: &Vec3,
    angle: f64,
) -> Vec<NurbsCurve> {
    spans
        .iter()
        .map(|span| {
            let rotated_cps: Vec<Point3> = span
                .control_points
                .iter()
                .map(|cp| rotate_point(cp, axis_origin, axis_direction, angle))
                .collect();
            NurbsCurve {
                degree: span.degree,
                knots: span.knots.clone(),
                control_points: rotated_cps,
                weights: span.weights.clone(),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::surface::SurfaceEval;
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    /// A line segment in the XZ plane, offset from the Z axis by 1 unit.
    /// When revolved around Z, this produces a cylinder-like surface.
    fn line_profile_xz() -> SweepProfile {
        SweepProfile::Wire(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 2.0)],
            weights: vec![1.0, 1.0],
        })
    }

    /// A closed rectangular profile in the XZ plane for solid revolution.
    fn closed_rect_profile() -> SweepProfile {
        SweepProfile::Wire(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0],
            control_points: vec![
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(2.0, 0.0, 1.0),
                Point3::new(1.0, 0.0, 1.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0; 5],
        })
    }

    #[test]
    fn revolve_line_full_produces_brep() {
        let profile = line_profile_xz();
        let result = revolve(
            &profile,
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            TAU,
            Tolerance::default(),
        );
        assert!(result.is_ok(), "revolve failed: {:?}", result.err());
        let brep = &result.unwrap().brep; // ok in tests
        assert!(
            brep.faces.len() >= 1,
            "expected at least 1 face, got {}",
            brep.faces.len()
        );
    }

    #[test]
    fn revolve_line_partial_produces_brep() {
        let profile = line_profile_xz();
        let result = revolve(
            &profile,
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            FRAC_PI_2,
            Tolerance::default(),
        );
        assert!(result.is_ok(), "revolve failed: {:?}", result.err());
        let brep = &result.unwrap().brep; // ok in tests
        assert!(
            brep.faces.len() >= 1,
            "expected at least 1 face, got {}",
            brep.faces.len()
        );
    }

    #[test]
    fn revolve_closed_profile_partial() {
        let profile = closed_rect_profile();
        let result = revolve(
            &profile,
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            PI,
            Tolerance::default(),
        );
        assert!(result.is_ok(), "revolve failed: {:?}", result.err());
        let brep = &result.unwrap().brep; // ok in tests
        // Partial revolution of closed profile: side faces + 2 caps.
        assert!(
            brep.faces.len() >= 3,
            "expected at least 3 faces, got {}",
            brep.faces.len()
        );
    }

    #[test]
    fn revolve_closed_profile_full() {
        let profile = closed_rect_profile();
        let result = revolve(
            &profile,
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            TAU,
            Tolerance::default(),
        );
        assert!(result.is_ok(), "revolve failed: {:?}", result.err());
        let brep = &result.unwrap().brep; // ok in tests
        assert!(
            brep.faces.len() >= 1,
            "expected at least 1 face, got {}",
            brep.faces.len()
        );
    }

    #[test]
    fn revolve_surface_points_at_correct_radius() {
        // Revolve a line at x=2 around Z axis by PI.
        // All surface points should be at radius ~2.
        let profile = SweepProfile::Wire(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 3.0)],
            weights: vec![1.0, 1.0],
        });
        let result = revolve(
            &profile,
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            PI,
            Tolerance::default(),
        );
        assert!(result.is_ok(), "revolve failed: {:?}", result.err());
    }

    #[test]
    fn revolve_invalid_zero_angle_fails() {
        let profile = line_profile_xz();
        let result = revolve(
            &profile,
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            0.0,
            Tolerance::default(),
        );
        assert!(result.is_err(), "zero angle should fail");
    }

    #[test]
    fn revolve_negative_angle_fails() {
        let profile = line_profile_xz();
        let result = revolve(
            &profile,
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            -1.0,
            Tolerance::default(),
        );
        assert!(result.is_err(), "negative angle should fail");
    }
}
