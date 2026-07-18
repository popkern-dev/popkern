//! Cylindrical face draft: cylinder to cone transformation.
//!
//! When a cylindrical face is drafted, the cylinder becomes a cone whose
//! half-angle equals the draft angle. The apex is positioned along the
//! cylinder axis to keep the neutral line fixed.

use crate::error::DraftError;
use crate::neutral::NeutralLine;
use popkern_core::geometry::surface::SurfaceGeometry;
/// Draft a cylindrical face by converting it to a cone.
///
/// The cone's half-angle equals the draft angle, and the apex is positioned
/// along the cylinder axis at a distance `radius / tan(angle)` from the
/// neutral line origin, toward the shrinking end.
///
/// # Errors
///
/// Returns [`DraftError::InvalidFace`] if the surface is not a cylinder.
/// Returns [`DraftError::AngleTooLarge`] if the angle is too close to zero
/// (producing a degenerate cone apex at infinity).
pub fn draft_cylindrical(
    surface: &SurfaceGeometry,
    neutral: &NeutralLine,
    angle: f64,
) -> Result<SurfaceGeometry, DraftError> {
    let (origin, axis, radius) = match surface {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => (*origin, *axis, *radius),
        _ => {
            return Err(DraftError::InvalidFace(
                "draft_cylindrical called on non-cylindrical surface".into(),
            ));
        }
    };

    let abs_angle = angle.abs();
    if abs_angle < 1e-12 {
        // Zero draft angle: cylinder stays as-is.
        return Ok(surface.clone());
    }

    let tan_a = abs_angle.tan();
    if tan_a < 1e-12 {
        // Extremely small angle, effectively zero.
        return Ok(surface.clone());
    }

    let axis_norm = axis
        .normalized()
        .ok_or_else(|| DraftError::InvalidFace("cylinder axis is zero".into()))?;

    // Apex distance from the neutral station along the axis.
    let apex_distance = radius / tan_a;

    // Determine direction toward the shrinking end.
    // If angle > 0, the cone narrows in the positive axis direction.
    let apex_dir = if angle > 0.0 { axis_norm } else { -axis_norm };

    // The apex sits ON THE CYLINDER AXIS at the neutral line's axial
    // station — the cone must contain the neutral circle unchanged.
    // neutral.origin is a point on the WALL (the neutral line lies in
    // the face), so project it onto the axis first; offsetting from
    // the wall point directly produces an off-axis apex whose cone
    // misses the whole neutral circle (measured: apex x = −r).
    let axis_pt = origin + axis_norm * (neutral.origin - origin).dot(&axis_norm);
    let apex = axis_pt + apex_dir * apex_distance;

    // The cone axis points away from the apex into the body.
    let cone_axis = if angle > 0.0 { -axis_norm } else { axis_norm };

    Ok(SurfaceGeometry::Cone {
        apex,
        axis: cone_axis,
        half_angle: abs_angle,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::{Point3, Vec3};

    fn make_cylinder() -> SurfaceGeometry {
        SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        }
    }

    fn make_neutral() -> NeutralLine {
        NeutralLine {
            curve: popkern_core::geometry::NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
                weights: vec![1.0, 1.0],
            },
            direction: Vec3::new(1.0, 0.0, 0.0),
            origin: Point3::new(0.0, 0.0, 0.0),
        }
    }

    #[test]
    fn cylinder_to_cone_correct_half_angle() {
        let cyl = make_cylinder();
        let neutral = make_neutral();
        let angle = 10.0_f64.to_radians();

        let result = draft_cylindrical(&cyl, &neutral, angle).unwrap();
        if let SurfaceGeometry::Cone { half_angle, .. } = result {
            assert!(
                (half_angle - angle).abs() < 1e-10,
                "half_angle = {}, expected {}",
                half_angle,
                angle
            );
        } else {
            panic!("expected Cone, got {:?}", result);
        }
    }

    #[test]
    fn zero_angle_returns_cylinder() {
        let cyl = make_cylinder();
        let neutral = make_neutral();

        let result = draft_cylindrical(&cyl, &neutral, 0.0).unwrap();
        assert!(
            matches!(result, SurfaceGeometry::Cylinder { .. }),
            "zero angle should return cylinder"
        );
    }

    #[test]
    fn negative_angle_produces_cone() {
        let cyl = make_cylinder();
        let neutral = make_neutral();
        let angle = -10.0_f64.to_radians();

        let result = draft_cylindrical(&cyl, &neutral, angle).unwrap();
        if let SurfaceGeometry::Cone { half_angle, .. } = result {
            assert!(
                (half_angle - angle.abs()).abs() < 1e-10,
                "half_angle should be positive"
            );
        } else {
            panic!("expected Cone");
        }
    }

    #[test]
    fn non_cylinder_returns_error() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let neutral = make_neutral();
        let result = draft_cylindrical(&plane, &neutral, 0.1);
        assert!(result.is_err());
    }

    #[test]
    fn apex_position_correct() {
        let cyl = make_cylinder();
        let neutral = make_neutral();
        let angle = std::f64::consts::FRAC_PI_4; // 45 degrees

        let result = draft_cylindrical(&cyl, &neutral, angle).unwrap();
        if let SurfaceGeometry::Cone { apex, .. } = result {
            // radius=1, tan(45)=1, so apex should be 1.0 units along +Z from origin.
            assert!(
                (apex.z - 1.0).abs() < 1e-10,
                "apex z = {}, expected 1.0",
                apex.z
            );
        } else {
            panic!("expected Cone");
        }
    }
}
