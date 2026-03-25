//! Planar face draft via Rodrigues rotation.
//!
//! Rotates the face normal about the neutral line direction by the draft angle
//! using the Rodrigues rotation formula, producing a new tilted plane.

use crate::error::DraftError;
use crate::neutral::NeutralLine;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceGeometry;

/// Draft a planar face by rotating its normal about the neutral line.
///
/// Uses the Rodrigues rotation formula:
/// ```text
/// v_rot = v * cos(a) + (k x v) * sin(a) + k * (k . v) * (1 - cos(a))
/// ```
/// where `k` is the neutral line direction and `a` is the draft angle.
///
/// # Errors
///
/// Returns [`DraftError::InvalidFace`] if the surface is not a plane.
/// Returns [`DraftError::InvalidNeutral`] if the neutral direction is degenerate.
pub fn draft_planar(
    surface: &SurfaceGeometry,
    neutral: &NeutralLine,
    angle: f64,
) -> Result<SurfaceGeometry, DraftError> {
    let (_origin, normal, _u_axis, _v_axis) = match surface {
        SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        } => (*origin, *normal, *u_axis, *v_axis),
        _ => {
            return Err(DraftError::InvalidFace(
                "draft_planar called on non-planar surface".into(),
            ));
        }
    };

    let axis = neutral
        .direction
        .normalized()
        .ok_or_else(|| DraftError::InvalidNeutral("neutral direction is zero".into()))?;

    // Rodrigues rotation of the face normal.
    let rotated_normal = rodrigues_rotate(&normal, &axis, angle);

    let new_normal = rotated_normal
        .normalized()
        .ok_or_else(|| DraftError::InvalidNeutral("rotated normal is degenerate".into()))?;

    // Compute new u_axis and v_axis via Gram-Schmidt orthogonalization.
    let (new_u, new_v) = gram_schmidt_basis(&new_normal);

    Ok(SurfaceGeometry::Plane {
        origin: neutral.origin,
        normal: new_normal,
        u_axis: new_u,
        v_axis: new_v,
    })
}

/// Rodrigues rotation: rotate vector `v` about unit axis `k` by angle `a` (radians).
fn rodrigues_rotate(v: &Vec3, k: &Vec3, a: f64) -> Vec3 {
    let cos_a = a.cos();
    let sin_a = a.sin();
    let k_cross_v = k.cross(v);
    let k_dot_v = k.dot(v);
    *v * cos_a + k_cross_v * sin_a + *k * k_dot_v * (1.0 - cos_a)
}

/// Compute an orthonormal basis (u, v) given a normal vector using Gram-Schmidt.
fn gram_schmidt_basis(normal: &Vec3) -> (Vec3, Vec3) {
    // Pick a vector not parallel to normal.
    let candidate = if normal.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };

    // u = candidate - (candidate . normal) * normal, then normalize.
    let dot = candidate.dot(normal);
    let u_raw = candidate - *normal * dot;
    let u = u_raw.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));

    // v = normal x u
    let v = normal.cross(&u);
    let v = v.normalized().unwrap_or(Vec3::new(0.0, 1.0, 0.0));

    (u, v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;

    fn make_xy_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    fn make_neutral(direction: Vec3) -> NeutralLine {
        NeutralLine {
            curve: popkern_core::geometry::NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
                weights: vec![1.0, 1.0],
            },
            direction,
            origin: Point3::new(0.0, 0.0, 0.0),
        }
    }

    #[test]
    fn five_degree_draft_angle() {
        let plane = make_xy_plane();
        let neutral = make_neutral(Vec3::new(1.0, 0.0, 0.0));
        let angle = 5.0_f64.to_radians();

        let result = draft_planar(&plane, &neutral, angle).unwrap();
        if let SurfaceGeometry::Plane { normal, .. } = &result {
            let original_normal = Vec3::new(0.0, 0.0, 1.0);
            let angle_between = original_normal.angle_to(normal);
            assert!(
                (angle_between - angle).abs() < 0.01,
                "angle between normals = {} rad, expected {} rad",
                angle_between,
                angle
            );
        } else {
            panic!("expected Plane");
        }
    }

    #[test]
    fn zero_degree_draft_same_plane() {
        let plane = make_xy_plane();
        let neutral = make_neutral(Vec3::new(1.0, 0.0, 0.0));

        let result = draft_planar(&plane, &neutral, 0.0).unwrap();
        if let SurfaceGeometry::Plane { normal, .. } = &result {
            let original_normal = Vec3::new(0.0, 0.0, 1.0);
            let angle_between = original_normal.angle_to(normal);
            assert!(
                angle_between < 0.001,
                "zero draft should produce same normal, got angle = {}",
                angle_between
            );
        } else {
            panic!("expected Plane");
        }
    }

    #[test]
    fn forty_five_degree_draft() {
        let plane = make_xy_plane();
        let neutral = make_neutral(Vec3::new(1.0, 0.0, 0.0));
        let angle = 45.0_f64.to_radians();

        let result = draft_planar(&plane, &neutral, angle).unwrap();
        if let SurfaceGeometry::Plane { normal, .. } = &result {
            let original_normal = Vec3::new(0.0, 0.0, 1.0);
            let angle_between = original_normal.angle_to(normal);
            assert!(
                (angle_between - angle).abs() < 0.01,
                "angle between normals = {} rad, expected {} rad",
                angle_between,
                angle
            );
        } else {
            panic!("expected Plane");
        }
    }

    #[test]
    fn positive_vs_negative_angle() {
        let plane = make_xy_plane();
        let neutral = make_neutral(Vec3::new(1.0, 0.0, 0.0));

        let pos = draft_planar(&plane, &neutral, 0.1).unwrap();
        let neg = draft_planar(&plane, &neutral, -0.1).unwrap();

        if let (
            SurfaceGeometry::Plane { normal: n_pos, .. },
            SurfaceGeometry::Plane { normal: n_neg, .. },
        ) = (&pos, &neg)
        {
            // Positive and negative should tilt in opposite directions.
            // The y-components should have opposite signs (rotated about X axis).
            assert!(
                n_pos.y * n_neg.y <= 0.0 || (n_pos.y.abs() < 1e-10 && n_neg.y.abs() < 1e-10),
                "positive and negative drafts should tilt opposite: pos={:?}, neg={:?}",
                n_pos,
                n_neg
            );
        } else {
            panic!("expected Planes");
        }
    }

    #[test]
    fn non_planar_surface_returns_error() {
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let neutral = make_neutral(Vec3::new(1.0, 0.0, 0.0));
        let result = draft_planar(&cyl, &neutral, 0.1);
        assert!(result.is_err());
    }

    #[test]
    fn rodrigues_rotate_identity() {
        let v = Vec3::new(0.0, 0.0, 1.0);
        let k = Vec3::new(1.0, 0.0, 0.0);
        let rotated = rodrigues_rotate(&v, &k, 0.0);
        assert!((rotated.x - v.x).abs() < 1e-12);
        assert!((rotated.y - v.y).abs() < 1e-12);
        assert!((rotated.z - v.z).abs() < 1e-12);
    }

    #[test]
    fn gram_schmidt_produces_orthonormal() {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let (u, v) = gram_schmidt_basis(&normal);
        assert!((u.length() - 1.0).abs() < 1e-10);
        assert!((v.length() - 1.0).abs() < 1e-10);
        assert!(u.dot(&normal).abs() < 1e-10);
        assert!(v.dot(&normal).abs() < 1e-10);
        assert!(u.dot(&v).abs() < 1e-10);
    }
}
