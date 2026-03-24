//! 4x4 affine transformation matrix for 3D geometry.

use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// A 4x4 affine transformation matrix (row-major storage).
///
/// Combines rotation, translation, scaling, and shear into a single
/// homogeneous coordinate matrix.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Transform3 {
    /// Row-major 4x4 matrix elements.
    pub matrix: [[f64; 4]; 4],
}

impl Transform3 {
    /// The identity transformation (no change).
    pub fn identity() -> Self {
        Self {
            matrix: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// Pure translation by (dx, dy, dz).
    pub fn translation(dx: f64, dy: f64, dz: f64) -> Self {
        Self {
            matrix: [
                [1.0, 0.0, 0.0, dx],
                [0.0, 1.0, 0.0, dy],
                [0.0, 0.0, 1.0, dz],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// Apply this transformation to a 3D point.
    ///
    /// Treats the point as `[x, y, z, 1]` in homogeneous coordinates,
    /// multiplies by the matrix, and divides by w.
    pub fn apply_point(&self, p: &Point3) -> Point3 {
        let m = &self.matrix;
        let x = m[0][0] * p.x + m[0][1] * p.y + m[0][2] * p.z + m[0][3];
        let y = m[1][0] * p.x + m[1][1] * p.y + m[1][2] * p.z + m[1][3];
        let z = m[2][0] * p.x + m[2][1] * p.y + m[2][2] * p.z + m[2][3];
        let w = m[3][0] * p.x + m[3][1] * p.y + m[3][2] * p.z + m[3][3];
        if (w - 1.0).abs() < 1e-15 {
            Point3::new(x, y, z)
        } else {
            Point3::new(x / w, y / w, z / w)
        }
    }

    /// Apply this transformation to a 3D direction vector.
    ///
    /// Uses `[x, y, z, 0]` — ignores translation component.
    pub fn apply_direction(&self, v: &Vec3) -> Vec3 {
        let m = &self.matrix;
        Vec3::new(
            m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z,
        )
    }

    /// Rotation about the X axis by `angle` radians.
    pub fn rotation_x(angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        Self {
            matrix: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, c, -s, 0.0],
                [0.0, s, c, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// Rotation about the Y axis by `angle` radians.
    pub fn rotation_y(angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        Self {
            matrix: [
                [c, 0.0, s, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [-s, 0.0, c, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// Rotation about the Z axis by `angle` radians.
    pub fn rotation_z(angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        Self {
            matrix: [
                [c, -s, 0.0, 0.0],
                [s, c, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// Non-uniform scaling.
    pub fn scale(sx: f64, sy: f64, sz: f64) -> Self {
        Self {
            matrix: [
                [sx, 0.0, 0.0, 0.0],
                [0.0, sy, 0.0, 0.0],
                [0.0, 0.0, sz, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// Uniform scaling.
    pub fn scale_uniform(s: f64) -> Self {
        Self::scale(s, s, s)
    }

    /// Compose two transforms: `self` applied after `other`.
    ///
    /// `self.compose(&other).apply_point(p)` == `self.apply_point(&other.apply_point(p))`
    pub fn compose(&self, other: &Transform3) -> Transform3 {
        let a = &self.matrix;
        let b = &other.matrix;
        let mut result = [[0.0f64; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                result[i][j] =
                    a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j] + a[i][3] * b[3][j];
            }
        }
        Transform3 { matrix: result }
    }

    /// Compute the inverse transformation.
    ///
    /// Returns `CoreError::PreconditionFailed` if the matrix is singular.
    #[allow(clippy::needless_range_loop)]
    pub fn inverse(&self) -> Result<Transform3, crate::error::CoreError> {
        let mut augmented = [[0.0f64; 8]; 4];
        for i in 0..4 {
            for j in 0..4 {
                augmented[i][j] = self.matrix[i][j];
                augmented[i][j + 4] = if i == j { 1.0 } else { 0.0 };
            }
        }
        for col in 0..4 {
            let mut max_val = augmented[col][col].abs();
            let mut max_row = col;
            for row in (col + 1)..4 {
                let val = augmented[row][col].abs();
                if val > max_val {
                    max_val = val;
                    max_row = row;
                }
            }
            if max_val < 1e-15 {
                return Err(crate::error::CoreError::PreconditionFailed(
                    "singular matrix cannot be inverted".into(),
                ));
            }
            augmented.swap(col, max_row);
            let pivot = augmented[col][col];
            for j in 0..8 {
                augmented[col][j] /= pivot;
            }
            for row in 0..4 {
                if row == col {
                    continue;
                }
                let factor = augmented[row][col];
                for j in 0..8 {
                    augmented[row][j] -= factor * augmented[col][j];
                }
            }
        }
        let mut result = [[0.0f64; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                result[i][j] = augmented[i][j + 4];
            }
        }
        Ok(Transform3 { matrix: result })
    }

    /// Build a transform from a STEP-style axis placement.
    ///
    /// Maps the standard basis to the given coordinate system:
    /// - `origin` becomes the translation
    /// - `z_axis` (normalized) becomes the local Z direction
    /// - `x_axis` is orthogonalized against `z_axis` for the local X
    /// - Y is computed as Z × X
    ///
    /// Returns an error if either axis is degenerate (zero-length).
    pub fn from_axis_placement(
        origin: Point3,
        z_axis: Vec3,
        x_axis: Vec3,
    ) -> Result<Self, crate::error::CoreError> {
        let z = z_axis.normalized().ok_or_else(|| {
            crate::error::CoreError::PreconditionFailed(
                "z_axis is zero-length in axis placement".into(),
            )
        })?;
        let x_proj = x_axis - z * x_axis.dot(&z);
        let x = x_proj.normalized().ok_or_else(|| {
            crate::error::CoreError::PreconditionFailed(
                "x_axis is degenerate (parallel to z_axis or zero-length)".into(),
            )
        })?;
        let y = z.cross(&x);
        Ok(Self {
            matrix: [
                [x.x, y.x, z.x, origin.x],
                [x.y, y.y, z.y, origin.y],
                [x.z, y.z, z.z, origin.z],
                [0.0, 0.0, 0.0, 1.0],
            ],
        })
    }

    /// Check if two transforms are approximately equal.
    pub fn approx_eq(&self, other: &Self, tolerance: f64) -> bool {
        for i in 0..4 {
            for j in 0..4 {
                if (self.matrix[i][j] - other.matrix[i][j]).abs() > tolerance {
                    return false;
                }
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn identity_transform_point_unchanged() {
        let t = Transform3::identity();
        let p = Point3::new(1.0, 2.0, 3.0);
        let q = t.apply_point(&p);
        assert!((q.x - 1.0).abs() < 1e-12);
        assert!((q.y - 2.0).abs() < 1e-12);
        assert!((q.z - 3.0).abs() < 1e-12);
    }

    #[test]
    fn translation_moves_point() {
        let t = Transform3::translation(10.0, 20.0, 30.0);
        let p = Point3::new(1.0, 2.0, 3.0);
        let q = t.apply_point(&p);
        assert!((q.x - 11.0).abs() < 1e-12);
        assert!((q.y - 22.0).abs() < 1e-12);
        assert!((q.z - 33.0).abs() < 1e-12);
    }

    #[test]
    fn translation_does_not_affect_direction() {
        let t = Transform3::translation(10.0, 20.0, 30.0);
        let v = Vec3::new(1.0, 0.0, 0.0);
        let u = t.apply_direction(&v);
        assert!((u.x - 1.0).abs() < 1e-12);
        assert!((u.y).abs() < 1e-12);
        assert!((u.z).abs() < 1e-12);
    }

    #[test]
    fn rotation_z_90_degrees() {
        let t = Transform3::rotation_z(FRAC_PI_2);
        let p = Point3::new(1.0, 0.0, 0.0);
        let q = t.apply_point(&p);
        assert!((q.x).abs() < 1e-12);
        assert!((q.y - 1.0).abs() < 1e-12);
        assert!((q.z).abs() < 1e-12);
    }

    #[test]
    fn rotation_x_90_degrees() {
        let t = Transform3::rotation_x(FRAC_PI_2);
        let p = Point3::new(0.0, 1.0, 0.0);
        let q = t.apply_point(&p);
        assert!((q.x).abs() < 1e-12);
        assert!((q.y).abs() < 1e-12);
        assert!((q.z - 1.0).abs() < 1e-12);
    }

    #[test]
    fn rotation_y_90_degrees() {
        let t = Transform3::rotation_y(FRAC_PI_2);
        let p = Point3::new(0.0, 0.0, 1.0);
        let q = t.apply_point(&p);
        assert!((q.x - 1.0).abs() < 1e-12);
        assert!((q.y).abs() < 1e-12);
        assert!((q.z).abs() < 1e-12);
    }

    #[test]
    fn uniform_scale() {
        let t = Transform3::scale_uniform(2.0);
        let p = Point3::new(1.0, 2.0, 3.0);
        let q = t.apply_point(&p);
        assert!((q.x - 2.0).abs() < 1e-12);
        assert!((q.y - 4.0).abs() < 1e-12);
        assert!((q.z - 6.0).abs() < 1e-12);
    }

    #[test]
    fn non_uniform_scale() {
        let t = Transform3::scale(2.0, 3.0, 4.0);
        let p = Point3::new(1.0, 1.0, 1.0);
        let q = t.apply_point(&p);
        assert!((q.x - 2.0).abs() < 1e-12);
        assert!((q.y - 3.0).abs() < 1e-12);
        assert!((q.z - 4.0).abs() < 1e-12);
    }

    #[test]
    fn compose_translation_then_rotation() {
        let t = Transform3::translation(1.0, 0.0, 0.0);
        let r = Transform3::rotation_z(FRAC_PI_2);
        let combined = r.compose(&t);
        let p = Point3::new(1.0, 0.0, 0.0);
        let q = combined.apply_point(&p);
        assert!((q.x).abs() < 1e-12);
        assert!((q.y - 2.0).abs() < 1e-12);
    }

    #[test]
    fn inverse_undoes_transform() {
        let t = Transform3::translation(3.0, 4.0, 5.0);
        let inv = t.inverse().expect("translation is invertible");
        let combined = inv.compose(&t);
        let p = Point3::new(7.0, 8.0, 9.0);
        let q = combined.apply_point(&p);
        assert!((q.x - 7.0).abs() < 1e-12);
        assert!((q.y - 8.0).abs() < 1e-12);
        assert!((q.z - 9.0).abs() < 1e-12);
    }

    #[test]
    fn inverse_of_rotation() {
        let t = Transform3::rotation_z(1.23);
        let inv = t.inverse().expect("rotation is invertible");
        let p = Point3::new(1.0, 2.0, 3.0);
        let q = inv.apply_point(&t.apply_point(&p));
        assert!((q.x - 1.0).abs() < 1e-10);
        assert!((q.y - 2.0).abs() < 1e-10);
        assert!((q.z - 3.0).abs() < 1e-10);
    }

    #[test]
    fn singular_matrix_has_no_inverse() {
        let mut t = Transform3::identity();
        t.matrix[0] = [0.0; 4];
        assert!(t.inverse().is_err());
    }

    #[test]
    fn from_axis_placement_identity_axes() {
        let t = Transform3::from_axis_placement(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
        )
        .expect("valid axes");
        let p = Point3::new(1.0, 2.0, 3.0);
        let q = t.apply_point(&p);
        assert!((q.x - 1.0).abs() < 1e-12);
        assert!((q.y - 2.0).abs() < 1e-12);
        assert!((q.z - 3.0).abs() < 1e-12);
    }

    #[test]
    fn from_axis_placement_with_offset() {
        let t = Transform3::from_axis_placement(
            Point3::new(10.0, 20.0, 30.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
        )
        .expect("valid axes");
        let p = Point3::new(0.0, 0.0, 0.0);
        let q = t.apply_point(&p);
        assert!((q.x - 10.0).abs() < 1e-12);
        assert!((q.y - 20.0).abs() < 1e-12);
        assert!((q.z - 30.0).abs() < 1e-12);
    }

    #[test]
    fn approx_eq_identical() {
        let t = Transform3::rotation_z(1.0);
        assert!(t.approx_eq(&t, 1e-12));
    }

    #[test]
    fn approx_eq_different() {
        let a = Transform3::translation(1.0, 0.0, 0.0);
        let b = Transform3::translation(2.0, 0.0, 0.0);
        assert!(!a.approx_eq(&b, 1e-6));
    }
}
