//! 4x4 affine transformation matrix for 3D geometry.

use crate::geometry::{Point3, Vec3};
use popkern_math::PortableMath;
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
        let (s, c) = angle.psin_cos();
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
        let (s, c) = angle.psin_cos();
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
        let (s, c) = angle.psin_cos();
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

    /// The uniform scale factor of this transform, when it is a
    /// similarity (rigid motion × uniform scale, mirror allowed):
    /// `Some(s)` with `s > 0` iff the linear part `L` satisfies
    /// `LᵀL = s²·I` within a relative band of `1e-9`.
    ///
    /// Analytic surface fields (plane normals, cylinder/cone/torus
    /// axes, radii) only stay analytic under similarities — a
    /// non-uniform scale or shear turns a cylinder into an elliptic
    /// cylinder that `SurfaceGeometry` cannot represent, so
    /// [`crate::topology::BRep::transform`] uses this to decline
    /// such transforms with a typed error instead of silently
    /// corrupting radii and direction lengths (PK-231: a
    /// `scale_uniform(0.6)` notch cube shipped plane normals of
    /// length 0.6, displacing every downstream fillet offset by
    /// 0.6×).
    pub fn similarity_scale(&self) -> Option<f64> {
        let m = &self.matrix;
        // Columns of the linear part.
        let cols = [
            Vec3::new(m[0][0], m[1][0], m[2][0]),
            Vec3::new(m[0][1], m[1][1], m[2][1]),
            Vec3::new(m[0][2], m[1][2], m[2][2]),
        ];
        let s_sq = cols[0].dot(&cols[0]);
        if s_sq <= f64::EPSILON {
            return None;
        }
        let band = 1.0e-9 * s_sq;
        for i in 0..3 {
            for j in i..3 {
                let expected = if i == j { s_sq } else { 0.0 };
                if (cols[i].dot(&cols[j]) - expected).abs() > band {
                    return None;
                }
            }
        }
        // A projective bottom row would make apply_point non-affine
        // (and apply_direction ignores w, so points and directions
        // would disagree).
        if m[3][0].abs() > 1.0e-12
            || m[3][1].abs() > 1.0e-12
            || m[3][2].abs() > 1.0e-12
            || (m[3][3] - 1.0).abs() > 1.0e-12
        {
            return None;
        }
        let s = s_sq.sqrt();
        // Snap rigid transforms to EXACTLY 1.0: composed rotations
        // carry ~1e-16 column-norm noise, and `radius * 1.0000…02`
        // would perturb every radius by an ULP on every translate/
        // rotate — enough to flip marginal welds (the Bishop-frame
        // seed lesson). Bit-exactness at s = 1 is part of the
        // contract.
        if (s - 1.0).abs() <= 1.0e-12 {
            return Some(1.0);
        }
        Some(s)
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

    /// Create a rotation transform around an arbitrary axis through a point.
    ///
    /// Uses Rodrigues' rotation formula composed with translation to/from
    /// the axis origin.
    pub fn rotation_around_axis(
        origin: Point3,
        axis: Vec3,
        angle: f64,
    ) -> Result<Self, crate::error::CoreError> {
        let u = axis.normalized().ok_or_else(|| {
            crate::error::CoreError::PreconditionFailed(
                "axis is zero-length in rotation_around_axis".into(),
            )
        })?;
        let (ux, uy, uz) = (u.x, u.y, u.z);
        let (sin_a, cos_a) = angle.psin_cos();
        let c = 1.0 - cos_a;
        let r = [
            [
                cos_a + ux * ux * c,
                ux * uy * c - uz * sin_a,
                ux * uz * c + uy * sin_a,
            ],
            [
                uy * ux * c + uz * sin_a,
                cos_a + uy * uy * c,
                uy * uz * c - ux * sin_a,
            ],
            [
                uz * ux * c - uy * sin_a,
                uz * uy * c + ux * sin_a,
                cos_a + uz * uz * c,
            ],
        ];
        // Translation so that origin is a fixed point: t = origin - R * origin
        let tx = origin.x - (r[0][0] * origin.x + r[0][1] * origin.y + r[0][2] * origin.z);
        let ty = origin.y - (r[1][0] * origin.x + r[1][1] * origin.y + r[1][2] * origin.z);
        let tz = origin.z - (r[2][0] * origin.x + r[2][1] * origin.y + r[2][2] * origin.z);
        Ok(Self {
            matrix: [
                [r[0][0], r[0][1], r[0][2], tx],
                [r[1][0], r[1][1], r[1][2], ty],
                [r[2][0], r[2][1], r[2][2], tz],
                [0.0, 0.0, 0.0, 1.0],
            ],
        })
    }

    /// Create a reflection (mirror) transform across a plane defined by
    /// origin and normal.
    ///
    /// Uses the Householder reflection: `I - 2 * n * nᵀ`, composed with
    /// translation so the plane origin is a fixed point.
    pub fn reflection(origin: Point3, normal: Vec3) -> Result<Self, crate::error::CoreError> {
        let n = normal.normalized().ok_or_else(|| {
            crate::error::CoreError::PreconditionFailed(
                "normal is zero-length in reflection".into(),
            )
        })?;
        let (nx, ny, nz) = (n.x, n.y, n.z);
        let h = [
            [1.0 - 2.0 * nx * nx, -2.0 * nx * ny, -2.0 * nx * nz],
            [-2.0 * ny * nx, 1.0 - 2.0 * ny * ny, -2.0 * ny * nz],
            [-2.0 * nz * nx, -2.0 * nz * ny, 1.0 - 2.0 * nz * nz],
        ];
        // Translation so that origin is a fixed point: t = origin - H * origin
        let tx = origin.x - (h[0][0] * origin.x + h[0][1] * origin.y + h[0][2] * origin.z);
        let ty = origin.y - (h[1][0] * origin.x + h[1][1] * origin.y + h[1][2] * origin.z);
        let tz = origin.z - (h[2][0] * origin.x + h[2][1] * origin.y + h[2][2] * origin.z);
        Ok(Self {
            matrix: [
                [h[0][0], h[0][1], h[0][2], tx],
                [h[1][0], h[1][1], h[1][2], ty],
                [h[2][0], h[2][1], h[2][2], tz],
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
#[path = "transform_tests.rs"]
mod tests;
