//! Rigid placements: where an instance sits in its parent.
//!
//! A [`Placement`] is a rotation plus a translation, validated at
//! every entry: orthonormal columns to [`RIGID_TOLERANCE`] and a
//! positive determinant. Scale, shear and mirror are refused, since a
//! mirror inverts a solid and scale changes what the part *is*, both
//! definition-level operations. A matrix within [`RIGID_TOLERANCE`]
//! but not exactly orthonormal (round-off from a STEP file, a long
//! chain of compositions) is re-orthonormalised from its z and x
//! columns; an exactly rigid one is kept bit for bit.

use popkern_core::geometry::{Point3, Transform3, Vec3};
use serde::{Deserialize, Serialize};

use crate::error::AssemblyError;

/// Largest `|RᵀR − I|` entry accepted as a rigid rotation.
pub const RIGID_TOLERANCE: f64 = 1e-9;

/// A rigid motion: rotation then translation, as a row-major 4×4.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "[[f64; 4]; 4]", into = "[[f64; 4]; 4]")]
pub struct Placement(Transform3);

impl Placement {
    /// No motion.
    pub fn identity() -> Placement {
        Placement(Transform3::identity())
    }

    /// A pure translation.
    pub fn translation(x: f64, y: f64, z: f64) -> Placement {
        Placement(Transform3::translation(x, y, z))
    }

    /// A rotation by `degrees` about the axis through `origin` along
    /// `axis` (right-handed).
    ///
    /// # Errors
    ///
    /// [`AssemblyError::NonRigidPlacement`] for a zero-length axis or
    /// a non-finite angle.
    pub fn rotation(origin: Point3, axis: Vec3, degrees: f64) -> Result<Placement, AssemblyError> {
        if !degrees.is_finite() {
            return Err(non_rigid("placement", f64::INFINITY, false));
        }
        Transform3::rotation_around_axis(origin, axis, degrees.to_radians())
            .map(Placement)
            .map_err(|_| non_rigid("placement", f64::INFINITY, false))
    }

    /// The frame with origin `origin`, local z along `z_dir` and local
    /// x along `x_dir` (orthogonalised against z): the STEP
    /// `AXIS2_PLACEMENT_3D` form.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::NonRigidPlacement`] when the axes are zero or
    /// parallel.
    pub fn frame(origin: Point3, z_dir: Vec3, x_dir: Vec3) -> Result<Placement, AssemblyError> {
        Transform3::from_axis_placement(origin, z_dir, x_dir)
            .map(Placement)
            .map_err(|_| non_rigid("placement", f64::INFINITY, false))
    }

    /// Validate an arbitrary matrix as a placement of `instance`.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::NonRigidPlacement`] when the matrix is not a
    /// rigid motion to [`RIGID_TOLERANCE`].
    pub fn from_matrix(matrix: Transform3, instance: &str) -> Result<Placement, AssemblyError> {
        let m = &matrix.matrix;
        let finite = m.iter().flatten().all(|v| v.is_finite());
        let affine = m[3] == [0.0, 0.0, 0.0, 1.0];
        let defect = if finite {
            orthonormality_defect(&matrix)
        } else {
            f64::INFINITY
        };
        let det = determinant(&matrix);
        if !finite || !affine || defect > RIGID_TOLERANCE || det <= 0.0 {
            return Err(non_rigid(instance, defect, det <= 0.0));
        }
        if defect == 0.0 {
            return Ok(Placement(matrix));
        }
        let placement = Placement(matrix);
        Placement::frame(placement.origin(), placement.z_axis(), placement.x_axis())
            .map_err(|_| non_rigid(instance, defect, false))
    }

    /// Wrap a matrix the caller built from rigid factors (no check);
    /// anything leaving the crate goes through [`Placement::from_matrix`].
    pub(crate) fn from_raw(matrix: Transform3) -> Placement {
        Placement(matrix)
    }

    /// The matrix.
    pub fn transform(&self) -> &Transform3 {
        &self.0
    }

    /// `self` followed by `after` (the motion `after ∘ self`).
    pub fn then(&self, after: &Placement) -> Placement {
        Placement(after.0.compose(&self.0))
    }

    /// The inverse motion.
    pub fn inverse(&self) -> Placement {
        // Rigid: Rᵀ and −Rᵀt, exact, no general inversion.
        let m = &self.0.matrix;
        let mut out = Transform3::identity();
        for (r, row) in out.matrix.iter_mut().enumerate().take(3) {
            for (c, value) in row.iter_mut().enumerate().take(3) {
                *value = m[c][r];
            }
            row[3] = -(m[0][r] * m[0][3] + m[1][r] * m[1][3] + m[2][r] * m[2][3]);
        }
        Placement(out)
    }

    /// Exactly the identity matrix.
    pub fn is_identity(&self) -> bool {
        self.0 == Transform3::identity()
    }

    /// No rotation (the 3×3 block is exactly the identity).
    pub fn is_translation(&self) -> bool {
        let m = &self.0.matrix;
        (0..3).all(|r| (0..3).all(|c| m[r][c] == if r == c { 1.0 } else { 0.0 }))
    }

    /// The translation (where the local origin lands).
    pub fn origin(&self) -> Point3 {
        let m = &self.0.matrix;
        Point3::new(m[0][3], m[1][3], m[2][3])
    }

    /// Where the local z axis points.
    pub fn z_axis(&self) -> Vec3 {
        let m = &self.0.matrix;
        Vec3::new(m[0][2], m[1][2], m[2][2])
    }

    /// Where the local x axis points.
    pub fn x_axis(&self) -> Vec3 {
        let m = &self.0.matrix;
        Vec3::new(m[0][0], m[1][0], m[2][0])
    }

    /// Apply to a point.
    pub fn apply_point(&self, p: &Point3) -> Point3 {
        self.0.apply_point(p)
    }

    /// The matrix as 16 column-major `f32`s (a GPU instance matrix).
    pub fn column_major_f32(&self) -> [f32; 16] {
        let m = &self.0.matrix;
        let mut out = [0.0f32; 16];
        for c in 0..4 {
            for r in 0..4 {
                out[c * 4 + r] = m[r][c] as f32;
            }
        }
        out
    }
}

impl Default for Placement {
    fn default() -> Self {
        Placement::identity()
    }
}

impl TryFrom<[[f64; 4]; 4]> for Placement {
    type Error = AssemblyError;
    fn try_from(matrix: [[f64; 4]; 4]) -> Result<Self, Self::Error> {
        Placement::from_matrix(Transform3 { matrix }, "placement")
    }
}

impl From<Placement> for [[f64; 4]; 4] {
    fn from(p: Placement) -> Self {
        p.0.matrix
    }
}

fn non_rigid(instance: &str, defect: f64, mirror: bool) -> AssemblyError {
    AssemblyError::NonRigidPlacement {
        instance: instance.to_string(),
        defect,
        mirror: if mirror { ", mirrored" } else { "" },
    }
}

/// Largest `|RᵀR − I|` entry of the 3×3 block.
fn orthonormality_defect(t: &Transform3) -> f64 {
    let m = &t.matrix;
    let mut worst: f64 = 0.0;
    for a in 0..3 {
        for b in 0..3 {
            let dot: f64 = (0..3).map(|r| m[r][a] * m[r][b]).sum();
            let expected = if a == b { 1.0 } else { 0.0 };
            worst = worst.max((dot - expected).abs());
        }
    }
    worst
}

fn determinant(t: &Transform3) -> f64 {
    let m = &t.matrix;
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

#[cfg(test)]
#[path = "placement_tests.rs"]
mod tests;
