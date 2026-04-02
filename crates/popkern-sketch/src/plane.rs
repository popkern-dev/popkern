//! 3D sketch plane for projecting 2D sketch geometry into world space.

use crate::error::SketchError;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

/// A plane in 3D space on which a 2D sketch is drawn.
///
/// Defined by an origin point and two orthonormal axes (u, v) that span
/// the plane, plus a normal vector perpendicular to both. The 2D sketch
/// coordinate `(x, y)` maps to `origin + x * u_axis + y * v_axis`.
#[derive(Debug, Clone, Copy)]
pub struct SketchPlane {
    /// Origin point of the plane (maps to sketch coordinate (0, 0)).
    pub origin: Point3,
    /// Unit vector along the sketch X axis.
    pub u_axis: Vec3,
    /// Unit vector along the sketch Y axis.
    pub v_axis: Vec3,
    /// Unit normal of the plane (`u_axis x v_axis`).
    pub normal: Vec3,
}

impl SketchPlane {
    /// Create a sketch plane with explicit origin and axis vectors.
    ///
    /// The caller is responsible for ensuring `u_axis` and `v_axis` are
    /// orthonormal. The normal is computed as `u_axis x v_axis`.
    ///
    /// # Errors
    ///
    /// Returns [`SketchError::DegenerateNormal`] if the cross product of
    /// the axes is zero-length.
    pub fn new(origin: Point3, u_axis: Vec3, v_axis: Vec3) -> Result<Self, SketchError> {
        let normal = u_axis
            .cross(&v_axis)
            .normalized()
            .ok_or(SketchError::DegenerateNormal)?;
        Ok(Self {
            origin,
            u_axis,
            v_axis,
            normal,
        })
    }

    /// The XY plane at the world origin (Z = 0).
    pub fn xy() -> Self {
        Self {
            origin: Point3::new(0.0, 0.0, 0.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        }
    }

    /// The XZ plane at the world origin (Y = 0).
    pub fn xz() -> Self {
        Self {
            origin: Point3::new(0.0, 0.0, 0.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 0.0, 1.0),
            normal: Vec3::new(0.0, -1.0, 0.0),
        }
    }

    /// The YZ plane at the world origin (X = 0).
    pub fn yz() -> Self {
        Self {
            origin: Point3::new(0.0, 0.0, 0.0),
            u_axis: Vec3::new(0.0, 1.0, 0.0),
            v_axis: Vec3::new(0.0, 0.0, 1.0),
            normal: Vec3::new(1.0, 0.0, 0.0),
        }
    }

    /// Create a plane from a point and a normal direction.
    ///
    /// The u and v axes are computed by crossing the normal with a
    /// candidate vector (Z-axis preferred, falls back to X-axis if the
    /// normal is nearly parallel to Z).
    ///
    /// # Errors
    ///
    /// Returns [`SketchError::DegenerateNormal`] if `normal` is
    /// zero-length.
    pub fn from_normal(origin: Point3, normal: Vec3) -> Result<Self, SketchError> {
        let normal = normal.normalized().ok_or(SketchError::DegenerateNormal)?;

        // Pick a candidate vector not parallel to normal.
        let candidate = if normal.cross(&Vec3::new(0.0, 0.0, 1.0)).length() > 1e-6 {
            Vec3::new(0.0, 0.0, 1.0)
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        };

        let u_axis = normal
            .cross(&candidate)
            .normalized()
            .ok_or(SketchError::DegenerateNormal)?;
        let v_axis = normal
            .cross(&u_axis)
            .normalized()
            .ok_or(SketchError::DegenerateNormal)?;

        // Ensure right-hand rule: u x v = normal
        // If u_axis cross v_axis points opposite to normal, flip v.
        let check = u_axis.cross(&v_axis);
        let v_axis = if check.dot(&normal) < 0.0 {
            -v_axis
        } else {
            v_axis
        };

        Ok(Self {
            origin,
            u_axis,
            v_axis,
            normal,
        })
    }

    /// Create a plane from three non-collinear points.
    ///
    /// `p1` is the origin, `p1 -> p2` defines the u-axis direction, and
    /// the normal is computed from `(p2 - p1) x (p3 - p1)`.
    ///
    /// # Errors
    ///
    /// Returns [`SketchError::CollinearPoints`] if the three points are
    /// collinear or any two coincide.
    pub fn from_points(p1: Point3, p2: Point3, p3: Point3) -> Result<Self, SketchError> {
        let raw_u = p2 - p1;
        let raw_v = p3 - p1;
        let raw_normal = raw_u.cross(&raw_v);

        let normal = raw_normal
            .normalized()
            .ok_or(SketchError::CollinearPoints)?;
        let u_axis = raw_u.normalized().ok_or(SketchError::CollinearPoints)?;

        // v = normal x u to ensure orthogonality
        let v_axis = normal
            .cross(&u_axis)
            .normalized()
            .ok_or(SketchError::CollinearPoints)?;

        Ok(Self {
            origin: p1,
            u_axis,
            v_axis,
            normal,
        })
    }

    /// Extract a sketch plane from a BRep face's planar surface.
    ///
    /// # Errors
    ///
    /// Returns [`SketchError::NonPlanarFace`] if the face has no surface
    /// or is not planar.
    pub fn from_face(brep: &BRep, face: FaceHandle) -> Result<Self, SketchError> {
        let face_data = brep.faces.get(face).ok_or(SketchError::NonPlanarFace)?;
        let sof = face_data
            .surface
            .as_ref()
            .ok_or(SketchError::NonPlanarFace)?;

        match &sof.geometry {
            SurfaceGeometry::Plane {
                origin,
                normal,
                u_axis,
                v_axis,
            } => {
                let (normal, u_axis, v_axis) = if sof.orientation {
                    (*normal, *u_axis, *v_axis)
                } else {
                    (-*normal, *v_axis, *u_axis)
                };
                Ok(Self {
                    origin: *origin,
                    u_axis,
                    v_axis,
                    normal,
                })
            }
            _ => Err(SketchError::NonPlanarFace),
        }
    }

    /// Extract a sketch plane from a BRep face, offset along the normal.
    ///
    /// Equivalent to [`from_face`](Self::from_face) but shifts the origin
    /// by `distance` along the face normal.
    ///
    /// # Errors
    ///
    /// Returns [`SketchError::NonPlanarFace`] if the face is not planar.
    pub fn from_face_offset(
        brep: &BRep,
        face: FaceHandle,
        distance: f64,
    ) -> Result<Self, SketchError> {
        let mut plane = Self::from_face(brep, face)?;
        plane.origin = plane.origin + plane.normal * distance;
        Ok(plane)
    }

    /// Reverse the plane orientation: negate the normal and swap u/v axes.
    pub fn flip(&mut self) {
        self.normal = -self.normal;
        std::mem::swap(&mut self.u_axis, &mut self.v_axis);
    }

    /// Convert a 2D sketch coordinate to a 3D world point.
    ///
    /// `origin + x * u_axis + y * v_axis`
    pub fn to_3d(&self, x: f64, y: f64) -> Point3 {
        self.origin + self.u_axis * x + self.v_axis * y
    }

    /// Project a 3D world point onto the plane, returning 2D coordinates.
    ///
    /// Uses dot products against the plane axes:
    /// `x = (point - origin) . u_axis`
    /// `y = (point - origin) . v_axis`
    pub fn to_2d(&self, point: Point3) -> (f64, f64) {
        let d = point - self.origin;
        (d.dot(&self.u_axis), d.dot(&self.v_axis))
    }
}

#[cfg(test)]
#[path = "plane_tests.rs"]
mod tests;
