//! Surface geometry types, evaluation trait, and topology binding.

use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// Evaluation interface for parametric surfaces.
pub trait SurfaceEval {
    /// Point on the surface at parameters (u, v).
    fn point_at(&self, u: f64, v: f64) -> Point3;
    /// Partial derivative with respect to u.
    fn partial_u(&self, u: f64, v: f64) -> Vec3;
    /// Partial derivative with respect to v.
    fn partial_v(&self, u: f64, v: f64) -> Vec3;
    /// Surface normal at (u, v). May not be unit length for degenerate points.
    fn normal_at(&self, u: f64, v: f64) -> Vec3;
    /// Valid parameter range ((u_min, u_max), (v_min, v_max)).
    fn parameter_range(&self) -> ((f64, f64), (f64, f64));
}

/// Geometric surface types.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SurfaceGeometry {
    /// An infinite plane.
    Plane {
        /// A point on the plane.
        origin: Point3,
        /// Unit normal of the plane.
        normal: Vec3,
        /// Unit vector along the u parameter direction.
        u_axis: Vec3,
        /// Unit vector along the v parameter direction.
        v_axis: Vec3,
    },
    /// A cylinder (infinite along axis).
    Cylinder {
        /// Point on the axis (origin of the parameterization).
        origin: Point3,
        /// Unit vector along the cylinder axis.
        axis: Vec3,
        /// Radius of the cylinder.
        radius: f64,
    },
    /// A cone.
    Cone {
        /// Apex (tip) of the cone.
        apex: Point3,
        /// Unit vector along the cone axis (pointing away from apex into the body).
        axis: Vec3,
        /// Half-angle between the axis and the surface generator, in radians.
        half_angle: f64,
    },
    /// A sphere.
    Sphere {
        /// Center of the sphere.
        center: Point3,
        /// Radius of the sphere.
        radius: f64,
    },
    /// A torus.
    Torus {
        /// Center of the torus.
        center: Point3,
        /// Unit vector along the torus symmetry axis.
        axis: Vec3,
        /// Distance from the torus center to the center of the tube.
        major_radius: f64,
        /// Radius of the tube.
        minor_radius: f64,
    },
}

/// Binding between a face and its geometric surface.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SurfaceOnFace {
    /// The geometric surface.
    pub geometry: SurfaceGeometry,
    /// True if surface normal agrees with face normal.
    pub orientation: bool,
}
