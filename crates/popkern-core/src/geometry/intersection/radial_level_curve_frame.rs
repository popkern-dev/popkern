//! Axisymmetric frame of the radial level-curve arm's analytic side
//! (cylinder or cone, PK-461/PK-799).

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

/// Axisymmetric frame of the analytic surface: unit axis with the same
/// chart basis `point_at` uses, plus the raw axis for the
/// `v = rel·axis/|axis|²` inversion. The surface is where the distance
/// to the axis equals `radius + slope·a`, `a` the axial coordinate from
/// `origin`: a cylinder (slope 0) or a cone (origin at the apex,
/// radius 0, slope tan of the half-angle; PK-799).
pub(super) struct AxisymmetricFrame {
    pub(super) origin: Point3,
    pub(super) axis_raw: Vec3,
    pub(super) axis: Vec3,
    pub(super) e1: Vec3,
    pub(super) e2: Vec3,
    pub(super) radius: f64,
    pub(super) slope: f64,
}

impl AxisymmetricFrame {
    /// The frame of a cylinder, or of a cone whose face stays on the
    /// `v > 0` nappe (the level field's zero set is that nappe alone, so
    /// a face reaching the apex or beyond declines to marching rather
    /// than tracing an incomplete answer). `None` for any other surface.
    pub(super) fn from_surface(cyl: &SurfaceGeometry, cyl_v_range: (f64, f64)) -> Option<Self> {
        let frame = match cyl {
            SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } => {
                let axis_unit = axis.normalized()?;
                let (e1, e2) = orthonormal_basis(axis_unit);
                AxisymmetricFrame {
                    origin: *origin,
                    axis_raw: *axis,
                    axis: axis_unit,
                    e1,
                    e2,
                    radius: *radius,
                    slope: 0.0,
                }
            }
            // A chamfer cone against a thread's flanks: the same radial
            // level-set model with the radius growing along the axis. The
            // marcher stepped past the flanks' clamped ends and paid a
            // multi-start projection per landing (an M5x12 screw's six
            // cone x flank pairs: 2.06 of its 2.5 SSI-seconds; PK-799).
            SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle,
            } => {
                // The chart's second basis vector is `axis × perp1`, of
                // length `|axis|`: a non-unit axis evaluates an elliptic
                // cone, which no radial level field models.
                if cyl_v_range.0.min(cyl_v_range.1) <= 0.0 || (axis.length() - 1.0).abs() > 1e-12 {
                    return None;
                }
                let axis_unit = axis.normalized()?;
                let (e1, e2) = orthonormal_basis(*axis);
                AxisymmetricFrame {
                    origin: *apex,
                    axis_raw: *axis,
                    axis: axis_unit,
                    e1,
                    e2,
                    radius: 0.0,
                    slope: half_angle.ptan(),
                }
            }
            _ => return None,
        };
        Some(frame)
    }
    pub(super) fn radial(&self, p: &Point3) -> f64 {
        let rel = *p - self.origin;
        (rel - self.axis * rel.dot(&self.axis)).length()
    }
    /// The level field whose zero set is the surface: distance to the
    /// axis less the surface's radius at `p`'s axial position.
    pub(super) fn level(&self, p: &Point3) -> f64 {
        let axial = (*p - self.origin).dot(&self.axis);
        self.radial(p) - (self.radius + self.slope * axial)
    }
    /// Distance from `p` to the surface near it (the level across the
    /// generator's slope).
    pub(super) fn distance(&self, p: &Point3) -> f64 {
        self.level(p).abs() / (1.0 + self.slope * self.slope).sqrt()
    }
    /// The surface point at chart `(theta, v)`.
    pub(super) fn chart_point(&self, theta: f64, v: f64) -> Point3 {
        let r = self.radius + self.slope * v * self.axis_raw.length();
        self.origin
            + self.e1 * (r * theta.pcos())
            + self.e2 * (r * theta.psin())
            + self.axis_raw * v
    }
    /// Principal-branch angle in (-π, π] and axial chart coordinate.
    pub(super) fn theta_v(&self, p: &Point3) -> (f64, f64) {
        let rel = *p - self.origin;
        let theta = rel.dot(&self.e2).patan2(rel.dot(&self.e1));
        let v = rel.dot(&self.axis_raw) / self.axis_raw.length_squared();
        (theta, v)
    }
}
