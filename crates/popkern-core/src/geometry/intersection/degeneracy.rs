//! Pluggable degeneracy handling for surface-surface intersection.

use crate::geometry::Point3;

/// Extension point for degeneracy handling during surface-surface intersection.
pub trait DegeneracyHandler {
    /// Called when surfaces are tangent at a starting point.
    /// Returns true to skip this starting point.
    fn on_tangent_starting_point(&self, point: &Point3, angle: f64) -> bool;
    /// Called when normals become parallel during marching.
    /// Returns true to stop marching.
    fn on_tangent_during_march(&self, point: &Point3, angle: f64) -> bool;
    /// Called when a surface normal is zero (singularity).
    /// Returns true to skip.
    fn on_singularity(&self, point: &Point3) -> bool;
}

/// Default: skip tangent starts, stop on tangent during march, skip singularities.
pub struct DefaultDegeneracyHandler;

impl DegeneracyHandler for DefaultDegeneracyHandler {
    fn on_tangent_starting_point(&self, _point: &Point3, _angle: f64) -> bool {
        true
    }
    fn on_tangent_during_march(&self, _point: &Point3, _angle: f64) -> bool {
        true
    }
    fn on_singularity(&self, _point: &Point3) -> bool {
        true
    }
}
