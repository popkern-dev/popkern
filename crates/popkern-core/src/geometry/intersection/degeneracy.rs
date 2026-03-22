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

/// Robust handler: attempts recovery for all degeneracies instead of skipping.
///
/// Unlike [`DefaultDegeneracyHandler`], this handler returns `false` for all
/// conditions, signalling to the caller that it should attempt recovery
/// (tangent tracing, normal estimation, adaptive stepping) rather than
/// skipping or stopping.
pub struct RobustDegeneracyHandler;

impl DegeneracyHandler for RobustDegeneracyHandler {
    fn on_tangent_starting_point(&self, _point: &Point3, _angle: f64) -> bool {
        false // don't skip — attempt tangent tracing
    }
    fn on_tangent_during_march(&self, _point: &Point3, _angle: f64) -> bool {
        false // don't stop — adapt
    }
    fn on_singularity(&self, _point: &Point3) -> bool {
        false // don't skip — attempt normal estimation
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_handler_skips_tangent_start() {
        let handler = DefaultDegeneracyHandler;
        assert!(handler.on_tangent_starting_point(&Point3::new(0.0, 0.0, 0.0), 0.1));
    }

    #[test]
    fn default_handler_stops_tangent_march() {
        let handler = DefaultDegeneracyHandler;
        assert!(handler.on_tangent_during_march(&Point3::new(1.0, 2.0, 3.0), 0.01));
    }

    #[test]
    fn default_handler_skips_singularity() {
        let handler = DefaultDegeneracyHandler;
        assert!(handler.on_singularity(&Point3::new(0.0, 0.0, 0.0)));
    }

    #[test]
    fn robust_handler_does_not_skip_tangent_start() {
        let handler = RobustDegeneracyHandler;
        assert!(!handler.on_tangent_starting_point(&Point3::new(0.0, 0.0, 0.0), 0.1));
    }

    #[test]
    fn robust_handler_does_not_stop_tangent_march() {
        let handler = RobustDegeneracyHandler;
        assert!(!handler.on_tangent_during_march(&Point3::new(1.0, 2.0, 3.0), 0.01));
    }

    #[test]
    fn robust_handler_does_not_skip_singularity() {
        let handler = RobustDegeneracyHandler;
        assert!(!handler.on_singularity(&Point3::new(0.0, 0.0, 0.0)));
    }
}
