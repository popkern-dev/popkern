//! Structured tolerance type for geometric comparisons.

use crate::geometry::Point3;
use serde::{Deserialize, Serialize};

/// Tolerance values for geometric comparisons.
///
/// Separates linear (distance), angular (radians), and parametric
/// (curve/surface parameter) tolerances. Stored per-BRep, never global.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Tolerance {
    /// Tolerance for distance comparisons (point coincidence, edge length).
    pub linear: f64,
    /// Tolerance for angle comparisons in radians.
    pub angular: f64,
    /// Tolerance for parameter-space comparisons.
    pub parametric: f64,
}

impl Default for Tolerance {
    fn default() -> Self {
        Self {
            linear: 1e-6,
            angular: 1e-4,
            parametric: 1e-8,
        }
    }
}

impl Tolerance {
    /// Set linear tolerance.
    pub fn with_linear(mut self, linear: f64) -> Self {
        self.linear = linear;
        self
    }

    /// Set angular tolerance.
    pub fn with_angular(mut self, angular: f64) -> Self {
        self.angular = angular;
        self
    }

    /// Set parametric tolerance.
    pub fn with_parametric(mut self, parametric: f64) -> Self {
        self.parametric = parametric;
        self
    }

    /// Returns true if two points are within linear tolerance.
    pub fn points_coincident(&self, a: &Point3, b: &Point3) -> bool {
        a.distance_to(b) < self.linear
    }

    /// Returns true if two angles (in radians) are within angular tolerance.
    pub fn angles_equal(&self, a: f64, b: f64) -> bool {
        (a - b).abs() < self.angular
    }

    /// Returns true if two parameter values are within parametric tolerance.
    pub fn params_equal(&self, t1: f64, t2: f64) -> bool {
        (t1 - t2).abs() < self.parametric
    }

    /// Returns true if two distances are within linear tolerance.
    pub fn distances_equal(&self, d1: f64, d2: f64) -> bool {
        (d1 - d2).abs() < self.linear
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_tolerance_values() {
        let tol = Tolerance::default();
        assert_eq!(tol.linear, 1e-6);
        assert_eq!(tol.angular, 1e-4);
        assert_eq!(tol.parametric, 1e-8);
    }

    #[test]
    fn builder_pattern() {
        let tol = Tolerance::default().with_linear(1e-3).with_angular(1e-2);
        assert_eq!(tol.linear, 1e-3);
        assert_eq!(tol.angular, 1e-2);
        assert_eq!(tol.parametric, 1e-8);
    }

    #[test]
    fn points_coincident_within_tolerance() {
        let tol = Tolerance::default();
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1e-7, 0.0, 0.0);
        assert!(tol.points_coincident(&a, &b));
    }

    #[test]
    fn points_not_coincident_outside_tolerance() {
        let tol = Tolerance::default();
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1e-5, 0.0, 0.0);
        assert!(!tol.points_coincident(&a, &b));
    }

    #[test]
    fn angles_equal_within_tolerance() {
        let tol = Tolerance::default();
        assert!(tol.angles_equal(1.0, 1.0 + 1e-5));
    }

    #[test]
    fn angles_not_equal_outside_tolerance() {
        let tol = Tolerance::default();
        assert!(!tol.angles_equal(1.0, 1.0 + 1e-3));
    }

    #[test]
    fn params_equal_within_tolerance() {
        let tol = Tolerance::default();
        assert!(tol.params_equal(0.5, 0.5 + 1e-9));
    }

    #[test]
    fn distances_equal_within_tolerance() {
        let tol = Tolerance::default();
        assert!(tol.distances_equal(5.0, 5.0 + 1e-7));
    }

    #[test]
    fn with_parametric_builder() {
        let tol = Tolerance::default().with_parametric(1e-3);
        assert_eq!(tol.parametric, 1e-3);
        // linear and angular should be unchanged
        assert_eq!(tol.linear, 1e-6);
        assert_eq!(tol.angular, 1e-4);
    }

    #[test]
    fn with_angular_builder() {
        let tol = Tolerance::default().with_angular(1e-1);
        assert_eq!(tol.angular, 1e-1);
        assert_eq!(tol.linear, 1e-6);
        assert_eq!(tol.parametric, 1e-8);
    }

    #[test]
    fn chained_builders() {
        let tol = Tolerance::default()
            .with_linear(1e-3)
            .with_angular(1e-2)
            .with_parametric(1e-1);
        assert_eq!(tol.linear, 1e-3);
        assert_eq!(tol.angular, 1e-2);
        assert_eq!(tol.parametric, 1e-1);
    }
}
