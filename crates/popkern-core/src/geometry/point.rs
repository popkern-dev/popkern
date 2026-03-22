//! 3D point type for positions in space.

use crate::geometry::Vec3;
use serde::{Deserialize, Serialize};
use std::ops::{Add, Sub};

/// A point in 3D space.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Point3 {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
    /// Z coordinate.
    pub z: f64,
}

impl Point3 {
    /// Creates a new point at the given coordinates.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Euclidean distance to another point.
    pub fn distance_to(&self, other: &Self) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Returns true if this point is within `tolerance` distance of `other`.
    pub fn approx_eq(&self, other: &Self, tolerance: f64) -> bool {
        self.distance_to(other) < tolerance
    }
}

impl Sub for Point3 {
    type Output = Vec3;
    fn sub(self, rhs: Self) -> Vec3 {
        Vec3::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Add<Vec3> for Point3 {
    type Output = Point3;
    fn add(self, rhs: Vec3) -> Point3 {
        Point3::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Sub<Vec3> for Point3 {
    type Output = Point3;
    fn sub(self, rhs: Vec3) -> Point3 {
        Point3::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point3_new() {
        let p = Point3::new(1.0, 2.0, 3.0);
        assert_eq!(p.x, 1.0);
        assert_eq!(p.y, 2.0);
        assert_eq!(p.z, 3.0);
    }

    #[test]
    fn point3_distance_to_self_is_zero() {
        let p = Point3::new(1.0, 2.0, 3.0);
        assert!((p.distance_to(&p)) < f64::EPSILON);
    }

    #[test]
    fn point3_distance_known_value() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        assert!((a.distance_to(&b) - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn point3_approx_eq_within_tolerance() {
        let a = Point3::new(1.0, 2.0, 3.0);
        let b = Point3::new(1.0 + 1e-8, 2.0, 3.0);
        assert!(a.approx_eq(&b, 1e-6));
    }

    #[test]
    fn point3_approx_eq_outside_tolerance() {
        let a = Point3::new(1.0, 2.0, 3.0);
        let b = Point3::new(1.1, 2.0, 3.0);
        assert!(!a.approx_eq(&b, 1e-6));
    }

    #[test]
    fn point_sub_point_gives_vec3() {
        let a = Point3::new(3.0, 4.0, 5.0);
        let b = Point3::new(1.0, 1.0, 1.0);
        let v = a - b;
        assert_eq!(v.x, 2.0);
        assert_eq!(v.y, 3.0);
        assert_eq!(v.z, 4.0);
    }

    #[test]
    fn point_add_vec3_gives_point() {
        let p = Point3::new(1.0, 2.0, 3.0);
        let v = Vec3::new(10.0, 20.0, 30.0);
        let q = p + v;
        assert_eq!(q.x, 11.0);
        assert_eq!(q.y, 22.0);
        assert_eq!(q.z, 33.0);
    }

    #[test]
    fn point_sub_vec3_gives_point() {
        let p = Point3::new(10.0, 20.0, 30.0);
        let v = Vec3::new(1.0, 2.0, 3.0);
        let q = p - v;
        assert_eq!(q.x, 9.0);
        assert_eq!(q.y, 18.0);
        assert_eq!(q.z, 27.0);
    }
}
