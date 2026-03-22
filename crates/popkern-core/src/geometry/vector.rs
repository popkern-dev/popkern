//! 3D vector type for directions, normals, and displacements.

use serde::{Deserialize, Serialize};
use std::ops::{Add, Mul, Neg, Sub};

/// A 3D vector (direction/displacement, distinct from Point3).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Vec3 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

impl Vec3 {
    /// Creates a new vector.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// The zero vector.
    pub fn zero() -> Self {
        Self { x: 0.0, y: 0.0, z: 0.0 }
    }

    /// Dot product.
    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Cross product.
    pub fn cross(&self, other: &Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    /// Squared length (avoids sqrt).
    pub fn length_squared(&self) -> f64 {
        self.dot(self)
    }

    /// Length (magnitude).
    pub fn length(&self) -> f64 {
        self.length_squared().sqrt()
    }

    /// Unit vector. Returns None if zero-length.
    pub fn normalized(&self) -> Option<Self> {
        let len = self.length();
        if len < f64::EPSILON {
            None
        } else {
            Some(Self { x: self.x / len, y: self.y / len, z: self.z / len })
        }
    }

    /// Angle between two vectors in radians.
    pub fn angle_to(&self, other: &Self) -> f64 {
        let cos = self.dot(other) / (self.length() * other.length());
        cos.clamp(-1.0, 1.0).acos()
    }
}

impl Add for Vec3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self { x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z }
    }
}

impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self { x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z }
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self { x: self.x * rhs, y: self.y * rhs, z: self.z * rhs }
    }
}

impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self {
        Self { x: -self.x, y: -self.y, z: -self.z }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn new_and_fields() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        assert_eq!(v.x, 1.0);
        assert_eq!(v.y, 2.0);
        assert_eq!(v.z, 3.0);
    }

    #[test]
    fn zero() {
        let v = Vec3::zero();
        assert_eq!(v.x, 0.0);
    }

    #[test]
    fn dot_product() {
        let a = Vec3::new(1.0, 0.0, 0.0);
        let b = Vec3::new(0.0, 1.0, 0.0);
        assert!((a.dot(&b)).abs() < f64::EPSILON);
        assert!((a.dot(&a) - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn cross_product() {
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 1.0, 0.0);
        let z = x.cross(&y);
        assert!((z.z - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn length_and_length_squared() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        assert!((v.length_squared() - 25.0).abs() < f64::EPSILON);
        assert!((v.length() - 5.0).abs() < f64::EPSILON);
    }

    #[test]
    fn normalized_unit_vector() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        let n = v.normalized().unwrap();
        assert!((n.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn normalized_zero_returns_none() {
        assert!(Vec3::zero().normalized().is_none());
    }

    #[test]
    fn angle_to_perpendicular() {
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 1.0, 0.0);
        assert!((x.angle_to(&y) - FRAC_PI_2).abs() < 1e-10);
    }

    #[test]
    fn add_sub_mul_neg() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        assert_eq!((a + b).x, 5.0);
        assert_eq!((a - b).x, -3.0);
        assert_eq!((a * 2.0).x, 2.0);
        assert_eq!((-a).x, -1.0);
    }
}
