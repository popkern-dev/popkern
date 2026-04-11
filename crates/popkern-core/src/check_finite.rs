//! Finite-value guards for geometric primitives.
//!
//! These functions reject NaN and Infinity early, before they propagate
//! through geometric computations and produce confusing downstream errors.

use crate::error::CoreError;
use crate::geometry::{Point3, Vec3};

/// Returns `Ok(())` if all components of `p` are finite (not NaN or Inf).
///
/// Returns `Err(CoreError::PreconditionFailed)` otherwise.
pub fn check_finite_point(p: &Point3) -> Result<(), CoreError> {
    if p.x.is_finite() && p.y.is_finite() && p.z.is_finite() {
        Ok(())
    } else {
        Err(CoreError::PreconditionFailed(format!(
            "non-finite point component: ({}, {}, {})",
            p.x, p.y, p.z
        )))
    }
}

/// Returns `Ok(())` if all components of `v` are finite (not NaN or Inf).
///
/// Returns `Err(CoreError::PreconditionFailed)` otherwise.
pub fn check_finite_vec(v: &Vec3) -> Result<(), CoreError> {
    if v.x.is_finite() && v.y.is_finite() && v.z.is_finite() {
        Ok(())
    } else {
        Err(CoreError::PreconditionFailed(format!(
            "non-finite vector component: ({}, {}, {})",
            v.x, v.y, v.z
        )))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_point_ok() {
        let p = Point3::new(1.0, 2.0, 3.0);
        assert!(check_finite_point(&p).is_ok());
    }

    #[test]
    fn nan_point_x() {
        let p = Point3::new(f64::NAN, 0.0, 0.0);
        assert!(check_finite_point(&p).is_err());
    }

    #[test]
    fn nan_point_y() {
        let p = Point3::new(0.0, f64::NAN, 0.0);
        assert!(check_finite_point(&p).is_err());
    }

    #[test]
    fn nan_point_z() {
        let p = Point3::new(0.0, 0.0, f64::NAN);
        assert!(check_finite_point(&p).is_err());
    }

    #[test]
    fn inf_point() {
        let p = Point3::new(f64::INFINITY, 0.0, 0.0);
        assert!(check_finite_point(&p).is_err());
    }

    #[test]
    fn neg_inf_point() {
        let p = Point3::new(0.0, f64::NEG_INFINITY, 0.0);
        assert!(check_finite_point(&p).is_err());
    }

    #[test]
    fn finite_vec_ok() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        assert!(check_finite_vec(&v).is_ok());
    }

    #[test]
    fn nan_vec_x() {
        let v = Vec3::new(f64::NAN, 0.0, 0.0);
        assert!(check_finite_vec(&v).is_err());
    }

    #[test]
    fn nan_vec_y() {
        let v = Vec3::new(0.0, f64::NAN, 0.0);
        assert!(check_finite_vec(&v).is_err());
    }

    #[test]
    fn nan_vec_z() {
        let v = Vec3::new(0.0, 0.0, f64::NAN);
        assert!(check_finite_vec(&v).is_err());
    }

    #[test]
    fn inf_vec() {
        let v = Vec3::new(0.0, 0.0, f64::INFINITY);
        assert!(check_finite_vec(&v).is_err());
    }

    #[test]
    fn neg_inf_vec() {
        let v = Vec3::new(f64::NEG_INFINITY, 0.0, 0.0);
        assert!(check_finite_vec(&v).is_err());
    }

    #[test]
    fn zero_point_ok() {
        let p = Point3::new(0.0, 0.0, 0.0);
        assert!(check_finite_point(&p).is_ok());
    }

    #[test]
    fn zero_vec_ok() {
        let v = Vec3::new(0.0, 0.0, 0.0);
        assert!(check_finite_vec(&v).is_ok());
    }

    #[test]
    fn error_message_contains_values() {
        let p = Point3::new(f64::NAN, 1.0, 2.0);
        let err = check_finite_point(&p).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("non-finite"));
    }
}
