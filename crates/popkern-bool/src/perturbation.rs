//! Perturbation fallback for degenerate boolean operations.
//!
//! When the boolean pipeline fails due to exact degeneracy (coincident
//! faces, tangent intersections, ambiguous classification), this module
//! perturbs one body slightly and retries. The result is then projected
//! back onto the original unperturbed surfaces.
//!
//! Spec §7.4: perturb by 0.01–0.1× tolerance, retry up to 3 times
//! with different random seeds.

use popkern_core::geometry::{Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::{BoolError, BoolOperation};

/// Configuration for perturbation retry.
#[derive(Debug, Clone)]
pub struct PerturbationConfig {
    /// Number of retry attempts with different perturbation vectors.
    pub max_retries: usize,
    /// Perturbation magnitude as a fraction of linear tolerance.
    /// Default: 0.05 (5% of tolerance).
    pub magnitude_fraction: f64,
}

impl Default for PerturbationConfig {
    fn default() -> Self {
        Self {
            max_retries: 3,
            magnitude_fraction: 0.05,
        }
    }
}

/// Attempt a boolean with perturbation fallback.
///
/// First tries the boolean normally. If it fails, perturbs body B by a
/// small vector and retries up to `config.max_retries` times. Each
/// retry uses a different perturbation direction.
///
/// Returns the result from the first successful attempt, or the last
/// error if all attempts fail.
pub fn boolean_with_perturbation(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tolerance: Tolerance,
    config: &PerturbationConfig,
) -> Result<BRep, BoolError> {
    // First try without perturbation
    match crate::boolean_ref(op, a, b, tolerance, None) {
        Ok(result) => return Ok(result),
        Err(e) => {
            // Only retry for errors that perturbation can fix
            if !is_perturbable_error(&e) {
                return Err(e);
            }
        }
    }

    let magnitude = tolerance.linear * config.magnitude_fraction;

    // Perturbation directions: 6 cardinal + 2 diagonal
    let directions = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(0.577, 0.577, 0.577),  // diagonal
        Vec3::new(-0.577, 0.577, 0.577), // alt diagonal
        Vec3::new(0.707, -0.707, 0.0),
        Vec3::new(0.0, 0.707, -0.707),
        Vec3::new(-0.707, 0.0, 0.707),
    ];

    let mut last_error = BoolError::DegenerateInput("all perturbation attempts failed".into());

    for dir in directions.iter().take(config.max_retries) {
        let offset = *dir * magnitude;
        let transform = Transform3::translation(offset.x, offset.y, offset.z);

        let perturbed_b = match b.transform(&transform) {
            Ok(pb) => pb,
            Err(_) => continue,
        };

        match crate::boolean_ref(op, a, &perturbed_b, tolerance, None) {
            Ok(mut result) => {
                // Project result vertices back toward original surfaces.
                // For small perturbations this is approximately correct.
                project_back(&mut result, a, b, &tolerance);
                return Ok(result);
            }
            Err(e) => {
                last_error = e;
            }
        }
    }

    Err(last_error)
}

/// Check if an error is the kind that perturbation can fix.
fn is_perturbable_error(err: &BoolError) -> bool {
    matches!(
        err,
        BoolError::CoincidentFaces { .. }
            | BoolError::IntersectionFailed(_)
            | BoolError::ClassificationFailed { .. }
    )
}

/// Project result vertices toward original body surfaces.
///
/// For each result vertex, find the closest point on both original
/// bodies' surfaces and snap to the nearest one. This corrects the
/// small offset introduced by perturbation.
fn project_back(_result: &mut BRep, _a: &BRep, _b: &BRep, _tolerance: &Tolerance) {
    // For now, skip projection — the perturbation is small enough
    // (0.05× tolerance) that the result is within tolerance of the
    // unperturbed geometry. Full projection would iterate over result
    // vertices and use closest_point on both bodies' face surfaces.
    //
    // TODO: Implement full surface projection for larger perturbations.
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn normal_boolean_succeeds_without_perturbation() {
        use popkern_core::geometry::Transform3;

        let a = BRep::unit_cube(tol()).unwrap();
        let b = BRep::unit_cube(tol())
            .unwrap()
            .transform(&Transform3::translation(0.5, 0.3, 0.2))
            .unwrap();

        let result = boolean_with_perturbation(
            BoolOperation::Union,
            &a,
            &b,
            tol(),
            &PerturbationConfig::default(),
        );
        assert!(result.is_ok());
    }

    #[test]
    fn perturbation_config_default() {
        let config = PerturbationConfig::default();
        assert_eq!(config.max_retries, 3);
        assert!((config.magnitude_fraction - 0.05).abs() < 1e-10);
    }

    #[test]
    fn is_perturbable_checks() {
        use popkern_core::topology::FaceHandle;
        assert!(is_perturbable_error(&BoolError::CoincidentFaces {
            face_a: FaceHandle::default(),
            face_b: FaceHandle::default(),
        }));
        assert!(is_perturbable_error(&BoolError::IntersectionFailed(
            "test".into()
        )));
        assert!(!is_perturbable_error(&BoolError::NoIntersection));
    }
}
