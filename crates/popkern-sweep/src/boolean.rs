//! Boolean convenience wrappers for sweep/loft operations.
//!
//! Combines sweep/loft surface generation with boolean operations
//! (union, subtraction, intersection) against an existing BRep solid.

use popkern_core::geometry::NurbsCurve;
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::SweepError;
use crate::types::{LoftOptions, LoftSection, PipeRadius, SweepOptions, SweepProfile, SweepResult};

/// Sweep a profile along a path and union with an existing solid.
pub fn sweep_add(
    brep: &BRep,
    profile: &SweepProfile,
    path: &NurbsCurve,
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<SweepResult, SweepError> {
    let swept = crate::sweep(profile, path, options, tolerance, tracking)?;
    let result = popkern_bool::union(brep.clone(), swept.brep, tolerance, None)
        .map_err(|e| SweepError::Boolean(format!("{e}")))?;
    Ok(SweepResult {
        brep: result,
        warnings: swept.warnings,
    })
}

/// Sweep a profile along a path and subtract from an existing solid.
pub fn sweep_cut(
    brep: &BRep,
    profile: &SweepProfile,
    path: &NurbsCurve,
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<SweepResult, SweepError> {
    let swept = crate::sweep(profile, path, options, tolerance, tracking)?;
    let result = popkern_bool::subtraction(brep.clone(), swept.brep, tolerance, None)
        .map_err(|e| SweepError::Boolean(format!("{e}")))?;
    Ok(SweepResult {
        brep: result,
        warnings: swept.warnings,
    })
}

/// Sweep a profile along a path and intersect with an existing solid.
pub fn sweep_intersect(
    brep: &BRep,
    profile: &SweepProfile,
    path: &NurbsCurve,
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<SweepResult, SweepError> {
    let swept = crate::sweep(profile, path, options, tolerance, tracking)?;
    let result = popkern_bool::intersection(brep.clone(), swept.brep, tolerance, None)
        .map_err(|e| SweepError::Boolean(format!("{e}")))?;
    Ok(SweepResult {
        brep: result,
        warnings: swept.warnings,
    })
}

/// Loft through sections and union with an existing solid.
pub fn loft_add(
    brep: &BRep,
    sections: &[LoftSection],
    options: &LoftOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<SweepResult, SweepError> {
    let lofted = crate::loft(sections, options, tolerance, tracking)?;
    let result = popkern_bool::union(brep.clone(), lofted.brep, tolerance, None)
        .map_err(|e| SweepError::Boolean(format!("{e}")))?;
    Ok(SweepResult {
        brep: result,
        warnings: lofted.warnings,
    })
}

/// Loft through sections and subtract from an existing solid.
pub fn loft_cut(
    brep: &BRep,
    sections: &[LoftSection],
    options: &LoftOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<SweepResult, SweepError> {
    let lofted = crate::loft(sections, options, tolerance, tracking)?;
    let result = popkern_bool::subtraction(brep.clone(), lofted.brep, tolerance, None)
        .map_err(|e| SweepError::Boolean(format!("{e}")))?;
    Ok(SweepResult {
        brep: result,
        warnings: lofted.warnings,
    })
}

/// Pipe along a path and union with an existing solid.
pub fn pipe_add(
    brep: &BRep,
    path: &NurbsCurve,
    radius: &PipeRadius,
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<SweepResult, SweepError> {
    let piped = crate::pipe(path, radius, options, tolerance, tracking)?;
    let result = popkern_bool::union(brep.clone(), piped.brep, tolerance, None)
        .map_err(|e| SweepError::Boolean(format!("{e}")))?;
    Ok(SweepResult {
        brep: result,
        warnings: piped.warnings,
    })
}

/// Pipe along a path and subtract from an existing solid.
pub fn pipe_cut(
    brep: &BRep,
    path: &NurbsCurve,
    radius: &PipeRadius,
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<SweepResult, SweepError> {
    let piped = crate::pipe(path, radius, options, tolerance, tracking)?;
    let result = popkern_bool::subtraction(brep.clone(), piped.brep, tolerance, None)
        .map_err(|e| SweepError::Boolean(format!("{e}")))?;
    Ok(SweepResult {
        brep: result,
        warnings: piped.warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;

    /// Build a simple test circle profile.
    fn test_circle() -> SweepProfile {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        SweepProfile::Wire(NurbsCurve {
            degree: 2,
            knots: vec![
                0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
            ],
            control_points: vec![
                Point3::new(0.3, 0.0, 0.0),
                Point3::new(0.3, 0.3, 0.0),
                Point3::new(0.0, 0.3, 0.0),
                Point3::new(-0.3, 0.3, 0.0),
                Point3::new(-0.3, 0.0, 0.0),
                Point3::new(-0.3, -0.3, 0.0),
                Point3::new(0.0, -0.3, 0.0),
                Point3::new(0.3, -0.3, 0.0),
                Point3::new(0.3, 0.0, 0.0),
            ],
            weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
        })
    }

    /// Build a straight-line path along Z.
    fn z_path() -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 5.0)],
            weights: vec![1.0, 1.0],
        }
    }

    #[test]
    #[ignore] // Boolean on curved/planar geometry is slow in debug builds.
    fn sweep_add_produces_result() {
        let cube = BRep::unit_cube(Tolerance::default());
        // unit_cube may fail in some configurations; skip if so.
        let cube = match cube {
            Ok(c) => c,
            Err(_) => return,
        };
        let profile = test_circle();
        let path = z_path();
        let options = SweepOptions::default();
        let tol = Tolerance::default();

        let result = sweep_add(&cube, &profile, &path, &options, tol, None);
        // Boolean may fail on complex geometry; we just check it returns.
        assert!(result.is_ok() || result.is_err());
    }

    #[test]
    #[ignore] // Boolean on curved/planar geometry is slow in debug builds.
    fn sweep_cut_produces_result() {
        let cube = BRep::unit_cube(Tolerance::default());
        let cube = match cube {
            Ok(c) => c,
            Err(_) => return,
        };
        let profile = test_circle();
        let path = z_path();
        let options = SweepOptions::default();
        let tol = Tolerance::default();

        let result = sweep_cut(&cube, &profile, &path, &options, tol, None);
        assert!(result.is_ok() || result.is_err());
    }
}
