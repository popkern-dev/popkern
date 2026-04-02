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
#[path = "boolean_tests.rs"]
mod tests;
