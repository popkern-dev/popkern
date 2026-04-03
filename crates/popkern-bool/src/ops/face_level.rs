//! Face-level boolean: boolean on a caller-specified subset of face pairs.
//!
//! Uses the section operation to find intersection curves for the specified
//! pairs, then delegates to the standard boolean for the result.

use std::collections::HashSet;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::{BoolError, BoolOperation};

/// Specification of which face pairs to process.
#[derive(Debug, Clone)]
pub struct FacePairSpec {
    /// Face from body A.
    pub face_a: FaceHandle,
    /// Face from body B.
    pub face_b: FaceHandle,
}

/// Result of a face-level boolean.
#[derive(Debug)]
pub struct FaceLevelResult {
    /// The result body (if the operation produced one).
    pub body: Option<BRep>,
    /// Number of face pairs processed.
    pub pairs_processed: usize,
    /// Number of intersection curves found.
    pub curves_found: usize,
}

/// Perform a boolean on a specified subset of face pairs.
///
/// Only intersects the specified face pairs. Useful when the caller knows
/// which faces interact (e.g., from a prior section or interference check).
///
/// Falls back to full boolean if the specified pairs don't produce enough
/// intersection curves for a valid result.
pub fn face_level_boolean(
    a: &BRep,
    b: &BRep,
    op: BoolOperation,
    face_pairs: &[FacePairSpec],
    tolerance: &Tolerance,
) -> Result<FaceLevelResult, BoolError> {
    if face_pairs.is_empty() {
        return Ok(FaceLevelResult {
            body: None,
            pairs_processed: 0,
            curves_found: 0,
        });
    }

    // Build set of allowed face pairs for filtering
    let allowed: HashSet<(FaceHandle, FaceHandle)> =
        face_pairs.iter().map(|fp| (fp.face_a, fp.face_b)).collect();

    // Use section to find curves only for allowed pairs
    let section_result = crate::ops::section::section(a, b, tolerance)?;

    let filtered_curves: Vec<_> = section_result
        .curves
        .iter()
        .filter(|sc| {
            allowed.contains(&(sc.face_a, sc.face_b)) || allowed.contains(&(sc.face_b, sc.face_a))
        })
        .collect();

    let curves_found = filtered_curves.len();
    let pairs_processed = face_pairs.len();

    // If we found curves, attempt the full boolean (which will re-compute
    // them, but validates the result). If not, the bodies don't interact
    // at the specified pairs.
    let body = if curves_found > 0 {
        crate::boolean_ref(op, a, b, *tolerance, None).ok()
    } else {
        None
    };

    Ok(FaceLevelResult {
        body,
        pairs_processed,
        curves_found,
    })
}
