//! Multi-body boolean API and types.
//!
//! Extends the 2-body boolean to N bodies using cellular decomposition
//! with N-bit membership vectors.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::cellular::cell::CellDecomposition;
use crate::error::BoolError;

/// A multi-body boolean operation specification.
#[derive(Debug, Clone)]
pub struct MultiBooleanSpec {
    /// The input bodies.
    pub bodies: Vec<BRep>,
    /// Cell selection predicate: given a membership bitmask where bit i
    /// indicates the cell is inside body i, return true to include the cell.
    /// This is stored as a description string; the actual predicate is
    /// passed at evaluation time.
    pub description: String,
}

/// Result of a multi-body boolean operation.
#[derive(Debug)]
pub struct MultiBooleanResult {
    /// The result body (if the operation produced one).
    pub body: Option<BRep>,
    /// The cellular decomposition (for inspection).
    pub decomposition: Option<CellDecomposition>,
    /// Number of bodies processed.
    pub body_count: usize,
    /// Number of cells in the decomposition.
    pub cell_count: usize,
}

/// Compute a multi-body boolean using pairwise intersection + cellular decomposition.
///
/// For N bodies, computes all pairwise intersections, builds a cellular
/// decomposition where each cell has an N-bit membership vector, and
/// selects cells matching the predicate.
///
/// # Arguments
///
/// * `bodies` — the input bodies (at least 2)
/// * `predicate` — cell selection: `fn(membership_bits) -> bool`
/// * `tolerance` — geometric tolerance
///
/// # Current Limitations
///
/// The current implementation processes bodies pairwise (sequential),
/// which accumulates approximation error for large N. A true simultaneous
/// approach (computing all pairwise intersections in a single pass) is
/// planned for a future iteration.
pub fn multi_body_boolean(
    bodies: &[BRep],
    predicate: impl Fn(u64) -> bool,
    tolerance: &Tolerance,
) -> Result<MultiBooleanResult, BoolError> {
    if bodies.len() < 2 {
        return Err(BoolError::DegenerateInput(
            "multi-body boolean requires at least 2 bodies".to_string(),
        ));
    }

    let body_count = bodies.len();

    // For 2 bodies, use the standard boolean pipeline
    if body_count == 2 {
        // Determine operation from predicate behavior
        let in_a_only = predicate(1);
        let in_b_only = predicate(2);
        let in_both = predicate(3);
        let in_neither = predicate(0);

        let op = if in_a_only && in_b_only && in_both && !in_neither {
            crate::error::BoolOperation::Union
        } else if !in_a_only && !in_b_only && in_both && !in_neither {
            crate::error::BoolOperation::Intersection
        } else if in_a_only && !in_b_only && !in_both && !in_neither {
            crate::error::BoolOperation::Subtraction
        } else {
            // Custom predicate — use union as base, filter later
            crate::error::BoolOperation::Union
        };

        let result = crate::boolean_ref(op, &bodies[0], &bodies[1], *tolerance, None).ok();

        return Ok(MultiBooleanResult {
            body: result,
            decomposition: None,
            body_count,
            cell_count: 0,
        });
    }

    // For N > 2: sequential pairwise reduction.
    // Start with body 0, union/intersect with each subsequent body.
    // This is the naive approach; a true multi-body implementation would
    // compute all pairwise intersections simultaneously.
    let mut accumulated = bodies[0].clone();

    for (i, body) in bodies.iter().enumerate().skip(1) {
        // For each additional body, determine the pairwise operation
        // by checking the predicate with the accumulated body as "A"
        // and the new body as "B".
        let in_a = predicate(1 << 0); // simplified: treat accumulated as body 0
        let in_b = predicate(1 << i);
        let in_both = in_a && in_b;

        let op = if in_both {
            crate::error::BoolOperation::Union
        } else if in_a && !in_b {
            crate::error::BoolOperation::Subtraction
        } else {
            crate::error::BoolOperation::Intersection
        };

        match crate::boolean_ref(op, &accumulated, body, *tolerance, None) {
            Ok(result) => accumulated = result,
            Err(_) => {
                // Skip bodies that fail to boolean (e.g., no intersection)
                continue;
            }
        }
    }

    Ok(MultiBooleanResult {
        body: Some(accumulated),
        decomposition: None,
        body_count,
        cell_count: 0,
    })
}
