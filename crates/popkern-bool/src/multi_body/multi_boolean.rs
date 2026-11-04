//! The boolean of N bodies by a predicate on region membership.
//!
//! One partition of all the bodies (`crate::partition`), its regions
//! kept by the predicate: exact for every predicate, where the former
//! sequential pairwise fold inferred one operation per step and
//! evaluated "exactly one" or "at least two" wrongly. A partition that
//! declines falls back to two-operand booleans over the predicate's
//! disjunctive form (recorded), for at most
//! [`crate::FALLBACK_COUNTING_LIMIT`] bodies.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::BoolError;
use crate::partition::BodySet;

/// The regions of `bodies` (at least two) that `predicate` keeps, as one
/// body: `predicate` sees the set of bodies containing each region.
///
/// # Errors
///
/// [`BoolError::DegenerateInput`] for fewer than two bodies; the
/// partition's errors otherwise (see [`crate::nary_boolean_where`]).
pub fn multi_body_boolean(
    bodies: &[BRep],
    predicate: impl Fn(&BodySet) -> bool,
    tolerance: &Tolerance,
) -> Result<BRep, BoolError> {
    if bodies.len() < 2 {
        return Err(BoolError::DegenerateInput(
            "multi-body boolean requires at least 2 bodies".to_string(),
        ));
    }
    let refs: Vec<&BRep> = bodies.iter().collect();
    crate::partition::nary_boolean_where(
        &refs,
        predicate,
        *tolerance,
        &crate::BoolOptions::default(),
    )
}
