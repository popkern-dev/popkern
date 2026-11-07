//! Check that every ring of more than two fins is linked in its radial
//! order about the edge (§5 check 3).
//!
//! The wedge check (`cells`) reads which faces are neighbours about an
//! edge off the ring as linked; for a ring of two any order is the
//! order, but a longer ring must be the geometric order
//! ([`BRep::radial_order`]) or a wedge would join faces with others
//! between them. An order the geometry cannot decide is its own finding
//! (coincident faces, faces crossing inside the edge), never a guess.
//! Ordering evaluates the faces, so the check runs in `validate_all`
//! only.

use crate::error::ValidationError;
use crate::topology::{BRep, HalfEdgeHandle};

/// Whether `a` and `b` are the same cycle read from the same start.
fn same_from_start(a: &[HalfEdgeHandle], b: &[HalfEdgeHandle]) -> bool {
    a == b
}

impl BRep {
    /// Every ring of more than two fins out of its radial order, or
    /// whose order the geometry does not decide, at the body's
    /// tolerance.
    pub fn validate_ring_order_all(&self) -> Vec<ValidationError> {
        let tolerance = self.tolerance;
        let mut errors = Vec::new();
        for edge in self.edges.keys() {
            // Counted without collecting: most rings are two.
            if self.ring_summary(edge).count < 3 {
                continue;
            }
            let stored = self.radial_fins(edge);
            match self.radial_order(edge, &tolerance) {
                Ok(radial) if same_from_start(&stored, &radial) => {}
                Ok(radial) => errors.push(ValidationError::RingOrder {
                    edge,
                    stored,
                    radial,
                }),
                Err(error) => errors.push(ValidationError::RingOrderUndecided {
                    edge,
                    error: Box::new(error),
                }),
            }
        }
        errors
    }
}

#[cfg(test)]
#[path = "ring_order_tests.rs"]
mod tests;
