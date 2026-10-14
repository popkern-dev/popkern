//! The serve postprocess's validity gate (PK-804).
//!
//! Each structural pass runs only on a body that validates. The gate
//! measured `validate()` afresh before every pass — seven whole-body
//! validations per boolean (~0.3 ms each on a 255-face body), while a
//! typical boolean's passes find nothing to do. The verdict is kept
//! until a pass reports that it edited the body (a count of applied
//! or declined edits — declines count, since a declined transactional
//! edit is not proven bit-identical to its snapshot), then re-measured.

use popkern_core::topology::BRep;

/// A `validate()` verdict, re-measured only after a reported edit.
pub(crate) struct ValidityGate(Option<bool>);

impl ValidityGate {
    /// No verdict yet: the first [`ValidityGate::valid`] measures.
    pub(crate) fn new() -> Self {
        Self(None)
    }

    /// Whether `brep` validates (measured when no verdict is kept).
    pub(crate) fn valid(&mut self, brep: &BRep) -> bool {
        *self.0.get_or_insert_with(|| brep.validate().is_ok())
    }

    /// Drop the kept verdict when a pass reports `edits > 0`.
    pub(crate) fn edited(&mut self, edits: usize) {
        if edits > 0 {
            self.0 = None;
        }
    }
}
