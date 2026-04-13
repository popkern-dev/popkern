//! Undo/redo integration for boolean operations (§15).
//!
//! [`BoolOperationRecord`] captures everything a higher-level
//! application needs to undo, redo, or re-execute a boolean
//! operation: the operands, the operation kind, the tolerance, the
//! result, and the diagnostics. Applications managing a history
//! list store a `Vec<BoolOperationRecord>` and call the
//! [`undo_state`](BoolOperationRecord::undo_state) /
//! [`redo_state`](BoolOperationRecord::redo_state) accessors to
//! navigate the history.
//!
//! ## Why this is the right abstraction
//!
//! popkern-bool's boolean entry points already take owned `BRep`
//! inputs and return a fresh `BRep` result — they never mutate
//! the caller's data. This means the only thing an undo system
//! needs from the kernel is a place to *snapshot* the inputs and
//! the result so a higher-level history list can re-bind them on
//! undo / redo.
//!
//! The kernel does not manage the history list itself: that's the
//! application's concern (each app has its own history-list
//! semantics, undo limits, persistence model, etc.). The kernel
//! just provides the snapshot type and an entry point that
//! produces one.
//!
//! ## Example
//!
//! ```ignore
//! let record = popkern_bool::history::boolean_recorded(
//!     popkern_bool::BoolOperation::Union,
//!     &cube_a,
//!     &cube_b,
//!     tolerance,
//! )?;
//!
//! // The record is the application's history-list entry.
//! history.push(record);
//!
//! // On undo, restore the inputs:
//! let (a, b) = history.last().unwrap().undo_state();
//! current_model = (a.clone(), b.clone());
//!
//! // On redo, restore the cached result (no recomputation):
//! current_model = history.last().unwrap().redo_state().clone();
//!
//! // Or re-execute deterministically (e.g. after the inputs were
//! // restored from disk):
//! let new_result = history.last().unwrap().re_execute()?;
//! ```

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::{BoolError, BoolOperation, PhaseMetric};

/// A record of one boolean operation, suitable for storage in an
/// application-level undo/redo history list.
///
/// Created by [`boolean_recorded`]. The struct is `Clone` so an
/// application can deep-copy a record into a snapshot store, and
/// it owns its `input_a` / `input_b` / `result` BReps so it does
/// not depend on the lifetime of the original inputs.
///
/// Only metrics are kept (not the full [`crate::error::BoolDiagnostics`])
/// because the record is created on the success path — there's no
/// error or warning data to capture, and `BoolDiagnostics` is not
/// `Clone` due to embedded `thiserror::Error` types.
#[derive(Debug, Clone)]
pub struct BoolOperationRecord {
    /// Which boolean operation was performed.
    pub op: BoolOperation,
    /// Snapshot of the first operand at execution time.
    pub input_a: BRep,
    /// Snapshot of the second operand at execution time.
    pub input_b: BRep,
    /// The tolerance the operation ran with.
    pub tolerance: Tolerance,
    /// The result solid produced by the operation.
    pub result: BRep,
    /// Per-phase timing and counter samples (§16 telemetry) — for
    /// observability when replaying a history.
    pub phase_metrics: Vec<PhaseMetric>,
}

impl BoolOperationRecord {
    /// Inputs that were passed to the boolean — i.e. the state the
    /// caller should restore on undo.
    ///
    /// Returns references; clone them if the application needs
    /// owned values.
    pub fn undo_state(&self) -> (&BRep, &BRep) {
        (&self.input_a, &self.input_b)
    }

    /// The cached result the boolean produced — i.e. the state the
    /// caller should restore on redo, without paying the cost of
    /// recomputation.
    pub fn redo_state(&self) -> &BRep {
        &self.result
    }

    /// Re-run the boolean operation deterministically using the
    /// stored inputs and tolerance.
    ///
    /// Useful when the application reloads a history list from
    /// disk and wants to verify that the cached `result` still
    /// matches what the engine produces today (e.g. after a
    /// kernel upgrade), or when the cached `result` BRep has
    /// drifted out of sync with the rest of the model.
    pub fn re_execute(&self) -> Result<BRep, BoolError> {
        crate::boolean(
            self.op,
            self.input_a.clone(),
            self.input_b.clone(),
            self.tolerance,
            None,
        )
    }
}

/// Run a boolean operation and package the inputs, result, and
/// per-phase metrics into a [`BoolOperationRecord`] suitable for
/// an application-level history list.
///
/// Wraps [`crate::boolean_with_diagnostics`] so the per-phase
/// timing samples are captured for replay observability. On
/// failure returns the same `BoolError` the wrapped call would.
pub fn boolean_recorded(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tolerance: Tolerance,
) -> Result<BoolOperationRecord, BoolError> {
    let (result, diagnostics) = crate::boolean_with_diagnostics(op, a, b, tolerance)?;
    Ok(BoolOperationRecord {
        op,
        input_a: a.clone(),
        input_b: b.clone(),
        tolerance,
        result,
        phase_metrics: diagnostics.phase_metrics,
    })
}

#[cfg(test)]
#[path = "history_tests.rs"]
mod tests;
