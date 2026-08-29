//! Undo/redo integration for fillet and chamfer operations (§17).
//!
//! [`FilletOperationRecord`] captures everything a higher-level
//! application needs to undo, redo, or re-execute a fillet
//! operation: the operands, the specs, the options, the tolerance,
//! and the result BRep. Applications managing a history list store
//! a `Vec<FilletOperationRecord>` and call the
//! [`undo_state`](FilletOperationRecord::undo_state) /
//! [`redo_state`](FilletOperationRecord::redo_state) accessors to
//! navigate the history.
//!
//! The kernel does not manage the history list itself — that's the
//! application's concern (each app has its own history-list
//! semantics, undo limits, persistence model, etc.). The kernel
//! just provides the snapshot type and an entry point that
//! produces one. This mirrors `popkern-bool::history`.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::BlendError;
use crate::types::{
    BlendResult, ChamferOptions, ChamferSpec, FilletOptions, FilletRadius, FilletSpec,
};

/// A record of one fillet operation, suitable for storage in an
/// application-level undo/redo history list.
///
/// Created by [`fillet_recorded`]. The struct is `Clone` so an
/// application can deep-copy a record into a snapshot store, and it
/// owns its `input` / `result` BReps so it does not depend on the
/// lifetime of the original inputs.
#[derive(Debug, Clone)]
pub struct FilletOperationRecord {
    /// Snapshot of the input BRep at the time the fillet ran.
    pub input: BRep,
    /// The fillet specs that were applied.
    pub specs: Vec<FilletSpec>,
    /// The options that were applied.
    pub options: FilletOptions,
    /// The tolerance the operation ran with.
    pub tolerance: Tolerance,
    /// The result BRep produced by the operation.
    pub result: BRep,
}

impl FilletOperationRecord {
    /// Input that was passed to the fillet — i.e. the state the
    /// caller should restore on undo.
    ///
    /// Returns a reference; clone it if the application needs an
    /// owned value.
    pub fn undo_state(&self) -> &BRep {
        &self.input
    }

    /// The cached result the fillet produced — i.e. the state the
    /// caller should restore on redo, without paying the cost of
    /// recomputation.
    pub fn redo_state(&self) -> &BRep {
        &self.result
    }

    /// Re-run the fillet operation deterministically using the
    /// stored inputs, specs, options, and tolerance.
    ///
    /// Useful when the application reloads a history list from
    /// disk and wants to verify that the cached `result` still
    /// matches what the engine produces today (e.g. after a
    /// kernel upgrade), or when the cached `result` BRep has
    /// drifted out of sync with the rest of the model.
    ///
    /// # Errors
    /// Same as [`crate::fillet_edges`].
    pub fn re_execute(&self) -> Result<BRep, BlendError> {
        let result = crate::fillet_edges(
            &self.input,
            &self.specs,
            &self.options,
            self.tolerance,
            None,
        )?;
        Ok(result.brep)
    }

    /// Re-run the fillet on the cached input with replacement
    /// specs, returning a fresh record. The original record is
    /// unchanged — callers managing an undo/redo list typically
    /// push the new record onto the history stack.
    ///
    /// This is the record-based path for §14.7 in-place
    /// modification: the caller substitutes the desired specs (e.g.
    /// to change the radius of one or more edges) and gets back a
    /// new record whose `result` reflects those changes. Edge
    /// handles in `new_specs` must reference the original input
    /// BRep (`self.input`) — they're stable across the operation
    /// because the input is what's being re-filleted.
    ///
    /// # Errors
    /// Same as [`crate::fillet_edges`].
    pub fn re_execute_with_specs(
        &self,
        new_specs: &[FilletSpec],
    ) -> Result<FilletOperationRecord, BlendError> {
        let result =
            crate::fillet_edges(&self.input, new_specs, &self.options, self.tolerance, None)?;
        Ok(FilletOperationRecord {
            input: self.input.clone(),
            specs: new_specs.to_vec(),
            options: self.options.clone(),
            tolerance: self.tolerance,
            result: result.brep,
        })
    }

    /// Convenience: re-run the fillet with all specs replaced by a
    /// single uniform constant radius. Edge handles are preserved
    /// from the original specs. This is the simplest §14.7 entry
    /// point — bind a UI slider's `on_change` callback to it.
    ///
    /// # Errors
    /// Same as [`re_execute_with_specs`](Self::re_execute_with_specs).
    pub fn with_uniform_radius(
        &self,
        new_radius: f64,
    ) -> Result<FilletOperationRecord, BlendError> {
        let new_specs: Vec<FilletSpec> = self
            .specs
            .iter()
            .map(|s| FilletSpec {
                edge: s.edge,
                radius: FilletRadius::Constant(new_radius),
            })
            .collect();
        self.re_execute_with_specs(&new_specs)
    }
}

/// Run a fillet operation and package the input, specs, options,
/// tolerance, and result into a [`FilletOperationRecord`] suitable
/// for an application-level history list.
///
/// Wraps [`crate::fillet_edges`] so the existing validation, phase
/// pipeline, and error taxonomy all apply unchanged. On failure
/// returns the same [`BlendError`] the wrapped call would.
pub fn fillet_recorded(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
) -> Result<FilletOperationRecord, BlendError> {
    let result = crate::fillet_edges(brep, specs, options, tolerance, None)?;
    Ok(FilletOperationRecord {
        input: brep.clone(),
        specs: specs.to_vec(),
        options: options.clone(),
        tolerance,
        result: result.brep,
    })
}

/// A record of one chamfer operation. Same shape as
/// [`FilletOperationRecord`] but carrying chamfer specs.
#[derive(Debug, Clone)]
pub struct ChamferOperationRecord {
    /// Snapshot of the input BRep at the time the chamfer ran.
    pub input: BRep,
    /// The chamfer specs that were applied.
    pub specs: Vec<ChamferSpec>,
    /// The options that were applied.
    pub options: ChamferOptions,
    /// The tolerance the operation ran with.
    pub tolerance: Tolerance,
    /// The result BRep produced by the operation.
    pub result: BRep,
}

impl ChamferOperationRecord {
    /// Input that was passed to the chamfer.
    pub fn undo_state(&self) -> &BRep {
        &self.input
    }

    /// The cached result.
    pub fn redo_state(&self) -> &BRep {
        &self.result
    }

    /// Re-run the chamfer operation deterministically.
    pub fn re_execute(&self) -> Result<BRep, BlendError> {
        let result = crate::chamfer_edges(
            &self.input,
            &self.specs,
            &self.options,
            self.tolerance,
            None,
        )?;
        Ok(result.brep)
    }
}

/// Run a chamfer operation and package it into a
/// [`ChamferOperationRecord`].
pub fn chamfer_recorded(
    brep: &BRep,
    specs: &[ChamferSpec],
    options: &ChamferOptions,
    tolerance: Tolerance,
) -> Result<ChamferOperationRecord, BlendError> {
    let result = crate::chamfer_edges(brep, specs, options, tolerance, None)?;
    Ok(ChamferOperationRecord {
        input: brep.clone(),
        specs: specs.to_vec(),
        options: options.clone(),
        tolerance,
        result: result.brep,
    })
}

// Keep the compiler from complaining that `BlendResult` is unused
// when we build this module in isolation.
#[allow(dead_code)]
fn _uses_blend_result(_: BlendResult) {}

#[cfg(test)]
#[path = "history_tests.rs"]
mod tests;
