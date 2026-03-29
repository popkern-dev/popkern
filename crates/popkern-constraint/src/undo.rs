//! Undo/redo support via sketch state snapshots.
//!
//! Each mutation that should be undoable captures a snapshot of the sketch
//! state *before* the mutation. The undo log is a linear sequence of these
//! snapshots with a cursor indicating the current position. Undo restores
//! the snapshot; redo re-applies the mutation by restoring the *next*
//! snapshot (or the post-mutation state).
//!
//! Because [`SlotMap`](slotmap::SlotMap) keys are managed internally and
//! cannot be re-used after removal, a command-based approach cannot
//! correctly restore handles. Snapshot-based undo avoids this limitation.

use crate::constraint::ConstraintEntry;
use crate::elements::{Arc, Ellipse, Hyperbola, Parabola, Point, Spline};
use crate::pattern::{PatternDef, PatternHandle};
use crate::{
    ArcHandle, ConstraintHandle, EllipseHandle, HyperbolaHandle, Line, LineHandle, ParabolaHandle,
    PointHandle, SplineHandle,
};
use slotmap::SlotMap;

/// A complete snapshot of the mutable sketch state, captured before each
/// mutation for undo support.
#[derive(Debug, Clone)]
pub(crate) struct Snapshot {
    /// The parameter vector.
    pub params: Vec<f64>,
    /// All points.
    pub points: SlotMap<PointHandle, Point>,
    /// All lines.
    pub lines: SlotMap<LineHandle, Line>,
    /// All arcs.
    pub arcs: SlotMap<ArcHandle, Arc>,
    /// All ellipses.
    pub ellipses: SlotMap<EllipseHandle, Ellipse>,
    /// All parabolas.
    pub parabolas: SlotMap<ParabolaHandle, Parabola>,
    /// All hyperbolas.
    pub hyperbolas: SlotMap<HyperbolaHandle, Hyperbola>,
    /// All splines.
    pub splines: SlotMap<SplineHandle, Spline>,
    /// All constraints.
    pub constraints: SlotMap<ConstraintHandle, ConstraintEntry>,
    /// All patterns.
    pub patterns: SlotMap<PatternHandle, PatternDef>,
}

#[cfg(test)]
#[path = "undo_tests.rs"]
mod tests;
