//! Configuration data model for named sketch states.
//!
//! A [`Configuration`] captures a snapshot of each constraint's suppressed
//! state and dimensional value, plus each group's enabled state, so the
//! sketch can be quickly toggled between different design variants.

use crate::constraint::{Constraint, ConstraintHandle};
use crate::groups::GroupHandle;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

slotmap::new_key_type! {
    /// Opaque handle to a [`Configuration`] inside a [`Sketch`](crate::Sketch).
    pub struct ConfigHandle;
}

/// Snapshot of a single constraint's state for a named configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConstraintSnapshot {
    /// When `true` the constraint is treated as [`ConstraintKind::Driven`]
    /// (display-only) in this configuration.
    pub suppressed: bool,
    /// Dimensional value, if this constraint carries one.
    pub value: Option<f64>,
}

/// A named design configuration: a set of per-constraint state overrides and
/// per-group enabled flags, plus an optional warm-start parameter vector.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Configuration {
    /// Human-readable name for this configuration.
    pub name: String,
    /// Per-constraint state snapshots, keyed by handle.
    pub constraint_states: HashMap<ConstraintHandle, ConstraintSnapshot>,
    /// Per-group enabled flags, keyed by handle.
    pub group_states: HashMap<GroupHandle, bool>,
    /// Solved parameter vector from when the configuration was saved, used as
    /// a warm-start when activating.
    pub solved_params: Option<Vec<f64>>,
}

/// Extract the dimensional value from a constraint variant, if it has one.
pub(crate) fn get_constraint_value(constraint: &Constraint) -> Option<f64> {
    match constraint {
        Constraint::Distance(_, _, v) => Some(*v),
        Constraint::PointLineDistance(_, _, v) => Some(*v),
        Constraint::Angle(_, _, v) => Some(*v),
        Constraint::FocalDistance(_, v) => Some(*v),
        Constraint::EllipseEccentricity(_, v) => Some(*v),
        Constraint::HyperbolaEccentricity(_, v) => Some(*v),
        Constraint::OffsetLine(_, _, v) => Some(*v),
        Constraint::OffsetArc(_, _, v) => Some(*v),
        Constraint::OrdinateX(_, v) => Some(*v),
        Constraint::OrdinateY(_, v) => Some(*v),
        Constraint::SplineLength(_, v) => Some(*v),
        Constraint::RelativePosition(_, _, dx, _) => Some(*dx),
        _ => None,
    }
}

/// Set the dimensional value on a constraint variant in place.
///
/// Returns `true` if the value was written, `false` if the variant does not
/// carry a scalar dimension.
pub(crate) fn set_constraint_value_on(constraint: &mut Constraint, value: f64) -> bool {
    match constraint {
        Constraint::Distance(_, _, v) => {
            *v = value;
            true
        }
        Constraint::PointLineDistance(_, _, v) => {
            *v = value;
            true
        }
        Constraint::Angle(_, _, v) => {
            *v = value;
            true
        }
        Constraint::FocalDistance(_, v) => {
            *v = value;
            true
        }
        Constraint::EllipseEccentricity(_, v) => {
            *v = value;
            true
        }
        Constraint::HyperbolaEccentricity(_, v) => {
            *v = value;
            true
        }
        Constraint::OffsetLine(_, _, v) => {
            *v = value;
            true
        }
        Constraint::OffsetArc(_, _, v) => {
            *v = value;
            true
        }
        Constraint::OrdinateX(_, v) => {
            *v = value;
            true
        }
        Constraint::OrdinateY(_, v) => {
            *v = value;
            true
        }
        Constraint::SplineLength(_, v) => {
            *v = value;
            true
        }
        Constraint::RelativePosition(_, _, dx, _) => {
            *dx = value;
            true
        }
        _ => false,
    }
}
