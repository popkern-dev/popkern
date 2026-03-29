//! Pattern data model for linear, circular, and mirror patterns.
//!
//! A pattern replicates a set of seed elements (and their internal
//! constraints) across multiple instances with geometric offset constraints
//! tying each instance back to the seed.

use crate::constraint::ConstraintHandle;
use crate::elements::{
    ArcHandle, ElementHandle, EllipseHandle, HyperbolaHandle, LineHandle, ParabolaHandle,
    PointHandle, SplineHandle,
};
use serde::{Deserialize, Serialize};

slotmap::new_key_type! {
    /// Opaque handle to a [`PatternDef`] inside a [`Sketch`](crate::Sketch).
    pub struct PatternHandle;
}

/// The type of geometric replication a pattern performs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum PatternKind {
    /// Replicate along a direction vector with uniform spacing.
    Linear {
        /// Unit (or non-unit) direction vector `(dx, dy)`.
        direction: (f64, f64),
        /// Distance between consecutive instances along `direction`.
        spacing: f64,
    },
    /// Replicate around a centre point at uniform angular intervals.
    Circular {
        /// Centre of rotation.
        center: PointHandle,
        /// Total angle span (radians) distributed across instances.
        angle: f64,
    },
    /// Reflect across a line (always exactly one instance).
    Mirror {
        /// The mirror line.
        line: LineHandle,
    },
}

/// A complete pattern definition stored in the sketch.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PatternDef {
    /// What kind of pattern this is.
    pub kind: PatternKind,
    /// The original seed elements that are replicated.
    pub seed_elements: Vec<ElementHandle>,
    /// Constraints that are entirely internal to the seed set.
    pub seed_internal_constraints: Vec<ConstraintHandle>,
    /// The replicated instances (one per copy beyond the seed).
    pub instances: Vec<PatternInstance>,
}

/// One replicated copy of the seed elements within a pattern.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PatternInstance {
    /// Cloned elements corresponding to seed_elements (same order).
    pub elements: Vec<ElementHandle>,
    /// Phantom points created for arc/conic centre pinning.
    pub phantom_points: Vec<PointHandle>,
    /// Cloned copies of the seed's internal constraints.
    pub internal_constraints: Vec<ConstraintHandle>,
    /// Offset constraints tying this instance to the seed
    /// (RelativePosition, RotatedPosition, ReflectedPosition).
    pub offset_constraints: Vec<ConstraintHandle>,
    /// PinToCenter constraints for arc/conic phantom points.
    pub pin_constraints: Vec<ConstraintHandle>,
    /// Angle offset/mirror constraints for arcs and conics.
    pub angle_constraints: Vec<ConstraintHandle>,
    /// Whether this instance is suppressed (constraints set to Driven).
    pub suppressed: bool,
}

/// Remapping tables produced by element cloning.
///
/// Used to translate handles from seed space into instance space when
/// cloning constraints.
#[derive(Debug, Clone, Default)]
pub(crate) struct HandleRemap {
    /// Point handle remapping (seed -> instance).
    pub points: std::collections::HashMap<PointHandle, PointHandle>,
    /// Line handle remapping.
    pub lines: std::collections::HashMap<LineHandle, LineHandle>,
    /// Arc handle remapping.
    pub arcs: std::collections::HashMap<ArcHandle, ArcHandle>,
    /// Ellipse handle remapping.
    pub ellipses: std::collections::HashMap<EllipseHandle, EllipseHandle>,
    /// Parabola handle remapping.
    pub parabolas: std::collections::HashMap<ParabolaHandle, ParabolaHandle>,
    /// Hyperbola handle remapping.
    pub hyperbolas: std::collections::HashMap<HyperbolaHandle, HyperbolaHandle>,
    /// Spline handle remapping.
    pub splines: std::collections::HashMap<SplineHandle, SplineHandle>,
}
