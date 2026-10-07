//! Parametric hinges for the popkern kernel (PK-693, design:
//! `docs/architecture/hinges.md`): pinned butt hinges (a piano hinge is
//! one with many knuckles), print-in-place hinges with the pin integral
//! to one leaf, and living hinges.
//!
//! Everything builds in the **hinge frame**: pin axis +Z from `z = 0`
//! to `z = length`; leaf A extends toward −X and leaf B toward +X at an
//! open angle of 180° (flat); both mounting faces lie on the plane
//! `y = −R` (`R` the knuckle radius) and the front faces at
//! `y = −R + thickness`. An open angle θ turns leaf B about +Z by
//! (180° − θ), so θ = 0 folds the front faces together.
//!
//! Construction is transversal by design: the knuckle and its leaf are
//! one line–arc–line profile extruded along the axis (the leaf face
//! meets the knuckle cylinder along an edge the sketch builds, never a
//! tangency the boolean must find); the other leaf's knuckle zones go
//! in ONE comb-shaped cut whatever the knuckle count; bores and pins
//! are coaxial cylinders whose ends never sit flush with a knuckle
//! face. Clearances are explicit parameters, never tolerances.

mod construct;
mod leaf;
mod params;

#[cfg(test)]
#[path = "leaf_tests.rs"]
mod leaf_tests;
#[cfg(test)]
#[path = "params_tests.rs"]
mod params_tests;

pub use leaf::hinge_component;
pub use params::{
    HingeComponent, HingeHoles, HingeKind, HingeLayout, HingeParams, HolePosition, hinge_layout,
};

/// Why a hinge cannot be specified or built.
///
/// Specification variants carry the numbers a dialog needs to point at
/// the field and its limit before anything is built; construction
/// variants carry the kernel diagnostics whole, with the hinge stage
/// named (`"leaf A profile"`, `"comb notch"`, `"bore"`, …).
#[derive(Debug, thiserror::Error)]
pub enum HingeError {
    /// A number is out of range.
    #[error("invalid hinge spec: {field} = {value} must be {requirement}")]
    InvalidSpec {
        /// Which field (`"length"`, `"knuckles"`, …).
        field: &'static str,
        /// The offending value.
        value: f64,
        /// The requirement it violated.
        requirement: &'static str,
    },
    /// The knuckles and the gaps between them do not fit the length.
    #[error(
        "{knuckles} knuckles with {gap} gaps leave {segment} per knuckle on a {length} hinge: \
         use fewer knuckles, smaller gaps or a longer hinge"
    )]
    KnucklesDoNotFit {
        /// The knuckle count.
        knuckles: u32,
        /// The axial gap between knuckles.
        gap: f64,
        /// The hinge length.
        length: f64,
        /// The knuckle length that would result (≤ 0 or no longer than
        /// a gap).
        segment: f64,
    },
    /// The knuckle has no wall around the pin bore.
    #[error(
        "the knuckle (radius {knuckle_radius}) has no wall around the {bore_radius} bore: \
         raise the knuckle diameter or lower the pin diameter or clearance"
    )]
    NoKnuckleWall {
        /// The knuckle radius.
        knuckle_radius: f64,
        /// The bore radius (pin radius + clearance).
        bore_radius: f64,
    },
    /// The leaf ends before it clears the other leaf's knuckles.
    #[error(
        "leaf width {leaf_width} leaves no plate beyond the knuckles: it needs more than {needed}"
    )]
    PlateTooNarrow {
        /// The leaf width (axis to outer edge).
        leaf_width: f64,
        /// The width the knuckle and its clearance take.
        needed: f64,
    },
    /// A screw hole does not fit the leaf.
    #[error("screw hole {hole} needs {needed} of {dimension}, the leaf has {available}")]
    HoleDoesNotFit {
        /// The hole (1-based along the leaf).
        hole: u32,
        /// What does not fit (`"plate width"`, `"spacing along the
        /// leaf"`, `"thickness for the countersink"`).
        dimension: &'static str,
        /// What the hole needs.
        needed: f64,
        /// What the leaf has.
        available: f64,
    },
    /// A print-in-place hinge needs real gaps or it prints fused.
    #[error(
        "a print-in-place hinge needs {field} > 0 (the leaves would print fused); it is {value}"
    )]
    PrintClearance {
        /// `"clearance"` or `"gap"`.
        field: &'static str,
        /// The value given.
        value: f64,
    },
    /// A boolean combination declined (the input body is unchanged).
    #[error("hinge {stage}: boolean failed: {source}")]
    Boolean {
        /// The hinge stage being built.
        stage: String,
        /// The kernel error.
        source: Box<popkern_bool::BoolError>,
    },
    /// A constructed body failed kernel validation.
    #[error("hinge {stage}: produced an invalid body: {source}")]
    Validation {
        /// The hinge stage being built.
        stage: String,
        /// The kernel validation error.
        source: popkern_core::error::ValidationError,
    },
    /// A construction step outside the kernel error types failed
    /// (profile compile, extrude, placement, empty body).
    #[error("hinge {stage}: {detail}")]
    Construction {
        /// The hinge stage being built.
        stage: String,
        /// What went wrong.
        detail: String,
    },
    /// A screw hole cut declined.
    #[error("hinge screw hole: {0}")]
    Fastener(Box<popkern_fastener::FastenerError>),
}

impl From<popkern_fastener::FastenerError> for HingeError {
    fn from(e: popkern_fastener::FastenerError) -> HingeError {
        HingeError::Fastener(Box::new(e))
    }
}
