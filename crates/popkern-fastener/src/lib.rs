//! Parametric fastener generation for the popkern kernel (PK-460,
//! design: `docs/decisions/pk460-fastener-subsystem.md`): thread
//! specifications, standard dimension tables, and brep-level builders
//! for threaded rods, hex bolts, socket-cap and countersunk screws,
//! hex nuts, tapped/clearance holes, and external threads on existing
//! shanks.
//!
//! One thread mechanism serves every fastener: a **transversal
//! helical groove cut**. External threads subtract a thread-space
//! groove tool from a major-diameter shank; internal threads subtract
//! a clearance-grown thread-ridge tool outward from a minor-diameter
//! bore. Both tools cross the wall strictly transversally, so no
//! tangent or coincident surface pairs arise by construction.
//!
//! This crate sits in the operations tier (peers: blend, sweep) so
//! both the code-first facade (`popkern-part`, whose
//! `fastener` module wraps these builders in `Part` recipes) and the
//! parametric feature tree (`popkern-history` replay) can call it.
//! Everything here takes and returns bare validated [`BRep`]s in the
//! **insertion frame**: fastener axis +Z, origin at the bearing
//! plane, threaded end pointing +Z; receiving-hole cuts take an
//! arbitrary axis whose direction points into the material, so a
//! fastener built on the same axis mates by construction.
//!
//! Dimensions are plain numbers in the caller's model units; the
//! standard-size tables ([`standards`]) are nominal millimetres.
//!
//! [`BRep`]: popkern_core::topology::BRep

pub mod builders;
mod builders_shank;
mod closed_side_census_probe;
#[cfg(test)]
mod deep_plunge_probe_tests;
mod engine;
mod entry_ring;
pub mod holes;
#[cfg(test)]
mod serve_cost_probe_tests;
mod spec;
mod spec_standards;
pub mod standards;
mod thread_form;

#[cfg(test)]
#[path = "builders_tests.rs"]
mod builders_tests;
#[cfg(test)]
#[path = "coil_pierce_tests.rs"]
mod coil_pierce_tests;
#[cfg(test)]
#[path = "engine_tests.rs"]
mod engine_tests;
#[cfg(test)]
#[path = "helicoid_fidelity_tests.rs"]
mod helicoid_fidelity_tests;
#[cfg(test)]
#[path = "spec_tests.rs"]
mod spec_tests;
#[cfg(test)]
#[path = "thread_edge_fidelity_tests.rs"]
mod thread_edge_fidelity_tests;
#[cfg(test)]
#[path = "thread_form_tests.rs"]
mod thread_form_tests;

pub use holes::{ClearanceHoleOptions, HoleDepth, TapOptions};
pub use spec::{Handedness, ThreadSpec};
pub use standards::ClearanceFit;

/// Errors from fastener specification and construction.
///
/// Specification-level variants carry structured numbers a GUI dialog
/// can validate against before building anything; construction
/// variants carry the kernel diagnostics whole with the fastener
/// stage named (`"external thread"`, `"hex nut"`, `"tapped hole"`,
/// ...).
#[derive(Debug, thiserror::Error)]
pub enum FastenerError {
    /// A thread designation string was not understood.
    #[error("unknown thread designation '{designation}': {detail}")]
    UnknownDesignation {
        /// The string as given.
        designation: String,
        /// What failed and what grammars are accepted.
        detail: String,
    },
    /// A spec number is out of range.
    #[error("invalid thread spec: {field} = {value} must be {requirement}")]
    InvalidSpec {
        /// Which field ("major_diameter", "pitch", ...).
        field: &'static str,
        /// The offending value.
        value: f64,
        /// The requirement it violated.
        requirement: &'static str,
    },
    /// The thread profile would overlap itself or consume the crest
    /// flat: adjacent turns of the cut tool (or of the thread
    /// material) would merge. Reports the measured gap so a caller
    /// can shrink clearance, flank angle, or thread depth.
    #[error(
        "thread profile overlap: {detail} (computed gap {gap:.6}, \
         pitch {pitch:.6}, lead {lead:.6})"
    )]
    ProfileOverlap {
        /// Which invariant failed.
        detail: String,
        /// The signed gap that must stay positive.
        gap: f64,
        /// The thread pitch.
        pitch: f64,
        /// The thread lead (pitch × starts).
        lead: f64,
    },
    /// A sweep (tool body, shank, head) declined.
    #[error("fastener {stage}: sweep failed: {source}")]
    Sweep {
        /// The fastener stage being built.
        stage: String,
        /// The kernel error.
        source: popkern_sweep::SweepError,
    },
    /// A boolean combination declined (transactional: the input body
    /// is unchanged).
    #[error("fastener {stage}: boolean failed: {source}")]
    Boolean {
        /// The fastener stage being built.
        stage: String,
        /// The kernel error.
        source: popkern_bool::BoolError,
    },
    /// An edge chamfer declined.
    #[error("fastener {stage}: chamfer failed: {source}")]
    Blend {
        /// The fastener stage being built.
        stage: String,
        /// The kernel error.
        source: popkern_blend::BlendError,
    },
    /// A constructed body failed kernel validation.
    #[error("fastener {stage}: produced an invalid body: {source}")]
    Validation {
        /// The fastener stage being built.
        stage: String,
        /// The kernel validation error.
        source: popkern_core::error::ValidationError,
    },
    /// A construction step outside the kernel error types failed
    /// (profile compile, degenerate selection, empty body).
    #[error("fastener {stage}: {detail}")]
    Construction {
        /// The fastener stage being built.
        stage: String,
        /// What went wrong.
        detail: String,
    },
    /// The entry face is too steep for a coaxial 45° entry cone: past
    /// the bound the cone and the face approach tangency along the
    /// uphill rim (a near-tangent surface-intersection class) and the
    /// cut is no longer a lead-in chamfer. Disable the entry chamfer
    /// or countersink, or reorient the axis.
    #[error(
        "fastener {stage}: entry face tilt {tilt_deg:.1}° exceeds the \
         {limit_deg:.1}° bound for a coaxial 45° entry cone (cone and \
         face go near-tangent along the uphill rim); disable the entry \
         chamfer/countersink or reorient the axis"
    )]
    EntryTooSteep {
        /// The fastener stage being built.
        stage: String,
        /// The measured entry-face tilt, degrees off perpendicular.
        tilt_deg: f64,
        /// The supported bound, degrees.
        limit_deg: f64,
    },
}
