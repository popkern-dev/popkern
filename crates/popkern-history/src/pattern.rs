//! Types for feature patterns (linear, circular, curve-driven, mirror).
//!
//! A pattern takes a seed feature and replicates it across a set of
//! positions/orientations.  Each instance can be individually suppressed or
//! have parameter values overridden.

use crate::feature::Vec3Desc;
use crate::parameter::ParameterId;
use popkern_core::naming::EntityName;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Describes one direction axis of a linear pattern.
///
/// Stores the direction vector and the two quantities that jointly control
/// instance spacing: total count and total distance from the seed to the last
/// instance.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PatternDirection {
    /// The axis along which instances are placed.
    pub direction: Vec3Desc,
    /// Number of instances (including the seed) along this axis.
    pub count: ParameterId,
    /// Total distance from the seed to the furthest instance along this axis.
    pub total_distance: ParameterId,
}

/// Per-instance customisation applied on top of the seed feature's parameters.
///
/// If `suppressed` is `true` the instance is skipped entirely during replay.
/// Any entries in `param_overrides` replace the corresponding seed parameter
/// value for that instance only.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InstanceOverride {
    /// Whether this instance is suppressed (skipped) during replay.
    pub suppressed: bool,
    /// Parameter values that differ from the seed for this instance.
    pub param_overrides: BTreeMap<ParameterId, f64>,
}

/// Result of replaying a single pattern instance.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum InstanceOutcome {
    /// The instance was replayed successfully via the normal parametric path.
    Replayed,
    /// The parametric replay failed; a geometry fallback (e.g. copy-transform)
    /// was used instead.
    GeometryFallback,
    /// The instance was suppressed and intentionally skipped.
    Skipped,
    /// The instance could not be produced and no fallback was available.
    Failed {
        /// Human-readable description of the failure.
        error: String,
    },
}

/// Aggregate result of replaying all instances in a pattern.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum PatternResult {
    /// Every instance succeeded.
    AllSucceeded {
        /// Total number of instances that were produced.
        count: usize,
        /// Non-fatal warnings emitted during transform computation.
        #[serde(default)]
        warnings: Vec<String>,
    },
    /// Some instances succeeded and some failed.
    PartialSuccess {
        /// Instances that were produced, paired with their outcome.
        succeeded: Vec<(usize, InstanceOutcome)>,
        /// Instances that failed, paired with an error message.
        failed: Vec<(usize, String)>,
        /// Non-fatal warnings emitted during transform computation.
        #[serde(default)]
        warnings: Vec<String>,
    },
    /// Every instance failed; no geometry was produced.
    AllFailed {
        /// Per-instance error messages.
        errors: Vec<(usize, String)>,
    },
}

/// Describes the mirror plane used by a mirror-pattern feature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MirrorPlaneDesc {
    /// The global XY plane (Z = 0, normal = +Z).
    XY,
    /// The global XZ plane (Y = 0, normal = +Y).
    XZ,
    /// The global YZ plane (X = 0, normal = +X).
    YZ,
    /// A planar face in the model, referenced by stable name.
    OnFace(EntityName),
    /// An arbitrary plane defined by an explicit origin and unit normal.
    Custom {
        /// Plane origin in world space.
        origin: [f64; 3],
        /// Plane unit normal in world space.
        normal: [f64; 3],
    },
}
