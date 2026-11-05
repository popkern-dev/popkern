//! Structured feature-failure payloads (scripting surface spec §10).
//!
//! A kernel error crossing the replay boundary is projected into a
//! [`FeatureFailure`]: a serde-able, `Clone`-able class taxonomy that
//! carries every value a dialog or a script needs (blend kind, per-edge
//! failures with remediation, boolean stage, ...). The kernel error
//! types themselves stay in their crates — they own `Box<BRep>` partial
//! results and handles that must not be persisted — and the projection
//! is the one contract every skin (Rust, Python, Rhai, MCP) presents.
//!
//! `Display` is redundant with the fields: it is the kernel's own
//! message, kept verbatim so existing text-matching diagnostics are
//! unchanged.

use popkern_blend::error::BlendError;
use popkern_bool::BoolError;
use popkern_core::topology::BRep;
use serde::{Deserialize, Serialize};

pub use crate::failure_blend::{
    BlendDeclineKind, BlendDeclined, BlendFailureReason, BlendFailureRecord, RemediationRecord,
};

/// Why a feature failed to replay — one class per spec §10.2 row.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind")]
pub enum FeatureFailure {
    /// A fillet or chamfer declined transactionally.
    BlendDeclined(BlendDeclined),
    /// A boolean failed.
    BooleanFailed {
        /// `union` / `subtract` / `intersect`.
        operation: String,
        /// `BoolError` variant name.
        bool_kind: String,
        /// The kernel's message.
        detail: String,
        /// Merge sub-stage, when the kernel reported one.
        stage: Option<String>,
    },
    /// An extrude, revolve, sweep, loft, or pipe failed.
    SweepFailed {
        /// Operation identifier.
        feature: String,
        /// The kernel's message.
        detail: String,
    },
    /// A shell, thicken, or offset-solid failed.
    ShellFailed {
        /// Operation identifier.
        feature: String,
        /// The kernel's message.
        detail: String,
    },
    /// A draft failed.
    DraftFailed {
        /// The kernel's message.
        detail: String,
    },
    /// The result failed the serving-boundary validation.
    ValidationFailed {
        /// Operation identifier.
        feature: String,
        /// The validator's message.
        detail: String,
    },
    /// A parameter was rejected before any geometry work ran.
    InvalidParameter {
        /// Kernel spec path (e.g. `specs[0].radius`); skins translate
        /// it to the surface name.
        parameter: String,
        /// The rendered value.
        value: String,
        /// Why it was rejected.
        reason: String,
    },
    /// A generative design declined (PK-956): the solver's or the
    /// faceting's typed kind.
    GenerativeDeclined {
        /// The decline's kind (`KeptVolumeExceedsTarget`,
        /// `CellBudgetExceeded`, `DesignBody`, ...).
        generative_kind: String,
        /// The message.
        detail: String,
    },
    /// A piece rule selected no piece of the body (design 9): the
    /// feature is opened for editing rather than serving an empty body.
    PieceRuleMatchedNothing {
        /// Operation identifier (`split`).
        feature: String,
        /// The rule, as recorded.
        rule: String,
        /// How many pieces the body had.
        pieces: usize,
    },
    /// A cut removed no material (bench §4.5): no face of the result
    /// descends from the tool and the body kept every piece, so the
    /// tool missed the body (an `on_face` sketch placed off the face,
    /// a cut extruded away from the material) and the feature is
    /// opened for editing rather than served as a silent no-op.
    CutRemovedNothing {
        /// How many faces the tool body had.
        tool_faces: usize,
    },
    /// A cut sketched on a free plane, ending on no face, reaches more
    /// than one body: which it cuts is not chosen for the user. Sketch
    /// it on a face of the body, or end it on one.
    CutTargetAmbiguous {
        /// How many bodies the tool's box overlaps.
        candidates: usize,
    },
    /// An extrude or revolve up to (or offset from) a face could not
    /// end on it (brep-bench §4.4): the face lies `near`..`far` along
    /// the feature's travel (a distance for an extrude, an angle in
    /// radians for a revolve) across the section, and `reason` says why
    /// that end cannot be served (`curved_end_needs_boss`: only a boss
    /// into the face's own body follows a face the section does not
    /// meet all at once;
    /// `thin_behind_face`: that body is thinner behind the face than
    /// the end's depth; `face_out_of_reach`: a ray from the section
    /// misses the face's surface).
    UpToFaceDeclined {
        /// Why, as a stable kind (see above).
        reason: String,
        /// The nearest distance along the extrude at which the face
        /// lies under the section.
        near: f64,
        /// The farthest.
        far: f64,
    },
    /// The operation was cancelled by the caller.
    Cancelled,
    /// A failure with no richer class (reference resolution details,
    /// missing snapshots, replay preconditions).
    Other {
        /// The message.
        detail: String,
    },
}

impl FeatureFailure {
    /// Project a blend decline. `brep` is the body the blend ran on
    /// (for stable edge names in `failures[]`).
    pub fn blend(error: &BlendError, brep: Option<&BRep>) -> Self {
        FeatureFailure::BlendDeclined(BlendDeclined::project(error, brep))
    }

    /// Project a boolean failure.
    pub fn boolean(operation: &str, error: &BoolError) -> Self {
        let (bool_kind, stage) = match error {
            BoolError::NoIntersection => ("NoIntersection", None),
            BoolError::IntersectionFailed(_) => ("IntersectionFailed", None),
            BoolError::SplitFailed { .. } => ("SplitFailed", None),
            BoolError::ClassificationFailed { .. } => ("ClassificationFailed", None),
            BoolError::MergeFailed { stage, .. } => ("MergeFailed", Some(stage.clone())),
            BoolError::ValidationFailed(_) => ("ValidationFailed", None),
            BoolError::DegenerateInput(_) => ("DegenerateInput", None),
            BoolError::CoincidentFaces { .. } => ("CoincidentFaces", None),
            BoolError::SelfTrimFailed(_) => ("SelfTrimFailed", None),
            BoolError::CoreError(_) => ("CoreError", None),
            BoolError::Cancelled => ("Cancelled", None),
            BoolError::OpenShellNotSupported { .. } => ("OpenShellNotSupported", None),
            BoolError::InterferenceFailed(_) => ("InterferenceFailed", None),
            BoolError::SystemError(_) => ("SystemError", None),
            BoolError::NonManifoldResult { .. } => ("NonManifoldResult", None),
            BoolError::UnsupportedBodyKind { .. } => ("UnsupportedBodyKind", None),
            BoolError::TooManyBodies { .. } => ("TooManyBodies", None),
        };
        FeatureFailure::BooleanFailed {
            operation: operation.to_string(),
            bool_kind: bool_kind.to_string(),
            detail: error.to_string(),
            stage,
        }
    }

    /// A sweep-family failure (extrude, revolve, sweep, loft, pipe).
    pub fn sweep(feature: &str, error: &dyn std::fmt::Display) -> Self {
        FeatureFailure::SweepFailed {
            feature: feature.to_string(),
            detail: error.to_string(),
        }
    }

    /// A shell-family failure (shell, thicken, offset_solid).
    pub fn shell(feature: &str, error: &dyn std::fmt::Display) -> Self {
        FeatureFailure::ShellFailed {
            feature: feature.to_string(),
            detail: error.to_string(),
        }
    }

    /// A draft failure.
    pub fn draft(error: &dyn std::fmt::Display) -> Self {
        FeatureFailure::DraftFailed {
            detail: error.to_string(),
        }
    }

    /// A serving-boundary validation failure.
    pub fn validation(feature: &str, error: &dyn std::fmt::Display) -> Self {
        FeatureFailure::ValidationFailed {
            feature: feature.to_string(),
            detail: error.to_string(),
        }
    }

    /// A failure with only a message.
    pub fn other(detail: impl Into<String>) -> Self {
        FeatureFailure::Other {
            detail: detail.into(),
        }
    }

    /// The spec §10 class name (`kind`).
    pub fn kind(&self) -> &'static str {
        match self {
            FeatureFailure::BlendDeclined(_) => "BlendDeclined",
            FeatureFailure::BooleanFailed { .. } => "BooleanFailed",
            FeatureFailure::SweepFailed { .. } => "SweepFailed",
            FeatureFailure::ShellFailed { .. } => "ShellFailed",
            FeatureFailure::DraftFailed { .. } => "DraftFailed",
            FeatureFailure::ValidationFailed { .. } => "ValidationFailed",
            FeatureFailure::InvalidParameter { .. } => "InvalidParameter",
            FeatureFailure::GenerativeDeclined { .. } => "GenerativeDeclined",
            FeatureFailure::PieceRuleMatchedNothing { .. } => "PieceRuleMatchedNothing",
            FeatureFailure::CutRemovedNothing { .. } => "CutRemovedNothing",
            FeatureFailure::CutTargetAmbiguous { .. } => "CutTargetAmbiguous",
            FeatureFailure::UpToFaceDeclined { .. } => "UpToFaceDeclined",
            FeatureFailure::Cancelled => "Cancelled",
            FeatureFailure::Other { .. } => "Other",
        }
    }

    /// The blend payload, when this is a blend decline.
    pub fn as_blend(&self) -> Option<&BlendDeclined> {
        match self {
            FeatureFailure::BlendDeclined(b) => Some(b),
            _ => None,
        }
    }

    /// Whether the message contains `needle` (text-matching helper for
    /// diagnostics and tests).
    pub fn contains(&self, needle: &str) -> bool {
        self.to_string().contains(needle)
    }
}

impl std::fmt::Display for FeatureFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FeatureFailure::BlendDeclined(b) => f.write_str(&b.message),
            FeatureFailure::BooleanFailed { detail, .. }
            | FeatureFailure::SweepFailed { detail, .. }
            | FeatureFailure::ShellFailed { detail, .. }
            | FeatureFailure::DraftFailed { detail }
            | FeatureFailure::ValidationFailed { detail, .. }
            | FeatureFailure::GenerativeDeclined { detail, .. }
            | FeatureFailure::Other { detail } => f.write_str(detail),
            FeatureFailure::InvalidParameter {
                parameter,
                value,
                reason,
            } => write!(f, "invalid parameter `{parameter}` = {value}: {reason}"),
            FeatureFailure::PieceRuleMatchedNothing {
                feature,
                rule,
                pieces,
            } => write!(
                f,
                "{feature}: the rule {rule} selects none of the body's {pieces} pieces"
            ),
            FeatureFailure::CutRemovedNothing { tool_faces } => write!(
                f,
                "the cut removes no material: none of the tool's {tool_faces} faces reaches the body (check the sketch placement and the extrude direction)"
            ),
            FeatureFailure::CutTargetAmbiguous { candidates } => write!(
                f,
                "the cut reaches {candidates} bodies and names none: sketch it on a face of the body to cut, or end it on one"
            ),
            FeatureFailure::UpToFaceDeclined { reason, near, far } => match reason.as_str() {
                "curved_end_needs_boss" => write!(
                    f,
                    "the section meets the end face from {near} to {far} along its travel, not all at once: only a boss into the face's own body (offset 0) follows it"
                ),
                "thin_behind_face" => write!(
                    f,
                    "the section meets the end face from {near} to {far} along its travel, and the face's body is thinner behind it than that: the boss would come through"
                ),
                _ => write!(
                    f,
                    "the end face is out of reach of the section ({reason}, {near}..{far})"
                ),
            },
            FeatureFailure::Cancelled => f.write_str("operation cancelled by caller"),
        }
    }
}

impl From<String> for FeatureFailure {
    fn from(detail: String) -> Self {
        FeatureFailure::Other { detail }
    }
}

impl From<String> for Box<FeatureFailure> {
    fn from(detail: String) -> Self {
        Box::new(FeatureFailure::Other { detail })
    }
}

impl From<&str> for Box<FeatureFailure> {
    fn from(detail: &str) -> Self {
        Box::new(FeatureFailure::from(detail))
    }
}

impl From<&str> for FeatureFailure {
    fn from(detail: &str) -> Self {
        FeatureFailure::Other {
            detail: detail.to_string(),
        }
    }
}

impl From<crate::error::ReplayError> for FeatureFailure {
    /// The structured payload of a feature failure; other replay
    /// errors (resolution, core, tree) keep their message.
    fn from(e: crate::error::ReplayError) -> Self {
        match e {
            crate::error::ReplayError::FeatureFailed { error, .. } => *error,
            other => FeatureFailure::Other {
                detail: other.to_string(),
            },
        }
    }
}

#[cfg(test)]
#[path = "failure_tests.rs"]
mod failure_tests;
