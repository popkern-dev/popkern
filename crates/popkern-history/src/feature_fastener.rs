//! Parameter types for the fastener features (PK-460):
//! [`Feature::Fastener`](crate::feature::Feature::Fastener) bodies
//! and the
//! [`TappedHole`](crate::feature::Feature::TappedHole) /
//! [`ClearanceHole`](crate::feature::Feature::ClearanceHole) cuts.
//!
//! Threads are stored as plain resolved numbers (major diameter +
//! pitch through the parameter table), not designation strings: a
//! designation like `"M6"` is a skin-level convenience that binds
//! these parameters at record time, so the thread stays sliderable
//! afterwards.

use crate::parameter::ParameterId;
use serde::{Deserialize, Serialize};

/// Which fastener body a `Feature::Fastener` builds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FastenerKind {
    /// Hex-head bolt (ISO 4017 proportions by default).
    HexBolt,
    /// Socket-head cap screw (ISO 4762 proportions by default).
    SocketCapScrew,
    /// Countersunk flat-head screw (ISO 10642, 90°).
    CountersunkScrew,
    /// Plain threaded rod (no head).
    ThreadedRod,
    /// Hex nut with a modeled internal thread (ISO 4032 proportions).
    HexNut,
}

impl FastenerKind {
    /// The skin spelling of this kind (`hex | socket | countersunk |
    /// rod | nut`).
    pub fn spelling(self) -> &'static str {
        match self {
            FastenerKind::HexBolt => "hex",
            FastenerKind::SocketCapScrew => "socket",
            FastenerKind::CountersunkScrew => "countersunk",
            FastenerKind::ThreadedRod => "rod",
            FastenerKind::HexNut => "nut",
        }
    }
}

/// A thread specification as stored in the feature tree: resolved
/// dimensions (parameter table), winding, and start count.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ThreadDesc {
    /// Thread major (nominal) diameter.
    pub major_diameter: ParameterId,
    /// Thread pitch (axial crest-to-crest of adjacent turns).
    pub pitch: ParameterId,
    /// Winding sense. `false` = right-hand.
    #[serde(default)]
    pub left_handed: bool,
    /// Number of thread starts (≥ 1; the lead is `pitch × starts`).
    #[serde(default = "one_start")]
    pub starts: u32,
}

fn one_start() -> u32 {
    1
}

/// serde default for on-by-default flags.
pub(crate) fn default_true() -> bool {
    true
}

/// How a clearance hole is sized: an ISO 273 series or an explicit
/// diameter parameter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ClearanceFitDesc {
    /// ISO 273 close series.
    Close,
    /// ISO 273 normal (medium) series.
    Normal,
    /// ISO 273 loose (coarse) series.
    Loose,
    /// An explicit hole diameter.
    Custom(ParameterId),
}
