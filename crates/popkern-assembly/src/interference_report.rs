//! What an interference check is asked and what it reports (split
//! from `interference.rs` under the 500-line budget).

use serde::Serialize;

/// What to check.
#[derive(Debug, Clone)]
pub struct InterferenceOptions {
    /// Report pairs whose distance may be below this (0 reports only
    /// touching and interfering pairs).
    pub clearance: f64,
    /// Distances at or below this are contact; `None` takes the
    /// largest build tolerance among the parts.
    pub tolerance: Option<f64>,
    /// Tessellation chord height as a fraction of each part's size.
    pub chord_fraction: f64,
    /// Check only occurrences whose path is one of these or lies
    /// below one (`"rotor_L"` takes `rotor_L/shaft`); empty checks
    /// every occurrence. A drive-train sweep checks the parts that
    /// move against each other, not the whole machine at every step.
    pub among: Vec<String>,
    /// Build the overlap of interfering pairs to report its volume.
    /// Without it, a pair proven overlapping from its meshes alone
    /// (`mesh_bvh::certain_overlap`) is reported interfering with no
    /// volume — a swing needs the contact, not its size.
    pub overlap_volume: bool,
}

impl Default for InterferenceOptions {
    fn default() -> Self {
        InterferenceOptions {
            clearance: 0.0,
            tolerance: None,
            chord_fraction: 1e-3,
            among: Vec::new(),
            overlap_volume: true,
        }
    }
}

/// How a pair of occurrences relate.
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(tag = "status", rename_all = "snake_case")]
pub enum PairStatus {
    /// Apart: `distance` (upper bound when `refined`) and a certified
    /// `lower_bound` above the contact tolerance.
    Clear {
        /// The reported distance.
        distance: f64,
        /// Certified lower bound on the true distance.
        lower_bound: f64,
        /// Whether `distance` was refined onto the true faces.
        refined: bool,
    },
    /// In contact (distance within the tolerance), no overlap volume.
    Touching {
        /// The reported distance.
        distance: f64,
    },
    /// Overlapping solids.
    Interfering {
        /// Volume of the overlap, when the boolean built it.
        volume: Option<f64>,
        /// Why the boolean declined, when it did.
        decline: Option<String>,
    },
}

/// One reported pair.
#[derive(Debug, Clone, Serialize)]
pub struct PairReport {
    /// First occurrence path.
    pub a: String,
    /// Second occurrence path.
    pub b: String,
    /// The relation.
    #[serde(flatten)]
    pub status: PairStatus,
    /// The closest points (world), when apart or touching.
    pub closest: Option<[[f64; 3]; 2]>,
    /// The faces those points lie on, as entity paths.
    pub faces: Option<[String; 2]>,
}

/// Every reported pair, worst first (interfering, touching, then by
/// distance).
#[derive(Debug, Clone, Serialize)]
pub struct InterferenceReport {
    /// The pairs within the clearance.
    pub pairs: Vec<PairReport>,
    /// Part occurrences checked.
    pub occurrences: usize,
    /// Pairs whose bounds came within the clearance (each measured).
    pub candidate_pairs: usize,
    /// The contact tolerance used.
    pub tolerance: f64,
    /// The clearance asked for.
    pub clearance: f64,
}
