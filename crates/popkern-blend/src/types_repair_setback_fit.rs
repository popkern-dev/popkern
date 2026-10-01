//! The variational corner fit's record (PK-574): why the KP20 cap
//! declined at a setback corner, each fit attempt with its certified
//! seam bounds, and the level served.

use popkern_core::topology::VertexHandle;

use crate::types::ContinuityLevel;

/// The variational fit's mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SetbackFitMode {
    /// G1 and G2 rim rows at the first refinement level; served only
    /// when every cap seam certifies G2 under the commercial gates.
    CurvatureWeighted,
    /// G1 rim rows only; served when every cap seam certifies G1.
    Tangent,
    /// G1 rim rows with dual-jet corner cells (each rim owns its corner
    /// 2-jet, PK-595): the tangent fit's rescue where its shared corner
    /// control point cannot meet both strips (the setback corner
    /// theorem's G1 case on mixed-radius circular strips, PK-615).
    /// Tried only after [`Self::Tangent`] fails; served when every cap
    /// seam certifies G1.
    TangentCornerCells,
}

impl SetbackFitMode {
    /// The level every cap seam must certify for the mode to serve.
    pub fn level(self) -> ContinuityLevel {
        match self {
            SetbackFitMode::CurvatureWeighted => ContinuityLevel::G2,
            SetbackFitMode::Tangent | SetbackFitMode::TangentCornerCells => ContinuityLevel::G1,
        }
    }
}

/// Certified outcome of one fit attempt that reached the certifier.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SetbackFitCertified {
    /// Highest level every cap seam certified under the commercial gates
    /// (`G0` when even G1 failed).
    pub level: ContinuityLevel,
    /// Worst certified upper bounds `(g0, g1, g2)` over the cap↔neighbour
    /// (rim) seams.
    pub rim: [f64; 3],
    /// Worst certified upper bounds `(g0, g1, g2)` over the cap↔cap
    /// (spoke) seams.
    pub spoke: [f64; 3],
    /// Knot-insertion rounds the fit used.
    pub refinement_levels: usize,
    /// Free coefficients in the fit's last solve.
    pub unknowns: usize,
    /// The fit's smallest sampled `|S_u × S_v|` relative to its mean
    /// (the core fold gate declines below 0.02).
    pub min_regularity: f64,
    /// Where the fit split the sides: midpoints, or the spokes' feet
    /// nearest the fitted centre after the midpoint layout bent past the
    /// fairness bound (PK-618).
    pub layout: popkern_core::geometry::variational_cap::CapLayout,
}

/// How one fit attempt ended.
#[derive(Debug, Clone, PartialEq)]
pub enum SetbackFitOutcome {
    /// Certified at the mode's level and committed.
    Served(SetbackFitCertified),
    /// Reached the certifier but was not committed: certified below the
    /// mode's level, or a later postcondition (self-intersection, model
    /// validation) failed. `reason` names which.
    Rejected {
        /// The certified bounds of the rejected cap.
        certified: SetbackFitCertified,
        /// Why it was rejected.
        reason: String,
    },
    /// Declined before certification (a core-fit error such as a fold,
    /// or a failed topology step).
    Declined(String),
}

impl SetbackFitOutcome {
    /// The certified bounds, when the attempt reached the certifier.
    pub fn certified(&self) -> Option<&SetbackFitCertified> {
        match self {
            SetbackFitOutcome::Served(c) | SetbackFitOutcome::Rejected { certified: c, .. } => {
                Some(c)
            }
            SetbackFitOutcome::Declined(_) => None,
        }
    }
}

/// One fit attempt: its mode and how it ended.
#[derive(Debug, Clone, PartialEq)]
pub struct SetbackFitAttempt {
    /// Mode attempted.
    pub mode: SetbackFitMode,
    /// Served, rejected on its certificate, or declined earlier.
    pub outcome: SetbackFitOutcome,
    /// The fitted cap's interior fairness (largest sampled principal
    /// curvature times the corner's smallest radius, PK-614), when the
    /// attempt built a valid cap.
    pub interior: Option<f64>,
}

/// The fallback-rung record for one setback corner.
#[derive(Debug, Clone, PartialEq)]
pub struct SetbackFitDiagnostic {
    /// The apex vertex (valid on the input body).
    pub vertex: VertexHandle,
    /// Why the KP20 cap did not serve (the cap-stage decline).
    pub kp20_decline: String,
    /// Fit attempts in ladder order.
    pub attempts: Vec<SetbackFitAttempt>,
    /// Level served by the fit, `None` when every attempt declined and
    /// the corner fell through to the vertex-blend path.
    pub served: Option<ContinuityLevel>,
    /// The KP20 cap's interior fairness reading, when it exceeded the
    /// bound and the fit was tried for a fairer cap (PK-612).
    pub fairness: Option<SetbackCapFairness>,
}

/// Interior fairness of a corner's candidate caps (PK-612): the largest
/// sampled principal curvature over the cap's interior, times the
/// corner's smallest fillet radius.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SetbackCapFairness {
    /// The KP20 cap's reading.
    pub kp20: f64,
    /// The served fit's reading, when a fit served.
    pub fit: Option<f64>,
    /// The bound the KP20 cap exceeded.
    pub bound: f64,
}
