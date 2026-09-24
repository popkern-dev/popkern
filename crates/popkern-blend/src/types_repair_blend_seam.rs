//! PK-152 blend×blend corner-seam repair report types.
//!
//! When two blend strips on adjacent edges meet at a shared body
//! corner (the third edge staying sharp), the independent per-edge
//! constructions leave a stale corner tangle: a phantom flap face on
//! the support plane plus stub edges running through the original
//! corner at `(√2−1)·r` off both strips. The blend-seam-corner
//! settlement trims both strips to their mutual intersection curve
//! and retires the flap; these types record what it did and why it
//! declined. Design:
//! `docs/architecture/blend-blend-intersection-trimming.md`.

use popkern_core::topology::{FaceHandle, VertexHandle};

/// Structured outcome of the blend-seam-corner settlement pass.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BlendSeamCornerReport {
    /// Pair corners closed with an exact strip×strip seam curve
    /// (flap retired, stub chain dissolved, seam edge re-authored).
    pub seams_built: usize,
    /// Recognized pair corners the pass declined to repair, with the
    /// structural reason. A decline leaves the stale tangle in place
    /// — it never damages the body.
    pub declines: Vec<BlendSeamCornerDecline>,
}

/// One declined blend-seam-corner candidate.
#[derive(Debug, Clone, PartialEq)]
pub struct BlendSeamCornerDecline {
    /// The two blend strip faces meeting at the corner.
    pub strips: (FaceHandle, FaceHandle),
    /// The stale corner vertex at the tangle's centre.
    pub corner: VertexHandle,
    /// Why the repair did not apply.
    pub reason: BlendSeamCornerDeclineReason,
}

/// Structural reasons a recognized pair corner is not repaired.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlendSeamCornerDeclineReason {
    /// The strips carry different radii and the marching-intersector
    /// preconditions did not hold (retained as the structured
    /// boundary for corner shapes the marching path cannot serve).
    MismatchedRadii,
    /// The near anchor sits on one strip only, but its mismatch stub
    /// edge (toward the larger strip's terminal arc) or its rail
    /// edge could not be identified structurally.
    StubUnresolved,
    /// The marching intersector could not carry the seam from the
    /// far anchor to the support-plane exit (separating strips,
    /// missed exit, or step exhaustion).
    MarchUnresolved,
    /// The strip axes do not intersect within tolerance; the
    /// intersection leaves the analytic bisector-plane family
    /// (stage 2 territory).
    OffsetAxes,
    /// The strip axes are parallel — the parallel-interference
    /// class (PK-140 stage 3b), not a corner seam.
    ParallelAxes,
    /// Neither (or both) bisector branch planes contain the two
    /// trusted seam anchors, so the seam plane cannot be resolved.
    BranchUnresolved,
    /// A trusted seam anchor does not lie on the analytic
    /// intersection ellipse within tolerance.
    AnchorsOffSeam,
    /// The seam sweep between the anchors could not be resolved to
    /// the short branch nearest the stale corner.
    SweepUnresolved,
    /// The seam pcurve did not close against the neighbouring chart
    /// anchors or its sampled trace left the seam locus.
    PcurveClosureFailed,
    /// An Euler operator rejected the flap/chain surgery; the body
    /// was restored from the pre-corner snapshot.
    EulerRejected,
    /// The repaired body failed validation; the body was restored
    /// from the pre-corner snapshot.
    ValidationFailed,
}

/// Structured outcome of the mitre-corner settlement pass (PK-522,
/// PK-528): corners of three or more strips closed by trimming every
/// adjacent strip pair to its mutual intersection, all seams meeting
/// at the strips' common mitre point. Equal-radius cylinders get exact
/// ellipse seams; any other strips get marched ones.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MitreCornerReport {
    /// Corners closed with strip×strip seams meeting at one vertex
    /// (support-plane flaps retired).
    pub corners_built: usize,
    /// Unresolved corners the pass declined, with the reason. A
    /// decline restores the corner's pre-pass state; the
    /// corner-resolution gate then decides whether the operation
    /// serves.
    pub declines: Vec<MitreCornerDecline>,
}

/// One declined mitre corner.
#[derive(Debug, Clone, PartialEq)]
pub struct MitreCornerDecline {
    /// The corner vertex.
    pub corner: VertexHandle,
    /// Why the mitre was not built.
    pub reason: MitreCornerDeclineReason,
}

/// Structural reasons an unresolved corner is not mitred.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MitreCornerDeclineReason {
    /// The face-trim residue at the corner does not reduce to a
    /// polygon of rail crossings with one strip across each side
    /// (merged flaps bordered by anything but strip runs).
    UnrecognizedResidue,
    /// No point near the corner lies on every strip (for example an
    /// irregular corner of four or more strips): there is no mitre
    /// point for the seams to meet at.
    NoCommonPoint,
    /// A strip-pair seam from its rail crossing to the mitre point
    /// could not be authored: exact frame unresolved, a marched
    /// sample off either strip, or a fit off its surfaces.
    SeamUnresolved(BlendSeamCornerDeclineReason),
    /// An Euler operator rejected the flap surgery; the corner was
    /// restored.
    EulerRejected,
    /// The mitred body failed a seam or validation postcondition;
    /// the corner was restored.
    ValidationFailed,
}

impl std::fmt::Display for MitreCornerDeclineReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnrecognizedResidue => write!(
                f,
                "corner residue does not reduce to a rail-crossing polygon"
            ),
            Self::NoCommonPoint => write!(f, "no point lies on every strip"),
            Self::SeamUnresolved(reason) => write!(f, "strip-pair seam unresolved ({reason:?})"),
            Self::EulerRejected => write!(f, "flap surgery rejected"),
            Self::ValidationFailed => write!(f, "mitred corner failed its postconditions"),
        }
    }
}
