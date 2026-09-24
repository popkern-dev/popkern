//! PK-148 terminal-arc-trim repair report types.
//!
//! At an open fillet chain end the strip's terminal boundary must be
//! the cross-section arc of the blend surface. The pre-terminal-arc split
//! pipeline instead retained the two support-edge stubs (straight
//! line curves through the original corner vertex) while the strip's
//! pcurves already described the true arc — two disagreeing
//! geometric authorities on the same edges, measured as a
//! `(√2−1)·r` gap between the corner and the arc. The terminal arc
//! trim re-authors the stale side by construction; these types
//! record what it did and why it declined.

use popkern_core::topology::{FaceHandle, VertexHandle};

/// Structured outcome of the terminal-arc-trim settlement pass.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TerminalArcTrimReport {
    /// Chain-end corners whose vertex was repositioned onto the
    /// terminal arc and whose stub edges received arc curves and
    /// coherent cap-side pcurves.
    pub repaired_corners: usize,
    /// Candidate corners the pass recognized but declined to repair,
    /// with the structural reason. A decline leaves the stale
    /// geometry in place (pre-PK-148 behaviour) — it never damages
    /// the body.
    pub declines: Vec<TerminalArcTrimDecline>,
}

/// One declined terminal-arc-trim candidate.
#[derive(Debug, Clone, PartialEq)]
pub struct TerminalArcTrimDecline {
    /// Strip face whose loop carried the candidate corner.
    pub strip_face: FaceHandle,
    /// The stale corner vertex.
    pub vertex: VertexHandle,
    /// Why the repair did not apply.
    pub reason: TerminalArcTrimDeclineReason,
}

/// Structural reasons a recognized candidate is not repaired.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TerminalArcTrimDeclineReason {
    /// The strip surface is not an analytic cylinder; the terminal
    /// cross-section is not a circular arc and needs the sampled
    /// iso-curve construction (PK-139 stage 2 territory).
    NonCylinderStrip,
    /// The stub pair's twin half-edges do not bound one shared
    /// planar cap face, so there is no single end face to trim.
    NoCommonPlanarCap,
    /// The corner vertex is used by more edges than the two stubs;
    /// repositioning it would drag unrelated topology.
    VertexNotIsolated,
    /// The anchor vertices (the stubs' non-corner endpoints) are not
    /// on the strip cylinder's section by the end plane (a circle for
    /// a perpendicular end, an ellipse for an oblique one, PK-530).
    /// The trusted geometric inputs the repair is built from are then
    /// themselves unreliable, so nothing is authored.
    AnchorsOffSurface,
    /// The arc sweep could not be branch-resolved against the strip
    /// pcurve midpoint witness within tolerance.
    ArcBranchUnresolved,
    /// One anchor sits off the cylinder because an adjacent face
    /// rises inside the fillet's terminal footprint (support face
    /// wraps the obstruction), and the obstructed-end surgery
    /// (`terminal_obstructed_end`) could not establish or commit its
    /// supported single-crossing configuration.
    ObstructedEnd,
    /// The strip ends into a fan of planar faces (a slant edge running
    /// into a pyramid apex), and the fan-end settlement could not
    /// split the corner or re-author its arcs, or the settled body
    /// failed validation. The body was restored (PK-530).
    FanEndRejected,
}
