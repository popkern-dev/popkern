//! Error types for popkern-core.

use crate::topology::{
    EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle, RegionHandle, ShellHandle, VertexHandle,
};
use thiserror::Error;

/// How a region breaks the region structure
/// ([`ValidationError::RegionStructure`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
pub enum RegionDefect {
    /// The body's outside is not a live void region.
    #[error("the outside is missing or not a void region")]
    OutsideNotVoid,
    /// The outside has a peripheral shell: nothing encloses it.
    #[error("the outside has an outer shell")]
    OutsideBounded,
    /// A region other than the outside has no live outer shell of its
    /// own.
    #[error("no live outer shell of its own")]
    NoOuterShell,
    /// A shell names this region but the region is not live.
    #[error("a shell names it but it is not live")]
    MissingRegion,
}

/// Errors from geometric/topological validation checks.
#[derive(Debug, Error)]
pub enum ValidationError {
    /// A topology entity references another that no longer exists or
    /// does not own it (loop→face, half-edge→loop).
    #[error("dangling topology reference: {what}")]
    DanglingReference {
        /// Description of the broken reference.
        what: String,
    },
    /// Euler-Poincaré invariant violated.
    #[error("Euler-Poincaré violation on solid {solid:?}: expected {expected}, got {actual}")]
    EulerPoincare {
        /// The solid that failed.
        solid: RegionHandle,
        /// Expected value.
        expected: i64,
        /// Actual computed value.
        actual: i64,
    },

    /// Half-edge has no valid twin.
    #[error("Half-edge {half_edge:?} has no valid twin")]
    BrokenTwin {
        /// The half-edge with a broken twin.
        half_edge: HalfEdgeHandle,
    },

    /// Half-edge cycle is broken (following next doesn't return to start).
    #[error("Half-edge cycle broken at {half_edge:?} in loop {loop_:?}")]
    BrokenCycle {
        /// The half-edge where the cycle breaks.
        half_edge: HalfEdgeHandle,
        /// The loop containing the broken cycle.
        loop_: LoopHandle,
    },

    /// A loop's half-edge ends at a different vertex than its
    /// successor starts at: following `next` returns to the start, so
    /// the cycle check passes, but the loop jumps between vertices
    /// (measured: graph-merge served a sphere face whose outer loop
    /// jumped at three of its six junctions, 0.4 apart, through every
    /// structural check; PK-940, PK-943).
    #[error(
        "Half-edge {half_edge:?} in loop {loop_:?} ends at {end:?} but its successor \
         {next:?} starts at {next_origin:?}"
    )]
    LoopJunctionMismatch {
        /// The half-edge whose end misses its successor.
        half_edge: HalfEdgeHandle,
        /// The successor (`next`) half-edge.
        next: HalfEdgeHandle,
        /// The loop holding both.
        loop_: LoopHandle,
        /// Where `half_edge` ends (its twin's origin).
        end: VertexHandle,
        /// Where `next` starts.
        next_origin: VertexHandle,
    },

    /// Vertex's outgoing half-edge doesn't originate from it.
    #[error("Vertex {vertex:?} outgoing half-edge {half_edge:?} does not originate from it")]
    InconsistentVertex {
        /// The vertex with inconsistent outgoing half-edge.
        vertex: VertexHandle,
        /// The half-edge that should originate from the vertex.
        half_edge: HalfEdgeHandle,
    },

    /// Edge is non-manifold (not exactly 2 half-edges).
    #[error("Edge {edge:?} is non-manifold")]
    NonManifoldEdge {
        /// The non-manifold edge.
        edge: EdgeHandle,
    },

    /// Vertex has valence less than 2.
    #[error("Vertex {vertex:?} has valence {valence} (expected >= 2)")]
    LowValence {
        /// The low-valence vertex.
        vertex: VertexHandle,
        /// The actual valence.
        valence: usize,
    },

    /// Face normal orientation is inconsistent.
    #[error("Face {face:?} has inconsistent normal orientation")]
    InconsistentOrientation {
        /// The face with inconsistent orientation.
        face: FaceHandle,
    },

    /// Loop boundary is not geometrically closed.
    #[error("Loop {loop_:?} is not closed (gap: {gap})")]
    OpenLoop {
        /// The open loop.
        loop_: LoopHandle,
        /// The geometric gap distance.
        gap: f64,
    },

    /// A closed edge's canonical half-edge traverses the stored
    /// curve backward (the PK-303 canonical-forward convention,
    /// audited by `audit_canonical_forward` — PK-478).
    #[error("Edge {edge:?}: canonical half-edge runs the stored curve backward")]
    CanonicalBackward {
        /// The violating closed edge.
        edge: EdgeHandle,
    },

    /// An edge's curve window leaves the angular span its Arc or
    /// Ellipse declares (PK-906). Evaluation is periodic, so the edge
    /// still traces correctly, but projection and parameter-range
    /// queries answer for the declared span.
    #[error(
        "Edge {edge:?}: window [{t_start}, {t_end}] outside its curve's span \
         [{span_start}, {span_end}]"
    )]
    CurveWindowOutsideSpan {
        /// The violating edge.
        edge: EdgeHandle,
        /// The edge's window start.
        t_start: f64,
        /// The edge's window end.
        t_end: f64,
        /// The curve's declared span start.
        span_start: f64,
        /// The curve's declared span end.
        span_end: f64,
    },

    /// Edge curve deviates from adjacent surface.
    #[error("Edge {edge:?} deviates from surface by {max_deviation}")]
    CurveSurfaceDeviation {
        /// The deviating edge.
        edge: EdgeHandle,
        /// Maximum deviation distance.
        max_deviation: f64,
    },

    /// A solid region's shell faces the wrong way: its faces' folded
    /// normals must point out of the material, so the outer shell
    /// encloses positive volume and every inner (cavity) shell
    /// negative.
    #[error(
        "solid region {region:?}: shell {shell:?} ({}) encloses signed volume {volume}, \
         its normals point into the material",
        if *outer { "outer" } else { "cavity" }
    )]
    RegionSense {
        /// The solid region.
        region: RegionHandle,
        /// The shell facing the wrong way.
        shell: ShellHandle,
        /// Whether `shell` is the region's outer shell.
        outer: bool,
        /// The shell's folded-normal volume (approximate; its sign is
        /// the finding).
        volume: f64,
    },

    /// A fin of an edge is missing from the ring walked from the
    /// edge's first fin, or the first fin is not one of the edge's.
    #[error("fin {half_edge:?} of edge {edge:?} is not on the edge's ring")]
    RingMembership {
        /// The edge.
        edge: EdgeHandle,
        /// The fin off the ring (or the edge's stray first fin).
        half_edge: HalfEdgeHandle,
    },

    /// The two face sides bounding a wedge about an edge lie in
    /// different shells: the space between two radially adjacent faces
    /// would belong to two regions at once.
    #[error(
        "wedge at edge {edge:?} from fin {from:?} (shell {from_shell:?}) to fin {to:?} \
         (shell {to_shell:?}) spans two shells"
    )]
    WedgeShellMismatch {
        /// The edge.
        edge: EdgeHandle,
        /// The fin the wedge starts after.
        from: HalfEdgeHandle,
        /// The fin the wedge ends at.
        to: HalfEdgeHandle,
        /// The shell of `from`'s side facing the wedge.
        from_shell: ShellHandle,
        /// The shell of `to`'s side facing the wedge.
        to_shell: ShellHandle,
    },

    /// A shell's face sides, wire edges and acorn fall into more than
    /// one connected piece (through wedges and shared vertices).
    #[error("shell {shell:?} falls into {pieces} disconnected pieces")]
    ShellDisconnected {
        /// The shell.
        shell: ShellHandle,
        /// How many pieces.
        pieces: usize,
    },

    /// A shell's face sides do not form a closed surface: counting an
    /// edge once per wedge of the shell and a vertex once per fan of
    /// the shell, V - E + 2F - L must be 2 - 2g for a genus g >= 0.
    #[error("shell {shell:?}: V - E + 2F - L = {euler}, not 2 - 2g for any genus g >= 0")]
    ShellEuler {
        /// The shell.
        shell: ShellHandle,
        /// The computed V - E + 2F - L.
        euler: i64,
    },

    /// A region breaks the region structure.
    #[error("region {region:?}: {defect}")]
    RegionStructure {
        /// The region.
        region: RegionHandle,
        /// What is wrong with it.
        defect: RegionDefect,
    },

    /// Gap between adjacent faces at a shared edge.
    #[error("Shell gap at edge {edge:?} (gap: {gap})")]
    ShellGap {
        /// The edge with a gap.
        edge: EdgeHandle,
        /// The gap distance.
        gap: f64,
    },

    /// Two non-adjacent faces intersect each other.
    #[error("Self-intersection between face {face_a:?} and face {face_b:?}")]
    SelfIntersection {
        /// First face involved in intersection.
        face_a: FaceHandle,
        /// Second face involved in intersection.
        face_b: FaceHandle,
    },

    /// A face's stored orientation flag contradicts its outer loop's
    /// UV winding (interior-left invariant, the adjacent-face
    /// consistency half of check 7). Consumers folding the flag —
    /// edge convexity, mesher triangle winding, STEP `same_sense`,
    /// measured volume sign — mirror against neighbours on such a
    /// face.
    #[error(
        "Winding-dishonest face {face:?}: outer loop signed UV area {signed_uv_area:.3e} \
         contradicts orientation flag {orientation}"
    )]
    WindingDishonest {
        /// The face whose flag lies about its winding.
        face: FaceHandle,
        /// Signed UV area of the outer loop in topological order
        /// (positive = counter-clockwise).
        signed_uv_area: f64,
        /// The stored orientation flag.
        orientation: bool,
    },

    /// A face with surface geometry has loop half-edges without
    /// pcurves. Trim domains, loop-winding audits, and tessellation
    /// all read pcurves; a pcurve-less face silently degrades every
    /// one of them (measured: thicken served pcurve-less bodies for
    /// an unknown period — the winding audit deferred to its
    /// tangency witness and the mesher fell back to untrimmed
    /// full-grid faces; PK-331).
    #[error(
        "Face {face:?} carries surface geometry but {missing} of its {total} loop \
         half-edges have no pcurve"
    )]
    MissingPcurves {
        /// The face whose trim is unreadable.
        face: FaceHandle,
        /// Loop half-edges without a pcurve.
        missing: usize,
        /// Total loop half-edges on the face.
        total: usize,
    },

    /// A half-edge's pcurve, mapped through its face's surface, strays
    /// from its edge (ends, sense or trace): every trim consumer reads
    /// the pcurve (PK-819, PK-820).
    #[error(
        "Half-edge {half_edge:?} on face {face:?}: pcurve strays {max_deviation:.3e} from its edge"
    )]
    PcurveDeviation {
        /// The half-edge whose pcurve disagrees with its edge.
        half_edge: HalfEdgeHandle,
        /// The face owning the pcurve.
        face: FaceHandle,
        /// Worst accepted deviation (see `HalfEdgePcurveDeviation`).
        max_deviation: f64,
    },

    /// A face's loops leave a net winding around a periodic chart with
    /// no collapsed row, so they bound no region of it: some edge runs
    /// the long way round the period (PK-792).
    #[error(
        "Face {face:?}: its loops wind {winding_u} times round the u period and \
         {winding_v} times round the v period instead of cancelling"
    )]
    LoopPeriodWinding {
        /// The face whose trim bounds no region.
        face: FaceHandle,
        /// Net winding of all the face's loops round the u period.
        winding_u: i64,
        /// Net winding of all the face's loops round the v period.
        winding_v: i64,
    },

    /// A face's inner loop winds the same way as its outer loop: a hole
    /// must keep the face's material on its left like every loop, so
    /// it runs clockwise in UV when the orientation flag is `true`.
    /// Twins still pair across such a loop, so the shell closes with
    /// the material on opposite sides of it (measured: a block union
    /// served half inside out, signed volume -606.7 for +617.9, joined
    /// through a boss bottom's hole wound like its outer loop; PK-916,
    /// PK-923).
    #[error(
        "Face {face:?}: inner loop {loop_:?} signed UV area {signed_uv_area:.3e} \
         winds with the outer loop under orientation flag {orientation}"
    )]
    InnerLoopWinding {
        /// The face holding the loop.
        face: FaceHandle,
        /// The inner loop wound like the outer loop.
        loop_: LoopHandle,
        /// Signed UV area of the inner loop in topological order
        /// (positive = counter-clockwise).
        signed_uv_area: f64,
        /// The face's stored orientation flag.
        orientation: bool,
    },

    /// A face's surface normal reverses direction (surface folds over itself).
    #[error("Twisted surface on face {face:?}")]
    TwistedSurface {
        /// The face with a twisted surface.
        face: FaceHandle,
    },

    /// Precision constraint violated.
    #[error("precision: {}", .diagnostic.remediation)]
    Precision {
        /// The diagnostic describing the precision violation.
        diagnostic: crate::precision::PrecisionDiagnostic,
    },
}

/// Errors from healing analysis and repair operations.
#[derive(Debug, Error)]
pub enum HealingError {
    /// Gap repair failed on a specific edge.
    #[error("gap repair failed on edge {edge:?}: {reason}")]
    GapRepairFailed {
        /// The edge where repair failed.
        edge: EdgeHandle,
        /// Why the repair failed.
        reason: String,
    },

    /// Tolerance hierarchy enforcement failed.
    #[error("tolerance enforcement failed: {0}")]
    ToleranceEnforcementFailed(String),

    /// Surface simplification failed on a specific face.
    #[error("surface simplification failed on face {face:?}: {reason}")]
    SimplificationFailed {
        /// The face where simplification failed.
        face: FaceHandle,
        /// Why the simplification failed.
        reason: String,
    },

    /// Repair pipeline failed at a specific step.
    #[error("repair pipeline failed at step {step}: {reason}")]
    PipelineFailed {
        /// The pipeline step that failed.
        step: String,
        /// Why the step failed.
        reason: String,
    },

    /// A repair would displace geometry beyond the configured limit.
    #[error("max displacement {displacement:.6} exceeded limit {limit:.6} on edge {edge:?}")]
    DisplacementExceeded {
        /// The edge where displacement was too large.
        edge: EdgeHandle,
        /// The actual displacement.
        displacement: f64,
        /// The configured limit.
        limit: f64,
    },

    /// Repair iteration made the model worse (more defects after than before).
    #[error("repair made model worse: {before} defects before, {after} after")]
    RepairRegression {
        /// Defect count before the repair pass.
        before: usize,
        /// Defect count after the repair pass.
        after: usize,
    },
}

/// Top-level error type for popkern-core operations.
#[derive(Debug, Error)]
pub enum CoreError {
    /// A validation check failed.
    #[error(transparent)]
    Validation(#[from] ValidationError),

    /// A healing operation failed.
    #[error(transparent)]
    Healing(#[from] HealingError),

    /// A handle does not exist in the BRep.
    #[error("invalid handle: {0}")]
    InvalidHandle(String),

    /// A topological precondition was not met.
    #[error("precondition failed: {0}")]
    PreconditionFailed(String),

    /// A geometric operation violates tolerance.
    #[error("tolerance violation: {0}")]
    ToleranceViolation(String),

    /// Curve is degenerate (zero length or collapsed).
    #[error("degenerate curve: {0}")]
    DegenerateCurve(String),

    /// Surface is degenerate (zero area or collapsed).
    #[error("degenerate surface: {0}")]
    DegenerateSurface(String),

    /// Invalid NURBS data (bad knot vector, wrong control point count, etc).
    #[error("invalid NURBS: {0}")]
    InvalidNurbs(String),

    /// An operation that needs a manifold edge met a fin whose edge does
    /// not have exactly two fins (a free edge, a non-manifold edge, or a
    /// fin not linked into a ring).
    #[error("edge {edge:?} is not a manifold edge: its radial ring does not hold exactly two fins")]
    NoMate {
        /// The edge of the fin asked for its mate.
        edge: crate::topology::EdgeHandle,
    },
}

#[cfg(test)]
#[path = "error_tests.rs"]
mod tests;
