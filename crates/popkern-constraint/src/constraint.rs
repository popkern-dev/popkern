use crate::elements::{
    ArcHandle, ElementHandle, EllipseHandle, HyperbolaHandle, LineHandle, ParabolaHandle,
    PointHandle, SplineHandle,
};
use crate::groups::GroupHandle;
use serde::{Deserialize, Serialize};

slotmap::new_key_type! {
    /// Opaque handle to a [`ConstraintEntry`] inside a
    /// [`Sketch`](crate::Sketch).
    pub struct ConstraintHandle;
}

/// All supported geometric constraints.
///
/// Constraints are grouped into two tiers by implementation priority:
///
/// **Tier 1** — core constraints required for basic parametric modelling.
/// **Tier 2** — derived constraints that build on Tier 1.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Constraint {
    // ── Tier 1 ──────────────────────────────────────────────────────────────
    /// Force two points to share the same position.
    Coincident(PointHandle, PointHandle),

    /// Force a point to lie on an infinite line.
    PointOnLine(PointHandle, LineHandle),

    /// Force a point to lie on an arc (at the arc's radius from its centre).
    PointOnArc(PointHandle, ArcHandle),

    /// Fix the distance between two points to `value`.
    Distance(PointHandle, PointHandle, f64),

    /// Fix the perpendicular distance from a point to a line to `value`.
    PointLineDistance(PointHandle, LineHandle, f64),

    /// Fix the angle between two lines to `value` (radians).
    Angle(LineHandle, LineHandle, f64),

    /// Force a line to be horizontal (zero angle from the x-axis).
    Horizontal(LineHandle),

    /// Force a line to be vertical (90° from the x-axis).
    Vertical(LineHandle),

    /// Lock a point at its current position.
    ///
    /// The reference values are captured at constraint creation time.
    Fixed(PointHandle),

    /// Lock a line at its current position.
    ///
    /// The reference values are captured at constraint creation time.
    FixedLine(LineHandle),

    /// Force a line to be tangent to an arc at their intersection.
    Tangent(LineHandle, ArcHandle),

    /// Force two arcs to be tangent to each other.
    ArcTangent(ArcHandle, ArcHandle),

    /// Force two lines to have the same length.
    EqualLength(LineHandle, LineHandle),

    /// Force two arcs to have the same radius.
    EqualRadius(ArcHandle, ArcHandle),

    /// Force two lines to be parallel.
    Parallel(LineHandle, LineHandle),

    /// Force two lines to be perpendicular.
    Perpendicular(LineHandle, LineHandle),

    // ── Tier 2 ──────────────────────────────────────────────────────────────
    /// Force two points to be symmetric about a line.
    Symmetric(PointHandle, PointHandle, LineHandle),

    /// Force a point to lie at the midpoint of a line.
    Midpoint(PointHandle, LineHandle),

    /// Force two arcs to share the same centre.
    Concentric(ArcHandle, ArcHandle),

    /// Force two lines to be collinear (same infinite line).
    Collinear(LineHandle, LineHandle),

    // ── Conic constraints ──────────────────────────────────────────────────
    /// Force a point to lie on an ellipse.
    PointOnEllipse(PointHandle, EllipseHandle),

    /// Force a point to lie on a parabola.
    PointOnParabola(PointHandle, ParabolaHandle),

    /// Force a point to lie on a hyperbola.
    PointOnHyperbola(PointHandle, HyperbolaHandle),

    /// Force a line to be tangent to an ellipse.
    TangentLineEllipse(LineHandle, EllipseHandle),

    /// Force a line to be tangent to a parabola.
    TangentLineParabola(LineHandle, ParabolaHandle),

    /// Force a line to be tangent to a hyperbola.
    TangentLineHyperbola(LineHandle, HyperbolaHandle),

    /// Force an arc to be tangent to an ellipse (uses one auxiliary param).
    TangentArcEllipse(ArcHandle, EllipseHandle),

    /// Force an arc to be tangent to a parabola (uses one auxiliary param).
    TangentArcParabola(ArcHandle, ParabolaHandle),

    /// Force an arc to be tangent to a hyperbola (uses one auxiliary param).
    TangentArcHyperbola(ArcHandle, HyperbolaHandle),

    /// Force two ellipses to share the same centre.
    ConcentricEllipse(EllipseHandle, EllipseHandle),

    /// Force an arc and an ellipse to share the same centre.
    ConcentricArcEllipse(ArcHandle, EllipseHandle),

    /// Force two ellipses to have equal semi-major axis.
    EqualSemiAxis(EllipseHandle, EllipseHandle),

    /// Fix the focal distance of an ellipse to `value`.
    FocalDistance(EllipseHandle, f64),

    /// Fix the eccentricity of an ellipse to `value`.
    EllipseEccentricity(EllipseHandle, f64),

    /// Fix the eccentricity of a hyperbola to `value`.
    HyperbolaEccentricity(HyperbolaHandle, f64),

    /// Lock an ellipse at its current position and shape.
    FixedEllipse(EllipseHandle),

    /// Lock a parabola at its current position and shape.
    FixedParabola(ParabolaHandle),

    /// Lock a hyperbola at its current position and shape.
    FixedHyperbola(HyperbolaHandle),

    // ── Convenience arc constraints ────────────────────────────────────────
    /// Force two arcs to share the same centre AND radius (Concentric +
    /// EqualRadius combined as a single 3-equation constraint).
    ///
    /// Equations: cx1-cx2 = 0, cy1-cy2 = 0, r1-r2 = 0.
    Coradial(ArcHandle, ArcHandle),

    // ── Offset constraints ─────────────────────────────────────────────────
    /// Fix the perpendicular distance between two parallel lines to `d`.
    ///
    /// Uses the squared distance from the midpoint of L1 to the infinite
    /// line through L2. Assumes the lines are already parallel (enforce
    /// that separately with [`Constraint::Parallel`]).
    OffsetLine(LineHandle, LineHandle, f64),

    /// Fix the radial distance between two concentric arcs to `d`.
    ///
    /// Equation: (r1-r2)^2 - d^2 = 0 (squared form).  Concentricity is
    /// not enforced here — add a [`Constraint::Concentric`] separately.
    OffsetArc(ArcHandle, ArcHandle, f64),

    // ── Ordinate constraints ────────────────────────────────────────────────
    /// Fix the x coordinate of a point to `value`.
    ///
    /// Equation: px - value = 0.
    OrdinateX(PointHandle, f64),

    /// Fix the y coordinate of a point to `value`.
    ///
    /// Equation: py - value = 0.
    OrdinateY(PointHandle, f64),

    // ── Spline constraints ─────────────────────────────────────────────────
    /// Force a point to lie on a spline curve.
    ///
    /// Uses one auxiliary parameter for the spline parameter `t`.
    PointOnSpline(PointHandle, SplineHandle),

    /// Constrain the end tangent of a spline at endpoint `end_idx` (0 or 1)
    /// to be parallel to a line.
    SplineEndTangent(SplineHandle, usize, LineHandle),

    /// Constrain the end tangent of a spline at endpoint `end_idx` (0 or 1)
    /// to be tangent to an arc.
    SplineEndTangentToArc(SplineHandle, usize, ArcHandle),

    /// Fix the end tangent of a spline at endpoint `end_idx` (0 or 1) to the
    /// direction `(dx, dy)` (not required to be unit length).
    SplineFixedEndTangent(SplineHandle, usize, f64, f64),

    /// Fix the end curvature of a spline at endpoint `end_idx` (0 or 1) to
    /// a given signed value.
    SplineFixedEndCurvature(SplineHandle, usize, f64),

    /// Force a spline to be tangent to a line at their intersection.
    ///
    /// Uses one auxiliary parameter for the contact parameter `t`.
    SplineTangentToLine(SplineHandle, LineHandle),

    /// Force a spline to be tangent to an arc at their intersection.
    ///
    /// Uses one auxiliary parameter for the contact parameter `t`.
    SplineTangentToArc(SplineHandle, ArcHandle),

    /// Force a spline to be tangent to an ellipse at their intersection.
    ///
    /// Uses two auxiliary parameters: spline parameter `t` and ellipse
    /// contact angle.
    SplineTangentToEllipse(SplineHandle, EllipseHandle),

    /// Force a spline to be tangent to a parabola at their intersection.
    ///
    /// Uses two auxiliary parameters: spline parameter `t` and parabola
    /// contact parameter.
    SplineTangentToParabola(SplineHandle, ParabolaHandle),

    /// Force a spline to be tangent to a hyperbola at their intersection.
    ///
    /// Uses two auxiliary parameters: spline parameter `t` and hyperbola
    /// contact parameter.
    SplineTangentToHyperbola(SplineHandle, HyperbolaHandle),

    /// G1 (tangent) continuity between two splines at the specified endpoints
    /// (`end_idx` 0 or 1 for each).
    SplineG1(SplineHandle, usize, SplineHandle, usize),

    /// G2 (curvature) continuity between two splines at the specified
    /// endpoints (`end_idx` 0 or 1 for each).
    SplineG2(SplineHandle, usize, SplineHandle, usize),

    /// Fix the arc length of a spline to `value`.
    SplineLength(SplineHandle, f64),

    // ── Pattern constraints ────────────────────────────────────────────
    /// Constrain P2 to be offset from P1 by a fixed `(dx, dy)` vector.
    ///
    /// Equations: p2x - p1x - dx = 0, p2y - p1y - dy = 0.
    RelativePosition(PointHandle, PointHandle, f64, f64),

    /// Constrain P2 to be the rotation of P1 about Center by angle `theta`.
    ///
    /// Applies 2D rotation matrix centred at Center.
    RotatedPosition(PointHandle, PointHandle, PointHandle, f64),

    /// Constrain P2 to be the reflection of P1 across the infinite line
    /// through the endpoints of Line.
    ReflectedPosition(PointHandle, PointHandle, LineHandle),

    /// Pin a point to the centre (or vertex) of a curved element.
    ///
    /// Dispatches on element type: Arc/Ellipse centre, Parabola vertex,
    /// Hyperbola centre.
    PinToCenter(PointHandle, ElementHandle),

    /// Offset the start and end angles of Arc2 from Arc1 by `theta`.
    ///
    /// Equations: arc2.start - arc1.start - theta = 0,
    ///            arc2.end   - arc1.end   - theta = 0.
    AngleOffsetArc(ArcHandle, ArcHandle, f64),

    /// Offset the rotation parameter of a conic (Ellipse, Parabola, or
    /// Hyperbola) E2 from E1 by `theta`.
    ///
    /// Equation: e2.rotation - e1.rotation - theta = 0.
    AngleOffsetConic(ElementHandle, ElementHandle, f64),

    /// Mirror the arc angles of Arc2 to be the reflection of Arc1's angles
    /// about the angle of a mirror line. Note that start/end swap to
    /// preserve winding.
    MirrorAnglesArc(ArcHandle, ArcHandle, LineHandle),

    /// Mirror the rotation of a conic E2 to be the reflection of E1's
    /// rotation about the angle of a mirror line.
    MirrorAnglesConic(ElementHandle, ElementHandle, LineHandle),
}

// References and aux_param_count methods are in the constraint_refs submodule,
// split out to stay within the 500-line file limit.
#[path = "constraint_refs.rs"]
mod refs;

/// Whether a constraint drives the geometry or is merely displayed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConstraintKind {
    /// A driving constraint: influences the solver and changes the geometry.
    Driving,
    /// A driven constraint: a read-only annotation that does not move points.
    Driven,
}

/// Runtime satisfaction state of a single constraint.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConstraintStatus {
    /// Constraint is satisfied within tolerance.
    Satisfied,
    /// Constraint is not satisfied.
    Violated,
    /// Constraint duplicates the effect of another constraint.
    Redundant,
    /// Constraint conflicts with another constraint.
    Conflicting,
    /// Status has not been evaluated since the last geometry change.
    Unknown,
}

/// A constraint together with its metadata stored in the sketch.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConstraintEntry {
    /// The constraint itself.
    pub constraint: Constraint,
    /// Whether this entry drives the solver or is merely displayed.
    pub kind: ConstraintKind,
    /// Cached satisfaction state (updated after each solve).
    pub status: ConstraintStatus,
    /// Parameter values captured at creation time (used by [`Constraint::Fixed`]
    /// and [`Constraint::FixedLine`]).
    pub fixed_values: Vec<f64>,
    /// Optional group membership. A constraint belongs to at most one group.
    ///
    /// When the group is disabled the constraint is excluded from the solver.
    #[serde(default)]
    pub group: Option<GroupHandle>,
}
