use crate::elements::{
    ArcHandle, EllipseHandle, HyperbolaHandle, LineHandle, ParabolaHandle, PointHandle,
};
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
}

impl Constraint {
    /// Return `true` if this constraint references the given point.
    pub fn references_point(&self, ph: PointHandle) -> bool {
        match self {
            Constraint::Coincident(a, b) | Constraint::Distance(a, b, _) => *a == ph || *b == ph,
            Constraint::Fixed(p) => *p == ph,
            Constraint::PointOnLine(p, _) | Constraint::PointLineDistance(p, _, _) => *p == ph,
            Constraint::PointOnArc(p, _)
            | Constraint::PointOnEllipse(p, _)
            | Constraint::PointOnParabola(p, _)
            | Constraint::PointOnHyperbola(p, _) => *p == ph,
            Constraint::Symmetric(a, b, _) => *a == ph || *b == ph,
            Constraint::Midpoint(p, _) => *p == ph,
            // Line-only, arc-only, and conic-only constraints do not reference
            // points directly.
            Constraint::Horizontal(_)
            | Constraint::Vertical(_)
            | Constraint::FixedLine(_)
            | Constraint::Angle(..)
            | Constraint::Parallel(..)
            | Constraint::Perpendicular(..)
            | Constraint::EqualLength(..)
            | Constraint::Tangent(..)
            | Constraint::ArcTangent(..)
            | Constraint::EqualRadius(..)
            | Constraint::Concentric(..)
            | Constraint::Collinear(..)
            | Constraint::TangentLineEllipse(..)
            | Constraint::TangentLineParabola(..)
            | Constraint::TangentLineHyperbola(..)
            | Constraint::TangentArcEllipse(..)
            | Constraint::TangentArcParabola(..)
            | Constraint::TangentArcHyperbola(..)
            | Constraint::ConcentricEllipse(..)
            | Constraint::ConcentricArcEllipse(..)
            | Constraint::EqualSemiAxis(..)
            | Constraint::FocalDistance(..)
            | Constraint::EllipseEccentricity(..)
            | Constraint::HyperbolaEccentricity(..)
            | Constraint::FixedEllipse(..)
            | Constraint::FixedParabola(..)
            | Constraint::FixedHyperbola(..) => false,
        }
    }

    /// Return `true` if this constraint references the given line.
    pub fn references_line(&self, lh: LineHandle) -> bool {
        match self {
            Constraint::PointOnLine(_, l)
            | Constraint::PointLineDistance(_, l, _)
            | Constraint::Horizontal(l)
            | Constraint::Vertical(l)
            | Constraint::FixedLine(l)
            | Constraint::Midpoint(_, l) => *l == lh,
            Constraint::Angle(a, b, _)
            | Constraint::Parallel(a, b)
            | Constraint::Perpendicular(a, b)
            | Constraint::EqualLength(a, b)
            | Constraint::Collinear(a, b) => *a == lh || *b == lh,
            Constraint::Tangent(l, _)
            | Constraint::TangentLineEllipse(l, _)
            | Constraint::TangentLineParabola(l, _)
            | Constraint::TangentLineHyperbola(l, _) => *l == lh,
            Constraint::Symmetric(_, _, l) => *l == lh,
            // Point-only, arc-only, and conic-only constraints do not
            // reference lines.
            Constraint::Coincident(..)
            | Constraint::Distance(..)
            | Constraint::Fixed(_)
            | Constraint::PointOnArc(..)
            | Constraint::PointOnEllipse(..)
            | Constraint::PointOnParabola(..)
            | Constraint::PointOnHyperbola(..)
            | Constraint::ArcTangent(..)
            | Constraint::EqualRadius(..)
            | Constraint::Concentric(..)
            | Constraint::TangentArcEllipse(..)
            | Constraint::TangentArcParabola(..)
            | Constraint::TangentArcHyperbola(..)
            | Constraint::ConcentricEllipse(..)
            | Constraint::ConcentricArcEllipse(..)
            | Constraint::EqualSemiAxis(..)
            | Constraint::FocalDistance(..)
            | Constraint::EllipseEccentricity(..)
            | Constraint::HyperbolaEccentricity(..)
            | Constraint::FixedEllipse(..)
            | Constraint::FixedParabola(..)
            | Constraint::FixedHyperbola(..) => false,
        }
    }

    /// Return `true` if this constraint references the given ellipse.
    pub fn references_ellipse(&self, eh: EllipseHandle) -> bool {
        match self {
            Constraint::PointOnEllipse(_, e)
            | Constraint::TangentLineEllipse(_, e)
            | Constraint::TangentArcEllipse(_, e)
            | Constraint::ConcentricArcEllipse(_, e)
            | Constraint::FocalDistance(e, _)
            | Constraint::EllipseEccentricity(e, _)
            | Constraint::FixedEllipse(e) => *e == eh,
            Constraint::ConcentricEllipse(a, b) | Constraint::EqualSemiAxis(a, b) => {
                *a == eh || *b == eh
            }
            _ => false,
        }
    }

    /// Return `true` if this constraint references the given parabola.
    pub fn references_parabola(&self, ph: ParabolaHandle) -> bool {
        match self {
            Constraint::PointOnParabola(_, p)
            | Constraint::TangentLineParabola(_, p)
            | Constraint::TangentArcParabola(_, p)
            | Constraint::FixedParabola(p) => *p == ph,
            _ => false,
        }
    }

    /// Return `true` if this constraint references the given hyperbola.
    pub fn references_hyperbola(&self, hh: HyperbolaHandle) -> bool {
        match self {
            Constraint::PointOnHyperbola(_, h)
            | Constraint::TangentLineHyperbola(_, h)
            | Constraint::TangentArcHyperbola(_, h)
            | Constraint::HyperbolaEccentricity(h, _)
            | Constraint::FixedHyperbola(h) => *h == hh,
            _ => false,
        }
    }

    /// Return the number of auxiliary parameters this constraint requires.
    ///
    /// Auxiliary parameters are solver-internal unknowns (e.g. a contact
    /// angle for arc-to-conic tangency) that do not correspond to any
    /// geometric element.
    pub fn aux_param_count(&self) -> usize {
        match self {
            Constraint::TangentArcEllipse(..)
            | Constraint::TangentArcParabola(..)
            | Constraint::TangentArcHyperbola(..) => 1,
            _ => 0,
        }
    }

    /// Return `true` if this constraint references the given arc.
    pub fn references_arc(&self, ah: ArcHandle) -> bool {
        match self {
            Constraint::PointOnArc(_, a) | Constraint::Tangent(_, a) => *a == ah,
            Constraint::ArcTangent(a, b)
            | Constraint::EqualRadius(a, b)
            | Constraint::Concentric(a, b) => *a == ah || *b == ah,
            Constraint::TangentArcEllipse(a, _)
            | Constraint::TangentArcParabola(a, _)
            | Constraint::TangentArcHyperbola(a, _) => *a == ah,
            Constraint::ConcentricArcEllipse(a, _) => *a == ah,
            // Everything else does not reference arcs.
            _ => false,
        }
    }
}

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
}
