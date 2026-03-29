//! Element-reference query methods for [`Constraint`].
//!
//! Split from `constraint.rs` to stay within the 500-line file limit.

use super::Constraint;
use crate::elements::{
    ArcHandle, EllipseHandle, HyperbolaHandle, LineHandle, ParabolaHandle, PointHandle,
    SplineHandle,
};

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
            Constraint::OrdinateX(p, _) | Constraint::OrdinateY(p, _) => *p == ph,
            Constraint::PointOnSpline(p, _) => *p == ph,
            // Line-only, arc-only, conic-only, and spline-only constraints do
            // not reference points directly.
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
            | Constraint::FixedHyperbola(..)
            | Constraint::Coradial(..)
            | Constraint::OffsetLine(..)
            | Constraint::OffsetArc(..)
            | Constraint::SplineEndTangent(..)
            | Constraint::SplineEndTangentToArc(..)
            | Constraint::SplineFixedEndTangent(..)
            | Constraint::SplineFixedEndCurvature(..)
            | Constraint::SplineTangentToLine(..)
            | Constraint::SplineTangentToArc(..)
            | Constraint::SplineTangentToEllipse(..)
            | Constraint::SplineTangentToParabola(..)
            | Constraint::SplineTangentToHyperbola(..)
            | Constraint::SplineG1(..)
            | Constraint::SplineG2(..)
            | Constraint::SplineLength(..) => false,
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
            | Constraint::Collinear(a, b)
            | Constraint::OffsetLine(a, b, _) => *a == lh || *b == lh,
            Constraint::Tangent(l, _)
            | Constraint::TangentLineEllipse(l, _)
            | Constraint::TangentLineParabola(l, _)
            | Constraint::TangentLineHyperbola(l, _) => *l == lh,
            Constraint::Symmetric(_, _, l) => *l == lh,
            Constraint::SplineEndTangent(_, _, l) | Constraint::SplineTangentToLine(_, l) => {
                *l == lh
            }
            // Point-only, arc-only, conic-only, and non-line spline
            // constraints do not reference lines.
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
            | Constraint::FixedHyperbola(..)
            | Constraint::Coradial(..)
            | Constraint::OffsetArc(..)
            | Constraint::OrdinateX(..)
            | Constraint::OrdinateY(..)
            | Constraint::PointOnSpline(..)
            | Constraint::SplineEndTangentToArc(..)
            | Constraint::SplineFixedEndTangent(..)
            | Constraint::SplineFixedEndCurvature(..)
            | Constraint::SplineTangentToArc(..)
            | Constraint::SplineTangentToEllipse(..)
            | Constraint::SplineTangentToParabola(..)
            | Constraint::SplineTangentToHyperbola(..)
            | Constraint::SplineG1(..)
            | Constraint::SplineG2(..)
            | Constraint::SplineLength(..) => false,
        }
    }

    /// Return `true` if this constraint references the given arc.
    pub fn references_arc(&self, ah: ArcHandle) -> bool {
        match self {
            Constraint::PointOnArc(_, a) | Constraint::Tangent(_, a) => *a == ah,
            Constraint::ArcTangent(a, b)
            | Constraint::EqualRadius(a, b)
            | Constraint::Concentric(a, b)
            | Constraint::Coradial(a, b)
            | Constraint::OffsetArc(a, b, _) => *a == ah || *b == ah,
            Constraint::TangentArcEllipse(a, _)
            | Constraint::TangentArcParabola(a, _)
            | Constraint::TangentArcHyperbola(a, _) => *a == ah,
            Constraint::ConcentricArcEllipse(a, _) => *a == ah,
            Constraint::SplineEndTangentToArc(_, _, a) | Constraint::SplineTangentToArc(_, a) => {
                *a == ah
            }
            // Everything else does not reference arcs.
            _ => false,
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
            | Constraint::FixedEllipse(e)
            | Constraint::SplineTangentToEllipse(_, e) => *e == eh,
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
            | Constraint::FixedParabola(p)
            | Constraint::SplineTangentToParabola(_, p) => *p == ph,
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
            | Constraint::FixedHyperbola(h)
            | Constraint::SplineTangentToHyperbola(_, h) => *h == hh,
            _ => false,
        }
    }

    /// Return `true` if this constraint references the given spline.
    pub fn references_spline(&self, sh: SplineHandle) -> bool {
        match self {
            Constraint::PointOnSpline(_, s)
            | Constraint::SplineEndTangent(s, _, _)
            | Constraint::SplineEndTangentToArc(s, _, _)
            | Constraint::SplineFixedEndTangent(s, _, _, _)
            | Constraint::SplineFixedEndCurvature(s, _, _)
            | Constraint::SplineTangentToLine(s, _)
            | Constraint::SplineTangentToArc(s, _)
            | Constraint::SplineTangentToEllipse(s, _)
            | Constraint::SplineTangentToParabola(s, _)
            | Constraint::SplineTangentToHyperbola(s, _)
            | Constraint::SplineLength(s, _) => *s == sh,
            Constraint::SplineG1(s1, _, s2, _) | Constraint::SplineG2(s1, _, s2, _) => {
                *s1 == sh || *s2 == sh
            }
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
            Constraint::PointOnSpline(..)
            | Constraint::SplineTangentToLine(..)
            | Constraint::SplineTangentToArc(..) => 1,
            Constraint::SplineTangentToEllipse(..)
            | Constraint::SplineTangentToParabola(..)
            | Constraint::SplineTangentToHyperbola(..) => 2,
            _ => 0,
        }
    }
}
