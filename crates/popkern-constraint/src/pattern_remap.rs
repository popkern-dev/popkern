//! Constraint remapping and element-reference analysis for pattern cloning.
//!
//! Split from `pattern_clone.rs` to stay within the 500-line file limit.

use crate::Sketch;
use crate::constraint::Constraint;
use crate::elements::{ElementHandle, LineHandle, PointHandle};
use crate::pattern::HandleRemap;
use std::collections::HashSet;

/// Check whether all elements referenced by a constraint are in the
/// given set.
pub(crate) fn constraint_elements_in_set(
    constraint: &Constraint,
    sketch: &Sketch,
    set: &HashSet<ElementHandle>,
) -> bool {
    let refs = referenced_elements(constraint, sketch);
    refs.iter().all(|e| set.contains(e))
}

/// Collect all element handles referenced by a constraint.
pub(crate) fn referenced_elements(constraint: &Constraint, _sketch: &Sketch) -> Vec<ElementHandle> {
    match constraint {
        Constraint::Coincident(a, b) | Constraint::Distance(a, b, _) => {
            vec![ElementHandle::Point(*a), ElementHandle::Point(*b)]
        }
        Constraint::Fixed(p) | Constraint::OrdinateX(p, _) | Constraint::OrdinateY(p, _) => {
            vec![ElementHandle::Point(*p)]
        }
        Constraint::PointOnLine(p, l) | Constraint::PointLineDistance(p, l, _) => {
            vec![ElementHandle::Point(*p), ElementHandle::Line(*l)]
        }
        Constraint::PointOnArc(p, a) => {
            vec![ElementHandle::Point(*p), ElementHandle::Arc(*a)]
        }
        Constraint::Horizontal(l) | Constraint::Vertical(l) | Constraint::FixedLine(l) => {
            vec![ElementHandle::Line(*l)]
        }
        Constraint::Angle(a, b, _)
        | Constraint::Parallel(a, b)
        | Constraint::Perpendicular(a, b)
        | Constraint::EqualLength(a, b)
        | Constraint::Collinear(a, b)
        | Constraint::OffsetLine(a, b, _) => {
            vec![ElementHandle::Line(*a), ElementHandle::Line(*b)]
        }
        Constraint::Tangent(l, a) => {
            vec![ElementHandle::Line(*l), ElementHandle::Arc(*a)]
        }
        Constraint::ArcTangent(a1, a2)
        | Constraint::EqualRadius(a1, a2)
        | Constraint::Concentric(a1, a2)
        | Constraint::Coradial(a1, a2)
        | Constraint::OffsetArc(a1, a2, _) => {
            vec![ElementHandle::Arc(*a1), ElementHandle::Arc(*a2)]
        }
        Constraint::Symmetric(p1, p2, l) => {
            vec![
                ElementHandle::Point(*p1),
                ElementHandle::Point(*p2),
                ElementHandle::Line(*l),
            ]
        }
        Constraint::Midpoint(p, l) => {
            vec![ElementHandle::Point(*p), ElementHandle::Line(*l)]
        }
        Constraint::PointOnEllipse(p, e) => {
            vec![ElementHandle::Point(*p), ElementHandle::Ellipse(*e)]
        }
        Constraint::PointOnParabola(p, pa) => {
            vec![ElementHandle::Point(*p), ElementHandle::Parabola(*pa)]
        }
        Constraint::PointOnHyperbola(p, h) => {
            vec![ElementHandle::Point(*p), ElementHandle::Hyperbola(*h)]
        }
        Constraint::TangentLineEllipse(l, e) => {
            vec![ElementHandle::Line(*l), ElementHandle::Ellipse(*e)]
        }
        Constraint::TangentLineParabola(l, p) => {
            vec![ElementHandle::Line(*l), ElementHandle::Parabola(*p)]
        }
        Constraint::TangentLineHyperbola(l, h) => {
            vec![ElementHandle::Line(*l), ElementHandle::Hyperbola(*h)]
        }
        Constraint::TangentArcEllipse(a, e) => {
            vec![ElementHandle::Arc(*a), ElementHandle::Ellipse(*e)]
        }
        Constraint::TangentArcParabola(a, p) => {
            vec![ElementHandle::Arc(*a), ElementHandle::Parabola(*p)]
        }
        Constraint::TangentArcHyperbola(a, h) => {
            vec![ElementHandle::Arc(*a), ElementHandle::Hyperbola(*h)]
        }
        Constraint::ConcentricEllipse(e1, e2) | Constraint::EqualSemiAxis(e1, e2) => {
            vec![ElementHandle::Ellipse(*e1), ElementHandle::Ellipse(*e2)]
        }
        Constraint::ConcentricArcEllipse(a, e) => {
            vec![ElementHandle::Arc(*a), ElementHandle::Ellipse(*e)]
        }
        Constraint::FocalDistance(e, _) | Constraint::EllipseEccentricity(e, _) => {
            vec![ElementHandle::Ellipse(*e)]
        }
        Constraint::HyperbolaEccentricity(h, _) => {
            vec![ElementHandle::Hyperbola(*h)]
        }
        Constraint::FixedEllipse(e) => vec![ElementHandle::Ellipse(*e)],
        Constraint::FixedParabola(p) => vec![ElementHandle::Parabola(*p)],
        Constraint::FixedHyperbola(h) => vec![ElementHandle::Hyperbola(*h)],
        Constraint::PointOnSpline(p, s) => {
            vec![ElementHandle::Point(*p), ElementHandle::Spline(*s)]
        }
        Constraint::SplineEndTangent(s, _, l) | Constraint::SplineTangentToLine(s, l) => {
            vec![ElementHandle::Spline(*s), ElementHandle::Line(*l)]
        }
        Constraint::SplineEndTangentToArc(s, _, a) | Constraint::SplineTangentToArc(s, a) => {
            vec![ElementHandle::Spline(*s), ElementHandle::Arc(*a)]
        }
        Constraint::SplineFixedEndTangent(s, _, _, _)
        | Constraint::SplineFixedEndCurvature(s, _, _)
        | Constraint::SplineLength(s, _) => {
            vec![ElementHandle::Spline(*s)]
        }
        Constraint::SplineTangentToEllipse(s, e) => {
            vec![ElementHandle::Spline(*s), ElementHandle::Ellipse(*e)]
        }
        Constraint::SplineTangentToParabola(s, p) => {
            vec![ElementHandle::Spline(*s), ElementHandle::Parabola(*p)]
        }
        Constraint::SplineTangentToHyperbola(s, h) => {
            vec![ElementHandle::Spline(*s), ElementHandle::Hyperbola(*h)]
        }
        Constraint::SplineG1(s1, _, s2, _) | Constraint::SplineG2(s1, _, s2, _) => {
            vec![ElementHandle::Spline(*s1), ElementHandle::Spline(*s2)]
        }
        Constraint::RelativePosition(a, b, _, _) => {
            vec![ElementHandle::Point(*a), ElementHandle::Point(*b)]
        }
        Constraint::RotatedPosition(a, b, c, _) => {
            vec![
                ElementHandle::Point(*a),
                ElementHandle::Point(*b),
                ElementHandle::Point(*c),
            ]
        }
        Constraint::ReflectedPosition(a, b, l) => {
            vec![
                ElementHandle::Point(*a),
                ElementHandle::Point(*b),
                ElementHandle::Line(*l),
            ]
        }
        Constraint::PinToCenter(p, elem) => {
            vec![ElementHandle::Point(*p), *elem]
        }
        Constraint::AngleOffsetArc(a1, a2, _) | Constraint::MirrorAnglesArc(a1, a2, _) => {
            vec![ElementHandle::Arc(*a1), ElementHandle::Arc(*a2)]
        }
        Constraint::AngleOffsetConic(e1, e2, _) | Constraint::MirrorAnglesConic(e1, e2, _) => {
            vec![*e1, *e2]
        }
    }
}

/// Remap a single constraint's handles. Returns `None` if any
/// required handle is missing from the remap.
pub(crate) fn remap_constraint(constraint: &Constraint, remap: &HandleRemap) -> Option<Constraint> {
    let rp = |p: &PointHandle| remap.points.get(p).copied();
    let rl = |l: &LineHandle| remap.lines.get(l).copied();
    let ra = |a: &crate::ArcHandle| remap.arcs.get(a).copied();
    let re = |e: &crate::EllipseHandle| remap.ellipses.get(e).copied();
    let rpa = |p: &crate::ParabolaHandle| remap.parabolas.get(p).copied();
    let rh = |h: &crate::HyperbolaHandle| remap.hyperbolas.get(h).copied();
    let rs = |s: &crate::SplineHandle| remap.splines.get(s).copied();

    let relem = |e: &ElementHandle| -> Option<ElementHandle> {
        match e {
            ElementHandle::Point(p) => rp(p).map(ElementHandle::Point),
            ElementHandle::Line(l) => rl(l).map(ElementHandle::Line),
            ElementHandle::Arc(a) => ra(a).map(ElementHandle::Arc),
            ElementHandle::Ellipse(el) => re(el).map(ElementHandle::Ellipse),
            ElementHandle::Parabola(p) => rpa(p).map(ElementHandle::Parabola),
            ElementHandle::Hyperbola(h) => rh(h).map(ElementHandle::Hyperbola),
            ElementHandle::Spline(s) => rs(s).map(ElementHandle::Spline),
        }
    };

    Some(match constraint {
        Constraint::Coincident(a, b) => Constraint::Coincident(rp(a)?, rp(b)?),
        Constraint::Distance(a, b, d) => Constraint::Distance(rp(a)?, rp(b)?, *d),
        Constraint::Fixed(p) => Constraint::Fixed(rp(p)?),
        Constraint::FixedLine(l) => Constraint::FixedLine(rl(l)?),
        Constraint::PointOnLine(p, l) => Constraint::PointOnLine(rp(p)?, rl(l)?),
        Constraint::PointLineDistance(p, l, d) => Constraint::PointLineDistance(rp(p)?, rl(l)?, *d),
        Constraint::PointOnArc(p, a) => Constraint::PointOnArc(rp(p)?, ra(a)?),
        Constraint::Horizontal(l) => Constraint::Horizontal(rl(l)?),
        Constraint::Vertical(l) => Constraint::Vertical(rl(l)?),
        Constraint::Angle(a, b, v) => Constraint::Angle(rl(a)?, rl(b)?, *v),
        Constraint::Tangent(l, a) => Constraint::Tangent(rl(l)?, ra(a)?),
        Constraint::ArcTangent(a1, a2) => Constraint::ArcTangent(ra(a1)?, ra(a2)?),
        Constraint::EqualLength(a, b) => Constraint::EqualLength(rl(a)?, rl(b)?),
        Constraint::EqualRadius(a, b) => Constraint::EqualRadius(ra(a)?, ra(b)?),
        Constraint::Parallel(a, b) => Constraint::Parallel(rl(a)?, rl(b)?),
        Constraint::Perpendicular(a, b) => Constraint::Perpendicular(rl(a)?, rl(b)?),
        Constraint::Symmetric(p1, p2, l) => Constraint::Symmetric(rp(p1)?, rp(p2)?, rl(l)?),
        Constraint::Midpoint(p, l) => Constraint::Midpoint(rp(p)?, rl(l)?),
        Constraint::Concentric(a1, a2) => Constraint::Concentric(ra(a1)?, ra(a2)?),
        Constraint::Collinear(a, b) => Constraint::Collinear(rl(a)?, rl(b)?),
        Constraint::Coradial(a1, a2) => Constraint::Coradial(ra(a1)?, ra(a2)?),
        Constraint::OffsetLine(a, b, d) => Constraint::OffsetLine(rl(a)?, rl(b)?, *d),
        Constraint::OffsetArc(a1, a2, d) => Constraint::OffsetArc(ra(a1)?, ra(a2)?, *d),
        Constraint::OrdinateX(p, v) => Constraint::OrdinateX(rp(p)?, *v),
        Constraint::OrdinateY(p, v) => Constraint::OrdinateY(rp(p)?, *v),
        Constraint::PointOnEllipse(p, e) => Constraint::PointOnEllipse(rp(p)?, re(e)?),
        Constraint::PointOnParabola(p, pa) => Constraint::PointOnParabola(rp(p)?, rpa(pa)?),
        Constraint::PointOnHyperbola(p, h) => Constraint::PointOnHyperbola(rp(p)?, rh(h)?),
        Constraint::TangentLineEllipse(l, e) => Constraint::TangentLineEllipse(rl(l)?, re(e)?),
        Constraint::TangentLineParabola(l, p) => Constraint::TangentLineParabola(rl(l)?, rpa(p)?),
        Constraint::TangentLineHyperbola(l, h) => Constraint::TangentLineHyperbola(rl(l)?, rh(h)?),
        Constraint::TangentArcEllipse(a, e) => Constraint::TangentArcEllipse(ra(a)?, re(e)?),
        Constraint::TangentArcParabola(a, p) => Constraint::TangentArcParabola(ra(a)?, rpa(p)?),
        Constraint::TangentArcHyperbola(a, h) => Constraint::TangentArcHyperbola(ra(a)?, rh(h)?),
        Constraint::ConcentricEllipse(e1, e2) => Constraint::ConcentricEllipse(re(e1)?, re(e2)?),
        Constraint::ConcentricArcEllipse(a, e) => Constraint::ConcentricArcEllipse(ra(a)?, re(e)?),
        Constraint::EqualSemiAxis(e1, e2) => Constraint::EqualSemiAxis(re(e1)?, re(e2)?),
        Constraint::FocalDistance(e, d) => Constraint::FocalDistance(re(e)?, *d),
        Constraint::EllipseEccentricity(e, v) => Constraint::EllipseEccentricity(re(e)?, *v),
        Constraint::HyperbolaEccentricity(h, v) => Constraint::HyperbolaEccentricity(rh(h)?, *v),
        Constraint::FixedEllipse(e) => Constraint::FixedEllipse(re(e)?),
        Constraint::FixedParabola(p) => Constraint::FixedParabola(rpa(p)?),
        Constraint::FixedHyperbola(h) => Constraint::FixedHyperbola(rh(h)?),
        Constraint::PointOnSpline(p, s) => Constraint::PointOnSpline(rp(p)?, rs(s)?),
        Constraint::SplineEndTangent(s, i, l) => Constraint::SplineEndTangent(rs(s)?, *i, rl(l)?),
        Constraint::SplineEndTangentToArc(s, i, a) => {
            Constraint::SplineEndTangentToArc(rs(s)?, *i, ra(a)?)
        }
        Constraint::SplineFixedEndTangent(s, i, tx, ty) => {
            Constraint::SplineFixedEndTangent(rs(s)?, *i, *tx, *ty)
        }
        Constraint::SplineFixedEndCurvature(s, i, v) => {
            Constraint::SplineFixedEndCurvature(rs(s)?, *i, *v)
        }
        Constraint::SplineTangentToLine(s, l) => Constraint::SplineTangentToLine(rs(s)?, rl(l)?),
        Constraint::SplineTangentToArc(s, a) => Constraint::SplineTangentToArc(rs(s)?, ra(a)?),
        Constraint::SplineTangentToEllipse(s, e) => {
            Constraint::SplineTangentToEllipse(rs(s)?, re(e)?)
        }
        Constraint::SplineTangentToParabola(s, p) => {
            Constraint::SplineTangentToParabola(rs(s)?, rpa(p)?)
        }
        Constraint::SplineTangentToHyperbola(s, h) => {
            Constraint::SplineTangentToHyperbola(rs(s)?, rh(h)?)
        }
        Constraint::SplineG1(s1, e1, s2, e2) => Constraint::SplineG1(rs(s1)?, *e1, rs(s2)?, *e2),
        Constraint::SplineG2(s1, e1, s2, e2) => Constraint::SplineG2(rs(s1)?, *e1, rs(s2)?, *e2),
        Constraint::SplineLength(s, v) => Constraint::SplineLength(rs(s)?, *v),
        Constraint::RelativePosition(a, b, dx, dy) => {
            Constraint::RelativePosition(rp(a)?, rp(b)?, *dx, *dy)
        }
        Constraint::RotatedPosition(a, b, c, t) => {
            Constraint::RotatedPosition(rp(a)?, rp(b)?, rp(c)?, *t)
        }
        Constraint::ReflectedPosition(a, b, l) => {
            Constraint::ReflectedPosition(rp(a)?, rp(b)?, rl(l)?)
        }
        Constraint::PinToCenter(p, e) => Constraint::PinToCenter(rp(p)?, relem(e)?),
        Constraint::AngleOffsetArc(a1, a2, t) => Constraint::AngleOffsetArc(ra(a1)?, ra(a2)?, *t),
        Constraint::AngleOffsetConic(e1, e2, t) => {
            Constraint::AngleOffsetConic(relem(e1)?, relem(e2)?, *t)
        }
        Constraint::MirrorAnglesArc(a1, a2, l) => {
            Constraint::MirrorAnglesArc(ra(a1)?, ra(a2)?, rl(l)?)
        }
        Constraint::MirrorAnglesConic(e1, e2, l) => {
            Constraint::MirrorAnglesConic(relem(e1)?, relem(e2)?, rl(l)?)
        }
    })
}
