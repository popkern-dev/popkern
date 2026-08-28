//! Dispatch logic for constraint evaluation and equation counting.
//!
//! Split from `equations.rs` to stay within the 500-line file limit.

use super::{
    EquationOutput, conics, conics_arc, conics_tangent, patterns, splines, splines_tangent,
    splines_tangent_conic, tier2,
};
use crate::Sketch;
use crate::constraint::{Constraint, ConstraintEntry, ConstraintHandle};

/// Return the number of scalar equations produced by a constraint.
pub(crate) fn equation_count(constraint: &Constraint) -> usize {
    match constraint {
        Constraint::Coincident(..) => 2,
        Constraint::Fixed(..) => 2,
        Constraint::FixedLine(..) => 4,
        Constraint::Symmetric(..) => 2,
        Constraint::Midpoint(..) => 2,
        Constraint::Concentric(..) => 2,
        Constraint::ConcentricEllipse(..) | Constraint::ConcentricArcEllipse(..) => 2,
        Constraint::FixedEllipse(..) | Constraint::FixedHyperbola(..) => 5,
        Constraint::FixedParabola(..) => 4,
        Constraint::TangentArcEllipse(..)
        | Constraint::TangentArcParabola(..)
        | Constraint::TangentArcHyperbola(..) => 2,
        Constraint::Coradial(..) => 3,
        Constraint::Horizontal(..)
        | Constraint::Vertical(..)
        | Constraint::Distance(..)
        | Constraint::PointOnLine(..)
        | Constraint::PointOnArc(..)
        | Constraint::PointLineDistance(..)
        | Constraint::Angle(..)
        | Constraint::Tangent(..)
        | Constraint::ArcTangent(..)
        | Constraint::EqualLength(..)
        | Constraint::EqualRadius(..)
        | Constraint::Parallel(..)
        | Constraint::Perpendicular(..)
        | Constraint::Collinear(..)
        | Constraint::PointOnEllipse(..)
        | Constraint::PointOnParabola(..)
        | Constraint::PointOnHyperbola(..)
        | Constraint::TangentLineEllipse(..)
        | Constraint::TangentLineParabola(..)
        | Constraint::TangentLineHyperbola(..)
        | Constraint::EqualSemiAxis(..)
        | Constraint::FocalDistance(..)
        | Constraint::EllipseEccentricity(..)
        | Constraint::HyperbolaEccentricity(..)
        | Constraint::OffsetLine(..)
        | Constraint::OffsetArc(..)
        | Constraint::OrdinateX(..)
        | Constraint::OrdinateY(..)
        | Constraint::SplineEndTangent(..)
        | Constraint::SplineEndTangentToArc(..)
        | Constraint::SplineFixedEndCurvature(..)
        | Constraint::SplineLength(..)
        | Constraint::AngleOffsetConic(..)
        | Constraint::MirrorAnglesConic(..) => 1,
        Constraint::PointOnSpline(..)
        | Constraint::SplineFixedEndTangent(..)
        | Constraint::SplineTangentToLine(..)
        | Constraint::SplineTangentToArc(..)
        | Constraint::RelativePosition(..)
        | Constraint::RotatedPosition(..)
        | Constraint::ReflectedPosition(..)
        | Constraint::PinToCenter(..)
        | Constraint::AngleOffsetArc(..)
        | Constraint::MirrorAnglesArc(..) => 2,
        Constraint::SplineG1(..)
        | Constraint::SplineTangentToEllipse(..)
        | Constraint::SplineTangentToParabola(..)
        | Constraint::SplineTangentToHyperbola(..) => 3,
        Constraint::SplineG2(..) => 4,
    }
}

/// Evaluate the residuals and Jacobian for a single constraint.
///
/// `handle` is needed so that constraints with auxiliary parameters
/// (e.g. arc-to-conic tangency) can look up their auxiliary parameter
/// indices via [`Sketch::aux_param_map`].
///
/// `row_offset` is added to every row index so that constraints can be
/// stacked into a single system.
pub(crate) fn evaluate_constraint(
    sketch: &Sketch,
    handle: Option<ConstraintHandle>,
    entry: &ConstraintEntry,
    row_offset: usize,
) -> EquationOutput {
    match &entry.constraint {
        Constraint::Coincident(p1, p2) => super::eval_coincident(sketch, *p1, *p2, row_offset),
        Constraint::Fixed(p) => super::eval_fixed(sketch, *p, &entry.fixed_values, row_offset),
        Constraint::FixedLine(l) => {
            super::eval_fixed_line(sketch, *l, &entry.fixed_values, row_offset)
        }
        Constraint::Horizontal(l) => super::eval_horizontal(sketch, *l, row_offset),
        Constraint::Vertical(l) => super::eval_vertical(sketch, *l, row_offset),
        Constraint::Distance(p1, p2, d) => super::eval_distance(sketch, *p1, *p2, *d, row_offset),
        Constraint::PointOnLine(p, l) => super::eval_point_on_line(sketch, *p, *l, row_offset),
        Constraint::Parallel(l1, l2) => super::eval_parallel(sketch, *l1, *l2, row_offset),
        Constraint::Perpendicular(l1, l2) => {
            super::eval_perpendicular(sketch, *l1, *l2, row_offset)
        }
        Constraint::EqualLength(l1, l2) => super::eval_equal_length(sketch, *l1, *l2, row_offset),
        // Arc-related and tier-2 constraints are in the tier2 submodule.
        Constraint::PointOnArc(p, a) => tier2::eval_point_on_arc(sketch, *p, *a, row_offset),
        Constraint::PointLineDistance(p, l, d) => {
            tier2::eval_point_line_distance(sketch, *p, *l, *d, row_offset)
        }
        Constraint::Angle(l1, l2, theta) => tier2::eval_angle(sketch, *l1, *l2, *theta, row_offset),
        Constraint::Tangent(l, a) => tier2::eval_tangent(sketch, *l, *a, row_offset),
        Constraint::ArcTangent(a1, a2) => tier2::eval_arc_tangent(sketch, *a1, *a2, row_offset),
        Constraint::EqualRadius(a1, a2) => tier2::eval_equal_radius(sketch, *a1, *a2, row_offset),
        Constraint::Symmetric(p1, p2, l) => tier2::eval_symmetric(sketch, *p1, *p2, *l, row_offset),
        Constraint::Midpoint(p, l) => tier2::eval_midpoint(sketch, *p, *l, row_offset),
        Constraint::Concentric(a1, a2) => tier2::eval_concentric(sketch, *a1, *a2, row_offset),
        Constraint::Collinear(l1, l2) => tier2::eval_collinear(sketch, *l1, *l2, row_offset),
        Constraint::Coradial(a1, a2) => tier2::eval_coradial(sketch, *a1, *a2, row_offset),
        Constraint::OffsetLine(l1, l2, d) => {
            tier2::eval_offset_line(sketch, *l1, *l2, *d, row_offset)
        }
        Constraint::OffsetArc(a1, a2, d) => {
            tier2::eval_offset_arc(sketch, *a1, *a2, *d, row_offset)
        }
        Constraint::OrdinateX(p, v) => tier2::eval_ordinate_x(sketch, *p, *v, row_offset),
        Constraint::OrdinateY(p, v) => tier2::eval_ordinate_y(sketch, *p, *v, row_offset),
        // ── Conic constraints ──────────────────────────────────────────
        Constraint::PointOnEllipse(p, e) => {
            conics::eval_point_on_ellipse(sketch, *p, *e, row_offset)
        }
        Constraint::PointOnParabola(p, pa) => {
            conics::eval_point_on_parabola(sketch, *p, *pa, row_offset)
        }
        Constraint::PointOnHyperbola(p, h) => {
            conics::eval_point_on_hyperbola(sketch, *p, *h, row_offset)
        }
        Constraint::FixedEllipse(e) => {
            conics::eval_fixed_ellipse(sketch, *e, &entry.fixed_values, row_offset)
        }
        Constraint::FixedParabola(p) => {
            conics::eval_fixed_parabola(sketch, *p, &entry.fixed_values, row_offset)
        }
        Constraint::FixedHyperbola(h) => {
            conics::eval_fixed_hyperbola(sketch, *h, &entry.fixed_values, row_offset)
        }
        Constraint::FocalDistance(e, d) => conics::eval_focal_distance(sketch, *e, *d, row_offset),
        Constraint::EllipseEccentricity(e, ecc) => {
            conics::eval_ellipse_eccentricity(sketch, *e, *ecc, row_offset)
        }
        Constraint::HyperbolaEccentricity(h, ecc) => {
            conics::eval_hyperbola_eccentricity(sketch, *h, *ecc, row_offset)
        }
        Constraint::EqualSemiAxis(e1, e2) => {
            conics::eval_equal_semi_axis(sketch, *e1, *e2, row_offset)
        }
        Constraint::ConcentricEllipse(e1, e2) => {
            conics::eval_concentric_ellipse(sketch, *e1, *e2, row_offset)
        }
        Constraint::ConcentricArcEllipse(a, e) => {
            conics::eval_concentric_arc_ellipse(sketch, *a, *e, row_offset)
        }
        Constraint::TangentLineEllipse(l, e) => {
            conics_tangent::eval_tangent_line_ellipse(sketch, *l, *e, row_offset)
        }
        Constraint::TangentLineParabola(l, p) => {
            conics_tangent::eval_tangent_line_parabola(sketch, *l, *p, row_offset)
        }
        Constraint::TangentLineHyperbola(l, h) => {
            conics_tangent::eval_tangent_line_hyperbola(sketch, *l, *h, row_offset)
        }
        // ── Arc-to-conic tangency (aux param required) ─────────────
        Constraint::TangentArcEllipse(a, e) => {
            let aux_idx = handle
                .and_then(|h| sketch.aux_param_map.get(&h).copied())
                .unwrap_or(0);
            conics_arc::eval_tangent_arc_ellipse(sketch, *a, *e, aux_idx, row_offset)
        }
        Constraint::TangentArcParabola(a, p) => {
            let aux_idx = handle
                .and_then(|h| sketch.aux_param_map.get(&h).copied())
                .unwrap_or(0);
            conics_arc::eval_tangent_arc_parabola(sketch, *a, *p, aux_idx, row_offset)
        }
        Constraint::TangentArcHyperbola(a, h) => {
            let aux_idx = handle
                .and_then(|h| sketch.aux_param_map.get(&h).copied())
                .unwrap_or(0);
            conics_arc::eval_tangent_arc_hyperbola(sketch, *a, *h, aux_idx, row_offset)
        }
        // ── Spline constraints ─────────────────────────────────────────
        Constraint::PointOnSpline(p, s) => {
            splines::eval_point_on_spline(sketch, *p, *s, handle, row_offset)
        }
        Constraint::SplineEndTangent(s, end_idx, l) => {
            splines::eval_spline_end_tangent(sketch, *s, *end_idx, *l, row_offset)
        }
        Constraint::SplineEndTangentToArc(s, end_idx, a) => {
            splines::eval_spline_end_tangent_to_arc(sketch, *s, *end_idx, *a, row_offset)
        }
        Constraint::SplineFixedEndTangent(s, end_idx, _tx, _ty) => {
            splines::eval_spline_fixed_end_tangent(
                sketch,
                *s,
                *end_idx,
                &entry.fixed_values,
                row_offset,
            )
        }
        Constraint::SplineFixedEndCurvature(s, end_idx, val) => {
            splines::eval_spline_fixed_end_curvature(sketch, *s, *end_idx, *val, row_offset)
        }
        Constraint::SplineTangentToLine(s, l) => {
            splines_tangent::eval_spline_tangent_to_line(sketch, *s, *l, handle, row_offset)
        }
        Constraint::SplineTangentToArc(s, a) => {
            splines_tangent::eval_spline_tangent_to_arc(sketch, *s, *a, handle, row_offset)
        }
        Constraint::SplineTangentToEllipse(s, e) => {
            splines_tangent_conic::eval_spline_tangent_to_ellipse(
                sketch, *s, *e, handle, row_offset,
            )
        }
        Constraint::SplineTangentToParabola(s, p) => {
            splines_tangent_conic::eval_spline_tangent_to_parabola(
                sketch, *s, *p, handle, row_offset,
            )
        }
        Constraint::SplineTangentToHyperbola(s, h) => {
            splines_tangent_conic::eval_spline_tangent_to_hyperbola(
                sketch, *s, *h, handle, row_offset,
            )
        }
        Constraint::SplineG1(s1, e1, s2, e2) => {
            splines::eval_spline_g1(sketch, *s1, *e1, *s2, *e2, row_offset)
        }
        Constraint::SplineG2(s1, e1, s2, e2) => {
            splines::eval_spline_g2(sketch, *s1, *e1, *s2, *e2, row_offset)
        }
        Constraint::SplineLength(s, len) => {
            splines::eval_spline_length(sketch, *s, *len, row_offset)
        }
        // ── Pattern constraints ────────────────────────────────────────
        Constraint::RelativePosition(p1, p2, dx, dy) => {
            patterns::eval_relative_position(sketch, *p1, *p2, *dx, *dy, row_offset)
        }
        Constraint::RotatedPosition(p1, p2, c, theta) => {
            patterns::eval_rotated_position(sketch, *p1, *p2, *c, *theta, row_offset)
        }
        Constraint::ReflectedPosition(p1, p2, l) => {
            patterns::eval_reflected_position(sketch, *p1, *p2, *l, row_offset)
        }
        Constraint::PinToCenter(p, elem) => {
            patterns::eval_pin_to_center(sketch, *p, *elem, row_offset)
        }
        Constraint::AngleOffsetArc(a1, a2, theta) => {
            patterns::eval_angle_offset_arc(sketch, *a1, *a2, *theta, row_offset)
        }
        Constraint::AngleOffsetConic(e1, e2, theta) => {
            patterns::eval_angle_offset_conic(sketch, *e1, *e2, *theta, row_offset)
        }
        Constraint::MirrorAnglesArc(a1, a2, l) => {
            patterns::eval_mirror_angles_arc(sketch, *a1, *a2, *l, row_offset)
        }
        Constraint::MirrorAnglesConic(e1, e2, l) => {
            patterns::eval_mirror_angles_conic(sketch, *e1, *e2, *l, row_offset)
        }
    }
}
