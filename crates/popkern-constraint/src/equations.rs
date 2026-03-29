//! Residual and Jacobian evaluation for geometric constraints.
//!
//! Each constraint produces a set of residual equations F(x) = 0 and
//! sparse Jacobian entries (row, col, value) for use by the Newton-Raphson
//! solver.

pub(crate) mod conics;
mod conics_tangent;
mod tier2;

#[cfg(test)]
#[path = "equations_tests.rs"]
mod tests;

use crate::constraint::{Constraint, ConstraintEntry};
use crate::elements::{ArcHandle, LineHandle, PointHandle};
use crate::Sketch;

/// Output of evaluating a single constraint: residual vector and sparse
/// Jacobian entries.
pub(crate) struct EquationOutput {
    /// Residual values; all should be zero when the constraint is satisfied.
    pub residuals: Vec<f64>,
    /// Sparse Jacobian entries as `(row_index, param_index, derivative)`.
    pub jacobian: Vec<(usize, usize, f64)>,
}

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
        | Constraint::HyperbolaEccentricity(..) => 1,
    }
}

/// Evaluate the residuals and Jacobian for a single constraint.
///
/// `row_offset` is added to every row index so that constraints can be
/// stacked into a single system.
pub(crate) fn evaluate_constraint(
    sketch: &Sketch,
    entry: &ConstraintEntry,
    row_offset: usize,
) -> EquationOutput {
    match &entry.constraint {
        Constraint::Coincident(p1, p2) => eval_coincident(sketch, *p1, *p2, row_offset),
        Constraint::Fixed(p) => eval_fixed(sketch, *p, &entry.fixed_values, row_offset),
        Constraint::FixedLine(l) => eval_fixed_line(sketch, *l, &entry.fixed_values, row_offset),
        Constraint::Horizontal(l) => eval_horizontal(sketch, *l, row_offset),
        Constraint::Vertical(l) => eval_vertical(sketch, *l, row_offset),
        Constraint::Distance(p1, p2, d) => eval_distance(sketch, *p1, *p2, *d, row_offset),
        Constraint::PointOnLine(p, l) => eval_point_on_line(sketch, *p, *l, row_offset),
        Constraint::Parallel(l1, l2) => eval_parallel(sketch, *l1, *l2, row_offset),
        Constraint::Perpendicular(l1, l2) => eval_perpendicular(sketch, *l1, *l2, row_offset),
        Constraint::EqualLength(l1, l2) => eval_equal_length(sketch, *l1, *l2, row_offset),
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
        // ── Arc-to-conic tangency (Task 5) ───────────────────────────
        Constraint::TangentArcEllipse(..) => todo!("Task 5: arc-ellipse tangency"),
        Constraint::TangentArcParabola(..) => todo!("Task 5: arc-parabola tangency"),
        Constraint::TangentArcHyperbola(..) => todo!("Task 5: arc-hyperbola tangency"),
    }
}

// ── Helpers ────────────────────────────────────────────────────────────────

/// Extract point coordinates and parameter indices from a sketch.
pub(crate) fn pt(sketch: &Sketch, h: PointHandle) -> (f64, f64, usize, usize) {
    let p = &sketch.points[h];
    let [ix, iy] = p.param_indices();
    (sketch.params[ix], sketch.params[iy], ix, iy)
}

/// Extract line endpoint coordinates and parameter indices from a sketch.
pub(crate) fn ln(
    sketch: &Sketch,
    h: LineHandle,
) -> (f64, f64, usize, usize, f64, f64, usize, usize) {
    let l = &sketch.lines[h];
    let (x1, y1, ix1, iy1) = pt(sketch, l.start);
    let (x2, y2, ix2, iy2) = pt(sketch, l.end);
    (x1, y1, ix1, iy1, x2, y2, ix2, iy2)
}

/// Extract arc centre, radius, and parameter indices from a sketch.
pub(crate) fn ar(sketch: &Sketch, h: ArcHandle) -> (f64, f64, f64, usize, usize, usize) {
    let a = &sketch.arcs[h];
    let [icx, icy, ir, _, _] = a.param_indices();
    (
        sketch.params[icx],
        sketch.params[icy],
        sketch.params[ir],
        icx,
        icy,
        ir,
    )
}

// ── Tier 1 constraint evaluators ───────────────────────────────────────────

/// Coincident(P1, P2): F0 = x1-x2, F1 = y1-y2
fn eval_coincident(
    sketch: &Sketch,
    p1: PointHandle,
    p2: PointHandle,
    row: usize,
) -> EquationOutput {
    let (x1, y1, ix1, iy1) = pt(sketch, p1);
    let (x2, y2, ix2, iy2) = pt(sketch, p2);
    EquationOutput {
        residuals: vec![x1 - x2, y1 - y2],
        jacobian: vec![
            (row, ix1, 1.0),
            (row, ix2, -1.0),
            (row + 1, iy1, 1.0),
            (row + 1, iy2, -1.0),
        ],
    }
}

/// Fixed(P): F0 = x-x0, F1 = y-y0
fn eval_fixed(sketch: &Sketch, p: PointHandle, fixed: &[f64], row: usize) -> EquationOutput {
    let (x, y, ix, iy) = pt(sketch, p);
    let x0 = fixed[0];
    let y0 = fixed[1];
    EquationOutput {
        residuals: vec![x - x0, y - y0],
        jacobian: vec![(row, ix, 1.0), (row + 1, iy, 1.0)],
    }
}

/// FixedLine(L): four equations fixing both endpoints.
fn eval_fixed_line(sketch: &Sketch, l: LineHandle, fixed: &[f64], row: usize) -> EquationOutput {
    let (x1, y1, ix1, iy1, x2, y2, ix2, iy2) = ln(sketch, l);
    let (x10, y10, x20, y20) = (fixed[0], fixed[1], fixed[2], fixed[3]);
    EquationOutput {
        residuals: vec![x1 - x10, y1 - y10, x2 - x20, y2 - y20],
        jacobian: vec![
            (row, ix1, 1.0),
            (row + 1, iy1, 1.0),
            (row + 2, ix2, 1.0),
            (row + 3, iy2, 1.0),
        ],
    }
}

/// Horizontal(L): F = y2 - y1
fn eval_horizontal(sketch: &Sketch, l: LineHandle, row: usize) -> EquationOutput {
    let (_, y1, _, iy1, _, y2, _, iy2) = ln(sketch, l);
    EquationOutput {
        residuals: vec![y2 - y1],
        jacobian: vec![(row, iy1, -1.0), (row, iy2, 1.0)],
    }
}

/// Vertical(L): F = x2 - x1
fn eval_vertical(sketch: &Sketch, l: LineHandle, row: usize) -> EquationOutput {
    let (x1, _, ix1, _, x2, _, ix2, _) = ln(sketch, l);
    EquationOutput {
        residuals: vec![x2 - x1],
        jacobian: vec![(row, ix1, -1.0), (row, ix2, 1.0)],
    }
}

/// Distance(P1, P2, d): F = (x2-x1)^2 + (y2-y1)^2 - d^2
fn eval_distance(
    sketch: &Sketch,
    p1: PointHandle,
    p2: PointHandle,
    d: f64,
    row: usize,
) -> EquationOutput {
    let (x1, y1, ix1, iy1) = pt(sketch, p1);
    let (x2, y2, ix2, iy2) = pt(sketch, p2);
    let dx = x2 - x1;
    let dy = y2 - y1;
    EquationOutput {
        residuals: vec![dx * dx + dy * dy - d * d],
        jacobian: vec![
            (row, ix1, -2.0 * dx),
            (row, iy1, -2.0 * dy),
            (row, ix2, 2.0 * dx),
            (row, iy2, 2.0 * dy),
        ],
    }
}

/// PointOnLine(P, L): F = (px-x1)(y2-y1) - (py-y1)(x2-x1)
fn eval_point_on_line(
    sketch: &Sketch,
    p: PointHandle,
    l: LineHandle,
    row: usize,
) -> EquationOutput {
    let (px, py, ipx, ipy) = pt(sketch, p);
    let (x1, y1, ix1, iy1, x2, y2, ix2, iy2) = ln(sketch, l);
    let dx = x2 - x1;
    let dy = y2 - y1;
    // F = (px - x1) * dy - (py - y1) * dx
    let f = (px - x1) * dy - (py - y1) * dx;
    EquationOutput {
        residuals: vec![f],
        jacobian: vec![
            (row, ipx, dy),         // dF/dpx = y2 - y1
            (row, ipy, -dx),        // dF/dpy = -(x2 - x1)
            (row, ix1, py - y2),    // dF/dx1 = -dy + (py-y1)*0 ... = py - y2
            (row, iy1, x2 - px),    // dF/dy1 = -(px-x1) + dx  ... = x2 - px
            (row, ix2, -(py - y1)), // dF/dx2 = -(py - y1)
            (row, iy2, px - x1),    // dF/dy2 = px - x1
        ],
    }
}

/// Parallel(L1, L2): F = dx1*dy2 - dy1*dx2
fn eval_parallel(sketch: &Sketch, l1: LineHandle, l2: LineHandle, row: usize) -> EquationOutput {
    let (x11, y11, ix11, iy11, x12, y12, ix12, iy12) = ln(sketch, l1);
    let (x21, y21, ix21, iy21, x22, y22, ix22, iy22) = ln(sketch, l2);
    let dx1 = x12 - x11;
    let dy1 = y12 - y11;
    let dx2 = x22 - x21;
    let dy2 = y22 - y21;
    // F = dx1*dy2 - dy1*dx2
    EquationOutput {
        residuals: vec![dx1 * dy2 - dy1 * dx2],
        jacobian: vec![
            (row, ix11, -dy2), // d/dx11: -1 * dy2
            (row, iy11, dx2),  // d/dy11: -1 * (-dx2) = dx2
            (row, ix12, dy2),  // d/dx12: 1 * dy2
            (row, iy12, -dx2), // d/dy12: 1 * (-dx2)
            (row, ix21, dy1),  // d/dx21: -1 * (-dy1) = dy1
            (row, iy21, -dx1), // d/dy21: -1 * dx1
            (row, ix22, -dy1), // d/dx22: 1 * (-dy1)
            (row, iy22, dx1),  // d/dy22: dx1
        ],
    }
}

/// Perpendicular(L1, L2): F = dx1*dx2 + dy1*dy2
fn eval_perpendicular(
    sketch: &Sketch,
    l1: LineHandle,
    l2: LineHandle,
    row: usize,
) -> EquationOutput {
    let (x11, y11, ix11, iy11, x12, y12, ix12, iy12) = ln(sketch, l1);
    let (x21, y21, ix21, iy21, x22, y22, ix22, iy22) = ln(sketch, l2);
    let dx1 = x12 - x11;
    let dy1 = y12 - y11;
    let dx2 = x22 - x21;
    let dy2 = y22 - y21;
    // F = dx1*dx2 + dy1*dy2
    EquationOutput {
        residuals: vec![dx1 * dx2 + dy1 * dy2],
        jacobian: vec![
            (row, ix11, -dx2), // d/dx11
            (row, iy11, -dy2), // d/dy11
            (row, ix12, dx2),  // d/dx12
            (row, iy12, dy2),  // d/dy12
            (row, ix21, -dx1), // d/dx21
            (row, iy21, -dy1), // d/dy21
            (row, ix22, dx1),  // d/dx22
            (row, iy22, dy1),  // d/dy22
        ],
    }
}

/// EqualLength(L1, L2): F = |L1|^2 - |L2|^2
fn eval_equal_length(
    sketch: &Sketch,
    l1: LineHandle,
    l2: LineHandle,
    row: usize,
) -> EquationOutput {
    let (x11, y11, ix11, iy11, x12, y12, ix12, iy12) = ln(sketch, l1);
    let (x21, y21, ix21, iy21, x22, y22, ix22, iy22) = ln(sketch, l2);
    let dx1 = x12 - x11;
    let dy1 = y12 - y11;
    let dx2 = x22 - x21;
    let dy2 = y22 - y21;
    EquationOutput {
        residuals: vec![dx1 * dx1 + dy1 * dy1 - dx2 * dx2 - dy2 * dy2],
        jacobian: vec![
            (row, ix11, -2.0 * dx1),
            (row, iy11, -2.0 * dy1),
            (row, ix12, 2.0 * dx1),
            (row, iy12, 2.0 * dy1),
            (row, ix21, 2.0 * dx2),
            (row, iy21, 2.0 * dy2),
            (row, ix22, -2.0 * dx2),
            (row, iy22, -2.0 * dy2),
        ],
    }
}
