//! Constraint evaluators for pattern operations.
//!
//! These constraints express positional and angular relationships used by
//! linear, circular, and mirror patterns.

use super::{ln, pt, EquationOutput};
use crate::elements::{ArcHandle, ElementHandle, LineHandle, PointHandle};
use crate::Sketch;

// ── Position primitives ────────────────────────────────────────────────────

/// RelativePosition(P1, P2, dx, dy):
/// F1 = p2x - p1x - dx = 0
/// F2 = p2y - p1y - dy = 0
pub(crate) fn eval_relative_position(
    sketch: &Sketch,
    p1: PointHandle,
    p2: PointHandle,
    dx: f64,
    dy: f64,
    row: usize,
) -> EquationOutput {
    let (p1x, p1y, ip1x, ip1y) = pt(sketch, p1);
    let (p2x, p2y, ip2x, ip2y) = pt(sketch, p2);
    EquationOutput {
        residuals: vec![p2x - p1x - dx, p2y - p1y - dy],
        jacobian: vec![
            (row, ip1x, -1.0),
            (row, ip2x, 1.0),
            (row + 1, ip1y, -1.0),
            (row + 1, ip2y, 1.0),
        ],
    }
}

/// RotatedPosition(P1, P2, Center, theta):
/// rx = cx + (p1x-cx)*cos(theta) - (p1y-cy)*sin(theta)
/// ry = cy + (p1x-cx)*sin(theta) + (p1y-cy)*cos(theta)
/// F1 = p2x - rx = 0
/// F2 = p2y - ry = 0
pub(crate) fn eval_rotated_position(
    sketch: &Sketch,
    p1: PointHandle,
    p2: PointHandle,
    center: PointHandle,
    theta: f64,
    row: usize,
) -> EquationOutput {
    let (p1x, p1y, ip1x, ip1y) = pt(sketch, p1);
    let (p2x, p2y, ip2x, ip2y) = pt(sketch, p2);
    let (cx, cy, icx, icy) = pt(sketch, center);

    let cos_t = theta.cos();
    let sin_t = theta.sin();
    let dx = p1x - cx;
    let dy = p1y - cy;

    let rx = cx + dx * cos_t - dy * sin_t;
    let ry = cy + dx * sin_t + dy * cos_t;

    let f1 = p2x - rx;
    let f2 = p2y - ry;

    // Jacobian entries:
    // dF1/dp1x = -cos_t,  dF1/dp1y = sin_t,   dF1/dp2x = 1
    // dF1/dcx  = -(1 - cos_t) = cos_t - 1,  (from d/dcx: 1 + (-1)*cos_t = 1 - cos_t; negated)
    //   Actually: rx = cx + (p1x-cx)*cos_t - (p1y-cy)*sin_t
    //   drx/dcx = 1 - cos_t + 0_from_sin_term... wait:
    //   drx/dcx = 1 + d/dcx[(p1x-cx)*cos_t] + d/dcx[-(p1y-cy)*sin_t]
    //           = 1 - cos_t + sin_t * d(cy)/dcx ... no, cy is independent of cx
    //           = 1 - cos_t
    //   So dF1/dcx = -(1 - cos_t) = cos_t - 1
    //   drx/dcy = 0 + 0 - (-1)*sin_t = sin_t
    //   So dF1/dcy = -sin_t
    //
    // dF2/dp1x = -sin_t,  dF2/dp1y = -cos_t,  dF2/dp2y = 1
    //   dry/dcx = 0 - sin_t + 0 = -sin_t  ... wait:
    //   dry/dcx = d/dcx[cy + (p1x-cx)*sin_t + (p1y-cy)*cos_t]
    //           = 0 - sin_t + 0 = -sin_t
    //   So dF2/dcx = sin_t
    //   dry/dcy = 1 + 0 - cos_t = 1 - cos_t
    //   So dF2/dcy = -(1 - cos_t) = cos_t - 1

    EquationOutput {
        residuals: vec![f1, f2],
        jacobian: vec![
            (row, ip1x, -cos_t),
            (row, ip1y, sin_t),
            (row, ip2x, 1.0),
            (row, icx, cos_t - 1.0),
            (row, icy, -sin_t),
            (row + 1, ip1x, -sin_t),
            (row + 1, ip1y, -cos_t),
            (row + 1, ip2y, 1.0),
            (row + 1, icx, sin_t),
            (row + 1, icy, cos_t - 1.0),
        ],
    }
}

/// ReflectedPosition(P1, P2, Line):
/// Reflect P1 across the infinite line through Line's endpoints.
///
/// Let d = (lx2-lx1, ly2-ly1), len2 = d.d
/// foot = line_start + dot(p1-line_start, d)/len2 * d
/// reflected = 2*foot - p1
/// F1 = p2x - reflected_x = 0
/// F2 = p2y - reflected_y = 0
pub(crate) fn eval_reflected_position(
    sketch: &Sketch,
    p1: PointHandle,
    p2: PointHandle,
    l: LineHandle,
    row: usize,
) -> EquationOutput {
    let (p1x, p1y, ip1x, ip1y) = pt(sketch, p1);
    let (p2x, p2y, ip2x, ip2y) = pt(sketch, p2);
    let (lx1, ly1, ilx1, ily1, lx2, ly2, ilx2, ily2) = ln(sketch, l);

    let dx = lx2 - lx1;
    let dy = ly2 - ly1;
    let len2 = dx * dx + dy * dy;

    // t = dot(p1 - line_start, d) / len2
    let vx = p1x - lx1;
    let vy = p1y - ly1;
    let t = (vx * dx + vy * dy) / len2;

    // foot = (lx1 + t*dx, ly1 + t*dy)
    let fx = lx1 + t * dx;
    let fy = ly1 + t * dy;

    // reflected = 2*foot - p1
    let rx = 2.0 * fx - p1x;
    let ry = 2.0 * fy - p1y;

    let f1 = p2x - rx;
    let f2 = p2y - ry;

    // Jacobian: we need derivatives of rx, ry w.r.t. all 8 params
    // (p1x, p1y, p2x, p2y, lx1, ly1, lx2, ly2).
    //
    // rx = 2*(lx1 + t*dx) - p1x  where t = (vx*dx + vy*dy) / len2
    // Use chain rule through t.
    //
    // dt/dp1x = dx / len2
    // dt/dp1y = dy / len2
    // dt/dlx1 = [(-dx)*dx + vx*(-dx) + (-dy)*dy + vy*(-dy) ... no, careful:
    //   t = ((p1x-lx1)*dx + (p1y-ly1)*dy) / len2
    //     = ((p1x-lx1)*(lx2-lx1) + (p1y-ly1)*(ly2-ly1)) / ((lx2-lx1)^2 + (ly2-ly1)^2)
    //
    // This is complex, so compute numerator N and denominator D:
    //   N = vx*dx + vy*dy
    //   D = dx^2 + dy^2 = len2
    //   t = N/D
    //
    // Then rx = 2*(lx1 + (N/D)*dx) - p1x
    //         = 2*lx1 + 2*N*dx/D - p1x
    //
    // drx/dvar = 2*dlx1/dvar + 2*(dN/dvar * dx * D + N * ddx/dvar * D - N * dx * dD/dvar) / D^2 - dp1x/dvar
    //
    // Let's just compute it piece by piece using the quotient rule for t.

    let inv_len2 = 1.0 / len2;
    let inv_len4 = inv_len2 * inv_len2;
    let nn = vx * dx + vy * dy; // numerator of t

    // Helper: drx/dvar = 2*(dlx1/dvar + dt/dvar * dx + t * ddx/dvar) - dp1x/dvar
    // Similarly for ry.
    // dt/dvar = (dN/dvar * D - N * dD/dvar) / D^2

    // For each variable, compute (dN, dD, dlx1, ddx, dp1x, dly1, ddy, dp1y):
    // p1x: dN=dx, dD=0, dlx1=0, ddx=0, dp1x=1, dly1=0, ddy=0, dp1y=0
    // p1y: dN=dy, dD=0, dlx1=0, ddx=0, dp1x=0, dly1=0, ddy=0, dp1y=1
    // lx1: dN = -dx + vx*(-1) + 0 = -dx - vx + vx = ... careful:
    //   N = (p1x-lx1)*(lx2-lx1) + (p1y-ly1)*(ly2-ly1)
    //   dN/dlx1 = (-1)*(lx2-lx1) + (p1x-lx1)*(-1) + 0 = -dx - vx
    //   dD/dlx1 = 2*(lx2-lx1)*(-1) = -2*dx
    //   dlx1/dlx1 = 1
    //   ddx/dlx1 = -1

    // Let me just compute dt/dvar for each, then drx, dry.

    // dt/dp1x = dx * inv_len2
    // dt/dp1y = dy * inv_len2
    // dt/dlx1 = ((-dx - vx) * len2 - nn * (-2.0*dx)) * inv_len4
    //         = (-dx - vx + 2*nn*dx*inv_len2) * inv_len2
    //         hmm let me keep it cleaner:

    // dN/dlx1 = -(dx + vx),  dD/dlx1 = -2*dx
    let dt_dlx1 = (-(dx + vx) * len2 - nn * (-2.0 * dx)) * inv_len4;
    // dN/dly1 = -(dy + vy),  dD/dly1 = -2*dy
    let dt_dly1 = (-(dy + vy) * len2 - nn * (-2.0 * dy)) * inv_len4;
    // dN/dlx2 = vx,  dD/dlx2 = 2*dx
    let dt_dlx2 = (vx * len2 - nn * 2.0 * dx) * inv_len4;
    // dN/dly2 = vy,  dD/dly2 = 2*dy
    let dt_dly2 = (vy * len2 - nn * 2.0 * dy) * inv_len4;

    let dt_dp1x = dx * inv_len2;
    let dt_dp1y = dy * inv_len2;

    // drx/dvar = 2*(dlx1/dvar + dt/dvar * dx + t * ddx/dvar) - dp1x/dvar
    // ddx/dlx1 = -1, ddx/dlx2 = 1, ddx/others = 0
    // dlx1/dlx1 = 1, dlx1/others = 0
    // dp1x/dp1x = 1, dp1x/others = 0

    let drx_dp1x = 2.0 * dt_dp1x * dx - 1.0;
    let drx_dp1y = 2.0 * dt_dp1y * dx;
    let drx_dlx1 = 2.0 * (1.0 + dt_dlx1 * dx - t);
    let drx_dly1 = 2.0 * dt_dly1 * dx;
    let drx_dlx2 = 2.0 * (dt_dlx2 * dx + t);
    let drx_dly2 = 2.0 * dt_dly2 * dx;

    // dry/dvar = 2*(dly1/dvar + dt/dvar * dy + t * ddy/dvar) - dp1y/dvar
    let dry_dp1x = 2.0 * dt_dp1x * dy;
    let dry_dp1y = 2.0 * dt_dp1y * dy - 1.0;
    let dry_dlx1 = 2.0 * dt_dlx1 * dy;
    let dry_dly1 = 2.0 * (1.0 + dt_dly1 * dy - t);
    let dry_dlx2 = 2.0 * dt_dlx2 * dy;
    let dry_dly2 = 2.0 * (dt_dly2 * dy + t);

    // dF1/dvar = -drx/dvar (except dF1/dp2x = 1)
    // dF2/dvar = -dry/dvar (except dF2/dp2y = 1)
    EquationOutput {
        residuals: vec![f1, f2],
        jacobian: vec![
            (row, ip1x, -drx_dp1x),
            (row, ip1y, -drx_dp1y),
            (row, ip2x, 1.0),
            (row, ilx1, -drx_dlx1),
            (row, ily1, -drx_dly1),
            (row, ilx2, -drx_dlx2),
            (row, ily2, -drx_dly2),
            (row + 1, ip1x, -dry_dp1x),
            (row + 1, ip1y, -dry_dp1y),
            (row + 1, ip2y, 1.0),
            (row + 1, ilx1, -dry_dlx1),
            (row + 1, ily1, -dry_dly1),
            (row + 1, ilx2, -dry_dlx2),
            (row + 1, ily2, -dry_dly2),
        ],
    }
}

// ── PinToCenter + angle constraints ────────────────────────────────────────

/// PinToCenter(Point, ElementHandle): pin a point to an element's centre.
///
/// Dispatches on element type to find centre/vertex params.
/// F1 = px - cx = 0,  F2 = py - cy = 0.
pub(crate) fn eval_pin_to_center(
    sketch: &Sketch,
    p: PointHandle,
    elem: ElementHandle,
    row: usize,
) -> EquationOutput {
    let (px, py, ipx, ipy) = pt(sketch, p);
    let (cx, cy, icx, icy) = element_center(sketch, elem);
    EquationOutput {
        residuals: vec![px - cx, py - cy],
        jacobian: vec![
            (row, ipx, 1.0),
            (row, icx, -1.0),
            (row + 1, ipy, 1.0),
            (row + 1, icy, -1.0),
        ],
    }
}

/// Extract the centre (or vertex) coordinates and parameter indices for a
/// curved element.
fn element_center(sketch: &Sketch, elem: ElementHandle) -> (f64, f64, usize, usize) {
    match elem {
        ElementHandle::Arc(ah) => {
            let a = &sketch.arcs[ah];
            let s = a.param_start;
            (sketch.params[s], sketch.params[s + 1], s, s + 1)
        }
        ElementHandle::Ellipse(eh) => {
            let e = &sketch.ellipses[eh];
            let s = e.param_start;
            (sketch.params[s], sketch.params[s + 1], s, s + 1)
        }
        ElementHandle::Parabola(ph) => {
            let p = &sketch.parabolas[ph];
            let s = p.param_start;
            (sketch.params[s], sketch.params[s + 1], s, s + 1)
        }
        ElementHandle::Hyperbola(hh) => {
            let h = &sketch.hyperbolas[hh];
            let s = h.param_start;
            (sketch.params[s], sketch.params[s + 1], s, s + 1)
        }
        // Points, lines, and splines don't have a centre in this context.
        // Return zero with index 0 as a fallback (should not be reached
        // if constraints are validated at creation time).
        _ => (0.0, 0.0, 0, 0),
    }
}

/// AngleOffsetArc(Arc1, Arc2, theta):
/// F1 = arc2.start_angle - arc1.start_angle - theta = 0
/// F2 = arc2.end_angle - arc1.end_angle - theta = 0
pub(crate) fn eval_angle_offset_arc(
    sketch: &Sketch,
    a1: ArcHandle,
    a2: ArcHandle,
    theta: f64,
    row: usize,
) -> EquationOutput {
    let arc1 = &sketch.arcs[a1];
    let arc2 = &sketch.arcs[a2];
    let is1 = arc1.param_start + 3; // start_angle
    let ie1 = arc1.param_start + 4; // end_angle
    let is2 = arc2.param_start + 3;
    let ie2 = arc2.param_start + 4;

    let sa1 = sketch.params[is1];
    let ea1 = sketch.params[ie1];
    let sa2 = sketch.params[is2];
    let ea2 = sketch.params[ie2];

    EquationOutput {
        residuals: vec![sa2 - sa1 - theta, ea2 - ea1 - theta],
        jacobian: vec![
            (row, is1, -1.0),
            (row, is2, 1.0),
            (row + 1, ie1, -1.0),
            (row + 1, ie2, 1.0),
        ],
    }
}

/// AngleOffsetConic(E1, E2, theta):
/// F = e2.rotation - e1.rotation - theta = 0
///
/// Rotation parameter index depends on element type:
/// - Ellipse/Hyperbola: param_start + 4
/// - Parabola: param_start + 3
pub(crate) fn eval_angle_offset_conic(
    sketch: &Sketch,
    e1: ElementHandle,
    e2: ElementHandle,
    theta: f64,
    row: usize,
) -> EquationOutput {
    let (rot1, ir1) = element_rotation(sketch, e1);
    let (rot2, ir2) = element_rotation(sketch, e2);
    EquationOutput {
        residuals: vec![rot2 - rot1 - theta],
        jacobian: vec![(row, ir1, -1.0), (row, ir2, 1.0)],
    }
}

/// MirrorAnglesArc(Arc1, Arc2, Line):
/// line_angle = atan2(ly2-ly1, lx2-lx1)
/// F1 = arc2.start - (2*line_angle - arc1.end) = 0
/// F2 = arc2.end - (2*line_angle - arc1.start) = 0
///
/// Note: start/end swap to preserve winding under reflection.
pub(crate) fn eval_mirror_angles_arc(
    sketch: &Sketch,
    a1: ArcHandle,
    a2: ArcHandle,
    l: LineHandle,
    row: usize,
) -> EquationOutput {
    let (lx1, ly1, ilx1, ily1, lx2, ly2, ilx2, ily2) = ln(sketch, l);
    let ldx = lx2 - lx1;
    let ldy = ly2 - ly1;
    let line_angle = ldy.atan2(ldx);

    let arc1 = &sketch.arcs[a1];
    let arc2 = &sketch.arcs[a2];
    let is1 = arc1.param_start + 3;
    let ie1 = arc1.param_start + 4;
    let is2 = arc2.param_start + 3;
    let ie2 = arc2.param_start + 4;

    let sa1 = sketch.params[is1];
    let ea1 = sketch.params[ie1];
    let sa2 = sketch.params[is2];
    let ea2 = sketch.params[ie2];

    let f1 = sa2 - (2.0 * line_angle - ea1);
    let f2 = ea2 - (2.0 * line_angle - sa1);

    // d(line_angle)/d(ldx,ldy) via atan2:
    // d(atan2(ldy,ldx))/dldx = -ldy / (ldx^2 + ldy^2)
    // d(atan2(ldy,ldx))/dldy = ldx / (ldx^2 + ldy^2)
    let denom = ldx * ldx + ldy * ldy;
    let da_dldx = -ldy / denom;
    let da_dldy = ldx / denom;

    // Chain to line endpoints:
    // dldx/dlx1 = -1, dldx/dlx2 = 1
    // dldy/dly1 = -1, dldy/dly2 = 1
    let da_dlx1 = -da_dldx;
    let da_dly1 = -da_dldy;
    let da_dlx2 = da_dldx;
    let da_dly2 = da_dldy;

    // dF1/dvar: dF1/dsa2 = 1, dF1/dea1 = 1 (from -(-ea1)), dF1/dline_angle = -2
    // dF2/dvar: dF2/dea2 = 1, dF2/dsa1 = 1 (from -(-sa1)), dF2/dline_angle = -2

    EquationOutput {
        residuals: vec![f1, f2],
        jacobian: vec![
            (row, is2, 1.0),
            (row, ie1, 1.0),
            (row, ilx1, -2.0 * da_dlx1),
            (row, ily1, -2.0 * da_dly1),
            (row, ilx2, -2.0 * da_dlx2),
            (row, ily2, -2.0 * da_dly2),
            (row + 1, ie2, 1.0),
            (row + 1, is1, 1.0),
            (row + 1, ilx1, -2.0 * da_dlx1),
            (row + 1, ily1, -2.0 * da_dly1),
            (row + 1, ilx2, -2.0 * da_dlx2),
            (row + 1, ily2, -2.0 * da_dly2),
        ],
    }
}

/// MirrorAnglesConic(E1, E2, Line):
/// F = e2.rotation - (2*line_angle - e1.rotation) = 0
pub(crate) fn eval_mirror_angles_conic(
    sketch: &Sketch,
    e1: ElementHandle,
    e2: ElementHandle,
    l: LineHandle,
    row: usize,
) -> EquationOutput {
    let (lx1, ly1, ilx1, ily1, lx2, ly2, ilx2, ily2) = ln(sketch, l);
    let ldx = lx2 - lx1;
    let ldy = ly2 - ly1;
    let line_angle = ldy.atan2(ldx);

    let (rot1, ir1) = element_rotation(sketch, e1);
    let (rot2, ir2) = element_rotation(sketch, e2);

    let f = rot2 - (2.0 * line_angle - rot1);

    let denom = ldx * ldx + ldy * ldy;
    let da_dldx = -ldy / denom;
    let da_dldy = ldx / denom;

    let da_dlx1 = -da_dldx;
    let da_dly1 = -da_dldy;
    let da_dlx2 = da_dldx;
    let da_dly2 = da_dldy;

    // dF/dr1 = 1, dF/dr2 = 1, dF/dline_angle = -2
    EquationOutput {
        residuals: vec![f],
        jacobian: vec![
            (row, ir1, 1.0),
            (row, ir2, 1.0),
            (row, ilx1, -2.0 * da_dlx1),
            (row, ily1, -2.0 * da_dly1),
            (row, ilx2, -2.0 * da_dlx2),
            (row, ily2, -2.0 * da_dly2),
        ],
    }
}

/// Extract the rotation parameter value and index for a conic element.
fn element_rotation(sketch: &Sketch, elem: ElementHandle) -> (f64, usize) {
    match elem {
        ElementHandle::Ellipse(eh) => {
            let e = &sketch.ellipses[eh];
            let i = e.param_start + 4;
            (sketch.params[i], i)
        }
        ElementHandle::Parabola(ph) => {
            let p = &sketch.parabolas[ph];
            let i = p.param_start + 3;
            (sketch.params[i], i)
        }
        ElementHandle::Hyperbola(hh) => {
            let h = &sketch.hyperbolas[hh];
            let i = h.param_start + 4;
            (sketch.params[i], i)
        }
        // Fallback for non-conic types (should not be reached with
        // validated constraints).
        _ => (0.0, 0),
    }
}
