//! Arc-related and tier-2 constraint evaluators.
//!
//! Split from the main `equations` module to stay within the 500-line limit.

use super::{EquationOutput, ar, ln, pt};
use crate::Sketch;
use crate::elements::{ArcHandle, LineHandle, PointHandle};

/// PointOnArc(P, A): F = (px-cx)^2 + (py-cy)^2 - r^2
pub(crate) fn eval_point_on_arc(
    sketch: &Sketch,
    p: PointHandle,
    a: ArcHandle,
    row: usize,
) -> EquationOutput {
    let (px, py, ipx, ipy) = pt(sketch, p);
    let (cx, cy, r, icx, icy, ir) = ar(sketch, a);
    let dx = px - cx;
    let dy = py - cy;
    EquationOutput {
        residuals: vec![dx * dx + dy * dy - r * r],
        jacobian: vec![
            (row, ipx, 2.0 * dx),
            (row, ipy, 2.0 * dy),
            (row, icx, -2.0 * dx),
            (row, icy, -2.0 * dy),
            (row, ir, -2.0 * r),
        ],
    }
}

/// PointLineDistance(P, L, d):
/// F = [(px-x1)(y2-y1) - (py-y1)(x2-x1)]^2 / [(x2-x1)^2 + (y2-y1)^2] - d^2
///
/// We use the quotient form: F = N^2 / D - d^2, where
///   N = (px-x1)(y2-y1) - (py-y1)(x2-x1)
///   D = (x2-x1)^2 + (y2-y1)^2
///
/// dF/dv = (2*N*dN/dv * D - N^2 * dD/dv) / D^2
pub(crate) fn eval_point_line_distance(
    sketch: &Sketch,
    p: PointHandle,
    l: LineHandle,
    d: f64,
    row: usize,
) -> EquationOutput {
    let (px, py, ipx, ipy) = pt(sketch, p);
    let (x1, y1, ix1, iy1, x2, y2, ix2, iy2) = ln(sketch, l);
    let ldx = x2 - x1;
    let ldy = y2 - y1;
    let n = (px - x1) * ldy - (py - y1) * ldx;
    let denom = ldx * ldx + ldy * ldy;
    let f = n * n / denom - d * d;
    // Partial derivatives of N:
    //   dN/dpx = ldy, dN/dpy = -ldx
    //   dN/dx1 = py - y2, dN/dy1 = x2 - px
    //   dN/dx2 = -(py - y1), dN/dy2 = px - x1
    // Partial derivatives of D:
    //   dD/dx1 = -2*ldx, dD/dy1 = -2*ldy
    //   dD/dx2 = 2*ldx, dD/dy2 = 2*ldy
    //   dD/dpx = 0, dD/dpy = 0
    let d2 = denom * denom;
    let two_n = 2.0 * n;

    let df = |dn: f64, dd: f64| -> f64 { (two_n * dn * denom - n * n * dd) / d2 };

    EquationOutput {
        residuals: vec![f],
        jacobian: vec![
            (row, ipx, df(ldy, 0.0)),
            (row, ipy, df(-ldx, 0.0)),
            (row, ix1, df(py - y2, -2.0 * ldx)),
            (row, iy1, df(x2 - px, -2.0 * ldy)),
            (row, ix2, df(-(py - y1), 2.0 * ldx)),
            (row, iy2, df(px - x1, 2.0 * ldy)),
        ],
    }
}

/// Angle(L1, L2, theta):
/// F = atan2(dx1*dy2 - dy1*dx2, dx1*dx2 + dy1*dy2) - theta
///
/// Let cross = dx1*dy2 - dy1*dx2, dot = dx1*dx2 + dy1*dy2.
/// F = atan2(cross, dot) - theta.
/// dF/dv = (dot * d(cross)/dv - cross * d(dot)/dv) / (cross^2 + dot^2)
pub(crate) fn eval_angle(
    sketch: &Sketch,
    l1: LineHandle,
    l2: LineHandle,
    theta: f64,
    row: usize,
) -> EquationOutput {
    let (x11, y11, ix11, iy11, x12, y12, ix12, iy12) = ln(sketch, l1);
    let (x21, y21, ix21, iy21, x22, y22, ix22, iy22) = ln(sketch, l2);
    let dx1 = x12 - x11;
    let dy1 = y12 - y11;
    let dx2 = x22 - x21;
    let dy2 = y22 - y21;
    let cross = dx1 * dy2 - dy1 * dx2;
    let dot = dx1 * dx2 + dy1 * dy2;
    let f = cross.atan2(dot) - theta;
    let denom = cross * cross + dot * dot;

    // d(cross)/dv and d(dot)/dv for each variable:
    // dc/dx11 = -dy2, dd/dx11 = -dx2
    // dc/dy11 = dx2,  dd/dy11 = -dy2
    // dc/dx12 = dy2,  dd/dx12 = dx2
    // dc/dy12 = -dx2, dd/dy12 = dy2
    // dc/dx21 = dy1,  dd/dx21 = -dx1
    // dc/dy21 = -dx1, dd/dy21 = -dy1
    // dc/dx22 = -dy1, dd/dx22 = dx1
    // dc/dy22 = dx1,  dd/dy22 = dy1
    let df = |dc: f64, dd_val: f64| -> f64 { (dot * dc - cross * dd_val) / denom };

    EquationOutput {
        residuals: vec![f],
        jacobian: vec![
            (row, ix11, df(-dy2, -dx2)),
            (row, iy11, df(dx2, -dy2)),
            (row, ix12, df(dy2, dx2)),
            (row, iy12, df(-dx2, dy2)),
            (row, ix21, df(dy1, -dx1)),
            (row, iy21, df(-dx1, -dy1)),
            (row, ix22, df(-dy1, dx1)),
            (row, iy22, df(dx1, dy1)),
        ],
    }
}

/// Tangent(L, A): line tangent to arc.
/// F = [(cx-x1)(y2-y1) - (cy-y1)(x2-x1)]^2 / [(x2-x1)^2 + (y2-y1)^2] - r^2
///
/// Same quotient form as PointLineDistance with point = arc centre.
pub(crate) fn eval_tangent(
    sketch: &Sketch,
    l: LineHandle,
    a: ArcHandle,
    row: usize,
) -> EquationOutput {
    let (cx, cy, r, icx, icy, ir) = ar(sketch, a);
    let (x1, y1, ix1, iy1, x2, y2, ix2, iy2) = ln(sketch, l);
    let ldx = x2 - x1;
    let ldy = y2 - y1;
    let n = (cx - x1) * ldy - (cy - y1) * ldx;
    let denom = ldx * ldx + ldy * ldy;
    let f = n * n / denom - r * r;
    let d2 = denom * denom;
    let two_n = 2.0 * n;

    let df = |dn: f64, dd: f64| -> f64 { (two_n * dn * denom - n * n * dd) / d2 };

    // dN/dcx = ldy, dN/dcy = -ldx
    // dN/dx1 = cy - y2, dN/dy1 = x2 - cx
    // dN/dx2 = -(cy - y1), dN/dy2 = cx - x1
    EquationOutput {
        residuals: vec![f],
        jacobian: vec![
            (row, icx, df(ldy, 0.0)),
            (row, icy, df(-ldx, 0.0)),
            (row, ir, -2.0 * r),
            (row, ix1, df(cy - y2, -2.0 * ldx)),
            (row, iy1, df(x2 - cx, -2.0 * ldy)),
            (row, ix2, df(-(cy - y1), 2.0 * ldx)),
            (row, iy2, df(cx - x1, 2.0 * ldy)),
        ],
    }
}

/// ArcTangent(A1, A2): external tangency.
/// F = (cx1-cx2)^2 + (cy1-cy2)^2 - (r1+r2)^2
pub(crate) fn eval_arc_tangent(
    sketch: &Sketch,
    a1: ArcHandle,
    a2: ArcHandle,
    row: usize,
) -> EquationOutput {
    let (cx1, cy1, r1, icx1, icy1, ir1) = ar(sketch, a1);
    let (cx2, cy2, r2, icx2, icy2, ir2) = ar(sketch, a2);
    let dx = cx1 - cx2;
    let dy = cy1 - cy2;
    let rsum = r1 + r2;
    EquationOutput {
        residuals: vec![dx * dx + dy * dy - rsum * rsum],
        jacobian: vec![
            (row, icx1, 2.0 * dx),
            (row, icy1, 2.0 * dy),
            (row, ir1, -2.0 * rsum),
            (row, icx2, -2.0 * dx),
            (row, icy2, -2.0 * dy),
            (row, ir2, -2.0 * rsum),
        ],
    }
}

/// EqualRadius(A1, A2): F = r1 - r2
pub(crate) fn eval_equal_radius(
    sketch: &Sketch,
    a1: ArcHandle,
    a2: ArcHandle,
    row: usize,
) -> EquationOutput {
    let (_, _, r1, _, _, ir1) = ar(sketch, a1);
    let (_, _, r2, _, _, ir2) = ar(sketch, a2);
    EquationOutput {
        residuals: vec![r1 - r2],
        jacobian: vec![(row, ir1, 1.0), (row, ir2, -1.0)],
    }
}

/// Symmetric(P1, P2, L): midpoint on line + perpendicular.
/// F0 = (mx-lx1)*dy - (my-ly1)*dx  (midpoint on line)
/// F1 = (x2-x1)*dx + (y2-y1)*dy    (perpendicular)
/// where mx=(px1+px2)/2, my=(py1+py2)/2, dx=lx2-lx1, dy=ly2-ly1
pub(crate) fn eval_symmetric(
    sketch: &Sketch,
    p1: PointHandle,
    p2: PointHandle,
    l: LineHandle,
    row: usize,
) -> EquationOutput {
    let (px1, py1, ipx1, ipy1) = pt(sketch, p1);
    let (px2, py2, ipx2, ipy2) = pt(sketch, p2);
    let (lx1, ly1, ilx1, ily1, lx2, ly2, ilx2, ily2) = ln(sketch, l);
    let mx = (px1 + px2) / 2.0;
    let my = (py1 + py2) / 2.0;
    let ldx = lx2 - lx1;
    let ldy = ly2 - ly1;

    // F0 = (mx - lx1) * ldy - (my - ly1) * ldx
    let f0 = (mx - lx1) * ldy - (my - ly1) * ldx;
    // F1 = (px2 - px1) * ldx + (py2 - py1) * ldy
    let pdx = px2 - px1;
    let pdy = py2 - py1;
    let f1 = pdx * ldx + pdy * ldy;

    // F0 partials:
    // dF0/dpx1 = 0.5 * ldy, dF0/dpy1 = -0.5 * ldx
    // dF0/dpx2 = 0.5 * ldy, dF0/dpy2 = -0.5 * ldx
    // dF0/dlx1 = -ldy + (my - ly1) = my - ly2
    //   Actually: dF0/dlx1 = d/dlx1[(mx-lx1)*ldy - (my-ly1)*ldx]
    //     = -1*ldy + (mx-lx1)*(-1)_from_ldy_no... wait.
    //   ldy = ly2-ly1, ldx = lx2-lx1
    //   dF0/dlx1 = (mx-lx1)*d(ldy)/dlx1 + d(mx-lx1)/dlx1 * ldy
    //            - [(my-ly1)*d(ldx)/dlx1 + d(my-ly1)/dlx1 * ldx]
    //   d(ldy)/dlx1 = 0, d(mx-lx1)/dlx1 = -1
    //   d(ldx)/dlx1 = -1, d(my-ly1)/dlx1 = 0
    //   = -ldy + ldx_no... = -ldy - (-(my-ly1)*(-1))
    //   = -ldy + (my - ly1)
    //   Hmm, let me be more careful.
    //   F0 = (mx-lx1)*(ly2-ly1) - (my-ly1)*(lx2-lx1)
    //   dF0/dlx1 = (-1)*(ly2-ly1) - (my-ly1)*(-1) = -ldy + (my-ly1) = my - ly2
    // dF0/dly1 = (mx-lx1)*(-1) - (-1)*(lx2-lx1) = -(mx-lx1) + ldx = lx2-mx
    // dF0/dlx2 = 0 - (my-ly1)*1 = -(my-ly1) = ly1-my
    // dF0/dly2 = (mx-lx1)*1 - 0 = mx-lx1

    // F1 partials:
    // F1 = (px2-px1)*(lx2-lx1) + (py2-py1)*(ly2-ly1)
    // dF1/dpx1 = -ldx, dF1/dpy1 = -ldy
    // dF1/dpx2 = ldx,  dF1/dpy2 = ldy
    // dF1/dlx1 = -(px2-px1) = -pdx, dF1/dly1 = -(py2-py1) = -pdy
    // dF1/dlx2 = pdx, dF1/dly2 = pdy

    EquationOutput {
        residuals: vec![f0, f1],
        jacobian: vec![
            // F0 partials
            (row, ipx1, 0.5 * ldy),
            (row, ipy1, -0.5 * ldx),
            (row, ipx2, 0.5 * ldy),
            (row, ipy2, -0.5 * ldx),
            (row, ilx1, my - ly2),
            (row, ily1, lx2 - mx),
            (row, ilx2, ly1 - my),
            (row, ily2, mx - lx1),
            // F1 partials
            (row + 1, ipx1, -ldx),
            (row + 1, ipy1, -ldy),
            (row + 1, ipx2, ldx),
            (row + 1, ipy2, ldy),
            (row + 1, ilx1, -pdx),
            (row + 1, ily1, -pdy),
            (row + 1, ilx2, pdx),
            (row + 1, ily2, pdy),
        ],
    }
}

/// Midpoint(P, L): F0 = px - (x1+x2)/2, F1 = py - (y1+y2)/2
pub(crate) fn eval_midpoint(
    sketch: &Sketch,
    p: PointHandle,
    l: LineHandle,
    row: usize,
) -> EquationOutput {
    let (px, py, ipx, ipy) = pt(sketch, p);
    let (x1, y1, ix1, iy1, x2, y2, ix2, iy2) = ln(sketch, l);
    EquationOutput {
        residuals: vec![px - (x1 + x2) / 2.0, py - (y1 + y2) / 2.0],
        jacobian: vec![
            (row, ipx, 1.0),
            (row, ix1, -0.5),
            (row, ix2, -0.5),
            (row + 1, ipy, 1.0),
            (row + 1, iy1, -0.5),
            (row + 1, iy2, -0.5),
        ],
    }
}

/// Concentric(A1, A2): F0 = cx1-cx2, F1 = cy1-cy2
pub(crate) fn eval_concentric(
    sketch: &Sketch,
    a1: ArcHandle,
    a2: ArcHandle,
    row: usize,
) -> EquationOutput {
    let (cx1, cy1, _, icx1, icy1, _) = ar(sketch, a1);
    let (cx2, cy2, _, icx2, icy2, _) = ar(sketch, a2);
    EquationOutput {
        residuals: vec![cx1 - cx2, cy1 - cy2],
        jacobian: vec![
            (row, icx1, 1.0),
            (row, icx2, -1.0),
            (row + 1, icy1, 1.0),
            (row + 1, icy2, -1.0),
        ],
    }
}

/// Coradial(A1, A2): F0 = cx1-cx2, F1 = cy1-cy2, F2 = r1-r2
pub(crate) fn eval_coradial(
    sketch: &Sketch,
    a1: ArcHandle,
    a2: ArcHandle,
    row: usize,
) -> EquationOutput {
    let (cx1, cy1, r1, icx1, icy1, ir1) = ar(sketch, a1);
    let (cx2, cy2, r2, icx2, icy2, ir2) = ar(sketch, a2);
    EquationOutput {
        residuals: vec![cx1 - cx2, cy1 - cy2, r1 - r2],
        jacobian: vec![
            (row, icx1, 1.0),
            (row, icx2, -1.0),
            (row + 1, icy1, 1.0),
            (row + 1, icy2, -1.0),
            (row + 2, ir1, 1.0),
            (row + 2, ir2, -1.0),
        ],
    }
}

/// OffsetLine(L1, L2, d): squared distance from midpoint of L1 to line L2
/// equals d^2.
///
/// F = N^2 / D - d^2, where
///   N = (mx-x21)(y22-y21) - (my-y21)(x22-x21),   mx=(x11+x12)/2, my=(y11+y12)/2
///   D = (x22-x21)^2 + (y22-y21)^2
pub(crate) fn eval_offset_line(
    sketch: &Sketch,
    l1: LineHandle,
    l2: LineHandle,
    d: f64,
    row: usize,
) -> EquationOutput {
    let (x11, y11, ix11, iy11, x12, y12, ix12, iy12) = ln(sketch, l1);
    let (x21, y21, ix21, iy21, x22, y22, ix22, iy22) = ln(sketch, l2);
    let mx = (x11 + x12) / 2.0;
    let my = (y11 + y12) / 2.0;
    let ldx = x22 - x21;
    let ldy = y22 - y21;
    let n = (mx - x21) * ldy - (my - y21) * ldx;
    let denom = ldx * ldx + ldy * ldy;
    let f = n * n / denom - d * d;

    let d2 = denom * denom;
    let two_n = 2.0 * n;
    let df = |dn: f64, dd: f64| -> f64 { (two_n * dn * denom - n * n * dd) / d2 };

    // Partial derivatives of N w.r.t. midpoint components:
    //   dN/dmx = ldy, dN/dmy = -ldx
    // Midpoint partials: dmx/dx11 = dmx/dx12 = 0.5, dmy/dy11 = dmy/dy12 = 0.5
    // So dN/dx11 = 0.5*ldy, dN/dx12 = 0.5*ldy
    //    dN/dy11 = -0.5*ldx, dN/dy12 = -0.5*ldx
    // Partials of N w.r.t. L2 endpoints (same as PointLineDistance with
    // point = midpoint):
    //   dN/dx21 = my - y22, dN/dy21 = x22 - mx
    //   dN/dx22 = -(my - y21), dN/dy22 = mx - x21
    // Partials of D (same as PointLineDistance):
    //   dD/dx21 = -2*ldx, dD/dy21 = -2*ldy
    //   dD/dx22 = 2*ldx, dD/dy22 = 2*ldy
    EquationOutput {
        residuals: vec![f],
        jacobian: vec![
            (row, ix11, df(0.5 * ldy, 0.0)),
            (row, iy11, df(-0.5 * ldx, 0.0)),
            (row, ix12, df(0.5 * ldy, 0.0)),
            (row, iy12, df(-0.5 * ldx, 0.0)),
            (row, ix21, df(my - y22, -2.0 * ldx)),
            (row, iy21, df(x22 - mx, -2.0 * ldy)),
            (row, ix22, df(-(my - y21), 2.0 * ldx)),
            (row, iy22, df(mx - x21, 2.0 * ldy)),
        ],
    }
}

/// OffsetArc(A1, A2, d): F = (r1-r2)^2 - d^2
pub(crate) fn eval_offset_arc(
    sketch: &Sketch,
    a1: ArcHandle,
    a2: ArcHandle,
    d: f64,
    row: usize,
) -> EquationOutput {
    let (_, _, r1, _, _, ir1) = ar(sketch, a1);
    let (_, _, r2, _, _, ir2) = ar(sketch, a2);
    let diff = r1 - r2;
    EquationOutput {
        residuals: vec![diff * diff - d * d],
        jacobian: vec![(row, ir1, 2.0 * diff), (row, ir2, -2.0 * diff)],
    }
}

/// OrdinateX(P, value): F = px - value
pub(crate) fn eval_ordinate_x(
    sketch: &Sketch,
    p: PointHandle,
    value: f64,
    row: usize,
) -> EquationOutput {
    let (px, _, ipx, _) = pt(sketch, p);
    EquationOutput {
        residuals: vec![px - value],
        jacobian: vec![(row, ipx, 1.0)],
    }
}

/// OrdinateY(P, value): F = py - value
pub(crate) fn eval_ordinate_y(
    sketch: &Sketch,
    p: PointHandle,
    value: f64,
    row: usize,
) -> EquationOutput {
    let (_, py, _, ipy) = pt(sketch, p);
    EquationOutput {
        residuals: vec![py - value],
        jacobian: vec![(row, ipy, 1.0)],
    }
}

/// Collinear(L1, L2): one endpoint of L2 lies on line through L1.
/// F = (x21-x11)(y12-y11) - (y21-y11)(x12-x11)
pub(crate) fn eval_collinear(
    sketch: &Sketch,
    l1: LineHandle,
    l2: LineHandle,
    row: usize,
) -> EquationOutput {
    let (x11, y11, ix11, iy11, x12, y12, ix12, iy12) = ln(sketch, l1);
    let (x21, y21, ix21, iy21, _x22, _y22, _ix22, _iy22) = ln(sketch, l2);
    let dx = x12 - x11;
    let dy = y12 - y11;
    // F = (x21 - x11) * dy - (y21 - y11) * dx
    // This is PointOnLine with point = L2.start on line L1.
    let f = (x21 - x11) * dy - (y21 - y11) * dx;

    // dF/dx21 = dy (= y12-y11)
    // dF/dy21 = -dx (= -(x12-x11))
    // dF/dx11 = -dy + (y21-y11) ... same derivation as PointOnLine
    //   = y21 - y12
    // dF/dy11 = dx - (x21-x11) = x12 - x21
    // dF/dx12 = -(y21-y11)
    // dF/dy12 = x21-x11
    EquationOutput {
        residuals: vec![f],
        jacobian: vec![
            (row, ix21, dy),
            (row, iy21, -dx),
            (row, ix11, y21 - y12),
            (row, iy11, x12 - x21),
            (row, ix12, -(y21 - y11)),
            (row, iy12, x21 - x11),
        ],
    }
}
