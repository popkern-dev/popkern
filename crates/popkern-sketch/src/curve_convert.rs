// Items are pub(crate) for use by the BRep builder module.
// Suppress dead_code until those consumers are added.
#![allow(dead_code)]

//! 2D-to-3D curve conversion for sketch elements.
//!
//! Converts [`DirectedEdge`] entries (from the loop tracer) into 3D
//! [`CurveOnEdge`] bindings and 2D pcurves (`NurbsCurve` with z=0)
//! suitable for attaching to BRep edges and half-edges.

#[path = "curve_convert_spline.rs"]
mod spline;

use crate::error::SketchError;
use crate::loop_tracer::DirectedEdge;
use crate::plane::SketchPlane;
use popkern_constraint::{ElementHandle, Sketch};
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_math::PortableMath;

/// Number of sample points used when approximating conics as NURBS.
const CONIC_SAMPLE_COUNT: usize = 20;

/// Convert a 2D sketch element to a 3D [`CurveOnEdge`].
///
/// Dispatches on the element type stored in `directed`. Lines and arcs
/// produce analytic curve geometry; conics and splines are promoted to
/// NURBS.
pub(crate) fn convert_to_3d(
    directed: &DirectedEdge,
    sketch: &Sketch,
    plane: &SketchPlane,
) -> Result<CurveOnEdge, SketchError> {
    match &directed.element {
        ElementHandle::Line(_) => convert_line_3d(directed, plane),
        ElementHandle::Arc(h) => convert_arc_3d(directed, sketch, plane, *h),
        ElementHandle::Ellipse(_) => convert_conic_3d(directed, sketch, plane, &directed.element),
        ElementHandle::Parabola(_) => convert_conic_3d(directed, sketch, plane, &directed.element),
        ElementHandle::Hyperbola(_) => convert_conic_3d(directed, sketch, plane, &directed.element),
        ElementHandle::Spline(h) => spline::convert_spline_3d(directed, sketch, plane, *h),
        ElementHandle::Point(_) => Err(SketchError::InvalidTopology(
            "point elements cannot be converted to curves".into(),
        )),
    }
}

/// Convert a 2D sketch element to a pcurve in the sketch plane's UV
/// chart, parameterised in the DIRECTED edge's traversal direction
/// (`t_start` at `directed.start`). Arcs are exact rational quadratics
/// whose parameter is the 3-D arc's angle (see [`convert_to_3d`]);
/// every other element's pcurve runs on its own `[0, 1]`.
pub(crate) fn convert_to_pcurve(
    directed: &DirectedEdge,
    sketch: &Sketch,
    plane: &SketchPlane,
) -> Result<Pcurve, SketchError> {
    match &directed.element {
        ElementHandle::Line(_) => convert_line_pcurve(directed).map(Pcurve::from_nurbs),
        ElementHandle::Arc(h) => convert_arc_pcurve(directed, sketch, plane, *h),
        ElementHandle::Ellipse(_) | ElementHandle::Parabola(_) | ElementHandle::Hyperbola(_) => {
            convert_conic_pcurve(directed, sketch, &directed.element).map(Pcurve::from_nurbs)
        }
        ElementHandle::Spline(h) => spline::convert_spline_pcurve(directed, sketch, *h),
        ElementHandle::Point(_) => Err(SketchError::InvalidTopology(
            "point elements cannot be converted to pcurves".into(),
        )),
    }
}

/// A directed arc edge as sketch-plane geometry: centre, radius, the
/// angle of `directed.start` (sketch-plane `atan2`), and the SIGNED
/// angular sweep to `directed.end` in the traversal direction.
///
/// Sketch arcs run counter-clockwise from their start angle to their
/// end angle (clockwise when the end angle is the smaller); a graph
/// edge follows the element unless `reversed`. The sweep is derived
/// from that sense and the endpoint angles — never from the raw
/// element angles (a split full circle shares them between halves)
/// and never by subtracting two `atan2` values (the ±π branch cut
/// flipped every arc crossing it onto the complementary run; a
/// constraint-sketch circle extruded as two identical half-arcs,
/// validate red with a `√2` loop gap).
struct DirectedArc {
    cx: f64,
    cy: f64,
    r: f64,
    start_angle: f64,
    sweep: f64,
}

fn directed_arc(
    directed: &DirectedEdge,
    sketch: &Sketch,
    handle: popkern_constraint::ArcHandle,
) -> Result<DirectedArc, SketchError> {
    use std::f64::consts::TAU;
    let (cx, cy, r, elem_start, elem_end) = sketch
        .arc_params(handle)
        .ok_or_else(|| SketchError::InvalidTopology("arc handle is invalid".into()))?;
    let ccw = (elem_end >= elem_start) != directed.reversed;
    let start_angle = (directed.start.1 - cy).patan2(directed.start.0 - cx);
    let end_angle = (directed.end.1 - cy).patan2(directed.end.0 - cx);
    let forward = (end_angle - start_angle).rem_euclid(TAU);
    let magnitude = if ccw { forward } else { TAU - forward };
    // Coincident endpoints on an arc element mean a full turn (a
    // split circle's halves never coincide).
    let magnitude = if !(1e-12..=TAU - 1e-12).contains(&magnitude) {
        TAU
    } else {
        magnitude
    };
    Ok(DirectedArc {
        cx,
        cy,
        r,
        start_angle,
        sweep: if ccw { magnitude } else { -magnitude },
    })
}

// ---------------------------------------------------------------------------
// Line
// ---------------------------------------------------------------------------

/// Convert a line directed edge to a 3D CurveOnEdge.
fn convert_line_3d(
    directed: &DirectedEdge,
    plane: &SketchPlane,
) -> Result<CurveOnEdge, SketchError> {
    let start_3d = plane.to_3d(directed.start.0, directed.start.1);
    let end_3d = plane.to_3d(directed.end.0, directed.end.1);
    let direction = end_3d - start_3d;

    let (t_start, t_end) = if directed.reversed {
        (1.0, 0.0)
    } else {
        (0.0, 1.0)
    };

    Ok(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: start_3d,
            direction,
        },
        t_start,
        t_end,
    })
}

/// Convert a line directed edge to a 2D pcurve.
fn convert_line_pcurve(directed: &DirectedEdge) -> Result<NurbsCurve, SketchError> {
    let (sx, sy) = directed.start;
    let (ex, ey) = directed.end;
    Ok(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(sx, sy, 0.0), Point3::new(ex, ey, 0.0)],
        weights: vec![1.0, 1.0],
    })
}

// ---------------------------------------------------------------------------
// Arc
// ---------------------------------------------------------------------------

/// Convert an arc directed edge to a 3D CurveOnEdge.
///
/// The kernel's `Arc` evaluator measures angles in the frame
/// [`circle_basis`](popkern_core::geometry::analytic_curves::circle_basis)
/// derives from the normal — NOT the sketch plane's `(u, v)` — so the
/// start angle is re-measured in that frame from the 3-D start point
/// and the signed sweep (which both right-handed frames share) is
/// added to it. `t_start`/`t_end` run in the traversal direction.
fn convert_arc_3d(
    directed: &DirectedEdge,
    sketch: &Sketch,
    plane: &SketchPlane,
    handle: popkern_constraint::ArcHandle,
) -> Result<CurveOnEdge, SketchError> {
    let arc = directed_arc(directed, sketch, handle)?;
    let center_3d = plane.to_3d(arc.cx, arc.cy);
    let start_angle = arc_frame_angle(plane, &arc, directed.start);
    let end_angle = start_angle + arc.sweep;
    Ok(CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center: center_3d,
            normal: plane.normal,
            radius: arc.r,
            start_angle,
            end_angle,
        },
        t_start: start_angle,
        t_end: end_angle,
    })
}

/// The angle of sketch point `p` about the arc centre, measured in the
/// frame the kernel's `Arc` evaluator uses for `plane.normal`.
fn arc_frame_angle(plane: &SketchPlane, arc: &DirectedArc, p: (f64, f64)) -> f64 {
    use popkern_core::geometry::analytic_curves::circle_basis;
    let (u, v) = circle_basis(plane.normal);
    let d = plane.to_3d(p.0, p.1) - plane.to_3d(arc.cx, arc.cy);
    d.dot(&v).patan2(d.dot(&u))
}

/// Convert an arc directed edge to an EXACT rational-quadratic pcurve
/// whose parameter is the 3-D arc's angle (the same frame offset as
/// [`convert_arc_3d`]), so pcurve and edge curve agree at every span
/// joint. `t_start` is at `directed.start`; a clockwise traversal has
/// `t_start > t_end` over a counter-clockwise-spelled NURBS.
fn convert_arc_pcurve(
    directed: &DirectedEdge,
    sketch: &Sketch,
    plane: &SketchPlane,
    handle: popkern_constraint::ArcHandle,
) -> Result<Pcurve, SketchError> {
    let arc = directed_arc(directed, sketch, handle)?;
    let (a_lo, a_hi) = if arc.sweep >= 0.0 {
        (arc.start_angle, arc.start_angle + arc.sweep)
    } else {
        (arc.start_angle + arc.sweep, arc.start_angle)
    };
    let mut spans = crate::exact_arc::arc_spans(arc.cx, arc.cy, arc.r, a_lo, a_hi, 0.0);
    // Knot labels carry the 3-D frame's angle: both frames are
    // right-handed about the plane normal, so the offset between the
    // sketch-plane angle and the evaluator-frame angle is a constant.
    let offset = arc_frame_angle(plane, &arc, directed.start) - arc.start_angle;
    for s in &mut spans {
        s.t0 += offset;
        s.t1 += offset;
    }
    let nurbs = crate::exact_arc::assemble_quadratic(&spans);
    let mut pcurve = Pcurve::from_nurbs(nurbs);
    if arc.sweep < 0.0 {
        std::mem::swap(&mut pcurve.t_start, &mut pcurve.t_end);
    }
    Ok(pcurve)
}

// ---------------------------------------------------------------------------
// Splines: see `curve_convert_spline`.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Conics (Ellipse, Parabola, Hyperbola) — sample and fit
// ---------------------------------------------------------------------------

/// Convert a conic element to 3D via sampling and NURBS fitting.
fn convert_conic_3d(
    directed: &DirectedEdge,
    sketch: &Sketch,
    plane: &SketchPlane,
    element: &ElementHandle,
) -> Result<CurveOnEdge, SketchError> {
    let samples_2d = sample_element_2d(directed, sketch, element)?;

    let points_3d: Vec<Point3> = samples_2d.iter().map(|&(x, y)| plane.to_3d(x, y)).collect();

    let n = points_3d.len();
    let params: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
    let num_cp = (n / 2).max(4);

    let nurbs = popkern_core::geometry::fitting::fit_curve(&points_3d, &params, 3, num_cp)
        .map_err(|e| SketchError::InvalidTopology(format!("conic curve fitting failed: {e}")))?;

    Ok(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(nurbs),
        t_start: 0.0,
        t_end: 1.0,
    })
}

/// Convert a conic element to a 2D pcurve via sampling.
fn convert_conic_pcurve(
    directed: &DirectedEdge,
    sketch: &Sketch,
    element: &ElementHandle,
) -> Result<NurbsCurve, SketchError> {
    let samples_2d = sample_element_2d(directed, sketch, element)?;

    let points: Vec<Point3> = samples_2d
        .iter()
        .map(|&(x, y)| Point3::new(x, y, 0.0))
        .collect();

    let n = points.len();
    let params: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
    let num_cp = (n / 2).max(4);

    popkern_core::geometry::fitting::fit_curve(&points, &params, 3, num_cp)
        .map_err(|e| SketchError::InvalidTopology(format!("conic pcurve fitting failed: {e}")))
}

/// Sample a conic element at evenly-spaced points in 2D.
///
/// Returns `CONIC_SAMPLE_COUNT + 1` points from the directed start to
/// end (interpolating parameter linearly between the endpoints).
fn sample_element_2d(
    directed: &DirectedEdge,
    sketch: &Sketch,
    element: &ElementHandle,
) -> Result<Vec<(f64, f64)>, SketchError> {
    let n = CONIC_SAMPLE_COUNT;

    match element {
        ElementHandle::Ellipse(h) => {
            let (cx, cy, semi_a, semi_b, rotation) = sketch
                .ellipse_params(*h)
                .ok_or_else(|| SketchError::InvalidTopology("ellipse handle is invalid".into()))?;
            let cos_r = rotation.pcos();
            let sin_r = rotation.psin();

            // Use angle as parameter for sampling.
            let start_angle =
                angle_from_ellipse_point(directed.start, cx, cy, semi_a, semi_b, cos_r, sin_r);
            let end_angle =
                angle_from_ellipse_point(directed.end, cx, cy, semi_a, semi_b, cos_r, sin_r);

            Ok((0..=n)
                .map(|i| {
                    let t = start_angle + (end_angle - start_angle) * (i as f64 / n as f64);
                    let lx = semi_a * t.pcos();
                    let ly = semi_b * t.psin();
                    (cx + lx * cos_r - ly * sin_r, cy + lx * sin_r + ly * cos_r)
                })
                .collect())
        }
        ElementHandle::Parabola(h) => {
            let (vx, vy, focal_dist, rotation) = sketch
                .parabola_params(*h)
                .ok_or_else(|| SketchError::InvalidTopology("parabola handle is invalid".into()))?;
            let cos_r = rotation.pcos();
            let sin_r = rotation.psin();

            // Parabola: local (f*t^2, 2*f*t), t in [-1, 1]
            let t_start = -1.0;
            let t_end = 1.0;

            Ok((0..=n)
                .map(|i| {
                    let t = t_start + (t_end - t_start) * (i as f64 / n as f64);
                    let lx = focal_dist * t * t;
                    let ly = 2.0 * focal_dist * t;
                    (vx + lx * cos_r - ly * sin_r, vy + lx * sin_r + ly * cos_r)
                })
                .collect())
        }
        ElementHandle::Hyperbola(h) => {
            let (cx, cy, semi_a, semi_b, rotation) =
                sketch.hyperbola_params(*h).ok_or_else(|| {
                    SketchError::InvalidTopology("hyperbola handle is invalid".into())
                })?;
            let cos_r = rotation.pcos();
            let sin_r = rotation.psin();

            // Hyperbola: local (semi_a*cosh(t), semi_b*sinh(t)), t in [-1, 1]
            let t_start: f64 = -1.0;
            let t_end: f64 = 1.0;

            Ok((0..=n)
                .map(|i| {
                    let t = t_start + (t_end - t_start) * (i as f64 / n as f64);
                    let lx = semi_a * t.pcosh();
                    let ly = semi_b * t.psinh();
                    (cx + lx * cos_r - ly * sin_r, cy + lx * sin_r + ly * cos_r)
                })
                .collect())
        }
        _ => Err(SketchError::InvalidTopology(
            "unsupported element type for conic sampling".into(),
        )),
    }
}

/// Compute the angle parameter for a point on an ellipse.
fn angle_from_ellipse_point(
    pt: (f64, f64),
    cx: f64,
    cy: f64,
    semi_a: f64,
    semi_b: f64,
    cos_r: f64,
    sin_r: f64,
) -> f64 {
    // Rotate point into ellipse-local frame.
    let dx = pt.0 - cx;
    let dy = pt.1 - cy;
    let lx = dx * cos_r + dy * sin_r;
    let ly = -dx * sin_r + dy * cos_r;
    // atan2(ly / semi_b, lx / semi_a)
    (ly / semi_b).patan2(lx / semi_a)
}

// ---------------------------------------------------------------------------
// Spline
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "curve_convert_tests.rs"]
mod tests;
