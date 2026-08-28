// Items are pub(crate) for use by the BRep builder module.
// Suppress dead_code until those consumers are added.
#![allow(dead_code)]

//! 2D-to-3D curve conversion for sketch elements.
//!
//! Converts [`DirectedEdge`] entries (from the loop tracer) into 3D
//! [`CurveOnEdge`] bindings and 2D pcurves (`NurbsCurve` with z=0)
//! suitable for attaching to BRep edges and half-edges.

use crate::error::SketchError;
use crate::loop_tracer::DirectedEdge;
use crate::plane::SketchPlane;
use popkern_constraint::{ElementHandle, Sketch};
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;

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
        ElementHandle::Spline(h) => convert_spline_3d(directed, sketch, plane, *h),
        ElementHandle::Point(_) => Err(SketchError::InvalidTopology(
            "point elements cannot be converted to curves".into(),
        )),
    }
}

/// Convert a 2D sketch element to a 2D pcurve (NurbsCurve with z=0).
///
/// The resulting NURBS curve lives in the sketch plane's UV parameter
/// space. Control point z-coordinates are always 0.
pub(crate) fn convert_to_pcurve(
    directed: &DirectedEdge,
    sketch: &Sketch,
) -> Result<NurbsCurve, SketchError> {
    match &directed.element {
        ElementHandle::Line(_) => convert_line_pcurve(directed),
        ElementHandle::Arc(h) => convert_arc_pcurve(directed, sketch, *h),
        ElementHandle::Ellipse(_) | ElementHandle::Parabola(_) | ElementHandle::Hyperbola(_) => {
            convert_conic_pcurve(directed, sketch, &directed.element)
        }
        ElementHandle::Spline(h) => convert_spline_pcurve(directed, sketch, *h),
        ElementHandle::Point(_) => Err(SketchError::InvalidTopology(
            "point elements cannot be converted to pcurves".into(),
        )),
    }
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
/// Derives the actual angular range from the directed edge's start/end
/// positions, which may differ from the original arc parameters when a
/// full circle has been split into two half-arcs.
fn convert_arc_3d(
    directed: &DirectedEdge,
    sketch: &Sketch,
    plane: &SketchPlane,
    handle: popkern_constraint::ArcHandle,
) -> Result<CurveOnEdge, SketchError> {
    let (cx, cy, r, _original_start, _original_end) = sketch
        .arc_params(handle)
        .ok_or_else(|| SketchError::InvalidTopology("arc handle is invalid".into()))?;

    // Derive angles from the directed edge's 2D positions rather than
    // the sketch arc's original parameters. This correctly handles arcs
    // that have been split (e.g., full circles split into semicircles).
    let start_angle = (directed.start.1 - cy).atan2(directed.start.0 - cx);
    let end_angle = (directed.end.1 - cy).atan2(directed.end.0 - cx);

    let center_3d = plane.to_3d(cx, cy);

    Ok(CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center: center_3d,
            normal: plane.normal,
            radius: r,
            start_angle,
            end_angle,
        },
        t_start: start_angle,
        t_end: end_angle,
    })
}

/// Convert an arc directed edge to a 2D pcurve (sampled NURBS).
///
/// Derives the actual angular range from the directed edge's start/end
/// positions, which may differ from the original arc parameters when a
/// full circle has been split into two half-arcs.
fn convert_arc_pcurve(
    directed: &DirectedEdge,
    sketch: &Sketch,
    handle: popkern_constraint::ArcHandle,
) -> Result<NurbsCurve, SketchError> {
    let (cx, cy, r, _original_start, _original_end) = sketch
        .arc_params(handle)
        .ok_or_else(|| SketchError::InvalidTopology("arc handle is invalid".into()))?;

    // Derive angles from the directed edge's 2D positions.
    let start_angle = (directed.start.1 - cy).atan2(directed.start.0 - cx);
    let end_angle = (directed.end.1 - cy).atan2(directed.end.0 - cx);

    let n = CONIC_SAMPLE_COUNT;
    let points: Vec<Point3> = (0..=n)
        .map(|i| {
            let t = start_angle + (end_angle - start_angle) * (i as f64 / n as f64);
            let x = cx + r * t.cos();
            let y = cy + r * t.sin();
            Point3::new(x, y, 0.0)
        })
        .collect();

    let params: Vec<f64> = (0..=n).map(|i| i as f64 / n as f64).collect();
    let num_cp = (n / 2).max(4);

    popkern_core::geometry::fitting::fit_curve(&points, &params, 3, num_cp)
        .map_err(|e| SketchError::InvalidTopology(format!("arc pcurve fitting failed: {e}")))
}

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
            let cos_r = rotation.cos();
            let sin_r = rotation.sin();

            // Use angle as parameter for sampling.
            let start_angle =
                angle_from_ellipse_point(directed.start, cx, cy, semi_a, semi_b, cos_r, sin_r);
            let end_angle =
                angle_from_ellipse_point(directed.end, cx, cy, semi_a, semi_b, cos_r, sin_r);

            Ok((0..=n)
                .map(|i| {
                    let t = start_angle + (end_angle - start_angle) * (i as f64 / n as f64);
                    let lx = semi_a * t.cos();
                    let ly = semi_b * t.sin();
                    (cx + lx * cos_r - ly * sin_r, cy + lx * sin_r + ly * cos_r)
                })
                .collect())
        }
        ElementHandle::Parabola(h) => {
            let (vx, vy, focal_dist, rotation) = sketch
                .parabola_params(*h)
                .ok_or_else(|| SketchError::InvalidTopology("parabola handle is invalid".into()))?;
            let cos_r = rotation.cos();
            let sin_r = rotation.sin();

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
            let cos_r = rotation.cos();
            let sin_r = rotation.sin();

            // Hyperbola: local (semi_a*cosh(t), semi_b*sinh(t)), t in [-1, 1]
            let t_start: f64 = -1.0;
            let t_end: f64 = 1.0;

            Ok((0..=n)
                .map(|i| {
                    let t = t_start + (t_end - t_start) * (i as f64 / n as f64);
                    let lx = semi_a * t.cosh();
                    let ly = semi_b * t.sinh();
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
    (ly / semi_b).atan2(lx / semi_a)
}

// ---------------------------------------------------------------------------
// Spline
// ---------------------------------------------------------------------------

/// Convert a spline directed edge to a 3D CurveOnEdge.
///
/// Reads the spline's control points, knots, weights, and degree from
/// the sketch, transforms control points to 3D, and builds a NURBS
/// curve directly.
fn convert_spline_3d(
    directed: &DirectedEdge,
    sketch: &Sketch,
    plane: &SketchPlane,
    handle: popkern_constraint::SplineHandle,
) -> Result<CurveOnEdge, SketchError> {
    let cp_handles = sketch
        .spline_control_points(handle)
        .ok_or_else(|| SketchError::InvalidTopology("spline handle is invalid".into()))?;

    let degree = sketch
        .spline_degree(handle)
        .ok_or_else(|| SketchError::InvalidTopology("spline handle is invalid".into()))?;

    let knots = sketch
        .spline_knots(handle)
        .ok_or_else(|| SketchError::InvalidTopology("spline handle is invalid".into()))?
        .to_vec();

    let weights = sketch
        .spline_weights(handle)
        .ok_or_else(|| SketchError::InvalidTopology("spline handle is invalid".into()))?
        .to_vec();

    let control_points_3d: Result<Vec<Point3>, SketchError> = cp_handles
        .iter()
        .map(|&ph| {
            let (x, y) = sketch.point_position(ph).ok_or_else(|| {
                SketchError::InvalidTopology("spline control point is invalid".into())
            })?;
            Ok(plane.to_3d(x, y))
        })
        .collect();

    let control_points_3d = control_points_3d?;

    let (t_start, t_end) = if directed.reversed {
        (1.0, 0.0)
    } else {
        (0.0, 1.0)
    };

    Ok(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(NurbsCurve {
            degree,
            knots,
            control_points: control_points_3d,
            weights,
        }),
        t_start,
        t_end,
    })
}

/// Convert a spline directed edge to a 2D pcurve.
fn convert_spline_pcurve(
    _directed: &DirectedEdge,
    sketch: &Sketch,
    handle: popkern_constraint::SplineHandle,
) -> Result<NurbsCurve, SketchError> {
    let cp_handles = sketch
        .spline_control_points(handle)
        .ok_or_else(|| SketchError::InvalidTopology("spline handle is invalid".into()))?;

    let degree = sketch
        .spline_degree(handle)
        .ok_or_else(|| SketchError::InvalidTopology("spline handle is invalid".into()))?;

    let knots = sketch
        .spline_knots(handle)
        .ok_or_else(|| SketchError::InvalidTopology("spline handle is invalid".into()))?
        .to_vec();

    let weights = sketch
        .spline_weights(handle)
        .ok_or_else(|| SketchError::InvalidTopology("spline handle is invalid".into()))?
        .to_vec();

    let control_points: Result<Vec<Point3>, SketchError> = cp_handles
        .iter()
        .map(|&ph| {
            let (x, y) = sketch.point_position(ph).ok_or_else(|| {
                SketchError::InvalidTopology("spline control point is invalid".into())
            })?;
            Ok(Point3::new(x, y, 0.0))
        })
        .collect();

    Ok(NurbsCurve {
        degree,
        knots,
        control_points: control_points?,
        weights,
    })
}

#[cfg(test)]
#[path = "curve_convert_tests.rs"]
mod tests;
