//! Spline arm of the sketch-to-BRep curve conversion (split from
//! `curve_convert.rs` for the 500-line ratchet).
//!
//! Open splines convert as-is. A periodic (closed) sketch spline is
//! spelled as the equivalent UNCLAMPED open NURBS the kernel can
//! evaluate, and each piece the edge graph split it into carries its
//! own parameter sub-range (inverted from the piece's endpoints) on
//! both the 3-D curve and the pcurve.

use crate::error::SketchError;
use crate::loop_tracer::DirectedEdge;
use crate::plane::SketchPlane;
use popkern_constraint::Sketch;
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::pcurve::Pcurve;

/// The spline element's NURBS in 2-D (`z = 0`).
///
/// A periodic (closed) sketch spline stores `n` control points and a
/// one-period knot vector, evaluated with wrapped control points over
/// an extended knot vector. The kernel's evaluator is not periodic, so
/// the same curve is spelled as the equivalent UNCLAMPED open NURBS:
/// the extended knots (`degree` knots of the previous period
/// prepended, `degree` of the next appended) over `n + degree` control
/// points wrapped modulo `n` — the constraint crate's own evaluation
/// rule, so the two agree to round-off. Its parameter range is the
/// original period.
fn spline_nurbs_2d(
    sketch: &Sketch,
    handle: popkern_constraint::SplineHandle,
) -> Result<NurbsCurve, SketchError> {
    let invalid = || SketchError::InvalidTopology("spline handle is invalid".into());
    let cp_handles = sketch.spline_control_points(handle).ok_or_else(invalid)?;
    let degree = sketch.spline_degree(handle).ok_or_else(invalid)?;
    let knots = sketch.spline_knots(handle).ok_or_else(invalid)?.to_vec();
    let weights = sketch.spline_weights(handle).ok_or_else(invalid)?.to_vec();
    let periodic = sketch.spline_periodic(handle).ok_or_else(invalid)?;
    let control_points = cp_handles
        .iter()
        .map(|&ph| {
            sketch
                .point_position(ph)
                .map(|(x, y)| Point3::new(x, y, 0.0))
                .ok_or_else(|| {
                    SketchError::InvalidTopology("spline control point is invalid".into())
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let n = control_points.len();
    if !periodic || n == 0 || knots.len() <= n || degree == 0 || n < degree {
        return Ok(NurbsCurve {
            degree,
            knots,
            control_points,
            weights,
        });
    }
    let period = knots[n] - knots[0];
    let mut ext = Vec::with_capacity(n + 2 * degree + 1);
    ext.extend(knots[(n - degree)..n].iter().map(|k| k - period));
    ext.extend_from_slice(&knots[..=n]);
    ext.extend(knots[1..=degree].iter().map(|k| k + period));
    let wrapped_cps = (0..n + degree).map(|j| control_points[j % n]).collect();
    let wrapped_ws = (0..n + degree)
        .map(|j| weights.get(j % n).copied().unwrap_or(1.0))
        .collect();
    Ok(NurbsCurve {
        degree,
        knots: ext,
        control_points: wrapped_cps,
        weights: wrapped_ws,
    })
}

/// Parameter of the point on `curve` nearest the 2-D point `p`: dense
/// sampling over the clamped range, then ternary refinement.
fn spline_param_of(curve: &NurbsCurve, p: (f64, f64)) -> f64 {
    use popkern_core::geometry::curve::CurveEval;
    let (t0, t1) = curve.parameter_range();
    let dist = |t: f64| {
        let q = curve.point_at(t);
        (q.x - p.0).powi(2) + (q.y - p.1).powi(2)
    };
    const N: usize = 256;
    let mut best = (dist(t0), t0);
    for i in 1..=N {
        let t = t0 + (t1 - t0) * i as f64 / N as f64;
        let d = dist(t);
        if d < best.0 {
            best = (d, t);
        }
    }
    let h = (t1 - t0) / N as f64;
    let (mut lo, mut hi) = ((best.1 - h).max(t0), (best.1 + h).min(t1));
    for _ in 0..60 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if dist(m1) < dist(m2) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    (lo + hi) / 2.0
}

/// The spline parameter sub-range a directed edge covers, in traversal
/// direction. A periodic spline is split into pieces that all
/// reference the whole element (see `edge_graph`), so each piece's
/// range is inverted from its endpoints; the piece that returns to
/// the curve start ends at the range's end, not its start.
fn spline_edge_range(curve: &NurbsCurve, directed: &DirectedEdge) -> (f64, f64) {
    use popkern_core::geometry::curve::CurveEval;
    let (t0, t1) = curve.parameter_range();
    let eps = 1e-9 * (t1 - t0).abs().max(1.0);
    // The element's own direction runs from its lower parameter.
    let (from, to) = if directed.reversed {
        (directed.end, directed.start)
    } else {
        (directed.start, directed.end)
    };
    let mut ta = spline_param_of(curve, from);
    let mut tb = spline_param_of(curve, to);
    if tb <= ta + eps {
        if (tb - t0).abs() <= eps {
            tb = t1;
        } else if (ta - t1).abs() <= eps {
            ta = t0;
        }
    }
    if directed.reversed {
        (tb, ta)
    } else {
        (ta, tb)
    }
}

/// Convert a spline directed edge to a 3D CurveOnEdge: the whole
/// spline lifted to the plane, with the edge's parameter sub-range.
pub(super) fn convert_spline_3d(
    directed: &DirectedEdge,
    sketch: &Sketch,
    plane: &SketchPlane,
    handle: popkern_constraint::SplineHandle,
) -> Result<CurveOnEdge, SketchError> {
    let curve_2d = spline_nurbs_2d(sketch, handle)?;
    let (t_start, t_end) = spline_edge_range(&curve_2d, directed);
    let mut curve = curve_2d;
    for cp in &mut curve.control_points {
        *cp = plane.to_3d(cp.x, cp.y);
    }
    Ok(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(curve),
        t_start,
        t_end,
    })
}

/// Convert a spline directed edge to a 2D pcurve: the same spline in
/// the plane chart, on the same parameter sub-range as the 3-D curve.
pub(super) fn convert_spline_pcurve(
    directed: &DirectedEdge,
    sketch: &Sketch,
    handle: popkern_constraint::SplineHandle,
) -> Result<Pcurve, SketchError> {
    let curve = spline_nurbs_2d(sketch, handle)?;
    let (t_start, t_end) = spline_edge_range(&curve, directed);
    let mut pcurve = Pcurve::from_nurbs(curve);
    pcurve.t_start = t_start;
    pcurve.t_end = t_end;
    Ok(pcurve)
}
