//! Recognise a curve that traces an exact circle or circular arc and
//! spell it as the analytic `Arc` it is.
//!
//! Boolean intersection curves reach the result body as NURBS even
//! when the solver that produced them was analytic (plane × cylinder
//! gives a rational circle; PK-417: a bore's rims could not be
//! selected as `circular()` and every downstream consumer lost the
//! analytic path). The recogniser samples the curve over its edge
//! range, fits a circle in the Newell plane (Kåsa least squares) and
//! accepts only when every sample lies on that circle AND in its
//! plane within `gate` — the caller's tolerance, so a marching
//! approximation of a circle is promoted only within the accuracy the
//! body already records for it.

use crate::geometry::analytic_curves::circle_basis;
use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use crate::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

/// A circle fitted to a curve's samples.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CircleFit {
    /// Circle centre.
    pub center: Point3,
    /// Unit normal of the circle's plane (right-hand sense of the
    /// curve's traversal).
    pub normal: Vec3,
    /// Radius.
    pub radius: f64,
    /// Largest sample deviation from the fitted circle (radial or
    /// out-of-plane), for diagnostics.
    pub max_deviation: f64,
    /// Whether the curve's ends coincide (a full circle).
    pub closed: bool,
}

const SAMPLES: usize = 33;

/// Recognise a circle or circular arc in `curve` over `range`, within
/// `gate` (an absolute distance). `None` when the samples are not
/// coplanar-and-equidistant to that accuracy, when the arc is too
/// short to determine a circle, or when the curve is degenerate.
pub fn recognize_circular(
    curve: &dyn CurveEval,
    range: (f64, f64),
    gate: f64,
) -> Option<CircleFit> {
    let (t0, t1) = range;
    if !(t1 - t0).is_finite() || (t1 - t0).abs() <= 0.0 {
        return None;
    }
    let pts: Vec<Point3> = (0..SAMPLES)
        .map(|i| curve.point_at(t0 + (t1 - t0) * i as f64 / (SAMPLES - 1) as f64))
        .collect();
    let closed = pts[0].distance_to(&pts[SAMPLES - 1]) <= gate;
    // Samples that are collinear within the gate carry no circle: the
    // least-squares fit then answers an enormous radius whose sagitta
    // over the sample span is below the gate (a bore's straight seam
    // line read as a circle of radius 8078 with a zero-degree sweep at
    // the 1e-3 tolerance, and the degenerate arc tore the mesh).
    if !closed && chord_sagitta(&pts) <= gate {
        return None;
    }
    let m = pts.len() as f64;
    let c0 = Point3::new(
        pts.iter().map(|p| p.x).sum::<f64>() / m,
        pts.iter().map(|p| p.y).sum::<f64>() / m,
        pts.iter().map(|p| p.z).sum::<f64>() / m,
    );
    // Newell normal of the sample polygon (closed by its chord for an
    // open arc — a circular segment has the arc's plane normal).
    let mut n = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..pts.len() {
        let a = pts[i];
        let b = pts[(i + 1) % pts.len()];
        n.x += (a.y - b.y) * (a.z + b.z);
        n.y += (a.z - b.z) * (a.x + b.x);
        n.z += (a.x - b.x) * (a.y + b.y);
    }
    let normal = n.normalized()?;
    let (e1, e2) = orthonormal_basis(normal);
    let xy: Vec<(f64, f64)> = pts
        .iter()
        .map(|p| {
            let d = *p - c0;
            (d.dot(&e1), d.dot(&e2))
        })
        .collect();
    let (mut sx, mut sy, mut sxx, mut syy, mut sxy, mut sxz, mut syz, mut sz) =
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for &(x, y) in &xy {
        let z = x * x + y * y;
        sx += x;
        sy += y;
        sxx += x * x;
        syy += y * y;
        sxy += x * y;
        sxz += x * z;
        syz += y * z;
        sz += z;
    }
    let mm = [[sxx, sxy, sx], [sxy, syy, sy], [sx, sy, m]];
    let rhs = [sxz, syz, sz];
    let det3 = |a: &[[f64; 3]; 3]| -> f64 {
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    };
    let d0 = det3(&mm);
    if d0.abs() < 1e-30 {
        return None;
    }
    let solve = |col: usize| -> f64 {
        let mut t = mm;
        for (row, trow) in t.iter_mut().enumerate() {
            trow[col] = rhs[row];
        }
        det3(&t) / d0
    };
    let (a, b, c2) = (solve(0), solve(1), solve(2));
    let (ccx, ccy) = (0.5 * a, 0.5 * b);
    let radius = (c2 + ccx * ccx + ccy * ccy).max(0.0).sqrt();
    if !(radius.is_finite() && radius > gate) {
        return None;
    }
    let center = c0 + e1 * ccx + e2 * ccy;
    let mut max_deviation: f64 = 0.0;
    for p in &pts {
        let d = *p - center;
        max_deviation = max_deviation
            .max((d.length() - radius).abs())
            .max(d.dot(&normal).abs());
    }
    if max_deviation > gate {
        return None;
    }
    // The arc must subtend enough to determine its circle: a long chord,
    // or a bulge well above the gate. The chord rule alone, in units of
    // the gate, refused every quarter rim at a coarse tolerance (PK-779:
    // a 0.717 rim's 0.55 chord against 100 × 1e-2, its sagitta 0.21
    // resolving the circle twenty times over), so the boolean served
    // fitted rims only at 1e-2.
    let chord = pts[0].distance_to(&pts[SAMPLES / 2]);
    if !closed && chord < gate * 100.0 && chord_sagitta(&pts) < gate * 10.0 {
        return None;
    }
    // Orient the normal to the traversal's right-hand sense.
    let d_a = pts[0] - center;
    let d_b = pts[SAMPLES / 4] - center;
    let normal = if d_a.cross(&d_b).dot(&normal) >= 0.0 {
        normal
    } else {
        -normal
    };
    Some(CircleFit {
        center,
        normal,
        radius,
        max_deviation,
        closed,
    })
}

/// Largest distance of any sample from the line through the first
/// sample and the sample farthest from it — zero for a straight run.
fn chord_sagitta(pts: &[Point3]) -> f64 {
    let origin = pts[0];
    let far = pts
        .iter()
        .max_by(|a, b| {
            a.distance_to(&origin)
                .partial_cmp(&b.distance_to(&origin))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .copied()
        .unwrap_or(origin);
    let Some(dir) = (far - origin).normalized() else {
        return 0.0;
    };
    pts.iter()
        .map(|p| {
            let d = *p - origin;
            (d - dir * d.dot(&dir)).length()
        })
        .fold(0.0, f64::max)
}

/// The `Arc` edge curve for a recognised circle: angle-parameterised
/// in the arc evaluator's own frame (`circle_basis(normal)`), starting
/// at the curve's start point and running in its traversal sense; a
/// closed curve becomes the full turn.
pub fn arc_on_edge(curve: &dyn CurveEval, range: (f64, f64), fit: &CircleFit) -> CurveOnEdge {
    arc_on_edge_oriented(curve, range, fit, fit.normal)
}

/// As [`arc_on_edge`] with an explicit normal: `normal` may be the
/// fit's (the arc runs the curve's way) or its negation (the arc runs
/// the other way round — the caller then re-anchors the edge's
/// forward half-edge). The start angle is measured in
/// `circle_basis(normal)`, the evaluator's frame.
pub fn arc_on_edge_oriented(
    curve: &dyn CurveEval,
    range: (f64, f64),
    fit: &CircleFit,
    normal: Vec3,
) -> CurveOnEdge {
    let same_way = normal.dot(&fit.normal) >= 0.0;
    let (start_point, end_point) = if same_way {
        (curve.point_at(range.0), curve.point_at(range.1))
    } else {
        (curve.point_at(range.1), curve.point_at(range.0))
    };
    arc_between(fit, normal, start_point, end_point)
}

/// The `Arc` on a recognised circle from `start_point` to `end_point`
/// (each measured by its angle about the circle, so a point a little
/// off it still reads its own place), running positively about
/// `normal` in `circle_basis(normal)`, the evaluator's frame. A closed
/// fit becomes the full turn from `start_point`. Callers that know
/// the edge's vertices pass them (PK-944): a sampled curve's ends may
/// sit 1e-9 along the circle from the vertex they meet, and an arc
/// that starts there gives every iso-circle pcurve on it an end the
/// neighbouring edge's pcurve does not share.
pub fn arc_between(
    fit: &CircleFit,
    normal: Vec3,
    start_point: Point3,
    end_point: Point3,
) -> CurveOnEdge {
    let (e1, e2) = circle_basis(normal);
    let angle = |p: Point3| {
        let d = p - fit.center;
        d.dot(&e2).patan2(d.dot(&e1))
    };
    // A start within angle noise of the frame origin is spelled AS
    // the origin: a full turn then reads `[0, 2π]`, the principal
    // form every consumer's NURBS conversion keeps exact (a start of
    // -1e-16 would clamp a quarter turn away).
    let mut start_angle = angle(start_point);
    if start_angle.abs() < 1.0e-9 {
        start_angle = 0.0;
    }
    let end_angle = if fit.closed {
        start_angle + std::f64::consts::TAU
    } else {
        let mut end = angle(end_point);
        while end <= start_angle {
            end += std::f64::consts::TAU;
        }
        end
    };
    CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center: fit.center,
            normal,
            radius: fit.radius,
            start_angle,
            end_angle,
        },
        t_start: start_angle,
        t_end: end_angle,
    }
}

#[cfg(test)]
#[path = "circle_recognition_tests.rs"]
mod tests;
