//! Tracing a tangent curve from one of its contacts (tangency design
//! §4.8, W3c): two surfaces tangent along a curve meet in no transversal
//! intersection the marcher can follow (an analytic arm is
//! tangency-silent there), and contact discovery finds points of it,
//! each classified `TangentCurve` with the curve's direction. The trace
//! steps along that direction and settles each step back onto the curve:
//! across it the signed gap has a critical point (its valley floor), the
//! root of its derivative along the first surface's transverse
//! direction, read by central differences (a flat minimum located by its
//! derivative's root, never by its value). It stops where the tangency
//! ends (the gap past the linear reach), where either window ends, or
//! back at its start (a closed curve).

use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::uv_projection::project_point_to_uv_near;
use crate::geometry::{Point3, Vec3};

use super::contact_classify::{ContactScale, SurfaceContactKind};
use super::contact_discovery::SurfaceContact;

/// A parameter window `((u0, u1), (v0, v1))`.
type Window = ((f64, f64), (f64, f64));

/// Settled parameters on both surfaces and the gap there.
type Settled = ((f64, f64), (f64, f64), f64);

/// A tangent curve as traced: points on the first surface, in order.
#[derive(Debug, Clone, PartialEq)]
pub struct TangentTrace {
    /// The points, on the first surface.
    pub points: Vec<Point3>,
    /// Whether the curve closes (its last point steps back onto its
    /// first; the first is not repeated).
    pub closed: bool,
}

/// Steps along a curve before the trace gives up.
const MAX_STEPS: usize = 4096;

/// The uv direction on `s` at `uv` whose image is the 3-D direction `d`
/// (least squares onto the tangent plane), or `None` where the chart is
/// singular.
fn uv_of(s: &SurfaceGeometry, uv: (f64, f64), d: Vec3) -> Option<(f64, f64)> {
    let (su, sv) = s.first_partials(uv.0, uv.1);
    let (e, f, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
    let det = e * g - f * f;
    if !det.is_finite() || det.abs() <= 1e-30 {
        return None;
    }
    let (p, q) = (d.dot(&su), d.dot(&sv));
    Some(((g * p - f * q) / det, (e * q - f * p) / det))
}

/// `uv` brought into window `w` through `s`'s periods, or `None` when
/// it lies outside.
fn into_window(s: &SurfaceGeometry, w: Window, uv: (f64, f64)) -> Option<(f64, f64)> {
    let (pu, pv) = crate::geometry::chart::ChartDescriptor::evaluation_periods_of(s);
    let fit = |x: f64, (lo, hi): (f64, f64), period: Option<f64>| -> Option<f64> {
        let slack = 1e-12 * (hi - lo).abs().max(1.0);
        match period {
            Some(t) if t > 0.0 => {
                let y = crate::geometry::chart::ops::wrap_into(x, lo - slack, t);
                (y <= hi + slack).then_some(y)
            }
            _ => (x >= lo - slack && x <= hi + slack).then_some(x),
        }
    };
    Some((fit(uv.0, w.0, pu)?, fit(uv.1, w.1, pv)?))
}

/// The signed gap from `s1` at `uv1` to `s2` (along `s2`'s normal at the
/// foot found near `hint`), with the foot's parameters.
fn gap(
    s1: &SurfaceGeometry,
    uv1: (f64, f64),
    s2: &SurfaceGeometry,
    hint: (f64, f64),
) -> Option<(f64, (f64, f64))> {
    let p = s1.point_at(uv1.0, uv1.1);
    let uv2 = project_point_to_uv_near(s2, &p, hint.0, hint.1)?;
    let q = s2.point_at(uv2.0, uv2.1);
    let n = s2.normal_at(uv2.0, uv2.1).normalized()?;
    Some(((p - q).dot(&n), uv2))
}

/// `uv1` moved across the curve (along the uv image of `across`, a unit
/// 3-D direction) to the root of the gap's derivative there, at most
/// `reach` in space: the valley floor. Returns the settled parameters on
/// both surfaces and the gap.
fn settle(
    (s1, s2): (&SurfaceGeometry, &SurfaceGeometry),
    (mut uv1, mut uv2): ((f64, f64), (f64, f64)),
    across: Vec3,
    (h, reach): (f64, f64),
) -> Option<Settled> {
    let mut g0 = 0.0;
    for _ in 0..12 {
        let tau = uv_of(s1, uv1, across)?;
        let at = |s: f64| (uv1.0 + tau.0 * s, uv1.1 + tau.1 * s);
        let (g, foot) = gap(s1, uv1, s2, uv2)?;
        let (gp, _) = gap(s1, at(h), s2, foot)?;
        let (gm, _) = gap(s1, at(-h), s2, foot)?;
        g0 = g;
        uv2 = foot;
        let d1 = (gp - gm) / (2.0 * h);
        let d2 = (gp - 2.0 * g + gm) / (h * h);
        if !d1.is_finite() || !d2.is_finite() || d2.abs() <= 1e-12 {
            break;
        }
        let s = (-d1 / d2).clamp(-reach, reach);
        uv1 = at(s);
        if s.abs() <= 1e-3 * h {
            break;
        }
    }
    Some((uv1, uv2, g0))
}

/// Trace the tangent curve through `seed` (a contact of `s1` over `w1`
/// with `s2` over `w2`, classified `TangentCurve`): both ways from it,
/// in steps of a sixty-fourth of `scale.extent`, until the tangency ends
/// (the gap past the linear reach), a window ends, or the curve closes.
/// `None` when the seed is not on a tangent curve.
pub fn trace_tangent_curve(
    (s1, w1): (&SurfaceGeometry, Window),
    (s2, w2): (&SurfaceGeometry, Window),
    seed: &SurfaceContact,
    scale: ContactScale,
) -> Option<TangentTrace> {
    let SurfaceContactKind::TangentCurve { direction, .. } = &seed.classification.kind else {
        return None;
    };
    let reach = scale.linear * (1.0 + 1e-6);
    let step = (scale.extent / 64.0).max(10.0 * scale.linear);
    let diff = (1e-3 * step).max(10.0 * scale.linear);
    let mut halves: Vec<Vec<Point3>> = Vec::new();
    for sign in [1.0, -1.0] {
        let mut out = Vec::new();
        let (mut uv1, mut uv2) = (seed.uv1, seed.uv2);
        let mut t = *direction * sign;
        let mut p = seed.point;
        for k in 0..MAX_STEPS {
            let d = uv_of(s1, uv1, t)?;
            let Some(next) = into_window(s1, w1, (uv1.0 + d.0 * step, uv1.1 + d.1 * step)) else {
                break;
            };
            let n = s1.normal_at(next.0, next.1).normalized()?;
            let Some(across) = n.cross(&t).normalized() else {
                break;
            };
            let Some((u1, u2, g)) = settle((s1, s2), (next, uv2), across, (diff, step)) else {
                break;
            };
            let (Some(u1), Some(u2)) = (into_window(s1, w1, u1), into_window(s2, w2, u2)) else {
                break;
            };
            if g.abs() > reach {
                break;
            }
            let q = s1.point_at(u1.0, u1.1);
            if k > 2 && q.distance_to(&seed.point) < 0.75 * step {
                halves.push(out);
                return Some(TangentTrace {
                    points: closed_points(seed.point, halves.remove(0)),
                    closed: true,
                });
            }
            let Some(next_t) = (q - p).normalized() else {
                break;
            };
            (uv1, uv2, t, p) = (u1, u2, next_t, q);
            out.push(q);
        }
        halves.push(out);
    }
    let (fwd, back) = (halves.remove(0), halves.remove(0));
    let mut points: Vec<Point3> = back.into_iter().rev().collect();
    points.push(seed.point);
    points.extend(fwd);
    (points.len() >= 2).then_some(TangentTrace {
        points,
        closed: false,
    })
}

/// A closed trace's points: the seed, then the forward walk.
fn closed_points(seed: Point3, forward: Vec<Point3>) -> Vec<Point3> {
    std::iter::once(seed).chain(forward).collect()
}

#[cfg(test)]
#[path = "contact_trace_tests.rs"]
mod tests;
