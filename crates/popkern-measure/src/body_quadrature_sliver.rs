//! Boundary slivers for the quadrature volume (PK-872).
//!
//! A face mesh's boundary triangles end on straight (u, v) chords, but
//! the face's domain ends on its trim pcurves: where a trim is curved in
//! the chart, each chord cuts a sliver off the true domain (or adds one),
//! so the integral depends on how densely the chart samples the boundary
//! (measured: the rolling-ball dodecahedron's identical sphere caps read
//! -1.13e-5 of Steiner's volume under one chart axis and +2.7e-6 under
//! another). Here every boundary chord also integrates the region
//! between it and its trim: a fan from the chord's start over pcurve
//! samples to its end, in the triangle's own chart branch and winding,
//! so an outward bulge adds and an inward one subtracts.
//!
//! A chord whose ends do not both lie on one trim (a mesh slit, a face
//! with no pcurves, a sample far from the chord) keeps the chord: the
//! correction declines to the straight boundary it replaces.

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, FaceHandle};

use super::oriented_term;

/// Pcurve samples a trim is located against before refinement.
const LOCATE_SAMPLES: usize = 64;
/// Newton steps refining a located trim parameter.
const LOCATE_STEPS: usize = 6;
/// Fan panels between a chord's ends. The sliver's sag is quadratic in
/// the chord, so the panelled sliver is within 1/PANELS² of the true one.
const PANELS: usize = 16;

/// One trim of a face: its pcurve, located by samples.
pub(super) struct Trim<'a> {
    pcurve: &'a Pcurve,
    lo: f64,
    hi: f64,
    closed: bool,
    samples: Vec<(f64, (f64, f64))>,
}

/// The face's trims (its loops' pcurves on this face), or empty when
/// any half-edge lacks one.
pub(super) fn face_trims(brep: &BRep, face: FaceHandle) -> Vec<Trim<'_>> {
    let Some(data) = brep.faces.get(face) else {
        return Vec::new();
    };
    let mut trims = Vec::new();
    for &lp in std::iter::once(&data.outer_loop).chain(&data.inner_loops) {
        let Some(start) = brep.loops.get(lp).map(|l| l.half_edge) else {
            return Vec::new();
        };
        let mut he = start;
        loop {
            let Some(h) = brep.half_edges.get(he) else {
                return Vec::new();
            };
            let Some(pcurve) = brep.face_use_pcurve(face, he).or(h.pcurve.as_ref()) else {
                return Vec::new();
            };
            let (t0, t1) = pcurve.parameter_range();
            let (lo, hi) = (t0.min(t1), t0.max(t1));
            let samples: Vec<_> = (0..=LOCATE_SAMPLES)
                .map(|k| {
                    let t = lo + (hi - lo) * k as f64 / LOCATE_SAMPLES as f64;
                    let p = pcurve.point_at(t);
                    (t, (p.x, p.y))
                })
                .collect();
            trims.push(Trim {
                pcurve,
                lo,
                hi,
                closed: brep
                    .half_edges
                    .get(h.twin)
                    .is_some_and(|t| t.origin == h.origin),
                samples,
            });
            he = h.next;
            if he == start || trims.len() > 100_000 {
                break;
            }
        }
    }
    trims
}

/// A chord end lies on a trim within this share of the chord's length
/// in space.
const ON_TRIM: f64 = 1e-3;

/// `x` shifted by whole periods to lie nearest `to`.
fn near(x: f64, to: f64, period: Option<f64>) -> f64 {
    match period {
        Some(p) if p > 0.0 => x - p * ((x - to) / p).round(),
        _ => x,
    }
}

fn aligned(q: (f64, f64), to: (f64, f64), periods: (Option<f64>, Option<f64>)) -> (f64, f64) {
    (near(q.0, to.0, periods.0), near(q.1, to.1, periods.1))
}

fn uv_distance(a: (f64, f64), b: (f64, f64)) -> f64 {
    let (du, dv) = (a.0 - b.0, a.1 - b.1);
    (du * du + dv * dv).sqrt()
}

/// The trim parameter nearest `q`, refined by Newton on `|pc(t) − q|²`,
/// with its distance.
fn locate(trim: &Trim<'_>, q: (f64, f64), periods: (Option<f64>, Option<f64>)) -> (f64, f64) {
    let mut best = (trim.lo, f64::INFINITY);
    for &(t, p) in &trim.samples {
        let d = uv_distance(aligned(p, q, periods), q);
        if d < best.1 {
            best = (t, d);
        }
    }
    let mut t = best.0;
    for _ in 0..LOCATE_STEPS {
        let p = trim.pcurve.point_at(t);
        let p = aligned((p.x, p.y), q, periods);
        let d = trim.pcurve.tangent_at(t);
        let dd = d.x * d.x + d.y * d.y;
        if dd <= 0.0 {
            break;
        }
        let step = ((p.0 - q.0) * d.x + (p.1 - q.1) * d.y) / dd;
        t = (t - step).clamp(trim.lo, trim.hi);
    }
    let p = trim.pcurve.point_at(t);
    (t, uv_distance(aligned((p.x, p.y), q, periods), q))
}

/// `⅓ ∬ r · (∂u × ∂v)` over the region between the chord `a → b` (in a
/// triangle's chart branch and winding) and the trim through both ends,
/// or `None` when no single trim carries the chord. A chord end is on a
/// trim when the trim's nearest point lies within a thousandth of the
/// chord of it both in (u, v) and in space (and within `gate`): in
/// (u, v) alone a chart's lengths are not comparable (a 1 mm ring's two rims are
/// 1e-3 apart in v, its chords 1.57 long in u), and in space alone a
/// pole vertex's arbitrary u passes (a cone's apex at u = π matched its
/// seam's end at u = 2π). The nearest such trim carries it.
pub(super) fn sliver_term(
    surface: &SurfaceGeometry,
    trims: &[Trim<'_>],
    a: (f64, f64),
    b: (f64, f64),
    gate: f64,
) -> Option<f64> {
    let chord = uv_distance(a, b);
    if chord.is_nan() || chord <= 0.0 {
        return None;
    }
    let periods = surface.parameter_periodicity();
    let (pa, pb) = (surface.point_at(a.0, a.1), surface.point_at(b.0, b.1));
    // Within a thousandth of the chord, and never past the caller's
    // gate: the trim across a narrow face (a fillet's other rim at a
    // coarse tolerance) is a strip width away, far more than that.
    let gate = gate.min(ON_TRIM * pa.distance_to(&pb));
    let off = |trim: &Trim<'_>, t: f64, q: (f64, f64), p: popkern_core::geometry::Point3| {
        let c = trim.pcurve.point_at(t);
        let c = aligned((c.x, c.y), q, periods);
        surface.point_at(c.0, c.1).distance_to(&p)
    };
    let mut best: Option<(&Trim<'_>, f64, f64, f64)> = None;
    for trim in trims {
        let (ta, uva) = locate(trim, a, periods);
        let da = off(trim, ta, a, pa);
        if da > gate || uva > ON_TRIM * chord {
            continue;
        }
        let (tb, uvb) = locate(trim, b, periods);
        let worst = da.max(off(trim, tb, b, pb));
        if worst <= gate && uvb <= ON_TRIM * chord && ta != tb && best.is_none_or(|x| worst < x.3) {
            best = Some((trim, ta, tb, worst));
        }
    }
    let (trim, ta, tb, _) = best?;
    // A closed trim's short way round may cross its parameter seam.
    let span = trim.hi - trim.lo;
    let mut delta = tb - ta;
    if trim.closed && delta.abs() > 0.5 * span {
        delta -= span * delta.signum();
    }
    let mut points = Vec::with_capacity(PANELS + 1);
    points.push(a);
    for k in 1..PANELS {
        let f = k as f64 / PANELS as f64;
        let mut t = ta + delta * f;
        if t > trim.hi {
            t -= span;
        } else if t < trim.lo {
            t += span;
        }
        let p = trim.pcurve.point_at(t);
        let along = (a.0 + (b.0 - a.0) * f, a.1 + (b.1 - a.1) * f);
        let q = aligned((p.x, p.y), along, periods);
        // A sample off the chord's neighbourhood is another branch of
        // the trim (or a trim the chord does not follow).
        if uv_distance(q, along) > chord {
            return None;
        }
        points.push(q);
    }
    points.push(b);
    let sum: f64 = points
        .windows(2)
        .skip(1)
        .map(|w| oriented_term(surface, [a, w[0], w[1]]))
        .sum();
    sum.is_finite().then_some(sum)
}
