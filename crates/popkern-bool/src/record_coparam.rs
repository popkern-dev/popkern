//! Record uv curves co-parameterized with their 3D curve (PK-1031).
//!
//! An intersection record carries one curve in 3D and one in each
//! face's chart, over one parameter domain. The intersectors build them
//! independently: a circle comes back as a rational quadratic in 3D and
//! a straight line in a cylinder's chart, the same point set at
//! different speeds. Every consumer that evaluates a uv curve at a 3D
//! parameter (clip window ends, closed-record rotation, the vertex
//! table's record-end chart points, the pcurves of the record's edges)
//! then lands off its 3D point along the curve: measured 3e-3 to 2e-2
//! on cone and sphere walls, where the edge-face hit at the same point
//! is exact.
//!
//! A record whose uv image strays from its 3D curve beyond ten linear
//! tolerances, at matched fractions of the two domains, gets its uv
//! curve rebuilt: each sample of the 3D curve is located on the
//! existing uv curve (its chart branch kept, the search marching
//! forward from the start), and the located uv points are interpolated
//! over the 3D parameters, piecewise between the 3D curve's C0 knots.
//! The sample count doubles until the rebuilt curve's image follows
//! the 3D curve at every span midpoint within a tenth of the tolerance
//! or twice the uv curve's own shape error (the farthest a located
//! sample lies from the 3D point, a floor no sampling moves), stopping
//! at the first doubling that fails to halve the gap; a curve that
//! does not converge keeps its original uv (no worse than before).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, fitting::fit_curve};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::intersect::IntersectionRecord;

/// Fraction of the linear tolerance the rebuilt uv image may stray.
const GATE_FRACTION: f64 = 0.1;
/// A uv image within this many linear tolerances of its 3D curve at
/// matched fractions is left alone: a record's uv curves are fits, and
/// a helicoid's sit about one tolerance off their 3D curve by shape
/// alone (measured on an M6 thread cut: 9.4e-7..1.1e-6; re-fitting the
/// 98 just past one tolerance cost 0.6 s of a 7.4 s cut and moved
/// nothing), which no re-parameterization moves. The speed mismatch
/// this pass repairs reads 8.8e-4..2e-2.
const TRIGGER_TOLERANCES: f64 = 10.0;
/// Matched fractions sampled by the misfit check.
const CHECK_SAMPLES: usize = 16;
/// Initial spans of the rebuild, doubled while it misses the gate.
const INITIAL_SPANS: usize = 16;
/// Doublings: 16 spans reach 512.
const MAX_DOUBLINGS: usize = 5;
/// Samples of the uv curve the forward search starts from.
const LOCATE_SAMPLES: usize = 512;
/// Golden-section steps refining a located parameter.
const REFINE_STEPS: usize = 48;
/// uv distance within which two pieces join at a breakpoint.
const JOIN_TOLERANCE: f64 = 1e-9;

/// Co-parameterize every record's uv curves with its 3D curve.
pub(crate) fn coparameterize_records(
    records: &mut [IntersectionRecord],
    a: &BRep,
    b: &BRep,
    tolerance: &Tolerance,
) {
    let gate = GATE_FRACTION * tolerance.linear();
    let trigger = TRIGGER_TOLERANCES * tolerance.linear();
    let trace = popkern_core::trace::trace_enabled("uv-scramble");
    for (i, r) in records.iter_mut().enumerate() {
        let c3 = r.curve.curve_3d.clone();
        for (label, brep, face, uv) in [
            ("uv1", a, r.face_a, &mut r.curve.curve_uv1),
            ("uv2", b, r.face_b, &mut r.curve.curve_uv2),
        ] {
            let Some(surface) = brep
                .faces
                .get(face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry)
            else {
                continue;
            };
            let before = misfit(&c3, uv, surface);
            if before <= trigger {
                continue;
            }
            let rebuilt = rebuild(&c3, uv, surface, gate);
            if trace {
                eprintln!(
                    "  [uv-scramble] rec {i} {label}: image {before:.3e} off its 3D curve -> {}",
                    if rebuilt.is_some() {
                        "co-parameterized"
                    } else {
                        "kept (rebuild did not converge)"
                    }
                );
            }
            if let Some(curve) = rebuilt {
                *uv = curve;
            }
        }
    }
}

/// Largest distance between the uv image and the 3D curve at matched
/// fractions of their domains.
fn misfit(c3: &NurbsCurve, uv: &NurbsCurve, surface: &SurfaceGeometry) -> f64 {
    let (t0, t1) = c3.parameter_range();
    let (u0, u1) = uv.parameter_range();
    (0..=CHECK_SAMPLES)
        .map(|k| {
            let f = k as f64 / CHECK_SAMPLES as f64;
            let q = uv.point_at(u0 + (u1 - u0) * f);
            surface
                .point_at(q.x, q.y)
                .distance_to(&c3.point_at(t0 + (t1 - t0) * f))
        })
        .fold(0.0, f64::max)
}

/// The uv curve's parameter nearest `p` (by surface image), searched
/// forward from `from` on the dense table `table` of (τ, image).
fn locate(
    uv: &NurbsCurve,
    surface: &SurfaceGeometry,
    table: &[(f64, Point3)],
    from: usize,
    p: &Point3,
) -> (usize, f64, f64) {
    // The nearest table sample at or after `from`, stopping once the
    // distance has grown well past the best (the march is local: a
    // closed curve's far end returns to the start point).
    let mut best = (f64::INFINITY, from);
    for (j, (_, q)) in table.iter().enumerate().skip(from) {
        let d = q.distance_to(p);
        if d < best.0 {
            best = (d, j);
        } else if d > 4.0 * best.0 + 1e-12 && j > best.1 + 8 {
            break;
        }
    }
    let j = best.1;
    let (mut lo, mut hi) = (
        table[j.saturating_sub(1)].0,
        table[(j + 1).min(table.len() - 1)].0,
    );
    let dist = |tau: f64| {
        let q = uv.point_at(tau);
        surface.point_at(q.x, q.y).distance_to(p)
    };
    for _ in 0..REFINE_STEPS {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if dist(m1) < dist(m2) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let tau = 0.5 * (lo + hi);
    (j, tau, dist(tau))
}

/// 3D parameters of the C0 knots of `c3` (multiplicity at least its
/// degree) strictly inside its domain.
fn breakpoints(c3: &NurbsCurve) -> Vec<f64> {
    let (t0, t1) = c3.parameter_range();
    let mut out = Vec::new();
    let mut k = 0;
    while k < c3.knots.len() {
        let knot = c3.knots[k];
        let m = c3.knots[k..].iter().take_while(|&&x| x == knot).count();
        if m >= c3.degree && knot > t0 + 1e-12 && knot < t1 - 1e-12 {
            out.push(knot);
        }
        k += m;
    }
    out
}

/// `uv` rebuilt over `c3`'s parameter, or `None` when the rebuild does
/// not reach `gate`.
fn rebuild(
    c3: &NurbsCurve,
    uv: &NurbsCurve,
    surface: &SurfaceGeometry,
    gate: f64,
) -> Option<NurbsCurve> {
    let (t0, t1) = c3.parameter_range();
    let (u0, u1) = uv.parameter_range();
    let table: Vec<(f64, Point3)> = (0..=LOCATE_SAMPLES)
        .map(|j| {
            let tau = u0 + (u1 - u0) * j as f64 / LOCATE_SAMPLES as f64;
            let q = uv.point_at(tau);
            (tau, surface.point_at(q.x, q.y))
        })
        .collect();
    let bounds: Vec<f64> = std::iter::once(t0)
        .chain(breakpoints(c3))
        .chain(std::iter::once(t1))
        .collect();
    let mut spans = INITIAL_SPANS;
    let mut previous = f64::INFINITY;
    for _ in 0..=MAX_DOUBLINGS {
        // Spans per segment by its share of the domain, three at least.
        let per: Vec<usize> = bounds
            .windows(2)
            .map(|w| ((spans as f64 * (w[1] - w[0]) / (t1 - t0)).ceil() as usize).max(3))
            .collect();
        let mut params = vec![t0];
        for (w, &n) in bounds.windows(2).zip(&per) {
            params.extend((1..=n).map(|i| w[0] + (w[1] - w[0]) * i as f64 / n as f64));
        }
        // The ends are the uv curve's own: its branch at both ends.
        let mut from = 0;
        let mut floor = 0.0_f64;
        let mut points = Vec::with_capacity(params.len());
        for (k, &t) in params.iter().enumerate() {
            let tau = if k == 0 {
                u0
            } else if k + 1 == params.len() {
                u1
            } else {
                let (j, tau, d) = locate(uv, surface, &table, from, &c3.point_at(t));
                from = j.saturating_sub(1);
                floor = floor.max(d);
                tau
            };
            points.push(uv.point_at(tau));
        }
        let curve = piecewise(&points, &params, &per)?;
        let worst = params
            .windows(2)
            .map(|w| {
                let t = 0.5 * (w[0] + w[1]);
                let q = curve.point_at(t);
                surface.point_at(q.x, q.y).distance_to(&c3.point_at(t))
            })
            .fold(0.0, f64::max);
        if worst <= gate.max(2.0 * floor) {
            return Some(curve);
        }
        if worst > 0.5 * previous {
            return None;
        }
        previous = worst;
        spans *= 2;
    }
    None
}

/// Interpolate `points` at `params` piecewise (`per[i]` spans in
/// segment `i`), joined C0.
fn piecewise(points: &[Point3], params: &[f64], per: &[usize]) -> Option<NurbsCurve> {
    let mut start = 0;
    let mut joined: Option<NurbsCurve> = None;
    for &n in per {
        let end = start + n;
        let (pts, ts) = (points.get(start..=end)?, params.get(start..=end)?);
        let piece = fit_curve(pts, ts, 3.min(pts.len() - 1), pts.len()).ok()?;
        joined = Some(match joined {
            None => piece,
            Some(left) => left.concatenate_clamped(&piece, JOIN_TOLERANCE).ok()?,
        });
        start = end;
    }
    joined
}

#[cfg(test)]
#[path = "record_coparam_tests.rs"]
mod tests;
