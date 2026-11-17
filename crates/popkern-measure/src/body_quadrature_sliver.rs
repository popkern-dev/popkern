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

use popkern_core::geometry::Point3;
use popkern_core::geometry::chart::ChartDescriptor;
use popkern_core::geometry::chart::ops::{aligned, singular_point_near};
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, FaceHandle};

use super::{Integrand, oriented_term};

/// Pcurve samples a trim is located against before refinement.
const LOCATE_SAMPLES: usize = 64;
/// Newton steps refining a located trim parameter.
const LOCATE_STEPS: usize = 6;
/// Fan panels between a chord's ends at half resolution: the sliver is
/// integrated at `2 * PANELS` and `PANELS` panels and Richardson-
/// extrapolated (the inscribed fan's error is quadratic in the panel).
const PANELS: usize = 16;

/// One trim of a face: its pcurve, located by samples.
pub(super) struct Trim<'a> {
    pcurve: &'a Pcurve,
    /// The edge's 3-D curve: the true trim (a fitted pcurve can sit
    /// 1e-5 off it, PK-826).
    edge: Option<&'a popkern_core::geometry::curve::CurveOnEdge>,
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
                edge: brep.edges.get(h.edge).and_then(|e| e.curve.as_ref()),
                lo,
                hi,
                closed: h
                    .mate()
                    .and_then(|m| brep.half_edges.get(m))
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
/// A located pcurve point lies within this share of the chord of a
/// chord end in (u, v) when the edge curve decides in space.
const UV_BRANCH: f64 = 0.05;

fn uv_distance(a: (f64, f64), b: (f64, f64)) -> f64 {
    let (du, dv) = (a.0 - b.0, a.1 - b.1);
    (du * du + dv * dv).sqrt()
}

/// The trim parameter nearest `q`, refined by Newton on `|pc(t) − q|²`,
/// with its distance.
fn locate(trim: &Trim<'_>, q: (f64, f64), descriptor: &ChartDescriptor) -> (f64, f64) {
    let mut best = (trim.lo, f64::INFINITY);
    for &(t, p) in &trim.samples {
        let d = uv_distance(aligned(descriptor, p, q), q);
        if d < best.1 {
            best = (t, d);
        }
    }
    let mut t = best.0;
    for _ in 0..LOCATE_STEPS {
        let p = trim.pcurve.point_at(t);
        let p = aligned(descriptor, (p.x, p.y), q);
        let d = trim.pcurve.tangent_at(t);
        let dd = d.x * d.x + d.y * d.y;
        if dd <= 0.0 {
            break;
        }
        let step = ((p.0 - q.0) * d.x + (p.1 - q.1) * d.y) / dd;
        // A closed trim's parameter is periodic: a chord end at its
        // parameter seam steps across it. Clamping pinned such an end
        // to the far extreme, the chord matched no trim, and its sliver
        // went missing (measured: a cube ∩ sphere under a rotation that
        // put two cap circles' parameter seams on mesh chords, 1.6e-6
        // of the volume; PK-933).
        t = if trim.closed && trim.hi > trim.lo {
            trim.lo + (t - step - trim.lo).rem_euclid(trim.hi - trim.lo)
        } else {
            (t - step).clamp(trim.lo, trim.hi)
        };
    }
    let p = trim.pcurve.point_at(t);
    (t, uv_distance(aligned(descriptor, (p.x, p.y), q), q))
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
    descriptor: &ChartDescriptor,
    trims: &[Trim<'_>],
    a: (f64, f64),
    b: (f64, f64),
    gate: f64,
    kind: Integrand,
) -> Option<f64> {
    let chord = uv_distance(a, b);
    if chord.is_nan() || chord <= 0.0 {
        return None;
    }
    let (pa, pb) = (surface.point_at(a.0, a.1), surface.point_at(b.0, b.1));
    // Within a thousandth of the chord, and never past the caller's
    // gate: the trim across a narrow face (a fillet's other rim at a
    // coarse tolerance) is a strip width away, far more than that.
    let caller_gate = gate;
    let gate = gate.min(ON_TRIM * pa.distance_to(&pb));
    // The located parameter itself must image at the chord end: the
    // edge's own closest point says the end is ON the trim, not that
    // the parameter found in (u, v) is that point. In a chart whose
    // axes differ in scale the uv band admits another branch of the
    // trim: a kilometre Viviani lobe's two branches mirror each other
    // in u, 0.083 rad apart (416 units in space) while its chords run
    // 4.8 along v, and every chord near the pinch fanned over the
    // mirror branch (the common read 2.3e-4 high, PK-1146).
    let branch = UV_BRANCH * pa.distance_to(&pb) + caller_gate;
    let at = |trim: &Trim<'_>, t: f64, q: (f64, f64)| {
        let c = trim.pcurve.point_at(t);
        let c = aligned(descriptor, (c.x, c.y), q);
        surface.point_at(c.0, c.1)
    };
    // Off the trim in space: from the EDGE curve when the trim has one
    // (exact; a fitted pcurve sits up to ~1e-4 off it, PK-826, and the
    // hub of the propeller declined every blade-root chord on the
    // pcurve's distance), else from the pcurve's image.
    // Also the foot's edge parameter, when the trim has an edge.
    let off = |trim: &Trim<'_>, t: f64, q: (f64, f64), p: Point3| -> (f64, Option<f64>) {
        if let Some(edge) = trim.edge {
            let (lo, hi) = (edge.t_start.min(edge.t_end), edge.t_start.max(edge.t_end));
            return match edge.geometry.closest_point(&p) {
                Some(te) => {
                    let te = te.clamp(lo, hi);
                    (edge.geometry.point_at(te).distance_to(&p), Some(te))
                }
                None => (f64::INFINITY, None),
            };
        }
        let c = trim.pcurve.point_at(t);
        let c = aligned(descriptor, (c.x, c.y), q);
        (surface.point_at(c.0, c.1).distance_to(&p), None)
    };
    // In (u, v) the pcurve only picks the branch: within a twentieth of
    // the chord with an edge to decide in space, a thousandth without.
    let uv_band = |trim: &Trim<'_>| {
        if trim.edge.is_some() {
            UV_BRANCH * chord
        } else {
            ON_TRIM * chord
        }
    };
    // A chord end at a chart pole stands under the chord's other end
    // (`face_quadrature`), not at the trim's own longitude there: it is
    // found on the trim in space, at the trim's end imaging at the pole,
    // and the fan runs along the pole row to the trim (a row of no
    // area). Located in (u, v) it lay a longitude step off the trim and
    // the sliver went missing: Viviani's cap with the ball's poles on
    // its rim read 3.4e-4 short (PK-1169).
    let (pole_a, pole_b) = (
        at_pole(surface, descriptor, a),
        at_pole(surface, descriptor, b),
    );
    let find = |trim: &Trim<'_>, q: (f64, f64), p, pole: bool| -> (f64, f64) {
        if !pole {
            return locate(trim, q, descriptor);
        }
        let t = [trim.lo, trim.hi]
            .into_iter()
            .min_by(|&x, &y| {
                at(trim, x, q)
                    .distance_to(&p)
                    .total_cmp(&at(trim, y, q).distance_to(&p))
            })
            .unwrap_or(trim.lo);
        (t, 0.0)
    };
    // The cheap (u, v) tests first: the edge's global foot search runs
    // only for trims the chord could lie on (the same verdicts).
    let mut best: Option<Carrier<'_>> = None;
    for trim in trims {
        let (ta, uva) = find(trim, a, pa, pole_a);
        if uva > uv_band(trim) || at(trim, ta, a).distance_to(&pa) > branch {
            continue;
        }
        let (da, fa) = off(trim, ta, a, pa);
        if da > gate {
            continue;
        }
        let (tb, uvb) = find(trim, b, pb, pole_b);
        if uvb > uv_band(trim) || at(trim, tb, b).distance_to(&pb) > branch || ta == tb {
            continue;
        }
        let (db, fb) = off(trim, tb, b, pb);
        let worst = da.max(db);
        if worst <= gate && best.is_none_or(|x| worst < x.3) {
            best = Some((trim, ta, tb, worst, [fa, fb]));
        }
    }
    let (trim, ta, tb, _, feet) = best?;
    // The trim's own chart point at a pole end, on the chord's branch.
    let pole_point = |t: f64, q: (f64, f64)| {
        let c = trim.pcurve.point_at(t);
        aligned(descriptor, (c.x, c.y), q)
    };
    let (via_a, via_b) = (
        pole_a.then(|| pole_point(ta, a)),
        pole_b.then(|| pole_point(tb, b)),
    );
    // A closed trim's short way round may cross its parameter seam.
    let span = trim.hi - trim.lo;
    let mut delta = tb - ta;
    if trim.closed && delta.abs() > 0.5 * span {
        delta -= span * delta.signum();
    }
    // The edge parameters of the chord's ends: each fan sample's foot on
    // the edge lies between them, so Newton finds it from there. A
    // global search per sample (five Newton starts over the whole edge)
    // held ~90% of a threaded part's volume, 150 s on an arm link, and
    // on a multi-turn edge could settle on another turn (PK-1223). Ends
    // either side of a closed edge's parameter seam keep the search.
    let end_feet = match (trim.edge, feet) {
        (Some(edge), [Some(fa), Some(fb)]) => {
            let span = (edge.t_end - edge.t_start).abs();
            (!trim.closed || (fb - fa).abs() <= 0.5 * span).then_some((fa, fb))
        }
        _ => None,
    };
    let panels = 2 * PANELS;
    let mut points = Vec::with_capacity(panels + 1);
    for k in 1..panels {
        let f = k as f64 / panels as f64;
        let mut t = ta + delta * f;
        if t > trim.hi {
            t -= span;
        } else if t < trim.lo {
            t += span;
        }
        let p = trim.pcurve.point_at(t);
        let along = (a.0 + (b.0 - a.0) * f, a.1 + (b.1 - a.1) * f);
        let seed = aligned(descriptor, (p.x, p.y), along);
        // The sliver ends on the EDGE, not the fitted pcurve: the seed's
        // point on the edge curve, projected back near the seed
        // (measured: the propeller's fitted blade pcurves biased the
        // volume by -9e-7 per blade face against OCC).
        let q = match trim.edge {
            Some(edge) => {
                let on_surface = surface.point_at(seed.0, seed.1);
                let (lo, hi) = (edge.t_start.min(edge.t_end), edge.t_start.max(edge.t_end));
                let near = end_feet.and_then(|(fa, fb)| {
                    local_foot(&edge.geometry, &on_surface, fa + (fb - fa) * f, lo, hi)
                });
                let te = match near {
                    Some(te) => te,
                    None => edge.geometry.closest_point(&on_surface)?.clamp(lo, hi),
                };
                let on_edge = edge.geometry.point_at(te);
                let uv = popkern_core::geometry::uv_projection::project_point_to_uv_near(
                    surface, &on_edge, seed.0, seed.1,
                )?;
                aligned(descriptor, uv, along)
            }
            None => seed,
        };
        // A sample off the chord's neighbourhood is another branch of
        // the trim (or a trim the chord does not follow).
        if uv_distance(q, along) > chord {
            return None;
        }
        points.push(q);
    }
    // The fan from `a` over the trim's samples to `b`, through the
    // trim's pole points where an end is at a pole.
    let fan = |samples: &[(f64, f64)]| -> f64 {
        let pts: Vec<(f64, f64)> = via_a
            .into_iter()
            .chain(samples.iter().copied())
            .chain(via_b)
            .chain(std::iter::once(b))
            .collect();
        pts.windows(2)
            .map(|w| oriented_term(surface, [a, w[0], w[1]], kind))
            .sum()
    };
    // The fan inscribes the trim: its shortfall falls as 1/panels², so
    // the fine and half-resolution fans extrapolate it away
    // (Richardson). At 16 panels alone the hub's blade-root slivers
    // (5e-6 summed) left ~2e-8 of the propeller's volume (PK-818).
    let fine = fan(&points);
    let coarse: Vec<(f64, f64)> = points.iter().skip(1).step_by(2).copied().collect();
    let sum = fine + (fine - fan(&coarse)) / 3.0;
    sum.is_finite().then_some(sum)
}

/// The trim carrying a chord: its ends' trim parameters, their worst
/// distance off it, and their feet on its edge.
type Carrier<'a> = (&'a Trim<'a>, f64, f64, f64, [Option<f64>; 2]);

/// Newton steps of [`local_foot`] before it gives way to the global search.
const LOCAL_STEPS: usize = 12;

/// The foot of `p` on `curve` by Newton on (C(t) − p) · C′(t) from
/// `seed`, within `[lo, hi]`; `None` where the step does not settle or
/// the distance is not convex there (the caller's global search decides).
fn local_foot(curve: &dyn CurveEval, p: &Point3, seed: f64, lo: f64, hi: f64) -> Option<f64> {
    let mut t = seed.clamp(lo, hi);
    for _ in 0..LOCAL_STEPS {
        let diff = curve.point_at(t) - *p;
        let d1 = curve.derivative_at(t, 1);
        let fp = d1.dot(&d1) + diff.dot(&curve.derivative_at(t, 2));
        if fp.is_nan() || fp <= 0.0 {
            return None;
        }
        let next = (t - diff.dot(&d1) / fp).clamp(lo, hi);
        if (next - t).abs() <= 1e-12 * (1.0 + t.abs()) {
            return Some(next);
        }
        t = next;
    }
    None
}

/// Whether the chart point `q` lies on a pole or apex of `surface`, asked
/// of the chart service in space: its image within the linear tolerance
/// of a singular row's point (the fan from a pole runs along its row).
fn at_pole(surface: &SurfaceGeometry, descriptor: &ChartDescriptor, q: (f64, f64)) -> bool {
    singular_point_near(descriptor, &surface.point_at(q.0, q.1), descriptor.linear).is_some()
}

#[cfg(test)]
#[path = "body_quadrature_sliver_tests.rs"]
mod tests;
