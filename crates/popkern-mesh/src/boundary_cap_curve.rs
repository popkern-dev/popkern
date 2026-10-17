//! The B-rep edge curve point a boundary chord splits at (see the
//! parent module): the chord's parameter midpoint, or the curve point
//! nearest a cap's apex.

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, FaceHandle};

/// Samples of each B-rep edge curve, with their bounding box, kept for
/// the whole pass (the B-rep does not change under it).
pub(super) type CurveSamples =
    HashMap<popkern_core::topology::EdgeHandle, (Vec<Point3>, Point3, Point3)>;

/// Samples per edge curve for the nearest-parameter search.
const CURVE_SAMPLES: usize = 256;

/// The B-rep edge curve point halfway (in curve parameter) between `a`
/// and `b`, searched over the edges of `faces`; `None` when no edge
/// curve passes through both within `gate`.
pub(super) fn edge_curve_midpoint(
    brep: &BRep,
    faces: &[FaceHandle],
    a: Point3,
    b: Point3,
    gate: f64,
    cache: &mut CurveSamples,
    toward: Option<Point3>,
) -> Option<Point3> {
    let chord = a.distance_to(&b);
    let mut best: Option<(f64, Point3)> = None;
    for &face in faces {
        let Some(f) = brep.faces.get(face) else {
            continue;
        };
        for lh in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
            let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
                continue;
            };
            let mut cur = start;
            for _ in 0..=brep.half_edges.len() {
                let Some(he) = brep.half_edges.get(cur) else {
                    break;
                };
                if let Some(c) = brep.edges.get(he.edge).and_then(|e| e.curve.as_ref()) {
                    let at = |s: f64| c.t_start + (c.t_end - c.t_start) * s;
                    let (samples, lo, hi) = cache.entry(he.edge).or_insert_with(|| {
                        let pts: Vec<Point3> = (0..=CURVE_SAMPLES)
                            .map(|k| c.geometry.point_at(at(k as f64 / CURVE_SAMPLES as f64)))
                            .collect();
                        let lo = pts
                            .iter()
                            .fold(Point3::new(f64::MAX, f64::MAX, f64::MAX), |m, p| {
                                Point3::new(m.x.min(p.x), m.y.min(p.y), m.z.min(p.z))
                            });
                        let hi = pts
                            .iter()
                            .fold(Point3::new(f64::MIN, f64::MIN, f64::MIN), |m, p| {
                                Point3::new(m.x.max(p.x), m.y.max(p.y), m.z.max(p.z))
                            });
                        (pts, lo, hi)
                    });
                    // Sample boxes under-cover a curve by its sample sag;
                    // pad by the sample spacing's scale.
                    let pad = gate
                        + (hi.x - lo.x).max(hi.y - lo.y).max(hi.z - lo.z) / CURVE_SAMPLES as f64;
                    let outside = |p: Point3| {
                        p.x < lo.x - pad
                            || p.x > hi.x + pad
                            || p.y < lo.y - pad
                            || p.y > hi.y + pad
                            || p.z < lo.z - pad
                            || p.z > hi.z + pad
                    };
                    if !(outside(a) || outside(b)) {
                        let nearest = |p: Point3| -> (f64, f64) {
                            let d = |s: f64| c.geometry.point_at(at(s)).distance_to(&p);
                            let k = (0..samples.len())
                                .min_by(|x, y| {
                                    samples[*x]
                                        .distance_to(&p)
                                        .total_cmp(&samples[*y].distance_to(&p))
                                })
                                .unwrap_or(0);
                            let (mut lo, mut hi) = (
                                (k as f64 - 1.0).max(0.0) / CURVE_SAMPLES as f64,
                                (k as f64 + 1.0).min(CURVE_SAMPLES as f64) / CURVE_SAMPLES as f64,
                            );
                            for _ in 0..48 {
                                let (m1, m2) = (lo + (hi - lo) * 0.382, lo + (hi - lo) * 0.618);
                                if d(m1) <= d(m2) { hi = m2 } else { lo = m1 }
                            }
                            let s = 0.5 * (lo + hi);
                            (s, d(s))
                        };
                        let ((sa, da), (sb, db)) = (nearest(a), nearest(b));
                        // Every face use of the edge must carry a copy of the
                        // chord: a split landing on fewer opens the weld
                        // (PK-873: a sphere corner's rim chords at a BRep
                        // vertex keyed apart from its neighbour's by 5e-16,
                        // split alone for 12 rounds, 1542 open edges).
                        if da.max(db) <= gate && faces.len() >= face_uses(brep, he.edge) {
                            // The split parameter: the chord's midpoint, or
                            // the curve point nearest `toward` (a cap's apex:
                            // its foot), kept to the middle 70% of the span.
                            let s_split = match toward {
                                None => 0.5 * (sa + sb),
                                Some(q) => {
                                    let (lo0, hi0) = (sa.min(sb), sa.max(sb));
                                    let pad = 0.15 * (hi0 - lo0);
                                    let (mut lo, mut hi) = (lo0 + pad, hi0 - pad);
                                    let d = |s: f64| c.geometry.point_at(at(s)).distance_to(&q);
                                    for _ in 0..48 {
                                        let (m1, m2) =
                                            (lo + (hi - lo) * 0.382, lo + (hi - lo) * 0.618);
                                        if d(m1) <= d(m2) { hi = m2 } else { lo = m1 }
                                    }
                                    0.5 * (lo + hi)
                                }
                            };
                            let m = c.geometry.point_at(at(s_split));
                            // The parameter midpoint must sit over the chord
                            // (a closed curve's wrap picks the other half).
                            let mid = Point3::new(
                                0.5 * (a.x + b.x),
                                0.5 * (a.y + b.y),
                                0.5 * (a.z + b.z),
                            );
                            if m.distance_to(&mid) <= 0.5 * chord
                                && best.is_none_or(|(e, _)| da.max(db) < e)
                            {
                                best = Some((da.max(db), m));
                            }
                        }
                    }
                }
                cur = he.next;
                if cur == start {
                    break;
                }
            }
        }
    }
    best.map(|(_, m)| m)
}

/// Face uses of `edge`: its half-edges that lie in a live loop (two for
/// an edge between faces or a seam, one for a sheet's free edge).
fn face_uses(brep: &BRep, edge: popkern_core::topology::EdgeHandle) -> usize {
    let Some(first) = brep.edges.get(edge).map(|e| e.half_edge) else {
        return 0;
    };
    let twin = brep.half_edges.get(first).map(|h| h.twin);
    [Some(first), twin]
        .into_iter()
        .flatten()
        .filter(|&h| {
            brep.half_edges
                .get(h)
                .is_some_and(|h| brep.loops.contains_key(h.loop_))
        })
        .count()
}

#[cfg(test)]
#[path = "boundary_cap_curve_tests.rs"]
mod tests;
