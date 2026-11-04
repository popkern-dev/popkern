//! Vertex-edge interference (plan §4.3): a candidate lying on an
//! operand edge's or a record's interior gets that incidence too.
//!
//! The interferences that produce candidates see their own entities
//! only: a record's end knows its record and faces, not that it lands
//! mid-span of another record (a dovetail rail's flank trace ending on
//! the joint-shoulder trace) or of an operand edge. This pass tests
//! every candidate against every curve whose box it falls in, by the
//! curve's own closest point; within the candidate's radius plus the
//! curve's tolerance, and clear of the curve's ends by the same reach
//! (an end there is a V-V interference, already a candidate), it adds
//! the interior incidence. Resolution then applies the separation
//! rules to it like any other, and an interior incidence in a group is
//! a T-junction: a pave on that curve, by construction. A point of
//! one record never lies on another record of the same face pair: the
//! branch rule (`branches`) keeps the two branches apart.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, closest_point_in_window};
use popkern_core::topology::BRep;

use super::candidates::{CandidateVertex, Incidence};
use super::region::Bounds;
use crate::intersect::IntersectionRecord;

/// Add every interior incidence of `cands` on the operands' edges and
/// on `records`. `linear` is a curve's tolerance when it has none.
pub(crate) fn add_interior_incidences(
    cands: &mut [CandidateVertex],
    operands: [&BRep; 2],
    records: &[IntersectionRecord],
    linear: f64,
) {
    let r_max = cands.iter().map(|c| c.radius).fold(0.0_f64, f64::max);
    // Candidates by x: each curve scans only its box's x slab.
    let mut by_x: Vec<(f64, usize)> = cands
        .iter()
        .enumerate()
        .map(|(i, c)| (c.position.x, i))
        .collect();
    by_x.sort_by(|p, q| p.0.total_cmp(&q.0));
    let slab = |b: &Bounds| -> Vec<usize> {
        let lo = by_x.partition_point(|(x, _)| *x < b.lo[0]);
        let hi = by_x.partition_point(|(x, _)| *x <= b.hi[0]);
        by_x[lo..hi].iter().map(|&(_, i)| i).collect()
    };
    for (ri, record) in records.iter().enumerate() {
        let curve = &record.curve.curve_3d;
        let Some(bounds) = Bounds::of(&curve.control_points, r_max + linear) else {
            continue;
        };
        let (s0, s1) = curve.parameter_range();
        let ends = [curve.point_at(s0), curve.point_at(s1)];
        for ci in slab(&bounds) {
            let c = &mut cands[ci];
            if !bounds.holds(&c.position)
                || c.incidence.iter().any(|i| {
                    matches!(i, Incidence::Record { record, .. }
                        if *record == ri || super::branches::same_pair(records, ri, *record))
                })
            {
                continue;
            }
            let reach = c.radius + linear;
            if ends.iter().any(|e| e.distance_to(&c.position) < reach) {
                continue;
            }
            let Some(s) = curve.closest_point(&c.position) else {
                continue;
            };
            if curve.point_at(s).distance_to(&c.position) < reach {
                c.incidence.push(Incidence::Record {
                    record: ri,
                    s,
                    end: None,
                });
            }
        }
    }
    for (operand, brep) in operands.iter().enumerate() {
        for (eh, e) in &brep.edges {
            let Some(c) = e.curve.as_ref() else {
                continue;
            };
            let ends = [c.geometry.point_at(c.t_start), c.geometry.point_at(c.t_end)];
            let mut hull = ends.to_vec();
            crate::face_box::conservative_curve_extent(c, &mut hull);
            let tol = e.tolerance.unwrap_or(linear);
            let Some(bounds) = Bounds::of(&hull, r_max + tol) else {
                continue;
            };
            for ci in slab(&bounds) {
                let cand = &mut cands[ci];
                if !bounds.holds(&cand.position)
                    || cand.incidence.iter().any(|i| {
                        matches!(i, Incidence::Edge { operand: o, edge, .. } if *o == operand && *edge == eh)
                    })
                {
                    continue;
                }
                let reach = cand.radius + tol;
                if ends.iter().any(|p| p.distance_to(&cand.position) < reach) {
                    continue;
                }
                if clear_of(&c.geometry, (c.t_start, c.t_end), &cand.position, reach) {
                    continue;
                }
                let t = closest_point_in_window(&c.geometry, &cand.position, c.t_start, c.t_end);
                if c.geometry.point_at(t).distance_to(&cand.position) < reach {
                    cand.incidence.push(Incidence::Edge {
                        operand,
                        edge: eh,
                        t,
                    });
                }
            }
        }
    }
}

/// Deepest halving of a NURBS window in [`clear_of`].
const CLEARANCE_DEPTH: usize = 8;

/// Whether every point of `curve` over `window` is at least `reach`
/// from `p`, proven: by the distance to the whole line or circle
/// (closed form), or for a NURBS by halving the window and bounding
/// each piece by its control points' box (a positively weighted
/// rational piece lies in their hull). A box test passes every point
/// near a large rim (a disc's edge circle boxes the whole disc), and
/// the window search it would otherwise run costs 10-80 µs a point
/// (measured on a disc's NURBS rims). `false` when unproven.
pub(super) fn clear_of(curve: &CurveGeometry, window: (f64, f64), p: &Point3, reach: f64) -> bool {
    match curve {
        CurveGeometry::Nurbs(_) => {
            let (lo, hi) = (window.0.min(window.1), window.0.max(window.1));
            nurbs_clear_of(curve, (lo, hi), p, reach, 0)
        }
        _ => distance_lower_bound(curve, p) >= reach,
    }
}

fn nurbs_clear_of(
    curve: &CurveGeometry,
    (lo, hi): (f64, f64),
    p: &Point3,
    reach: f64,
    depth: usize,
) -> bool {
    let piece = curve.to_nurbs_in_range((lo, hi));
    if piece.control_points.is_empty() || piece.weights.iter().any(|w| *w <= 0.0) {
        return false;
    }
    let (mut min, mut max) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for q in &piece.control_points {
        for (k, x) in [q.x, q.y, q.z].into_iter().enumerate() {
            min[k] = min[k].min(x);
            max[k] = max[k].max(x);
        }
    }
    let gap: f64 = [p.x, p.y, p.z]
        .into_iter()
        .enumerate()
        .map(|(k, x)| (min[k] - x).max(x - max[k]).max(0.0).powi(2))
        .sum::<f64>()
        .sqrt();
    if gap >= reach {
        return true;
    }
    if depth == CLEARANCE_DEPTH {
        return false;
    }
    let mid = 0.5 * (lo + hi);
    nurbs_clear_of(curve, (lo, mid), p, reach, depth + 1)
        && nurbs_clear_of(curve, (mid, hi), p, reach, depth + 1)
}

/// A lower bound on the distance from `p` to any piece of `curve`: the
/// distance to the whole line or circle (exact, closed form), zero for
/// other curves.
fn distance_lower_bound(curve: &CurveGeometry, p: &Point3) -> f64 {
    match curve {
        CurveGeometry::Line { origin, direction } => {
            let len = direction.length();
            if len <= 0.0 {
                return 0.0;
            }
            (*p - *origin).cross(direction).length() / len
        }
        CurveGeometry::Circle {
            center,
            normal,
            radius,
        }
        | CurveGeometry::Arc {
            center,
            normal,
            radius,
            ..
        } => {
            let len = normal.length();
            if len <= 0.0 {
                return 0.0;
            }
            let d = *p - *center;
            let h = d.dot(normal) / len;
            let rho = (d.dot(&d) - h * h).max(0.0).sqrt();
            (h * h + (rho - radius) * (rho - radius)).sqrt()
        }
        _ => 0.0,
    }
}

/// Whether `c` lies on some edge's or record's interior: a T-junction
/// when its group commits. An edge-face hit's first edge is the edge
/// it pierces along, not a junction; every other edge incidence is one
/// this pass added.
pub(crate) fn is_t_junction(c: &CandidateVertex) -> bool {
    let own_edges = usize::from(c.source == super::candidates::SourceKind::EdgeFaceHit);
    let edges = c
        .incidence
        .iter()
        .filter(|i| matches!(i, Incidence::Edge { .. }))
        .count();
    edges > own_edges
        || c.incidence
            .iter()
            .any(|i| matches!(i, Incidence::Record { end: None, .. }))
}

/// Elongate the zone of every candidate lying on an operand edge where
/// that edge crosses a face of the other operand: radius / sin θ, θ the
/// angle between the edge's tangent and the face's tangent plane there,
/// capped at `cap`. The crossing of a curve with a surface at angle θ
/// is fixed to the tolerance across the edge but only to tol / sin θ
/// along it (measured: a thread flank's hit and the record end on the
/// same crown edge 1.6e-5 apart along it, against a 2e-6 radius).
pub(crate) fn add_crossing_zones(cands: &mut [CandidateVertex], operands: [&BRep; 2], cap: f64) {
    use popkern_core::geometry::surface::SurfaceEval;
    for c in cands.iter_mut() {
        let mut zone = c.radius;
        for inc in &c.incidence {
            let Incidence::Edge { operand, edge, t } = *inc else {
                continue;
            };
            let Some(curve) = operands[operand]
                .edges
                .get(edge)
                .and_then(|e| e.curve.as_ref())
            else {
                continue;
            };
            let tangent = curve.geometry.tangent_at(t);
            let tl = tangent.length();
            if tl <= 0.0 {
                continue;
            }
            for face_inc in &c.incidence {
                let Incidence::Face {
                    operand: fo, face, ..
                } = *face_inc
                else {
                    continue;
                };
                if fo == operand {
                    continue;
                }
                let Some(surface) = operands[fo]
                    .faces
                    .get(face)
                    .and_then(|f| f.surface.as_ref())
                else {
                    continue;
                };
                let Some((u, v)) = surface.geometry.closest_point(&c.position) else {
                    continue;
                };
                let n = surface.geometry.normal_at(u, v);
                let nl = n.length();
                if nl <= 0.0 {
                    continue;
                }
                let sin = (tangent.dot(&n) / (tl * nl)).abs();
                let along = if sin * cap > c.radius {
                    c.radius / sin
                } else {
                    cap
                };
                zone = zone.max(along.min(cap));
            }
        }
        c.zone = zone;
    }
}
