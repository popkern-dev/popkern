//! Records of different body pairs crossing on a face they share
//! (design 7.1 step 3).
//!
//! Two bodies' records meet only at their ends or at T-junctions. With
//! three bodies, a record of (A, B) and a record of (A, C) both lie on
//! a face of A and cross where all three surfaces meet: an interior
//! point of both, which no two-body interference finds. Each crossing
//! becomes a node, and the graph's junction split then cuts both
//! records there.

use std::collections::BTreeMap;

use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::geometry::intersection::{CurveIntersection, intersect_curves};
use popkern_core::tolerance::Tolerance;

use crate::merge::graph::{EdgeSource, FaceRef, IntersectionGraph, Operand};

/// The record half-edge pairs registered to each face, by key (the
/// lower half-edge index), with the body pair that made them.
fn records_by_face(graph: &IntersectionGraph) -> BTreeMap<FaceRef, Vec<(usize, [Operand; 2])>> {
    let mut out: BTreeMap<FaceRef, Vec<(usize, [Operand; 2])>> = BTreeMap::new();
    for (face, hes) in &graph.face_edges {
        let mut keys: Vec<(usize, [Operand; 2])> = hes
            .iter()
            .filter_map(|&h| match graph.half_edges[h].source {
                EdgeSource::Shared { operands, .. } => {
                    Some((h.min(graph.half_edges[h].twin), operands))
                }
                _ => None,
            })
            .collect();
        keys.sort_by_key(|k| k.0);
        keys.dedup_by_key(|k| k.0);
        out.insert(*face, keys);
    }
    out
}

fn curve_of(graph: &IntersectionGraph, key: usize) -> Option<&CurveOnEdge> {
    let he = &graph.half_edges[key];
    he.curve
        .as_ref()
        .or_else(|| graph.half_edges[he.twin].curve.as_ref())
}

/// Split every record at its interior crossings with records of other
/// body pairs on a common face; returns the splits made. Each crossing
/// carries its parameter on both records, so a record is cut at its
/// crossings in curve order directly (the graph's general node splitter
/// restarts its scan after every split: quadratic in a dense bundle).
pub(super) fn split_at_record_crossings(
    graph: &mut IntersectionGraph,
    snap: f64,
    tolerance: &Tolerance,
) -> usize {
    let mut on_record: BTreeMap<usize, Vec<(f64, popkern_core::geometry::Point3)>> =
        BTreeMap::new();
    for keys in records_by_face(graph).values() {
        for (i, &(k1, p1)) in keys.iter().enumerate() {
            for &(k2, p2) in &keys[i + 1..] {
                if p1 == p2 {
                    continue;
                }
                let (Some(c1), Some(c2)) = (curve_of(graph, k1), curve_of(graph, k2)) else {
                    continue;
                };
                let r1 = (c1.t_start.min(c1.t_end), c1.t_start.max(c1.t_end));
                let r2 = (c2.t_start.min(c2.t_end), c2.t_start.max(c2.t_end));
                if r1.0 >= r1.1 || r2.0 >= r2.1 {
                    continue;
                }
                let Ok(found) = intersect_curves(&c1.geometry, r1, &c2.geometry, r2, tolerance)
                else {
                    continue;
                };
                let ends = |k: usize| {
                    let he = &graph.half_edges[k];
                    [
                        graph.nodes[he.origin].position,
                        graph.nodes[graph.half_edges[he.twin].origin].position,
                    ]
                };
                let hits: Vec<(f64, f64, popkern_core::geometry::Point3)> = found
                    .intersections
                    .iter()
                    .filter_map(|hit| match *hit {
                        CurveIntersection::Point { point, t1, t2 } => Some((t1, t2, point)),
                        _ => None,
                    })
                    .collect();
                let hits = one_per_pass((c1, r1), (c2, r2), hits, tolerance.linear());
                for (t1, t2, point) in hits {
                    let at_end = ends(k1)
                        .iter()
                        .chain(ends(k2).iter())
                        .any(|e| e.distance_to(&point) <= snap);
                    if !at_end {
                        on_record.entry(k1).or_default().push((t1, point));
                        on_record.entry(k2).or_default().push((t2, point));
                    }
                }
            }
        }
    }
    // One vertex per crossing place: a crossing adopts a node already
    // within the graph identity's adoption reach (20 tolerances, as a
    // record end adopts one), else an earlier crossing there. Triple
    // points of different pairs crowd where many bodies meet near one
    // point (measured: two 2.2e-6 apart in a six-rod star, minting a
    // micro edge the wall walked four times: a radial tie).
    let wide = tolerance.linear() * 20.0;
    let mut anchors: Vec<popkern_core::geometry::Point3> = Vec::new();
    for list in on_record.values_mut() {
        for (_, p) in list.iter_mut() {
            let at = match graph.find_node(p, wide) {
                Some(n) => graph.nodes[n].position,
                None => match anchors.iter().find(|a| a.distance_to(p) <= wide) {
                    Some(a) => *a,
                    None => {
                        anchors.push(*p);
                        *p
                    }
                },
            };
            *p = at;
        }
    }
    // In each record's own direction: its forward half-edge runs the
    // curve window from `t_start` to `t_end`.
    for (&key, list) in on_record.iter_mut() {
        let forward = curve_of(graph, key).is_none_or(|c| c.t_start <= c.t_end);
        list.sort_by(|a, b| {
            if forward {
                a.0.total_cmp(&b.0)
            } else {
                b.0.total_cmp(&a.0)
            }
        });
        list.dedup_by(|a, b| a.1.distance_to(&b.1) <= snap);
    }
    let records: Vec<(usize, &[(f64, popkern_core::geometry::Point3)])> =
        on_record.iter().map(|(&k, v)| (k, &v[..])).collect();
    crate::merge::graph_paves::split_records_at_paves(graph, &records, snap)
}

/// The crossings of curves `c1` (over `r1`) and `c2` (over `r2`), as
/// `(t1, t2, point)`, with every pair joined by a stretch of `c1` lying
/// within `linear` of `c2` taken as ONE meeting at the stretch's
/// middle: the curves pass within the tolerance there, so the two
/// crossings are one place (plan §5.4: a pass within ε is one vertex).
/// Two balls overlapping by half a tolerance on a box's top face have
/// equators crossing twice 2.4e-3 apart, at most 5e-7 apart between;
/// cut at both, the records bounded a lens the partition joined both
/// sides of (SidesJoined), and it fell back. A closed `c1` is also
/// tried the other way round its seam.
fn one_per_pass(
    (c1, r1): (&CurveOnEdge, (f64, f64)),
    (c2, r2): (&CurveOnEdge, (f64, f64)),
    mut hits: Vec<(f64, f64, popkern_core::geometry::Point3)>,
    linear: f64,
) -> Vec<(f64, f64, popkern_core::geometry::Point3)> {
    use popkern_core::geometry::curve::{CurveEval, closest_point_in_window};
    const STATIONS: usize = 16;
    let g1: &dyn CurveEval = &c1.geometry;
    let g2: &dyn CurveEval = &c2.geometry;
    let period = r1.1 - r1.0;
    let closed = g1.point_at(r1.0).distance_to(&g1.point_at(r1.1)) <= linear;
    // c1 at parameter `t` folded into its window (a closed curve wraps).
    let at = |t: f64| {
        let t = if closed && t > r1.1 { t - period } else { t };
        g1.point_at(t.clamp(r1.0, r1.1))
    };
    let foot = |p: &popkern_core::geometry::Point3| {
        let t = closest_point_in_window(g2, p, r2.0, r2.1);
        (t, g2.point_at(t))
    };
    // The stretch of c1 from `a` up to `b` (b > a) within `linear` of c2.
    let within = |a: f64, b: f64| {
        (1..STATIONS).all(|k| {
            let p = at(a + (b - a) * k as f64 / STATIONS as f64);
            foot(&p).1.distance_to(&p) <= linear
        })
    };
    hits.sort_by(|x, y| x.0.total_cmp(&y.0));
    if hits.len() < 2 {
        return hits;
    }
    // Clusters of consecutive crossings, each as its first and last
    // parameter on c1 (the last may pass the window's end on a closed
    // c1) and the hits it holds.
    let mut clusters: Vec<(f64, f64, usize)> = Vec::new();
    for h in &hits {
        match clusters.last_mut() {
            Some(c) if within(c.1, h.0) => {
                c.1 = h.0;
                c.2 += 1;
            }
            _ => clusters.push((h.0, h.0, 1)),
        }
    }
    if closed && clusters.len() > 1 {
        let (first, last) = (clusters[0], clusters[clusters.len() - 1]);
        if within(last.1, first.0 + period) {
            clusters[0] = (last.0, first.1 + period, first.2 + last.2);
            clusters.pop();
        }
    }
    let mut out = Vec::with_capacity(clusters.len());
    let mut k = 0;
    for &(a, b, n) in &clusters {
        if n == 1 {
            // A lone crossing keeps its exact point.
            if let Some(h) = hits.iter().find(|h| h.0 == a) {
                out.push(*h);
            }
            k += 1;
            continue;
        }
        // The middle of the stretch, between the curves.
        let tm = 0.5 * (a + b);
        let p = at(tm);
        let (t2, q) = foot(&p);
        let t1 = if closed && tm > r1.1 { tm - period } else { tm };
        let mid = popkern_core::geometry::Point3::new(
            0.5 * (p.x + q.x),
            0.5 * (p.y + q.y),
            0.5 * (p.z + q.z),
        );
        out.push((t1, t2, mid));
        k += n;
    }
    debug_assert!(k == hits.len());
    out
}
