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
                for hit in found.intersections {
                    let CurveIntersection::Point { point, t1, t2 } = hit else {
                        continue;
                    };
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
