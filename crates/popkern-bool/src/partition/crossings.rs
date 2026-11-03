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

use crate::merge::graph::{EdgeSource, FaceRef, IntersectionGraph, NodeSource, Operand};

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

/// Add a node at every interior crossing of two records of different
/// body pairs on a common face; returns the nodes added or found.
pub(super) fn add_record_crossings(
    graph: &mut IntersectionGraph,
    snap: f64,
    tolerance: &Tolerance,
) -> Vec<usize> {
    let mut points = Vec::new();
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
                    let CurveIntersection::Point { point, .. } = hit else {
                        continue;
                    };
                    let at_end = ends(k1)
                        .iter()
                        .chain(ends(k2).iter())
                        .any(|e| e.distance_to(&point) <= snap);
                    if !at_end {
                        points.push(point);
                    }
                }
            }
        }
    }
    points
        .into_iter()
        .map(|p| match graph.find_node(&p, snap) {
            Some(n) => n,
            None => graph.add_node(p, NodeSource::Intersection),
        })
        .collect()
}
