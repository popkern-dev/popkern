//! Seam weld + vertex unification for the V-free consolidation
//! (PK-153): dead boundary sides with coincident antiparallel spans
//! re-twin pairwise (duplicate edge records unify; curve authority
//! prefers the non-linear record), and vertex records duplicated
//! across welded seams merge by position.

use std::collections::BTreeMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};

use super::side_span;

/// Re-twin dead boundary sides whose spans coincide antiparallel,
/// unifying the duplicate edge records (curve authority prefers the
/// side that carries a non-linear curve).
pub(in crate::apex_construction) fn weld_dead_coincident_sides(
    brep: &mut BRep,
) -> (usize, Vec<VertexHandle>) {
    let tol = brep.tolerance.linear().max(1.0e-9) * 10.0;
    let dead: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, data)| !brep.half_edges.contains_key(data.twin))
        .map(|(handle, _)| handle)
        .collect();
    let mut welded = 0usize;
    let mut seam_vertices: Vec<VertexHandle> = Vec::new();
    let mut used: std::collections::BTreeSet<HalfEdgeHandle> = std::collections::BTreeSet::new();
    for i in 0..dead.len() {
        if used.contains(&dead[i]) {
            continue;
        }
        let Some((a_start, a_end)) = side_span(brep, dead[i]) else {
            continue;
        };
        for j in (i + 1)..dead.len() {
            if used.contains(&dead[j]) {
                continue;
            }
            let Some((b_start, b_end)) = side_span(brep, dead[j]) else {
                continue;
            };
            if a_start.distance_to(&b_end) > tol || a_end.distance_to(&b_start) > tol {
                continue;
            }
            let (side_a, side_b) = (dead[i], dead[j]);
            let edge_a = brep.half_edges[side_a].edge;
            let edge_b = brep.half_edges[side_b].edge;
            // Curve authority: keep the non-linear (blend arc) record.
            let is_curved = |edge: popkern_core::topology::EdgeHandle| {
                brep.edges
                    .get(edge)
                    .and_then(|e| e.curve.as_ref())
                    .is_some_and(|c| !matches!(c.geometry, CurveGeometry::Line { .. }))
            };
            let (keep, drop, drop_side) = if edge_a == edge_b {
                (edge_a, None, side_b)
            } else if is_curved(edge_b) && !is_curved(edge_a) {
                (edge_b, Some(edge_a), side_a)
            } else {
                (edge_a, Some(edge_b), side_b)
            };
            brep.half_edges[side_a].twin = side_b;
            brep.half_edges[side_b].twin = side_a;
            brep.half_edges[drop_side].edge = keep;
            brep.half_edges[drop_side].pcurve = None;
            if let Some(drop_edge) = drop {
                brep.edges.remove(drop_edge);
                brep.boundary_edges.remove(&drop_edge);
            }
            brep.boundary_edges.remove(&keep);
            if let Some(edge_data) = brep.edges.get_mut(keep) {
                edge_data.half_edge = side_a;
            }
            for side in [side_a, side_b] {
                if let Some(data) = brep.half_edges.get(side) {
                    seam_vertices.push(data.origin);
                    if let Some(next) = brep.half_edges.get(data.next) {
                        seam_vertices.push(next.origin);
                    }
                }
            }
            used.insert(side_a);
            used.insert(side_b);
            welded += 1;
            break;
        }
    }
    if crate::direct_apex_trace::enabled() {
        for side in &dead {
            if used.contains(side) || !brep.half_edges.contains_key(*side) {
                continue;
            }
            if let Some((start, end)) = side_span(brep, *side) {
                let face = brep
                    .loops
                    .get(brep.half_edges[*side].loop_)
                    .map(|l| format!("{:?}", l.face))
                    .unwrap_or_else(|| "?".into());
                eprintln!(
                    "[direct-apex-trace] consolidation: unwelded dead side {side:?} \
                     face={face} span=({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                    start.x, start.y, start.z, end.x, end.y, end.z
                );
            }
        }
    }
    (welded, seam_vertices)
}

/// Merge coincident vertex records referenced across welded seams.
/// Scoped to the weld-adjacent vertices — position-coincident distinct
/// vertices elsewhere (periodic seams) are by-design and untouched.
pub(in crate::apex_construction) fn merge_coincident_seam_vertices(
    brep: &mut BRep,
    seam: &[VertexHandle],
) -> usize {
    let tol = brep.tolerance.linear().max(1.0e-9) * 10.0;
    let mut unique: Vec<VertexHandle> = seam.to_vec();
    unique.sort();
    unique.dedup();
    let vertices: Vec<(VertexHandle, Point3)> = unique
        .into_iter()
        .filter_map(|handle| {
            brep.vertices
                .get(handle)
                .map(|data| (handle, data.position))
        })
        .collect();
    let mut canonical: BTreeMap<VertexHandle, VertexHandle> = BTreeMap::new();
    for i in 0..vertices.len() {
        for j in (i + 1)..vertices.len() {
            if vertices[i].1.distance_to(&vertices[j].1) <= tol {
                let target = *canonical.get(&vertices[i].0).unwrap_or(&vertices[i].0);
                canonical.entry(vertices[j].0).or_insert(target);
            }
        }
    }
    if canonical.is_empty() {
        return 0;
    }
    let handles: Vec<HalfEdgeHandle> = brep.half_edges.keys().collect();
    for handle in handles {
        let origin = brep.half_edges[handle].origin;
        if let Some(target) = canonical.get(&origin) {
            brep.half_edges[handle].origin = *target;
        }
    }
    let mut merged = 0usize;
    for (duplicate, _) in canonical {
        if brep.vertices.remove(duplicate).is_some() {
            merged += 1;
        }
    }
    merged
}
