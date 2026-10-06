//! Removal of operand edges duplicated by a shared intersection edge
//! (split from `graph_split.rs` under the file-size ratchet).
//!
//! After `split_edges_at_shared_nodes`, an operand edge segment may
//! connect the same pair of nodes as a shared intersection curve on
//! the same face — a coplanar overlap's record running along the
//! operand's own boundary. The walk cannot distinguish the two, so the
//! operand duplicate is removed.
//!
//! A duplicate is the SAME CURVE, not merely the same endpoints. Two
//! distinct curves can join one node pair: on a flush cap of two
//! overlapping cylinders, the cap's outer rim arc and the other
//! cylinder's arc inside the cap both run between the two wall
//! junctions. Removing the rim arc as a "duplicate" left the cap with
//! only its lens to walk — the crescent never became a fragment and
//! the two-cylinder union served a zero-volume body (measured).
//! Endpoint coincidence stays exact; curve coincidence is tested at
//! the operand segment's midpoint against the shared curve, whose
//! parameterization may differ (a fitted record against an exact
//! rational arc), so the test is a point-to-curve distance.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};

use super::graph::{EdgeSource, IntersectionGraph};

/// Endpoint coincidence: both curves end at the same graph nodes.
const ENDPOINT_TOL: f64 = 1e-10;

/// Curve coincidence: the operand segment's midpoint lies within this
/// fraction of the segment's extent from the shared curve. Generous
/// against record fit slop; distinct curves through one node pair
/// separate by a sizeable fraction of the chord.
const MIDPOINT_FRACTION: f64 = 0.05;

/// Remove non-shared half-edges that duplicate a shared half-edge on
/// the same face: same end nodes AND the same curve (see the module
/// docs).
pub(crate) fn remove_collinear_duplicates(graph: &mut IntersectionGraph) {
    let mut to_remove: Vec<usize> = Vec::new();

    for he_indices in graph.face_edges.values() {
        let ends = |idx: usize| -> (Point3, Point3) {
            let he = &graph.half_edges[idx];
            (
                graph.nodes[he.origin].position,
                graph.nodes[graph.half_edges[he.twin].origin].position,
            )
        };
        let shared: Vec<usize> = he_indices
            .iter()
            .copied()
            .filter(|&idx| matches!(graph.half_edges[idx].source, EdgeSource::Shared { .. }))
            .collect();
        if shared.is_empty() {
            continue;
        }

        for &idx in he_indices {
            if matches!(graph.half_edges[idx].source, EdgeSource::Shared { .. }) {
                continue;
            }
            let (p0, p1) = ends(idx);
            // Adaptive squared-distance predicate: the decision is
            // order-independent at the tolerance boundary.
            let close = |a: &Point3, b: &Point3| -> bool {
                popkern_core::adaptive_precision::adaptive_dist_sq_lt(
                    (a.x, a.y, a.z),
                    (b.x, b.y, b.z),
                    ENDPOINT_TOL * ENDPOINT_TOL,
                )
            };
            let is_dup = shared.iter().any(|&s| {
                let (sp0, sp1) = ends(s);
                let same_ends = (close(&sp0, &p0) && close(&sp1, &p1))
                    || (close(&sp0, &p1) && close(&sp1, &p0));
                same_ends && same_curve(curve_of(graph, idx), curve_of(graph, s), p0, p1)
            });
            if is_dup {
                to_remove.push(idx);
            }
        }
    }

    for &idx in &to_remove {
        if let Some(f) = graph.half_edges[idx].face
            && let Some(edges) = graph.face_edges.get_mut(&f)
        {
            edges.retain(|&e| e != idx);
        }
        let origin = graph.half_edges[idx].origin;
        graph.nodes[origin].outgoing.retain(|&e| e != idx);
        graph.half_edges[idx].face = None;
    }
}

/// The curve of half-edge `idx`, read from its twin when the pair
/// stores it on one side only (shared record pairs do).
fn curve_of(graph: &IntersectionGraph, idx: usize) -> Option<&CurveOnEdge> {
    let he = &graph.half_edges[idx];
    he.curve
        .as_ref()
        .or_else(|| graph.half_edges[he.twin].curve.as_ref())
}

/// Whether the operand segment `operand` (ending at `p0`, `p1`) lies
/// on the shared curve `shared`. Without a curve on either side the
/// endpoint match alone decides (the historical straight-segment
/// contract).
fn same_curve(
    operand: Option<&CurveOnEdge>,
    shared: Option<&CurveOnEdge>,
    p0: Point3,
    p1: Point3,
) -> bool {
    let (Some(operand), Some(shared)) = (operand, shared) else {
        return true;
    };
    let mid = operand
        .geometry
        .point_at(0.5 * (operand.t_start + operand.t_end));
    let extent = p0.distance_to(&p1).max(mid.distance_to(&p0));
    let Some(t) = shared.geometry.closest_point(&mid) else {
        return true;
    };
    shared.geometry.point_at(t).distance_to(&mid) <= MIDPOINT_FRACTION * extent
}
