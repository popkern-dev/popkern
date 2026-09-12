//! Insert original BRep topology into the intersection graph.
//!
//! Copies vertices, edges, and face-loop structure from one operand's BRep
//! into the shared [`IntersectionGraph`], deduplicating vertices by position.

// Used by later pipeline stages; not yet called from the main merge path.
#![allow(dead_code)]

use std::collections::HashMap;

use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};
use popkern_core::validation::CanonicalForwardState;

use super::graph::{EdgeSource, FaceRef, IntersectionGraph, NodeSource, Operand};

#[path = "closed_side_census.rs"]
pub(crate) mod closed_side_census;
#[path = "graph_build_closed_side.rs"]
mod graph_build_closed_side;
use graph_build_closed_side::{ClosedSideContext, closed_he_runs_curve_forward};

/// Insert all topology from a BRep operand into the graph.
///
/// Returns `(vertex_map, half_edge_map)` mapping original handles to graph
/// indices (node indices and half-edge indices respectively).
pub(crate) fn insert_operand_topology(
    graph: &mut IntersectionGraph,
    brep: &BRep,
    operand: Operand,
    tolerance: f64,
) -> (HashMap<VertexHandle, usize>, HashMap<HalfEdgeHandle, usize>) {
    let mut vertex_map: HashMap<VertexHandle, usize> = HashMap::new();
    let mut he_map: HashMap<HalfEdgeHandle, usize> = HashMap::new();

    // 1. Insert all vertices as nodes (dedup by position).
    for (vh, vertex) in &brep.vertices {
        let node_idx = if let Some(existing) = graph.find_node(&vertex.position, tolerance) {
            existing
        } else {
            let source = match operand {
                Operand::A => NodeSource::OriginalA(vh),
                Operand::B => NodeSource::OriginalB(vh),
            };
            graph.add_node(vertex.position, source)
        };
        vertex_map.insert(vh, node_idx);
    }

    // 2. Insert all edges as half-edge pairs.
    //
    // PK-478 stage 4: the canonical-forward audit is the
    // loop-COVARIANT authority on closed-edge traversal (verdicts
    // are reversal-stable because reversal swaps canonical bits and
    // pcurve directions together), so a determinate verdict
    // registers blindly and the geometric vote witness — whose
    // ray-parity votes are reversal-INVARIANT and therefore provably
    // unable to register a reversed body loop-consistently — is
    // scoped to the audit-indeterminate remainder only.
    let canonical_audit: HashMap<popkern_core::topology::EdgeHandle, CanonicalForwardState> = brep
        .audit_canonical_forward()
        .into_iter()
        .map(|v| (v.edge, v.state))
        .collect();
    let mut closed_side_ctx = ClosedSideContext::new(brep);
    for (eh, edge) in &brep.edges {
        let he_h = edge.half_edge;
        let he = &brep.half_edges[he_h];
        let twin_h = he.twin;

        // Skip if we already inserted this edge (both half-edges mapped).
        if he_map.contains_key(&he_h) {
            continue;
        }

        let twin_he = &brep.half_edges[twin_h];

        let start_node = vertex_map[&he.origin];
        let end_node = vertex_map[&twin_he.origin];

        let source = match operand {
            Operand::A => EdgeSource::OriginalA(eh),
            Operand::B => EdgeSource::OriginalB(eh),
        };

        // CLOSED-EDGE SIDE FIDELITY (audit-first since PK-478 stage
        // 4 — the vote witness below now handles only the edges the
        // canonical-forward audit cannot measure): the chart
        // embedding orients a
        // run by the graph pair's curve window (fwd side = curve
        // increasing-t), and its endpoint-node authority is void on
        // a closed edge (both ends are one node) — so for closed
        // edges the he_map side assignment IS the traversal-direction
        // registration. Mapping the canonical half-edge to the
        // forward side blindly registers its face on the wrong side
        // whenever that face's loop actually traverses the stored
        // curve backward. Measured on the angled tap's bore wall
        // (PK-465): the bottom rim ran backward in the wall's chart,
        // both rims wound the SAME way, the boundary shoelace read
        // the mirror sign, σ flipped, and every wedge match on the
        // wall mirrored into a spur-weave. The witness is pcurve-free
        // (closed-edge pcurves are direction-untrustworthy — legacy
        // producers store edge-aligned pcurves): at the rim midpoint,
        // each face's material must lie LEFT of its traversal
        // (flag-folded normal × tangent points toward that face's
        // boundary centroid), and the two faces must vote opposite
        // ways or there is no verdict and the historical mapping
        // stands. SEAM edges (both half-edges on one face) are
        // exempt: their sides are chart branches of ONE face, and
        // swapping them scrambles the tuned pcurve-branch handling
        // (measured: a cube−sphere seam deviated 0.11 from its
        // face surface under a swap).
        let is_seam = {
            let face_of = |lh| brep.loops.get(lh).map(|l| l.face);
            face_of(he.loop_) == face_of(twin_he.loop_)
        };
        let canonical_forward = if start_node == end_node && !is_seam {
            use CanonicalForwardState as S;
            let (verdict, authority) = match canonical_audit.get(&eh) {
                Some(S::Ok) | Some(S::OkByChain { .. }) => (Some(true), "audit"),
                // A violated canonical bit forces a side SWAP only
                // when the canonical side's own pcurve lies forward:
                // there the chart embedding falls back to the pair's
                // curve-window direction and the side mapping is the
                // only direction authority (the PK-465 bore wall).
                // With honest pcurves the embedding already runs the
                // loop truth and the historical canonical→forward
                // side identity is load-bearing convention — weld
                // partners key on side identity (measured twice:
                // swapping honest-pcurve rims broke the closed-rim
                // shell trim welds, and audit-swapping every Violated
                // edge broke the offset campaign gates).
                Some(S::ViolatedByChain { lying_half_edge }) if *lying_half_edge == he_h => {
                    (Some(false), "audit")
                }
                Some(S::Violated) | Some(S::ViolatedByChain { .. }) => (Some(true), "audit"),
                _ => (
                    closed_he_runs_curve_forward(&mut closed_side_ctx, he_h),
                    "vote-witness",
                ),
            };
            if popkern_core::trace::trace_enabled("boolean") {
                eprintln!(
                    "  [closed-side] {operand:?} edge {eh:?} canonical-forward {verdict:?} ({authority})"
                );
            }
            verdict.unwrap_or(true)
        } else {
            true
        };

        let (pcurve_fwd, pcurve_twin) = if canonical_forward {
            (he.pcurve.clone(), twin_he.pcurve.clone())
        } else {
            (twin_he.pcurve.clone(), he.pcurve.clone())
        };

        let (fwd_idx, twin_idx) = graph.add_half_edge_pair_with_pcurves(
            start_node,
            end_node,
            source,
            edge.curve.clone(),
            pcurve_fwd,
            pcurve_twin,
        );

        if canonical_forward {
            he_map.insert(he_h, fwd_idx);
            he_map.insert(twin_h, twin_idx);
        } else {
            he_map.insert(he_h, twin_idx);
            he_map.insert(twin_h, fwd_idx);
        }
    }

    // 3. Assign faces to half-edges by walking loops.
    for (fh, face) in &brep.faces {
        let face_ref = FaceRef { operand, face: fh };

        // Walk outer loop.
        walk_loop_and_assign(graph, brep, face.outer_loop, face_ref, &he_map);

        // Walk inner loops.
        for &inner_loop_h in &face.inner_loops {
            walk_loop_and_assign(graph, brep, inner_loop_h, face_ref, &he_map);
        }
    }

    (vertex_map, he_map)
}

/// Walk a loop in the BRep and assign the face to each half-edge in the graph.
fn walk_loop_and_assign(
    graph: &mut IntersectionGraph,
    brep: &BRep,
    loop_h: popkern_core::topology::LoopHandle,
    face_ref: FaceRef,
    he_map: &HashMap<HalfEdgeHandle, usize>,
) {
    let start_he = brep.loops[loop_h].half_edge;
    let mut current = start_he;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        if let Some(&graph_he_idx) = he_map.get(&current) {
            graph.set_half_edge_face(graph_he_idx, face_ref);
        }
        current = brep.half_edges[current].next;
        if current == start_he {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;

    #[test]
    fn insert_cube_topology() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let mut graph = IntersectionGraph::new();

        let (vertex_map, he_map) = insert_operand_topology(&mut graph, &cube, Operand::A, 1e-6);

        // A cube has 8 vertices, 12 edges (24 half-edges), 6 faces.
        assert_eq!(graph.nodes.len(), 8, "Expected 8 nodes for a cube");
        assert_eq!(vertex_map.len(), 8);

        assert_eq!(
            graph.half_edges.len(),
            24,
            "Expected 24 half-edges (12 edge pairs) for a cube"
        );
        assert_eq!(he_map.len(), 24);

        assert_eq!(graph.face_edges.len(), 6, "Expected 6 faces in face_edges");

        // Each face of a cube has 4 half-edges.
        for (face_ref, he_indices) in &graph.face_edges {
            assert_eq!(
                he_indices.len(),
                4,
                "Face {:?} should have 4 half-edges, got {}",
                face_ref,
                he_indices.len()
            );
        }
    }

    #[test]
    fn insert_cube_vertex_dedup() {
        // Insert the same cube twice as A and B — vertices should be deduplicated.
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let mut graph = IntersectionGraph::new();

        let (vm_a, _) = insert_operand_topology(&mut graph, &cube, Operand::A, 1e-6);
        let (vm_b, _) = insert_operand_topology(&mut graph, &cube, Operand::B, 1e-6);

        // Both operands share the same 8 positions, so we should still have 8 nodes.
        assert_eq!(graph.nodes.len(), 8);
        assert_eq!(vm_a.len(), 8);
        assert_eq!(vm_b.len(), 8);
    }

    #[test]
    fn inserted_edges_have_curves() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let mut graph = IntersectionGraph::new();
        let (_, he_map) = insert_operand_topology(&mut graph, &cube, Operand::A, 1e-6);

        // Each edge in the cube should have geometry.
        // The forward half-edge of each pair carries the curve.
        let edges_with_curves = graph
            .half_edges
            .iter()
            .filter(|he| he.curve.is_some())
            .count();

        // 12 edges, each pair has 1 curve on the forward half-edge.
        assert_eq!(
            edges_with_curves, 12,
            "Expected 12 half-edges with curves (one per edge pair)"
        );

        // Verify he_map covers all BRep half-edges.
        assert_eq!(he_map.len(), cube.half_edges.len());
    }
}
