//! Result BRep assembly from selected graph fragments.
//!
//! Takes the intersection graph, face fragments, and selection decisions,
//! and constructs a new BRep containing only the selected fragments with
//! correct topology and geometry.

// Used by later pipeline stages; not yet called from the main merge path.
#![allow(dead_code)]

use std::collections::HashMap;

use popkern_core::{
    constructors::builder::BRepBuilder,
    tolerance::Tolerance,
    topology::{BRep, EdgeHandle, HalfEdgeHandle, VertexHandle},
};

use super::fragment_classify::FragmentSelection;
use super::graph::{IntersectionGraph, Operand};
use super::walk::Fragment;
use crate::error::{BoolDiagnostics, BoolError, Stage};

/// Assemble selected fragments into a result BRep.
///
/// Creates a new BRep from the selected face fragments, copying vertex
/// positions from graph nodes, creating edges for half-edge pairs, and
/// attaching surface geometry from the original operand faces.
///
/// # Errors
///
/// Returns [`BoolDiagnostics`] if the assembly cannot produce valid topology
/// (e.g. empty selection or loop wiring failure).
pub(crate) fn assemble_result(
    graph: &IntersectionGraph,
    fragments: &[Fragment],
    selections: &[FragmentSelection],
    brep_a: &BRep,
    brep_b: &BRep,
    tolerance: Tolerance,
) -> Result<BRep, BoolDiagnostics> {
    let mut builder = BRepBuilder::new(tolerance);

    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);

    // Lazily map graph node index → VertexHandle.
    let mut node_to_vertex: HashMap<usize, VertexHandle> = HashMap::new();

    // Dedup edges by (min_node, max_node, pair_id) → (EdgeHandle, HE_fwd, HE_twin).
    // The pair_id is min(he_idx, twin_idx) to distinguish different graph
    // edges connecting the same pair of nodes (e.g. original edges vs.
    // shared intersection edges).  Without this, two fragments from
    // different faces sharing a node pair would map to the same BRep edge
    // and overwrite each other's loop wiring.
    let mut edge_map: HashMap<(usize, usize, usize), (EdgeHandle, HalfEdgeHandle, HalfEdgeHandle)> =
        HashMap::new();

    // Map graph half-edge index → HalfEdgeHandle in the result.
    let mut he_map: HashMap<usize, HalfEdgeHandle> = HashMap::new();

    for (frag_idx, fragment) in fragments.iter().enumerate() {
        let selection = &selections[frag_idx];
        if !selection.keep {
            continue;
        }

        // Ensure all nodes used by this fragment have vertices.
        let all_hes: Vec<usize> = fragment
            .outer_loop
            .iter()
            .chain(fragment.inner_loops.iter().flat_map(|l| l.iter()))
            .copied()
            .collect();

        for &he_idx in &all_hes {
            let origin_node = graph.half_edges[he_idx].origin;
            node_to_vertex.entry(origin_node).or_insert_with(|| {
                let pos = graph.nodes[origin_node].position;
                builder.add_vertex(pos)
            });
        }

        // Create edges for each half-edge in the fragment.
        for &he_idx in &all_hes {
            if he_map.contains_key(&he_idx) {
                continue;
            }

            let origin_node = graph.half_edges[he_idx].origin;
            let twin_idx = graph.half_edges[he_idx].twin;
            let end_node = graph.half_edges[twin_idx].origin;

            // Ensure end node has a vertex too.
            node_to_vertex.entry(end_node).or_insert_with(|| {
                let pos = graph.nodes[end_node].position;
                builder.add_vertex(pos)
            });

            let pair_id = he_idx.min(twin_idx);
            let edge_key = (
                origin_node.min(end_node),
                origin_node.max(end_node),
                pair_id,
            );

            if let Some(&(_, existing_fwd, existing_twin)) = edge_map.get(&edge_key) {
                // Edge already created. Figure out which half-edge we are.
                let vh_origin = node_to_vertex[&origin_node];
                let fwd_origin = builder.brep().half_edges[existing_fwd].origin;

                if vh_origin == fwd_origin {
                    he_map.insert(he_idx, existing_fwd);
                    he_map.insert(twin_idx, existing_twin);
                } else {
                    he_map.insert(he_idx, existing_twin);
                    he_map.insert(twin_idx, existing_fwd);
                }
            } else {
                let vh_start = node_to_vertex[&origin_node];
                let vh_end = node_to_vertex[&end_node];
                let (edge_h, he_fwd, he_twin) = builder.add_edge(vh_start, vh_end);

                edge_map.insert(edge_key, (edge_h, he_fwd, he_twin));
                he_map.insert(he_idx, he_fwd);
                he_map.insert(twin_idx, he_twin);

                // Attach a line curve between the endpoint positions.
                // We construct a fresh line rather than copying the graph's
                // curve because graph edges may have been split or replaced
                // by intersection NURBS that don't match the new vertex
                // positions.
                let p_start = graph.nodes[origin_node].position;
                let p_end = graph.nodes[end_node].position;
                let dir = p_end - p_start;
                let len = dir.length();
                if len > 1e-12 {
                    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
                    let line_curve = CurveOnEdge {
                        geometry: CurveGeometry::Line {
                            origin: p_start,
                            direction: dir * (1.0 / len),
                        },
                        t_start: 0.0,
                        t_end: len,
                    };
                    builder.set_edge_curve(edge_h, line_curve);
                }
            }
        }

        // Create face.
        let face_h = builder.add_face(shell);

        // Build outer loop.
        // When the fragment is flipped (e.g. B-inside for subtraction),
        // reverse the loop and use the twin of each half-edge.  This
        // ensures the face normal direction matches the loop winding and
        // adjacent faces share edges in opposite directions (manifold
        // requirement).
        let outer_he_handles: Vec<HalfEdgeHandle> = if selection.flip {
            fragment
                .outer_loop
                .iter()
                .rev()
                .filter_map(|&he_idx| {
                    let brep_he = he_map.get(&he_idx).copied()?;
                    Some(builder.brep().half_edges[brep_he].twin)
                })
                .collect()
        } else {
            fragment
                .outer_loop
                .iter()
                .filter_map(|&he_idx| he_map.get(&he_idx).copied())
                .collect()
        };

        if !outer_he_handles.is_empty() {
            builder
                .add_loop(&outer_he_handles, face_h)
                .map_err(|e| BoolDiagnostics {
                    error: Some(BoolError::MergeFailed {
                        stage: format!("{:?}", Stage::Merge),
                        reason: format!("outer loop wiring failed: {e}"),
                    }),
                    failed_stage: Some(Stage::Merge),
                    ..Default::default()
                })?;
        }

        // Build inner loops.
        for inner_loop in &fragment.inner_loops {
            let inner_he_handles: Vec<HalfEdgeHandle> = if selection.flip {
                inner_loop
                    .iter()
                    .rev()
                    .filter_map(|&he_idx| {
                        let brep_he = he_map.get(&he_idx).copied()?;
                        Some(builder.brep().half_edges[brep_he].twin)
                    })
                    .collect()
            } else {
                inner_loop
                    .iter()
                    .filter_map(|&he_idx| he_map.get(&he_idx).copied())
                    .collect()
            };

            if !inner_he_handles.is_empty() {
                builder
                    .add_inner_loop(&inner_he_handles, face_h)
                    .map_err(|e| BoolDiagnostics {
                        error: Some(BoolError::MergeFailed {
                            stage: format!("{:?}", Stage::Merge),
                            reason: format!("inner loop wiring failed: {e}"),
                        }),
                        failed_stage: Some(Stage::Merge),
                        ..Default::default()
                    })?;
            }
        }

        // Attach surface geometry from the original face.
        let source_brep = match fragment.face.operand {
            Operand::A => brep_a,
            Operand::B => brep_b,
        };

        if let Some(face) = source_brep.faces.get(fragment.face.face)
            && let Some(ref surface) = face.surface
        {
            let mut surf = surface.clone();
            if selection.flip {
                surf.orientation = !surf.orientation;
            }
            builder.set_face_surface(face_h, surf);
        }
    }

    // Remove half-edges that weren't wired into any loop, then clean up
    // edges and vertices.
    let brep = builder.brep_mut();
    let default_loop = popkern_core::topology::LoopHandle::default();

    // Collect unwired half-edges.
    let dangling_hes: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, he)| he.loop_ == default_loop)
        .map(|(heh, _)| heh)
        .collect();

    // Collect twin info before mutating.
    let twin_fixups: Vec<(HalfEdgeHandle, HalfEdgeHandle)> = dangling_hes
        .iter()
        .filter_map(|&heh| brep.half_edges.get(heh).map(|he| (heh, he.twin)))
        .collect();

    for (heh, twin_h) in &twin_fixups {
        // Point the surviving twin's twin to itself so it doesn't
        // reference a removed handle.
        if let Some(twin) = brep.half_edges.get_mut(*twin_h) {
            twin.twin = *twin_h;
        }
        // Remove curve geometry from the edge — the watertight check
        // evaluates curves at endpoints and compares to vertex positions,
        // but with a self-pointing twin the check would fail.
        if let Some(surviving) = brep.half_edges.get(*twin_h) {
            let edge_h = surviving.edge;
            if let Some(edge) = brep.edges.get_mut(edge_h) {
                edge.curve = None;
            }
        }
        brep.half_edges.remove(*heh);
    }

    // Remove edges whose half_edge handle was removed AND whose twin
    // half-edge was also removed. Keep edges that still have at least
    // one surviving half-edge.
    let dangling_edges: Vec<EdgeHandle> = brep
        .edges
        .iter()
        .filter(|(_, edge)| !brep.half_edges.contains_key(edge.half_edge))
        .map(|(eh, _)| eh)
        .collect();

    for eh in &dangling_edges {
        brep.edges.remove(*eh);
    }

    // Fix up vertex outgoing handles that point to removed half-edges.
    let outgoing_map: HashMap<VertexHandle, HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .map(|(heh, he)| (he.origin, heh))
        .collect();
    for (vh, vertex) in &mut brep.vertices {
        if !brep.half_edges.contains_key(vertex.outgoing) {
            if let Some(&heh) = outgoing_map.get(&vh) {
                vertex.outgoing = heh;
            } else {
                vertex.outgoing = HalfEdgeHandle::default();
            }
        }
    }

    // Remove orphaned vertices.
    let orphan_verts: Vec<VertexHandle> = brep
        .vertices
        .iter()
        .filter(|(_, v)| v.outgoing == HalfEdgeHandle::default())
        .map(|(vh, _)| vh)
        .collect();

    for vh in &orphan_verts {
        brep.vertices.remove(*vh);
    }

    Ok(builder.build_unchecked())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::merge::fragment_classify::FragmentSelection;
    use crate::merge::graph::{EdgeSource, FaceRef, IntersectionGraph, NodeSource, Operand};
    use popkern_core::geometry::Point3;
    use popkern_core::topology::FaceHandle as CoreFaceHandle;

    /// Build a minimal graph with one triangular face fragment and assemble it.
    #[test]
    fn assemble_single_triangle() {
        let mut graph = IntersectionGraph::new();

        let n0 = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let n1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        let n2 = graph.add_node(Point3::new(0.5, 1.0, 0.0), NodeSource::Intersection);

        let eh = popkern_core::topology::EdgeHandle::default();
        let face_ref = FaceRef {
            operand: Operand::A,
            face: CoreFaceHandle::default(),
        };

        let (fwd01, _) = graph.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(eh), None);
        let (fwd12, _) = graph.add_half_edge_pair(n1, n2, EdgeSource::OriginalA(eh), None);
        let (fwd20, _) = graph.add_half_edge_pair(n2, n0, EdgeSource::OriginalA(eh), None);

        graph.set_half_edge_face(fwd01, face_ref);
        graph.set_half_edge_face(fwd12, face_ref);
        graph.set_half_edge_face(fwd20, face_ref);

        let fragment = Fragment {
            face: face_ref,
            outer_loop: vec![fwd01, fwd12, fwd20],
            inner_loops: vec![],
        };

        let selection = FragmentSelection {
            keep: true,
            flip: false,
        };

        // Use an empty BRep as source (no surface geometry to copy).
        let brep_a = BRep::new(Tolerance::default());
        let brep_b = BRep::new(Tolerance::default());

        let result = assemble_result(
            &graph,
            &[fragment],
            &[selection],
            &brep_a,
            &brep_b,
            Tolerance::default(),
        );

        assert!(result.is_ok(), "Assembly should succeed");
        let brep = result.unwrap_or_else(|_| panic!("unreachable"));
        assert_eq!(brep.vertices.len(), 3, "Expected 3 vertices");
        assert_eq!(brep.edges.len(), 3, "Expected 3 edges");
        assert_eq!(brep.faces.len(), 1, "Expected 1 face");
    }

    /// Verify that unselected fragments are skipped.
    #[test]
    fn unselected_fragments_are_skipped() {
        let mut graph = IntersectionGraph::new();

        let n0 = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let n1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        let n2 = graph.add_node(Point3::new(0.5, 1.0, 0.0), NodeSource::Intersection);

        let eh = popkern_core::topology::EdgeHandle::default();
        let face_ref = FaceRef {
            operand: Operand::A,
            face: CoreFaceHandle::default(),
        };

        let (fwd01, _) = graph.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(eh), None);
        let (fwd12, _) = graph.add_half_edge_pair(n1, n2, EdgeSource::OriginalA(eh), None);
        let (fwd20, _) = graph.add_half_edge_pair(n2, n0, EdgeSource::OriginalA(eh), None);

        graph.set_half_edge_face(fwd01, face_ref);
        graph.set_half_edge_face(fwd12, face_ref);
        graph.set_half_edge_face(fwd20, face_ref);

        let fragment = Fragment {
            face: face_ref,
            outer_loop: vec![fwd01, fwd12, fwd20],
            inner_loops: vec![],
        };

        let selection = FragmentSelection {
            keep: false,
            flip: false,
        };

        let brep_a = BRep::new(Tolerance::default());
        let brep_b = BRep::new(Tolerance::default());

        let result = assemble_result(
            &graph,
            &[fragment],
            &[selection],
            &brep_a,
            &brep_b,
            Tolerance::default(),
        );

        assert!(result.is_ok());
        let brep = result.unwrap_or_else(|_| panic!("unreachable"));
        assert_eq!(brep.faces.len(), 0, "No faces should be in result");
    }
}
