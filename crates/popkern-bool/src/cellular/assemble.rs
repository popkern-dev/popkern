//! Order-independent BRep assembly from selected graph fragments.
//!
//! Unlike `merge/assemble.rs`, this assembly pre-creates all edges from
//! a canonical set (sorted by node indices) before wiring loops. This
//! makes the output independent of fragment iteration order, enabling
//! BVH-based face-pair filtering.

use std::collections::HashMap;

use popkern_core::{
    constructors::builder::BRepBuilder,
    geometry::curve::{CurveGeometry, CurveOnEdge},
    tolerance::Tolerance,
    topology::{BRep, HalfEdgeHandle, VertexHandle},
};

use crate::error::{BoolDiagnostics, BoolError, Stage};
use crate::merge::fragment_classify::FragmentSelection;
use crate::merge::graph::{IntersectionGraph, Operand};
use crate::merge::walk::Fragment;

/// Assemble a BRep from selected fragments using order-independent edge creation.
///
/// 1. Collect all vertices needed by kept fragments
/// 2. Collect all unique edges (canonical: min_node, max_node)
/// 3. Create all vertices and edges upfront
/// 4. Wire fragment loops using the pre-created edges
/// 5. Attach surface geometry
pub(crate) fn cellular_assemble(
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

    // Step 1: Collect all nodes needed by kept fragments.
    let mut needed_nodes = std::collections::BTreeSet::new();
    let kept: Vec<usize> = (0..fragments.len())
        .filter(|&i| selections[i].keep)
        .collect();

    for &fi in &kept {
        let sel = &selections[fi];
        let frag = &fragments[fi];
        for &ghe in frag
            .outer_loop
            .iter()
            .chain(frag.inner_loops.iter().flat_map(|l| l.iter()))
        {
            let ghe = if sel.flip {
                graph.half_edges[ghe].twin
            } else {
                ghe
            };
            needed_nodes.insert(graph.half_edges[ghe].origin);
            let twin = graph.half_edges[ghe].twin;
            needed_nodes.insert(graph.half_edges[twin].origin);
        }
    }

    // Step 2: Create vertices (sorted by node index for determinism).
    let mut node_to_vertex: HashMap<usize, VertexHandle> = HashMap::new();
    for &node in &needed_nodes {
        let vh = builder.add_vertex(graph.nodes[node].position);
        node_to_vertex.insert(node, vh);
    }

    // Step 3: Collect all unique edges (canonical: sorted node pair).
    // Key: (min_node, max_node), Value: graph half-edge index (for curve data).
    let mut canonical_edges: std::collections::BTreeMap<(usize, usize), usize> =
        std::collections::BTreeMap::new();

    for &fi in &kept {
        let sel = &selections[fi];
        let frag = &fragments[fi];
        for &raw_ghe in frag
            .outer_loop
            .iter()
            .chain(frag.inner_loops.iter().flat_map(|l| l.iter()))
        {
            let ghe = if sel.flip {
                graph.half_edges[raw_ghe].twin
            } else {
                raw_ghe
            };
            let origin = graph.half_edges[ghe].origin;
            let twin = graph.half_edges[ghe].twin;
            let end = graph.half_edges[twin].origin;
            let key = if origin <= end {
                (origin, end)
            } else {
                (end, origin)
            };
            canonical_edges.entry(key).or_insert(ghe);
        }
    }

    // Step 4: Create all edges upfront (deterministic order from BTreeMap).
    // Map: (node_a, node_b) → (he_fwd, he_twin) where fwd goes a→b.
    let mut edge_he_map: HashMap<(usize, usize), (HalfEdgeHandle, HalfEdgeHandle)> = HashMap::new();

    for (&(na, nb), &sample_ghe) in &canonical_edges {
        let vh_a = node_to_vertex[&na];
        let vh_b = node_to_vertex[&nb];
        let (edge_h, he_fwd, he_twin) = builder.add_edge(vh_a, vh_b);

        edge_he_map.insert((na, nb), (he_fwd, he_twin));

        // Attach curve geometry.
        // For shared (intersection) edges: use the graph's actual intersection
        // curve. Both incident faces share the exact same CurveOnEdge, enforcing
        // the shared geometry principle (spec §4.5, §9.2) for watertightness.
        // For original/split edges and self-loops: use the graph's stored curve
        // if available, otherwise synthesize a line from endpoints.
        let is_shared = matches!(
            graph.half_edges[sample_ghe].source,
            crate::merge::graph::EdgeSource::Shared { .. }
        );
        let graph_curve = graph.half_edges[sample_ghe].curve.as_ref().or_else(|| {
            let twin = graph.half_edges[sample_ghe].twin;
            graph.half_edges[twin].curve.as_ref()
        });

        if is_shared || na == nb {
            // Shared edge or self-loop: use the graph's curve (shared between
            // both incident faces via the single CurveOnEdge on the Edge).
            if let Some(curve) = graph_curve {
                builder.set_edge_curve(edge_h, curve.clone());
            }
        } else {
            // Original/split edge: synthesize line from endpoints.
            let p_start = graph.nodes[na].position;
            let p_end = graph.nodes[nb].position;
            let dir = p_end - p_start;
            let len = dir.length();
            if len > 1e-12 {
                let line = CurveOnEdge {
                    geometry: CurveGeometry::Line {
                        origin: p_start,
                        direction: dir * (1.0 / len),
                    },
                    t_start: 0.0,
                    t_end: len,
                };
                builder.set_edge_curve(edge_h, line);
            }
        }
    }

    // Step 5: Wire fragment loops.
    // Track which half-edges have been used (to handle shared edges).
    // When both fwd and twin are used, create a duplicate edge.
    let mut used_hes: std::collections::HashSet<HalfEdgeHandle> = std::collections::HashSet::new();

    for &fi in &kept {
        let sel = &selections[fi];
        let frag = &fragments[fi];
        let face_h = builder.add_face(shell);

        // Resolve all graph half-edges to BRep half-edges upfront,
        // before calling builder methods that borrow mutably.
        let mut resolve_loop = |raw_hes: &[usize], flip: bool| -> Vec<HalfEdgeHandle> {
            let mut result = Vec::new();
            for &raw_ghe in raw_hes {
                let ghe = if flip {
                    graph.half_edges[raw_ghe].twin
                } else {
                    raw_ghe
                };
                let origin = graph.half_edges[ghe].origin;
                let twin_idx = graph.half_edges[ghe].twin;
                let end = graph.half_edges[twin_idx].origin;

                let (na, nb) = if origin <= end {
                    (origin, end)
                } else {
                    (end, origin)
                };

                let Some(&(he_fwd, he_twin)) = edge_he_map.get(&(na, nb)) else {
                    continue;
                };

                // fwd goes na→nb. We need origin→end.
                let brep_he = if origin == na { he_fwd } else { he_twin };

                if !used_hes.contains(&brep_he) {
                    used_hes.insert(brep_he);
                    result.push(brep_he);
                } else {
                    // Try twin
                    let twin_brep = if brep_he == he_fwd { he_twin } else { he_fwd };
                    if !used_hes.contains(&twin_brep) {
                        used_hes.insert(twin_brep);
                        result.push(twin_brep);
                    }
                }
            }
            result
        };

        let raw_outer: Vec<usize> = if sel.flip {
            frag.outer_loop.iter().rev().copied().collect()
        } else {
            frag.outer_loop.clone()
        };
        let outer_hes = resolve_loop(&raw_outer, sel.flip);

        if !outer_hes.is_empty() {
            builder
                .add_loop(&outer_hes, face_h)
                .map_err(|e| BoolDiagnostics {
                    error: Some(BoolError::MergeFailed {
                        stage: format!("{:?}", Stage::Merge),
                        reason: format!("cellular outer loop: {e}"),
                    }),
                    failed_stage: Some(Stage::Merge),
                    ..Default::default()
                })?;
        }

        for inner in &frag.inner_loops {
            let raw_inner: Vec<usize> = if sel.flip {
                inner.iter().rev().copied().collect()
            } else {
                inner.clone()
            };
            let inner_hes = resolve_loop(&raw_inner, sel.flip);
            if !inner_hes.is_empty() {
                let _ = builder.add_inner_loop(&inner_hes, face_h);
            }
        }

        // Attach surface geometry
        let source_brep = match frag.face.operand {
            Operand::A => brep_a,
            Operand::B => brep_b,
        };
        if let Some(face) = source_brep.faces.get(frag.face.face)
            && let Some(ref surface) = face.surface
        {
            let mut surf = surface.clone();
            if sel.flip {
                surf.orientation = !surf.orientation;
            }
            builder.set_face_surface(face_h, surf);
        }
    }

    // Cleanup: remove unwired half-edges, dangling edges, orphan vertices.
    cleanup_unwired(&mut builder);

    Ok(builder.build_unchecked())
}

/// Remove half-edges not wired into any loop, and their orphaned edges/vertices.
fn cleanup_unwired(builder: &mut BRepBuilder) {
    let brep = builder.brep_mut();
    let default_loop = popkern_core::topology::LoopHandle::default();

    let dangling: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, he)| he.loop_ == default_loop)
        .map(|(h, _)| h)
        .collect();

    let twins: Vec<(HalfEdgeHandle, HalfEdgeHandle)> = dangling
        .iter()
        .filter_map(|&h| brep.half_edges.get(h).map(|he| (h, he.twin)))
        .collect();

    for (h, twin_h) in &twins {
        if let Some(twin) = brep.half_edges.get_mut(*twin_h) {
            twin.twin = *twin_h;
        }
        if let Some(surviving) = brep.half_edges.get(*twin_h) {
            let edge_h = surviving.edge;
            if let Some(edge) = brep.edges.get_mut(edge_h) {
                edge.curve = None;
            }
        }
        brep.half_edges.remove(*h);
    }

    // Remove fully-orphaned edges
    let dead_edges: Vec<_> = brep
        .edges
        .iter()
        .filter(|(_, e)| !brep.half_edges.contains_key(e.half_edge))
        .filter(|(eh, _)| !brep.half_edges.iter().any(|(_, he)| he.edge == *eh))
        .map(|(eh, _)| eh)
        .collect();
    for eh in dead_edges {
        brep.edges.remove(eh);
    }

    // Fix vertex outgoing handles
    let outgoing: HashMap<VertexHandle, HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .map(|(h, he)| (he.origin, h))
        .collect();
    for (vh, v) in &mut brep.vertices {
        if !brep.half_edges.contains_key(v.outgoing) {
            v.outgoing = outgoing
                .get(&vh)
                .copied()
                .unwrap_or(HalfEdgeHandle::default());
        }
    }

    // Remove orphan vertices
    let orphans: Vec<VertexHandle> = brep
        .vertices
        .iter()
        .filter(|(_, v)| v.outgoing == HalfEdgeHandle::default())
        .map(|(h, _)| h)
        .collect();
    for vh in orphans {
        brep.vertices.remove(vh);
    }
}
