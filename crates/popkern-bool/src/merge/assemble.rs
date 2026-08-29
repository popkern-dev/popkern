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
    topology::{BRep, HalfEdgeHandle, VertexHandle},
};

use super::fragment_classify::FragmentSelection;
#[path = "assemble_curve_adopt.rs"]
mod curve_adopt;
pub(crate) use curve_adopt::subrange_healed_curve;
#[path = "assemble_strand_gate.rs"]
mod strand_gate;
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

    // Map graph half-edge index → BRep HalfEdgeHandle (global across
    // fragments).  When an edge is created for graph half-edge `ghe`,
    // we map BOTH `ghe → he_fwd` and `twin(ghe) → he_rev`.  When a
    // second fragment needs the twin direction, it finds the mapping
    // already present and reuses the opposite BRep half-edge — one BRep
    // edge per shared intersection edge.
    let mut he_map: HashMap<usize, HalfEdgeHandle> = HashMap::new();

    // Process non-flipped fragments first so their half-edge mappings
    // are established before flipped fragments resolve to twins.  This
    // ensures shared intersection edges are used by A-fragments first,
    // with B-fragments getting the twin direction via flip.
    let mut ordered_indices: Vec<usize> = (0..fragments.len())
        .filter(|&i| selections[i].keep)
        .collect();
    ordered_indices.sort_by_key(|&i| (!selections[i].flip) as u8);

    for &frag_idx in &ordered_indices {
        let fragment = &fragments[frag_idx];
        let selection = &selections[frag_idx];
        if popkern_core::trace::trace_enabled("boolean") {
            eprintln!(
                "  [gm-assemble] kept {:?} outer {} inner {} flip={} ghes={:?}",
                fragment.face,
                fragment.outer_loop.len(),
                fragment.inner_loops.len(),
                selection.flip,
                fragment.outer_loop
            );
        }

        // Collect all graph half-edge indices this fragment needs,
        // resolving flips (graph twin) where necessary.  For shared
        // intersection edges, the flip ensures that operand A's
        // fragment uses one direction and operand B's fragment uses
        // the opposite direction — no aliasing.
        let resolve_he = |he_idx: usize| -> usize {
            if selection.flip {
                graph.half_edges[he_idx].twin
            } else {
                he_idx
            }
        };

        let all_needed: Vec<usize> = fragment
            .outer_loop
            .iter()
            .chain(fragment.inner_loops.iter().flat_map(|l| l.iter()))
            .map(|&he| resolve_he(he))
            .collect();

        // Ensure all nodes used by this fragment have vertices.
        for &ghe_idx in &all_needed {
            let origin_node = graph.half_edges[ghe_idx].origin;
            node_to_vertex.entry(origin_node).or_insert_with(|| {
                let pos = graph.nodes[origin_node].position;
                builder.add_vertex(pos)
            });
        }

        // Create BRep edges for each graph half-edge this fragment needs.
        // If the graph half-edge (or its twin) is already mapped, reuse
        // the existing BRep edge — no duplicates.
        for &ghe_idx in &all_needed {
            if he_map.contains_key(&ghe_idx) {
                // Already mapped (either directly or via twin). Skip.
                continue;
            }

            let origin_node = graph.half_edges[ghe_idx].origin;
            let twin_ghe = graph.half_edges[ghe_idx].twin;
            let end_node = graph.half_edges[twin_ghe].origin;

            // Ensure end node has a vertex too.
            node_to_vertex.entry(end_node).or_insert_with(|| {
                let pos = graph.nodes[end_node].position;
                builder.add_vertex(pos)
            });

            let vh_start = node_to_vertex[&origin_node];
            let vh_end = node_to_vertex[&end_node];
            let (edge_h, he_fwd, he_twin) = builder.add_edge(vh_start, vh_end);

            // Map BOTH directions so that whichever graph half-edge a
            // later fragment needs, it finds the correct BRep half-edge
            // without creating a duplicate edge.
            he_map.insert(ghe_idx, he_fwd);
            he_map.insert(twin_ghe, he_twin);

            // Attach curve geometry.
            // For shared (intersection) edges and self-loops: use the graph's
            // actual curve for watertightness (shared geometry principle §4.5).
            // For original/split edges: synthesize a line from endpoints.
            let is_shared = matches!(
                graph.half_edges[ghe_idx].source,
                super::graph::EdgeSource::Shared { .. }
            );

            // Determine which graph half-edge owns the curve.
            let (curve_he, actual_curve) = if graph.half_edges[ghe_idx].curve.is_some() {
                (ghe_idx, graph.half_edges[ghe_idx].curve.as_ref())
            } else {
                (twin_ghe, graph.half_edges[twin_ghe].curve.as_ref())
            };

            if let Some(curve_ref) = actual_curve.filter(|_| is_shared || origin_node == end_node) {
                let mut curve = curve_ref.clone();
                // Align parameterization to match BRep edge direction.
                let curve_origin = graph.half_edges[curve_he].origin;
                if curve_origin != origin_node && origin_node != end_node {
                    std::mem::swap(&mut curve.t_start, &mut curve.t_end);
                }
                builder.set_edge_curve(edge_h, curve);
            } else {
                let p_start = graph.nodes[origin_node].position;
                let p_end = graph.nodes[end_node].position;
                // Original/split edge: ADOPT the stored curve when its
                // endpoints coincide with this edge's nodes — the
                // unconditional line synthesis degraded every curved
                // ORIGINAL boundary edge to its chord (measured on the
                // first torus union: eight grid arcs flattened, each a
                // CurveSurfaceDeviation up to 0.30 on the served body,
                // and OCC read the STEP as a NULL shape). Same rule as
                // cellular assemble: fragments of a SPLIT curved edge
                // whose stored curve spans the whole original fail the
                // endpoint match and keep the line fallback.
                let adopted = actual_curve.and_then(|c| {
                    use popkern_core::geometry::curve::CurveEval;
                    let c0 = c.geometry.point_at(c.t_start);
                    let c1 = c.geometry.point_at(c.t_end);
                    let tol = 1.0e-6;
                    if c0.distance_to(&p_start) < tol && c1.distance_to(&p_end) < tol {
                        Some(c.clone())
                    } else if c0.distance_to(&p_end) < tol && c1.distance_to(&p_start) < tol {
                        let mut rev = c.clone();
                        std::mem::swap(&mut rev.t_start, &mut rev.t_end);
                        Some(rev)
                    } else {
                        // SPLIT piece with a full-span stored curve:
                        // subrange + end-heal (see
                        // `subrange_healed_curve`).
                        subrange_healed_curve(c, &p_start, &p_end)
                    }
                });
                if let Some(curve) = adopted {
                    builder.set_edge_curve(edge_h, curve);
                } else {
                    if popkern_core::trace::trace_enabled("boolean") {
                        eprintln!(
                            "  [line-fallback] ghe {ghe_idx} src {:?} stored={} ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                            graph.half_edges[ghe_idx].source,
                            actual_curve.is_some(),
                            p_start.x,
                            p_start.y,
                            p_start.z,
                            p_end.x,
                            p_end.y,
                            p_end.z
                        );
                    }
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

            // Wire pcurves from graph half-edges onto result half-edges.
            if let Some(ref pcurve) = graph.half_edges[ghe_idx].pcurve {
                builder.set_half_edge_pcurve(he_fwd, pcurve.clone());
            }
            if let Some(ref pcurve) = graph.half_edges[twin_ghe].pcurve {
                builder.set_half_edge_pcurve(he_twin, pcurve.clone());
            }
        }

        // Create face.
        let face_h = builder.add_face(shell);
        let default_loop = popkern_core::topology::LoopHandle::default();

        // For shared intersection edges, the pcurve on each graph
        // half-edge corresponds to a specific face's UV space (forward →
        // face A, twin → face B).  When a fragment uses the "wrong"
        // direction (e.g. face A's inner loop walks the twin direction),
        // the BRep half-edge would get face B's pcurve, which is in the
        // wrong UV space.  This helper returns the correct pcurve for
        // `ghe` given the fragment's `face_ref`: if the graph half-edge's
        // face matches, use its pcurve; otherwise use the twin's pcurve
        // (which belongs to the other face and hence the correct UV space).
        let correct_pcurve_for_face = |ghe: usize,
                                       face_ref: &super::graph::FaceRef|
         -> Option<popkern_core::geometry::pcurve::Pcurve> {
            let ghe_face = &graph.half_edges[ghe].face;
            let twin = graph.half_edges[ghe].twin;
            let twin_face = &graph.half_edges[twin].face;

            // If this graph half-edge's face matches the fragment's
            // face, its pcurve is in the correct UV space.
            if ghe_face.as_ref() == Some(face_ref) {
                return graph.half_edges[ghe].pcurve.clone();
            }
            // If the twin's face matches, the twin's pcurve is in the
            // correct UV space but parameterized for the OPPOSITE walk
            // direction — reverse it. Verbatim borrowing silently
            // flipped closed rims (both rims of a boolean-rebuilt
            // cylinder wall wound the same way; endpoint matching
            // cannot catch direction on a closed curve).
            if twin_face.as_ref() == Some(face_ref) {
                return graph.half_edges[twin].pcurve.as_ref().map(|p| p.reversed());
            }
            // Neither matches (original edges, not shared) — use the
            // graph half-edge's own pcurve as-is.
            graph.half_edges[ghe].pcurve.clone()
        };

        // Resolve graph half-edges to BRep half-edges for a loop.
        // When both directions of an edge are already wired, creates a
        // new BRep edge pair and updates he_map so that later fragments
        // can find the duplicate's twin — preventing dangling half-edges.
        let resolve_loop_hes = |raw_hes: &[usize],
                                flip: bool,
                                hmap: &mut HashMap<usize, HalfEdgeHandle>,
                                nmap: &HashMap<usize, VertexHandle>,
                                bld: &mut BRepBuilder|
         -> Vec<HalfEdgeHandle> {
            let mut handles = Vec::new();
            for &raw_ghe in raw_hes {
                let ghe = if flip {
                    graph.half_edges[raw_ghe].twin
                } else {
                    raw_ghe
                };
                let graph_twin = graph.half_edges[ghe].twin;

                // Try primary mapping.
                if let Some(&brep_he) = hmap.get(&ghe)
                    && bld.brep().half_edges[brep_he].loop_ == default_loop
                {
                    // Fix pcurve: ensure it matches this fragment's face.
                    if let Some(pc) = correct_pcurve_for_face(ghe, &fragment.face) {
                        bld.set_half_edge_pcurve(brep_he, pc);
                    }
                    handles.push(brep_he);
                    continue;
                }
                // Try graph twin's mapping.
                if let Some(&alt) = hmap.get(&graph_twin)
                    && bld.brep().half_edges[alt].loop_ == default_loop
                {
                    // Fix pcurve: the alt half-edge may have the twin's
                    // pcurve which belongs to a different face's UV space.
                    if let Some(pc) = correct_pcurve_for_face(ghe, &fragment.face) {
                        bld.set_half_edge_pcurve(alt, pc);
                    }
                    handles.push(alt);
                    continue;
                }
                // Both wired — create a new edge pair and update the map
                // so subsequent fragments can use the twin.
                if popkern_core::trace::trace_enabled("boolean") {
                    eprintln!(
                        "  [assemble-dup] ghe {ghe} (twin {graph_twin}) for {:?} src {:?}",
                        fragment.face, graph.half_edges[ghe].source
                    );
                }
                let origin_node = graph.half_edges[ghe].origin;
                let end_node = graph.half_edges[graph_twin].origin;
                let Some(&vh_start) = nmap.get(&origin_node) else {
                    continue;
                };
                let Some(&vh_end) = nmap.get(&end_node) else {
                    continue;
                };
                let (_edge_h, he_fwd, he_twin) = bld.add_edge(vh_start, vh_end);
                // Wire correct pcurves: forward gets this fragment's face
                // pcurve, twin gets the other face's pcurve.
                if let Some(pc) = correct_pcurve_for_face(ghe, &fragment.face) {
                    bld.set_half_edge_pcurve(he_fwd, pc);
                }
                if let Some(ref pcurve) = graph.half_edges[graph_twin].pcurve {
                    bld.set_half_edge_pcurve(he_twin, pcurve.clone());
                }
                // Overwrite mappings so later fragments pick up the new
                // twin instead of the already-wired originals.
                hmap.insert(ghe, he_fwd);
                hmap.insert(graph_twin, he_twin);
                handles.push(he_fwd);
            }
            handles
        };

        // Build outer loop.  For flipped fragments, reverse the iteration
        // order to flip winding.
        let raw_outer: Vec<usize> = if selection.flip {
            fragment.outer_loop.iter().rev().copied().collect()
        } else {
            fragment.outer_loop.clone()
        };
        let outer_he_handles = resolve_loop_hes(
            &raw_outer,
            selection.flip,
            &mut he_map,
            &node_to_vertex,
            &mut builder,
        );

        if !outer_he_handles.is_empty() {
            builder.add_loop(&outer_he_handles, face_h).map_err(|e| {
                BoolDiagnostics::from_error(
                    BoolError::MergeFailed {
                        stage: format!("{:?}", Stage::Merge),
                        reason: format!("outer loop wiring failed: {e}"),
                    },
                    Stage::Merge,
                )
            })?;
        }

        // Build inner loops.
        for inner_loop in &fragment.inner_loops {
            let raw_inner: Vec<usize> = if selection.flip {
                inner_loop.iter().rev().copied().collect()
            } else {
                inner_loop.clone()
            };
            let inner_he_handles = resolve_loop_hes(
                &raw_inner,
                selection.flip,
                &mut he_map,
                &node_to_vertex,
                &mut builder,
            );

            if !inner_he_handles.is_empty() {
                builder
                    .add_inner_loop(&inner_he_handles, face_h)
                    .map_err(|e| {
                        BoolDiagnostics::from_error(
                            BoolError::MergeFailed {
                                stage: format!("{:?}", Stage::Merge),
                                reason: format!("inner loop wiring failed: {e}"),
                            },
                            Stage::Merge,
                        )
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

    // HONEST-DECLINE strand gate (PK-304/305) — mechanics and the
    // calibrated benign classes live in `assemble_strand_gate.rs`.
    strand_gate::decline_on_lethal_strands(
        brep,
        graph,
        &he_map,
        &dangling_hes,
        brep_a,
        brep_b,
        &tolerance,
    )?;

    // Mask remaining benign dangling halves (self-twin + curve strip)
    // and clean up edges/vertices — `assemble_strand_gate.rs`.
    strand_gate::mask_dangling_halves(brep, graph, &he_map, &dangling_hes);

    Ok(builder.build_unchecked())
}

#[cfg(test)]
#[path = "assemble_tests.rs"]
mod tests;
