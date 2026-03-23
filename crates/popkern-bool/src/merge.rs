//! Stage 5: Merge selected face fragments into a single result BRep.
//!
//! Takes faces selected by Stage 4 from two split operands and combines them
//! into one valid closed BRep by copying topology, merging coincident vertices,
//! and stitching boundary edges.

use std::collections::HashMap;

use popkern_core::{
    tolerance::Tolerance,
    topology::{
        BRep, Edge, EdgeHandle, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle,
        Shell, ShellHandle, Solid, Vertex, VertexHandle,
    },
};

use crate::error::BoolDiagnostics;
use crate::select::SelectionResult;
use crate::split::SplitResult;

/// Maps old handles from a source BRep to new handles in the result BRep.
#[derive(Debug, Default)]
struct HandleMap {
    /// Old vertex handle to new vertex handle.
    vertices: HashMap<VertexHandle, VertexHandle>,
    /// Old half-edge handle to new half-edge handle.
    half_edges: HashMap<HalfEdgeHandle, HalfEdgeHandle>,
    /// Old edge handle to new edge handle.
    edges: HashMap<EdgeHandle, EdgeHandle>,
    /// Old loop handle to new loop handle.
    loops: HashMap<LoopHandle, LoopHandle>,
    /// Old face handle to new face handle.
    faces: HashMap<FaceHandle, FaceHandle>,
}

/// Build a result BRep from selected face fragments of two split operands.
///
/// Copies selected faces from both operands into a new BRep, merges coincident
/// vertices (within tolerance), stitches boundary half-edges as twins, and
/// optionally flips surface normals as directed by the selection result.
///
/// # Errors
///
/// Returns [`BoolDiagnostics`] if the merge fails (e.g. topology cannot be
/// stitched into a closed manifold).
pub(crate) fn merge_faces(
    split: &SplitResult,
    selection: &SelectionResult,
    tolerance: &Tolerance,
) -> Result<BRep, BoolDiagnostics> {
    let mut result = BRep::new(*tolerance);

    // Create a solid and shell in the result.
    let solid_handle = result.solids.insert(Solid {
        shell: ShellHandle::default(),
    });
    let shell_handle = result.shells.insert(Shell {
        face: FaceHandle::default(),
        solid: solid_handle,
    });
    result.solids[solid_handle].shell = shell_handle;

    // Maps from old handles (in operands) to new handles (in result).
    let mut map_a = HandleMap::default();
    let mut map_b = HandleMap::default();

    // Copy selected faces from operand A.
    let mut first_face = None;
    for face_h in split.a.faces.keys() {
        if selection.keep_a.contains(&face_h) {
            let new_fh = copy_face(
                &split.a,
                face_h,
                &mut result,
                shell_handle,
                &mut map_a,
                tolerance,
            )?;
            if first_face.is_none() {
                first_face = Some(new_fh);
            }
            // Flip surface if selection says to.
            if selection.flip_a.contains(&face_h) {
                flip_surface(&mut result, new_fh);
            }
        }
    }

    // Copy selected faces from operand B.
    for face_h in split.b.faces.keys() {
        if selection.keep_b.contains(&face_h) {
            let new_fh = copy_face(
                &split.b,
                face_h,
                &mut result,
                shell_handle,
                &mut map_b,
                tolerance,
            )?;
            if first_face.is_none() {
                first_face = Some(new_fh);
            }
            if selection.flip_b.contains(&face_h) {
                flip_surface(&mut result, new_fh);
            }
        }
    }

    // Point the shell at its first face.
    if let Some(fh) = first_face {
        result.shells[shell_handle].face = fh;
    }

    // Stitch boundary edges: pair up unpaired half-edges sharing the same
    // two vertices in opposite order, and merge their Edge entities.
    stitch_boundary_edges(&mut result);

    // Remove orphaned edges (edges whose half-edge no longer exists or
    // whose half-edge points to a different edge).
    remove_orphan_edges(&mut result);

    Ok(result)
}

/// Copy a single face (and all its topology) from `source` into `result`.
///
/// Vertices are deduplicated by position within tolerance against existing
/// result vertices.
fn copy_face(
    source: &BRep,
    face_h: FaceHandle,
    result: &mut BRep,
    shell: ShellHandle,
    map: &mut HandleMap,
    tolerance: &Tolerance,
) -> Result<FaceHandle, BoolDiagnostics> {
    let face = &source.faces[face_h];

    // Create the face and outer loop in the result (with placeholder handles).
    let new_loop_h = result.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: FaceHandle::default(),
    });
    map.loops.insert(face.outer_loop, new_loop_h);

    let new_face_h = result.faces.insert(Face {
        outer_loop: new_loop_h,
        inner_loops: Vec::new(),
        shell,
        surface: face.surface.clone(),
    });
    map.faces.insert(face_h, new_face_h);
    result.loops[new_loop_h].face = new_face_h;

    // Copy the outer loop's half-edges.
    copy_loop_half_edges(source, face.outer_loop, result, new_loop_h, map, tolerance)?;

    // Copy inner loops (holes).
    for &inner_loop_h in &face.inner_loops {
        let new_inner_h = result.loops.insert(Loop {
            half_edge: HalfEdgeHandle::default(),
            face: new_face_h,
        });
        map.loops.insert(inner_loop_h, new_inner_h);
        copy_loop_half_edges(source, inner_loop_h, result, new_inner_h, map, tolerance)?;
        result.faces[new_face_h].inner_loops.push(new_inner_h);
    }

    Ok(new_face_h)
}

/// Copy all half-edges of a loop from `source` into `result`.
///
/// Walks around the loop following `next` pointers. Vertices are deduplicated
/// by position. Edges and half-edges are created fresh; twin pointers are set
/// to self (unpaired) and fixed up later by `stitch_boundary_edges`.
fn copy_loop_half_edges(
    source: &BRep,
    loop_h: LoopHandle,
    result: &mut BRep,
    new_loop_h: LoopHandle,
    map: &mut HandleMap,
    tolerance: &Tolerance,
) -> Result<(), BoolDiagnostics> {
    let loop_ = &source.loops[loop_h];
    let start_he = loop_.half_edge;

    // Collect half-edge handles around the loop.
    let mut he_handles = Vec::new();
    let mut current = start_he;
    loop {
        he_handles.push(current);
        current = source.half_edges[current].next;
        if current == start_he {
            break;
        }
    }

    // First pass: create all half-edges, edges, and vertices.
    let mut new_he_handles = Vec::with_capacity(he_handles.len());
    for &he_h in &he_handles {
        let he = &source.half_edges[he_h];

        // Map or create the origin vertex.
        let new_origin = map_or_create_vertex(source, he.origin, result, map, tolerance);

        // Reuse existing edge if already mapped (intra-operand twin sharing),
        // otherwise create a new edge.
        let new_edge_h = if let Some(&existing) = map.edges.get(&he.edge) {
            existing
        } else {
            let edge = &source.edges[he.edge];
            let new_eh = result.edges.insert(Edge {
                half_edge: HalfEdgeHandle::default(),
                curve: edge.curve.clone(),
            });
            map.edges.insert(he.edge, new_eh);
            new_eh
        };

        // Create the half-edge with self-twin (unpaired for now).
        let new_he_h = result.half_edges.insert(HalfEdge {
            origin: new_origin,
            twin: HalfEdgeHandle::default(), // placeholder
            next: HalfEdgeHandle::default(), // placeholder
            prev: HalfEdgeHandle::default(), // placeholder
            edge: new_edge_h,
            loop_: new_loop_h,
        });
        // Point twin at self until stitching.
        result.half_edges[new_he_h].twin = new_he_h;
        map.half_edges.insert(he_h, new_he_h);

        // Point the edge at this half-edge if not yet set.
        if result.edges[new_edge_h].half_edge == HalfEdgeHandle::default() {
            result.edges[new_edge_h].half_edge = new_he_h;
        }

        new_he_handles.push(new_he_h);
    }

    // Second pass: link next/prev pointers around the loop.
    let n = new_he_handles.len();
    for i in 0..n {
        let next_i = (i + 1) % n;
        let prev_i = (i + n - 1) % n;
        result.half_edges[new_he_handles[i]].next = new_he_handles[next_i];
        result.half_edges[new_he_handles[i]].prev = new_he_handles[prev_i];
    }

    // Point the loop at its first half-edge.
    if let Some(&first) = new_he_handles.first() {
        result.loops[new_loop_h].half_edge = first;
    }

    // Set vertex outgoing pointers.
    for &new_he_h in &new_he_handles {
        let origin = result.half_edges[new_he_h].origin;
        result.vertices[origin].outgoing = new_he_h;
    }

    Ok(())
}

/// Find or create a vertex in `result` that matches the position of `old_v`
/// from `source`, deduplicating by position within tolerance.
fn map_or_create_vertex(
    source: &BRep,
    old_v: VertexHandle,
    result: &mut BRep,
    map: &mut HandleMap,
    tolerance: &Tolerance,
) -> VertexHandle {
    // Already mapped?
    if let Some(&existing) = map.vertices.get(&old_v) {
        return existing;
    }

    let position = source.vertices[old_v].position;

    // Check if any existing result vertex is coincident.
    for (vh, vertex) in &result.vertices {
        if tolerance.points_coincident(&position, &vertex.position) {
            map.vertices.insert(old_v, vh);
            return vh;
        }
    }

    // Create a new vertex.
    let new_vh = result.vertices.insert(Vertex {
        position,
        outgoing: HalfEdgeHandle::default(),
    });
    map.vertices.insert(old_v, new_vh);
    new_vh
}

/// Pair up unpaired half-edges that share the same two vertices in opposite
/// order. After this pass, manifold edges will have proper twin linkage and
/// shared Edge entities.
fn stitch_boundary_edges(result: &mut BRep) {
    // Collect all unpaired half-edges (those whose twin points at themselves).
    let unpaired: Vec<HalfEdgeHandle> = result
        .half_edges
        .keys()
        .filter(|&heh| result.half_edges[heh].twin == heh)
        .collect();

    // Build a map from (origin, dest) -> half-edge handle for fast lookup.
    // Destination is the origin of the next half-edge.
    let mut edge_map: HashMap<(VertexHandle, VertexHandle), HalfEdgeHandle> =
        HashMap::with_capacity(unpaired.len());

    for &heh in &unpaired {
        let origin = result.half_edges[heh].origin;
        let next = result.half_edges[heh].next;
        let dest = result.half_edges[next].origin;
        edge_map.insert((origin, dest), heh);
    }

    // For each unpaired half-edge (A->B), look for one going (B->A).
    for &heh in &unpaired {
        // Skip if already stitched in this pass.
        if result.half_edges[heh].twin != heh {
            continue;
        }

        let origin = result.half_edges[heh].origin;
        let next = result.half_edges[heh].next;
        let dest = result.half_edges[next].origin;

        if let Some(&twin_heh) = edge_map.get(&(dest, origin))
            && twin_heh != heh
            && result.half_edges[twin_heh].twin == twin_heh
        {
            // Pair as twins.
            result.half_edges[heh].twin = twin_heh;
            result.half_edges[twin_heh].twin = heh;

            // Merge edge handles: both half-edges must share the same Edge.
            // Keep the edge from `heh`, discard the one from `twin_heh`.
            let keep_edge = result.half_edges[heh].edge;
            let discard_edge = result.half_edges[twin_heh].edge;

            if keep_edge != discard_edge {
                // Point the twin's half-edge at the kept edge.
                result.half_edges[twin_heh].edge = keep_edge;

                // Update the kept edge to point at `heh` as canonical.
                result.edges[keep_edge].half_edge = heh;

                // If the kept edge has no curve but the discarded one does,
                // take the discarded edge's curve.
                if result.edges[keep_edge].curve.is_none() {
                    let curve = result.edges[discard_edge].curve.clone();
                    result.edges[keep_edge].curve = curve;
                }

                // Mark the discarded edge for removal (we do it in a
                // cleanup pass to avoid invalidating iterators).
                // For now, set its half_edge to default so it's detected
                // as orphaned.
                result.edges[discard_edge].half_edge = HalfEdgeHandle::default();
            }
        }
    }
}

/// Remove edges that are orphaned (their half_edge is default or points to a
/// half-edge that references a different edge).
fn remove_orphan_edges(result: &mut BRep) {
    let orphans: Vec<EdgeHandle> = result
        .edges
        .keys()
        .filter(|&eh| {
            let edge = &result.edges[eh];
            let heh = edge.half_edge;
            // Orphan if half_edge is default.
            if heh == HalfEdgeHandle::default() {
                return true;
            }
            // Orphan if the half-edge it points to references a different edge.
            if let Some(he) = result.half_edges.get(heh) {
                he.edge != eh
            } else {
                true
            }
        })
        .collect();

    for eh in orphans {
        result.edges.remove(eh);
    }
}

/// Flip the surface orientation on a face.
fn flip_surface(result: &mut BRep, face_h: FaceHandle) {
    if let Some(ref mut surf) = result.faces[face_h].surface {
        surf.orientation = !surf.orientation;
    }
}

#[cfg(test)]
#[path = "merge_tests.rs"]
mod merge_tests;
