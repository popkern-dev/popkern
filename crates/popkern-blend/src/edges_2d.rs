//! Edge-sequence 2D corner operations — the standard
//! edge fillet / chamfer family for 2D sheet-
//! body wires.
//!
//! All four functions here take a set of edges and apply a
//! corner operation at every **interior corner** of the chain
//! they imply. "Interior corner" means a vertex shared between
//! at least two edges from the set; chain endpoints (vertices
//! touched by only one edge from the set) are skipped.
//!
//! The constant-radius variants ([`fillet_edges_2d`],
//! [`chamfer_edges_2d`]) process corners in sorted vertex-handle
//! order — cheap to compute, no chain walk required.
//!
//! The variable-parameter variants
//! ([`fillet_edges_2d_variable`], [`chamfer_edges_2d_variable`])
//! walk the chain in topological order via
//! [`crate::corner_2d_helpers::walk_chain_corners`] and linearly
//! interpolate the radius / distance between a start and end
//! value. Open chains and closed loops are both accepted; for
//! closed loops the walk start is the smallest-slotmap-key
//! edge, so callers wanting a smooth closed profile should pin
//! `start == end`.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::chamfer_2d::{chamfer_corner_2d, chamfer_corners_2d};
use crate::corner_2d_helpers::topology_failure;
use crate::error::BlendError;
use crate::fillet_2d::{fillet_corner_2d, fillet_corners_2d};

/// Walk an edge-set chain and return its interior-corner
/// vertices **in walk order**. Open chains (two degree-1
/// endpoints) yield `N − 1` corners for `N` edges, skipping
/// both endpoints. Closed loops (no degree-1 vertices) yield
/// `N` corners; the walk starts at the smallest-slotmap-key
/// edge for determinism. Disconnected sets, branching junctions
/// (degree ≥ 3), and empty inputs each return a specific error.
fn walk_chain_corners(brep: &BRep, edges: &[EdgeHandle]) -> Result<Vec<VertexHandle>, String> {
    use std::collections::{HashMap, HashSet};

    if edges.is_empty() {
        return Ok(Vec::new());
    }

    let mut vertex_edges: HashMap<VertexHandle, Vec<EdgeHandle>> = HashMap::new();
    for &eh in edges {
        let he = brep.edges[eh].half_edge;
        let v_from = brep.half_edges[he].origin;
        let v_to = brep.half_edges[brep.half_edges[he].twin].origin;
        vertex_edges.entry(v_from).or_default().push(eh);
        vertex_edges.entry(v_to).or_default().push(eh);
    }

    for (v, es) in &vertex_edges {
        if es.len() > 2 {
            return Err(format!(
                "chain walker: vertex {v:?} is incident to {} edges in the set (expected ≤ 2)",
                es.len()
            ));
        }
    }

    let mut endpoints: Vec<VertexHandle> = vertex_edges
        .iter()
        .filter_map(|(v, es)| (es.len() == 1).then_some(*v))
        .collect();
    endpoints.sort();

    let (start_vertex, start_edge, is_closed) = match endpoints.len() {
        0 => {
            let mut sorted_edges: Vec<EdgeHandle> = edges.to_vec();
            sorted_edges.sort();
            let first_edge = sorted_edges[0];
            let he = brep.edges[first_edge].half_edge;
            let start_v = brep.half_edges[he].origin;
            (start_v, first_edge, true)
        }
        2 => {
            let start_v = endpoints[0];
            let start_e = vertex_edges[&start_v][0];
            (start_v, start_e, false)
        }
        _ => {
            return Err(format!(
                "chain walker: found {} chain endpoints (expected 0 for closed or 2 for open)",
                endpoints.len()
            ));
        }
    };

    let mut corners: Vec<VertexHandle> = Vec::new();
    let mut visited_edges: HashSet<EdgeHandle> = HashSet::new();
    let mut current_vertex = start_vertex;
    let mut current_edge = start_edge;

    if is_closed {
        corners.push(start_vertex);
    }

    loop {
        visited_edges.insert(current_edge);
        let he = brep.edges[current_edge].half_edge;
        let v_from = brep.half_edges[he].origin;
        let v_to = brep.half_edges[brep.half_edges[he].twin].origin;
        let next_vertex = if v_from == current_vertex {
            v_to
        } else {
            v_from
        };

        let incident = &vertex_edges[&next_vertex];
        let next_edge_opt = incident
            .iter()
            .find(|&&e| e != current_edge && !visited_edges.contains(&e))
            .copied();

        match next_edge_opt {
            Some(next_edge) => {
                corners.push(next_vertex);
                current_vertex = next_vertex;
                current_edge = next_edge;
            }
            None => break,
        }
    }

    if visited_edges.len() != edges.len() {
        return Err(format!(
            "chain walker: edge set is disconnected (walked {} of {} edges)",
            visited_edges.len(),
            edges.len()
        ));
    }

    Ok(corners)
}

/// Derive interior-corner vertices from an edge set: vertices
/// touched by ≥ 2 edges in the input. Chain endpoints (single-
/// edge incidents) are skipped. Result is sorted by vertex
/// handle for deterministic collision semantics in the batch
/// fillet / chamfer paths.
fn corner_vertices_for_edge_set(brep: &BRep, edges: &[EdgeHandle]) -> Vec<VertexHandle> {
    use std::collections::HashMap;
    let mut counts: HashMap<VertexHandle, usize> = HashMap::new();
    for &edge_h in edges {
        let he = brep.edges[edge_h].half_edge;
        let v_from = brep.half_edges[he].origin;
        let v_to = brep.half_edges[brep.half_edges[he].twin].origin;
        *counts.entry(v_from).or_insert(0) += 1;
        *counts.entry(v_to).or_insert(0) += 1;
    }
    let mut corners: Vec<VertexHandle> = counts
        .into_iter()
        .filter_map(|(v, c)| (c >= 2).then_some(v))
        .collect();
    corners.sort();
    corners
}

/// Walk the outer loop of a face and collect its edges in loop
/// order. Inner loops (holes) are not visited — callers who need
/// them should drive [`fillet_edges_2d`] directly with a custom
/// edge set.
fn face_outer_loop_edges(brep: &BRep, face: FaceHandle) -> Vec<EdgeHandle> {
    let loop_h = brep.faces[face].outer_loop;
    let start_he = brep.loops[loop_h].half_edge;
    let mut edges = Vec::new();
    edges.push(brep.half_edges[start_he].edge);
    let mut current = brep.half_edges[start_he].next;
    while current != start_he {
        edges.push(brep.half_edges[current].edge);
        current = brep.half_edges[current].next;
    }
    edges
}

fn check_face_handle(brep: &BRep, face: FaceHandle) -> Result<(), BlendError> {
    if brep.faces.get(face).is_none() {
        return Err(BlendError::InvalidParameter {
            parameter: "face".into(),
            value: format!("{face:?}"),
            reason: "face handle not present in BRep".into(),
        });
    }
    Ok(())
}

/// Constant-radius edge-sequence 2D fillet. Fillets every
/// interior corner of the chain implied by `edges` with the
/// same `radius`, in sorted vertex-handle order. Collision
/// semantics match [`fillet_corners_2d`]. An empty slice, a
/// single-edge slice, or a set whose edges don't share a common
/// vertex is a no-op.
///
/// # Errors
///
/// Returns the first per-corner error encountered, wrapped
/// unchanged.
pub fn fillet_edges_2d(brep: &BRep, edges: &[EdgeHandle], radius: f64) -> Result<BRep, BlendError> {
    let corners = corner_vertices_for_edge_set(brep, edges);
    let specs: Vec<(VertexHandle, f64)> = corners.into_iter().map(|v| (v, radius)).collect();
    fillet_corners_2d(brep, &specs)
}

/// Variable-radius edge-sequence fillet — the 2D analog of
/// a 2D profile fillet on a wire. Walks the
/// chain implied by `edges` in topological order and fillets
/// every interior corner with a radius linearly interpolated
/// between `r_start` (at the first corner in walk order) and
/// `r_end` (at the last). Degenerates to [`fillet_edges_2d`]
/// when `r_start == r_end`.
///
/// For closed loops the walk starts at the smallest-slotmap-key
/// edge and the interpolation still applies linearly, so
/// closed-loop callers who want a smooth profile should pin
/// `r_start == r_end`.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] if either radius is
///   non-finite or non-positive.
/// - [`BlendError::PartialFailure`] / `InvalidTopology` if the
///   edge set isn't a single connected chain, or if any
///   per-corner fillet fails.
pub fn fillet_edges_2d_variable(
    brep: &BRep,
    edges: &[EdgeHandle],
    r_start: f64,
    r_end: f64,
) -> Result<BRep, BlendError> {
    for (name, r) in [("r_start", r_start), ("r_end", r_end)] {
        if !r.is_finite() || r <= 0.0 {
            return Err(BlendError::InvalidParameter {
                parameter: name.into(),
                value: format!("{r}"),
                reason: "radius must be a positive finite value".into(),
            });
        }
    }
    let corners_ordered = walk_chain_corners(brep, edges)
        .map_err(|reason| topology_failure(brep, VertexHandle::default(), reason))?;
    let n = corners_ordered.len();
    if n == 0 {
        return Ok(brep.clone());
    }
    let mut current = brep.clone();
    for (i, v) in corners_ordered.into_iter().enumerate() {
        let t = if n == 1 {
            0.0
        } else {
            i as f64 / (n - 1) as f64
        };
        let r = r_start + t * (r_end - r_start);
        current = fillet_corner_2d(&current, v, r)?;
    }
    Ok(current)
}

/// Constant-distance edge-sequence 2D chamfer. Applies a
/// symmetric chamfer of size `dist` (both legs equal) to every
/// interior corner of the chain implied by `edges`. Collision
/// semantics match [`chamfer_corners_2d`]. An empty slice, a
/// single-edge slice, or a disconnected set is a no-op.
pub fn chamfer_edges_2d(brep: &BRep, edges: &[EdgeHandle], dist: f64) -> Result<BRep, BlendError> {
    let corners = corner_vertices_for_edge_set(brep, edges);
    let specs: Vec<(VertexHandle, f64, f64)> =
        corners.into_iter().map(|v| (v, dist, dist)).collect();
    chamfer_corners_2d(brep, &specs)
}

/// Variable-distance edge-sequence chamfer. Walks the chain
/// implied by `edges` in topological order and applies a
/// symmetric chamfer with a distance linearly interpolated
/// between `d_start` (at the first corner in walk order) and
/// `d_end` (at the last). Closed loops walk from the smallest-
/// slotmap-key edge.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] if either distance is
///   non-finite or non-positive.
/// - [`BlendError::PartialFailure`] / `InvalidTopology` if the
///   edge set isn't a single connected chain, or if any
///   per-corner chamfer fails.
pub fn chamfer_edges_2d_variable(
    brep: &BRep,
    edges: &[EdgeHandle],
    d_start: f64,
    d_end: f64,
) -> Result<BRep, BlendError> {
    for (name, d) in [("d_start", d_start), ("d_end", d_end)] {
        if !d.is_finite() || d <= 0.0 {
            return Err(BlendError::InvalidParameter {
                parameter: name.into(),
                value: format!("{d}"),
                reason: "distance must be a positive finite value".into(),
            });
        }
    }
    let corners_ordered = walk_chain_corners(brep, edges)
        .map_err(|reason| topology_failure(brep, VertexHandle::default(), reason))?;
    let n = corners_ordered.len();
    if n == 0 {
        return Ok(brep.clone());
    }
    let mut current = brep.clone();
    for (i, v) in corners_ordered.into_iter().enumerate() {
        let t = if n == 1 {
            0.0
        } else {
            i as f64 / (n - 1) as f64
        };
        let d = d_start + t * (d_end - d_start);
        current = chamfer_corner_2d(&current, v, d, d)?;
    }
    Ok(current)
}

/// Constant-radius fillet of every corner on a face's outer
/// loop. Thin wrapper over [`fillet_edges_2d`] that auto-derives
/// the edge chain from the face handle — useful for the common
/// "fillet the whole sheet" operation.
///
/// On a sheet body either the front or back face yields the
/// same result: the chain walker in [`fillet_edges_2d_variable`]
/// / this path starts at the smallest-slotmap-key edge, not the
/// order in which edges were collected.
pub fn fillet_face_corners_2d(
    brep: &BRep,
    face: FaceHandle,
    radius: f64,
) -> Result<BRep, BlendError> {
    check_face_handle(brep, face)?;
    let edges = face_outer_loop_edges(brep, face);
    fillet_edges_2d(brep, &edges, radius)
}

/// Variable-radius linear profile over every corner on a face's
/// outer loop. Thin wrapper over [`fillet_edges_2d_variable`].
pub fn fillet_face_corners_2d_variable(
    brep: &BRep,
    face: FaceHandle,
    r_start: f64,
    r_end: f64,
) -> Result<BRep, BlendError> {
    check_face_handle(brep, face)?;
    let edges = face_outer_loop_edges(brep, face);
    fillet_edges_2d_variable(brep, &edges, r_start, r_end)
}

/// Constant-distance symmetric chamfer of every corner on a
/// face's outer loop. Thin wrapper over [`chamfer_edges_2d`].
pub fn chamfer_face_corners_2d(
    brep: &BRep,
    face: FaceHandle,
    dist: f64,
) -> Result<BRep, BlendError> {
    check_face_handle(brep, face)?;
    let edges = face_outer_loop_edges(brep, face);
    chamfer_edges_2d(brep, &edges, dist)
}

/// Variable-distance linear profile over every corner on a
/// face's outer loop. Thin wrapper over
/// [`chamfer_edges_2d_variable`].
pub fn chamfer_face_corners_2d_variable(
    brep: &BRep,
    face: FaceHandle,
    d_start: f64,
    d_end: f64,
) -> Result<BRep, BlendError> {
    check_face_handle(brep, face)?;
    let edges = face_outer_loop_edges(brep, face);
    chamfer_edges_2d_variable(brep, &edges, d_start, d_end)
}

#[cfg(test)]
#[path = "edges_2d_tests.rs"]
mod tests;
