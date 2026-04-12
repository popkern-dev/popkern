//! Wound healing: fill holes in shell boundaries.
//!
//! When faces are deleted or imports are incomplete, holes remain as
//! loops of free boundary edges. This module provides strategies for
//! detecting and reporting holes to restore shell closure. Currently
//! supports detection and Extension-based fill analysis; actual surface
//! construction is deferred to a future iteration that implements
//! Coons-patch or energy-minimizing fill surfaces.

use std::collections::{HashMap, HashSet};

use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, VertexHandle};

/// Strategy for filling a hole.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HoleFillStrategy {
    /// Extend adjacent surfaces across the hole.
    Extension,
    /// Construct a fill surface minimizing bending energy (future).
    EnergyMinimizing,
}

/// Result of a hole fill operation.
#[derive(Debug, Clone)]
pub struct HoleFillResult {
    /// Number of holes detected.
    pub holes_detected: usize,
    /// Number of holes filled.
    pub holes_filled: usize,
    /// Strategy used.
    pub strategy: HoleFillStrategy,
}

/// A detected hole (loop of free boundary edges).
#[derive(Debug, Clone)]
pub struct DetectedHole {
    /// Ordered list of edges forming the hole boundary.
    pub boundary_edges: Vec<EdgeHandle>,
    /// Number of sides (edges) in the hole.
    pub num_sides: usize,
}

/// Identify free boundary edges in the BRep.
///
/// A free edge is one where both half-edges belong to the same face
/// (the edge borders only one face). Returns the set of such edges.
fn find_free_edges(brep: &BRep) -> HashSet<EdgeHandle> {
    let mut free = HashSet::new();

    for (edge_handle, edge_data) in &brep.edges {
        let he = match brep.half_edges.get(edge_data.half_edge) {
            Some(h) => h,
            None => continue,
        };
        let twin = match brep.half_edges.get(he.twin) {
            Some(t) => t,
            None => {
                // Missing twin means this is a boundary edge.
                free.insert(edge_handle);
                continue;
            }
        };

        // Both half-edges in the same loop means the edge is free.
        let loop_a = he.loop_;
        let loop_b = twin.loop_;

        let face_a = brep.loops.get(loop_a).map(|l| l.face);
        let face_b = brep.loops.get(loop_b).map(|l| l.face);

        if face_a == face_b {
            free.insert(edge_handle);
        }
    }

    free
}

/// Get the two endpoint vertex handles for an edge.
///
/// Returns `None` if the topology is broken.
fn edge_endpoints(brep: &BRep, edge: EdgeHandle) -> Option<(VertexHandle, VertexHandle)> {
    let edge_data = brep.edges.get(edge)?;
    let he = brep.half_edges.get(edge_data.half_edge)?;
    let twin = brep.half_edges.get(he.twin)?;
    Some((he.origin, twin.origin))
}

/// Build adjacency: map each vertex to the free edges incident on it.
fn build_vertex_edge_adjacency(
    brep: &BRep,
    free_edges: &HashSet<EdgeHandle>,
) -> HashMap<VertexHandle, Vec<EdgeHandle>> {
    let mut adjacency: HashMap<VertexHandle, Vec<EdgeHandle>> = HashMap::new();

    for &edge in free_edges {
        if let Some((v_start, v_end)) = edge_endpoints(brep, edge) {
            adjacency.entry(v_start).or_default().push(edge);
            adjacency.entry(v_end).or_default().push(edge);
        }
    }

    adjacency
}

/// Trace a chain of free edges starting from a given edge, following
/// vertex adjacency until the chain closes (returns to the starting
/// vertex) or dead-ends.
fn trace_hole_chain(
    brep: &BRep,
    start_edge: EdgeHandle,
    adjacency: &HashMap<VertexHandle, Vec<EdgeHandle>>,
    visited: &mut HashSet<EdgeHandle>,
) -> Option<Vec<EdgeHandle>> {
    let (v_start, v_next) = edge_endpoints(brep, start_edge)?;

    let mut chain = vec![start_edge];
    visited.insert(start_edge);
    let mut current_vertex = v_next;

    // Walk along free edges until we close the loop or get stuck.
    // Guard against infinite loops with a bound on chain length.
    let max_steps = adjacency.len().max(1);
    for _ in 0..max_steps {
        if current_vertex == v_start && chain.len() >= 3 {
            // Closed loop found.
            return Some(chain);
        }

        // Find the next unvisited free edge at this vertex.
        let candidates = adjacency.get(&current_vertex)?;
        let &next_edge = candidates.iter().find(|&&e| !visited.contains(&e))?;

        visited.insert(next_edge);
        chain.push(next_edge);

        // Advance to the other endpoint of next_edge.
        let (ea, eb) = edge_endpoints(brep, next_edge)?;
        current_vertex = if ea == current_vertex { eb } else { ea };
    }

    None // Chain exceeded maximum length without closing.
}

/// Detect holes in the BRep shell -- loops of free boundary edges.
///
/// A hole is a connected loop of edges where each edge borders only
/// one face. Returns detected holes ordered by size (smallest first).
pub fn detect_holes(brep: &BRep) -> Vec<DetectedHole> {
    let free_edges = find_free_edges(brep);
    if free_edges.is_empty() {
        return Vec::new();
    }

    let adjacency = build_vertex_edge_adjacency(brep, &free_edges);
    let mut visited = HashSet::new();
    let mut holes = Vec::new();

    // Try to trace a closed chain starting from each unvisited free edge.
    // Collect edge handles first to avoid borrow conflicts.
    let edge_list: Vec<EdgeHandle> = free_edges.iter().copied().collect();
    for edge in edge_list {
        if visited.contains(&edge) {
            continue;
        }
        if let Some(chain) = trace_hole_chain(brep, edge, &adjacency, &mut visited) {
            let num_sides = chain.len();
            holes.push(DetectedHole {
                boundary_edges: chain,
                num_sides,
            });
        }
    }

    // Sort by number of sides (smallest holes first -- easier to fill).
    holes.sort_by_key(|h| h.num_sides);
    holes
}

/// Attempt to fill detected holes.
///
/// Currently uses surface extension strategy: for each hole, finds
/// the adjacent face surfaces and evaluates whether extending them
/// could close the hole. This is a detection + reporting operation;
/// actual surface construction (Coons patch, energy-minimizing fill)
/// requires significant infrastructure and is deferred to a future
/// iteration.
///
/// Returns the number of holes detected and how many could
/// potentially be filled.
pub fn fill_holes(
    brep: &mut BRep,
    strategy: HoleFillStrategy,
    _tolerance: &Tolerance,
) -> HoleFillResult {
    let holes = detect_holes(brep);
    let holes_detected = holes.len();
    let mut holes_filled = 0usize;

    match strategy {
        HoleFillStrategy::Extension => {
            // For each hole, check if adjacent faces have surface geometry
            // that could plausibly be extended. Count fillable holes.
            for hole in &holes {
                let has_adjacent_surfaces = hole.boundary_edges.iter().any(|&edge| {
                    let edge_data = match brep.edges.get(edge) {
                        Some(e) => e,
                        None => return false,
                    };
                    let he = match brep.half_edges.get(edge_data.half_edge) {
                        Some(h) => h,
                        None => return false,
                    };
                    let loop_data = match brep.loops.get(he.loop_) {
                        Some(l) => l,
                        None => return false,
                    };
                    let face_data = match brep.faces.get(loop_data.face) {
                        Some(f) => f,
                        None => return false,
                    };
                    face_data.surface.is_some()
                });

                if has_adjacent_surfaces {
                    // Hole is fillable in principle -- actual fill is deferred.
                    holes_filled += 1;
                }
            }
        }
        HoleFillStrategy::EnergyMinimizing => {
            // Energy-minimizing fill is not yet implemented.
            // Report detection only.
        }
    }

    HoleFillResult {
        holes_detected,
        holes_filled,
        strategy,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn valid_cube_has_no_holes() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let holes = detect_holes(&brep);
        assert!(holes.is_empty(), "valid cube should have no boundary holes");
    }

    #[test]
    fn valid_cube_has_no_free_edges() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let free = find_free_edges(&brep);
        assert!(
            free.is_empty(),
            "valid cube should have no free boundary edges"
        );
    }

    #[test]
    fn fill_holes_on_valid_cube_is_noop() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let result = fill_holes(&mut brep, HoleFillStrategy::Extension, &tol);
        assert_eq!(result.holes_detected, 0);
        assert_eq!(result.holes_filled, 0);
        assert_eq!(result.strategy, HoleFillStrategy::Extension);
    }

    #[test]
    fn empty_brep_has_no_holes() {
        let brep = BRep::new(Tolerance::default());
        let holes = detect_holes(&brep);
        assert!(holes.is_empty());
    }

    #[test]
    fn energy_minimizing_reports_detection_only() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let result = fill_holes(&mut brep, HoleFillStrategy::EnergyMinimizing, &tol);
        assert_eq!(result.holes_detected, 0);
        assert_eq!(result.holes_filled, 0);
        assert_eq!(result.strategy, HoleFillStrategy::EnergyMinimizing);
    }
}
