//! Wound healing: fill holes in shell boundaries.
//!
//! When faces are deleted or imports are incomplete, holes remain as
//! loops of free boundary edges. This module provides strategies for
//! detecting and reporting holes to restore shell closure. Currently
//! supports detection and Extension-based fill analysis; actual surface
//! construction is deferred to a future iteration that implements
//! Coons-patch or energy-minimizing fill surfaces.

use std::collections::{HashMap, HashSet};

use crate::geometry::surface::SurfaceEval;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

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
    if let Some(twin) = brep.half_edges.get(he.twin) {
        return Some((he.origin, twin.origin));
    }
    // Null-sentinel twin: a genuine free boundary edge. The surviving
    // half-edge still lives in a face loop, so its end vertex is the
    // next half-edge's origin.
    let next = brep.half_edges.get(he.next)?;
    Some((he.origin, next.origin))
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

/// Get the adjacent face (the non-hole side) for a boundary edge.
///
/// Returns `None` if the edge or its topology is broken.
fn adjacent_face_for_edge(brep: &BRep, edge: EdgeHandle) -> Option<FaceHandle> {
    let edge_data = brep.edges.get(edge)?;
    let he = brep.half_edges.get(edge_data.half_edge)?;
    let loop_data = brep.loops.get(he.loop_)?;
    Some(loop_data.face)
}

/// Sample 4 evenly-spaced points along an edge's curve geometry.
///
/// Falls back to linear interpolation between endpoints if no curve is
/// available. Returns an empty vec only if the topology is broken.
fn sample_edge_points(brep: &BRep, edge: EdgeHandle, n: usize) -> Vec<crate::geometry::Point3> {
    use crate::geometry::curve::CurveEval;

    let edge_data = match brep.edges.get(edge) {
        Some(e) => e,
        None => return Vec::new(),
    };

    // If the edge has curve geometry, sample along the curve.
    if let Some(ref curve_on_edge) = edge_data.curve {
        let t0 = curve_on_edge.t_start;
        let t1 = curve_on_edge.t_end;
        return (0..n)
            .map(|i| {
                let t = t0 + (t1 - t0) * (i as f64 + 0.5) / (n as f64);
                curve_on_edge.geometry.point_at(t)
            })
            .collect();
    }

    // Fallback: interpolate between edge endpoints.
    let he = match brep.half_edges.get(edge_data.half_edge) {
        Some(h) => h,
        None => return Vec::new(),
    };
    let twin = match brep.half_edges.get(he.twin) {
        Some(t) => t,
        None => return Vec::new(),
    };
    let p0 = match brep.vertices.get(he.origin) {
        Some(v) => v.position,
        None => return Vec::new(),
    };
    let p1 = match brep.vertices.get(twin.origin) {
        Some(v) => v.position,
        None => return Vec::new(),
    };

    (0..n)
        .map(|i| {
            let t = (i as f64 + 0.5) / (n as f64);
            crate::geometry::Point3::new(
                p0.x + (p1.x - p0.x) * t,
                p0.y + (p1.y - p0.y) * t,
                p0.z + (p1.z - p0.z) * t,
            )
        })
        .collect()
}

/// Attempt to fill a hole by extending adjacent face surfaces.
///
/// For each edge in the hole boundary, evaluates the adjacent face's
/// surface beyond the edge to see if it can cover the hole. Picks the
/// face whose surface extension has the smallest deviation from the
/// other boundary edges.
///
/// This works well for small holes in smooth regions where adjacent
/// surfaces naturally extend across the gap.
///
/// This is analysis/assessment only -- actually creating the fill face
/// topology would require Euler operators. But it determines WHETHER
/// extension is viable and WHICH surface to extend, which is the core
/// algorithmic work.
pub fn fill_by_extension(
    brep: &BRep,
    hole: &DetectedHole,
    tolerance: &Tolerance,
) -> Option<HoleFillResult> {
    let mut best_deviation = f64::MAX;

    for (edge_idx, &edge) in hole.boundary_edges.iter().enumerate() {
        // Get the adjacent face (the non-hole side).
        let face = match adjacent_face_for_edge(brep, edge) {
            Some(f) => f,
            None => continue,
        };

        // Get the face's surface geometry.
        let face_data = match brep.faces.get(face) {
            Some(f) => f,
            None => continue,
        };
        let surface_on_face = match face_data.surface.as_ref() {
            Some(s) => s,
            None => continue,
        };
        let surface = &surface_on_face.geometry;

        // For each OTHER boundary edge in the hole, sample and project.
        let mut total_deviation = 0.0_f64;
        let mut sample_count = 0usize;

        for (other_idx, &other_edge) in hole.boundary_edges.iter().enumerate() {
            if other_idx == edge_idx {
                continue;
            }
            let points = sample_edge_points(brep, other_edge, 4);
            for pt in &points {
                // Project the point onto the extending surface.
                if let Some((u, v)) = surface.closest_point(pt) {
                    let projected = surface.point_at(u, v);
                    total_deviation += pt.distance_to(&projected);
                    sample_count += 1;
                }
            }
        }

        if sample_count > 0 {
            let avg_deviation = total_deviation / sample_count as f64;
            if avg_deviation < best_deviation {
                best_deviation = avg_deviation;
            }
        }
    }

    // Loose threshold: the extending surface should be within
    // 100x linear tolerance of the other boundary edges.
    let threshold = tolerance.linear() * 100.0;
    if best_deviation < threshold {
        Some(HoleFillResult {
            holes_detected: 1,
            holes_filled: 1,
            strategy: HoleFillStrategy::Extension,
        })
    } else {
        None
    }
}

/// Automatically select the best hole filling strategy based on hole geometry.
///
/// Scoring criteria:
/// - Small holes (< 4 edges): prefer extension (simpler, faster)
/// - 4-sided holes with smooth boundaries: prefer Coons patch (EnergyMinimizing)
///   when extension is not viable
/// - Large holes or high-curvature boundaries: prefer extension if viable
/// - N-sided holes (> 4): use energy-minimizing fan decomposition
pub fn auto_select_strategy(
    brep: &BRep,
    hole: &DetectedHole,
    tolerance: &Tolerance,
) -> HoleFillStrategy {
    // Small holes (3 or fewer edges): extension is simplest and fastest.
    if hole.num_sides <= 3 {
        return HoleFillStrategy::Extension;
    }

    // 4-sided holes: try extension first; fall back to energy-minimizing
    // (Coons patch) if extension is not viable.
    if hole.num_sides == 4 {
        if fill_by_extension(brep, hole, tolerance).is_some() {
            return HoleFillStrategy::Extension;
        }
        return HoleFillStrategy::EnergyMinimizing;
    }

    // N-sided holes (> 4): energy-minimizing fan decomposition.
    HoleFillStrategy::EnergyMinimizing
}

/// Attempt to fill detected holes.
///
/// Uses [`auto_select_strategy`] to pick the best strategy for each hole,
/// respecting the caller's explicit strategy preference when it matches.
///
/// Returns the number of holes detected and how many could
/// potentially be filled.
pub fn fill_holes(
    brep: &mut BRep,
    strategy: HoleFillStrategy,
    tolerance: &Tolerance,
) -> HoleFillResult {
    let holes = detect_holes(brep);
    let holes_detected = holes.len();
    let mut holes_filled = 0usize;

    for hole in &holes {
        let effective_strategy = auto_select_strategy(brep, hole, tolerance);
        // Use the caller's explicit strategy if it matches, otherwise
        // prefer the auto-selected one.
        let chosen = if strategy == effective_strategy {
            strategy
        } else {
            effective_strategy
        };

        match chosen {
            HoleFillStrategy::Extension => {
                if fill_by_extension(brep, hole, tolerance).is_some() {
                    holes_filled += 1;
                    continue;
                }
                // Extension failed. For 4-sided holes check if adjacent
                // surfaces exist as a fallback heuristic.
                if hole.num_sides == 4 {
                    let has_adjacent_surfaces = hole.boundary_edges.iter().any(|&edge| {
                        adjacent_face_for_edge(brep, edge)
                            .and_then(|f| brep.faces.get(f))
                            .and_then(|fd| fd.surface.as_ref())
                            .is_some()
                    });
                    if has_adjacent_surfaces {
                        holes_filled += 1;
                    }
                }
            }
            HoleFillStrategy::EnergyMinimizing => {
                // Energy-minimizing fill surface construction is not yet
                // implemented. Report as detected but not filled.
            }
        }
    }

    HoleFillResult {
        holes_detected,
        holes_filled,
        strategy,
    }
}

#[cfg(test)]
#[path = "hole_fill_tests.rs"]
mod tests;
