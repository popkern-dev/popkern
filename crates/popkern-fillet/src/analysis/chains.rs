//! Chain detection and vertex grouping for fillet edge analysis.
//!
//! Groups edges into G1-continuous chains and identifies vertices
//! where multiple fillet edges converge.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

use super::EdgeAnalysis;

/// Detect G1-continuous edge chains.
///
/// Groups edges into chains where consecutive edges share a vertex and
/// their tangent directions at the shared vertex are within `angular_tolerance`
/// radians. Returns a list of chains, each being a list of indices into
/// the `analyses` slice.
pub fn detect_chains(
    brep: &BRep,
    analyses: &[EdgeAnalysis],
    angular_tolerance: f64,
) -> Vec<Vec<usize>> {
    let n = analyses.len();
    if n == 0 {
        return Vec::new();
    }

    let mut visited = vec![false; n];
    let mut chains: Vec<Vec<usize>> = Vec::new();

    for start_idx in 0..n {
        if visited[start_idx] {
            continue;
        }
        visited[start_idx] = true;
        let mut chain = vec![start_idx];

        // Extend forward from the end vertex.
        extend_chain(
            brep,
            analyses,
            &mut chain,
            &mut visited,
            angular_tolerance,
            false,
        );
        // Extend backward from the start vertex.
        extend_chain(
            brep,
            analyses,
            &mut chain,
            &mut visited,
            angular_tolerance,
            true,
        );

        chains.push(chain);
    }

    chains
}

/// Extend a chain in one direction (forward from end vertex, or backward from start).
fn extend_chain(
    brep: &BRep,
    analyses: &[EdgeAnalysis],
    chain: &mut Vec<usize>,
    visited: &mut [bool],
    angular_tolerance: f64,
    prepend: bool,
) {
    loop {
        let current_idx = if prepend {
            chain[0]
        } else {
            chain[chain.len() - 1]
        };

        let shared_vertex = if prepend {
            analyses[current_idx].vertex_start
        } else {
            analyses[current_idx].vertex_end
        };

        // Find an unvisited neighbor sharing this vertex.
        let neighbor = find_tangent_neighbor(
            brep,
            analyses,
            current_idx,
            shared_vertex,
            visited,
            angular_tolerance,
        );

        match neighbor {
            Some(next_idx) => {
                visited[next_idx] = true;
                if prepend {
                    chain.insert(0, next_idx);
                } else {
                    chain.push(next_idx);
                }
            }
            None => break,
        }
    }
}

/// Find an unvisited edge sharing `vertex` with `current_idx` that is tangent-continuous.
fn find_tangent_neighbor(
    brep: &BRep,
    analyses: &[EdgeAnalysis],
    current_idx: usize,
    vertex: VertexHandle,
    visited: &[bool],
    angular_tolerance: f64,
) -> Option<usize> {
    let current = &analyses[current_idx];

    for (idx, analysis) in analyses.iter().enumerate() {
        if visited[idx] {
            continue;
        }

        // Check if this edge shares the vertex.
        let shares_start = analysis.vertex_start == vertex;
        let shares_end = analysis.vertex_end == vertex;
        if !shares_start && !shares_end {
            continue;
        }

        // Check tangent continuity at the shared vertex.
        if let Some(is_tangent) =
            check_tangent_continuity(brep, current, analysis, vertex, angular_tolerance)
            && is_tangent
        {
            return Some(idx);
        }
    }

    None
}

/// Check if two edges are tangent-continuous at a shared vertex.
///
/// Returns `Some(true)` if the angle between tangents is within tolerance,
/// `Some(false)` if not, or `None` if tangent evaluation fails.
fn check_tangent_continuity(
    brep: &BRep,
    edge_a: &EdgeAnalysis,
    edge_b: &EdgeAnalysis,
    vertex: VertexHandle,
    angular_tolerance: f64,
) -> Option<bool> {
    let tangent_a = edge_tangent_at_vertex(brep, edge_a.edge, vertex)?;
    let tangent_b = edge_tangent_at_vertex(brep, edge_b.edge, vertex)?;

    // One tangent points toward the vertex, the other away.
    // Compare by taking the absolute value of the dot product.
    let dot = tangent_a.dot(&tangent_b).abs();
    let angle = dot.clamp(-1.0, 1.0).acos();

    Some(angle < angular_tolerance)
}

/// Get the tangent vector of an edge at a given vertex.
fn edge_tangent_at_vertex(brep: &BRep, edge: EdgeHandle, vertex: VertexHandle) -> Option<Vec3> {
    let edge_data = brep.edges.get(edge)?;
    let curve = edge_data.curve.as_ref()?;

    let he = brep.half_edges.get(edge_data.half_edge)?;
    let t = if he.origin == vertex {
        curve.t_start
    } else {
        curve.t_end
    };

    let tangent = curve.geometry.tangent_at(t);
    if tangent.length() < 1e-15 {
        return None;
    }
    Some(tangent)
}

/// Vertex → indices of DISTINCT blend edges incident there, built
/// from each edge's `(vertex_start, vertex_end)` pair.
///
/// A closed edge (start == end — a circular hole rim, a boss seam
/// circle) contributes ONE incidence at its closure vertex. The
/// closure vertex is the closed loop's parameterization seam, not a
/// junction between two fillets: counting it twice manufactured a
/// phantom 2-edge corner there, and the smooth vertex blend then
/// declined against a corner that does not exist (PK-296 — every
/// closed fillet in Smooth mode shipped a spurious Smooth→Mitred
/// fallback warning). The closed-fillet surgery owns that vertex;
/// vertex-blend machinery engages only when a SECOND distinct
/// fillet edge genuinely ends there.
pub(crate) fn vertex_incidences(
    ends: impl Iterator<Item = (VertexHandle, VertexHandle)>,
) -> std::collections::BTreeMap<VertexHandle, Vec<usize>> {
    let mut vertex_edges: std::collections::BTreeMap<VertexHandle, Vec<usize>> =
        std::collections::BTreeMap::new();
    for (idx, (start, end)) in ends.enumerate() {
        vertex_edges.entry(start).or_default().push(idx);
        if end != start {
            vertex_edges.entry(end).or_default().push(idx);
        }
    }
    vertex_edges
}

/// Group vertices where multiple fillet edges meet.
///
/// Returns a list of (vertex, vec of edge indices) for each vertex
/// touched by two or more DISTINCT edges in `analyses` (a closed
/// edge counts once at its closure vertex — see
/// [`vertex_incidences`]).
pub fn group_vertices(analyses: &[EdgeAnalysis]) -> Vec<(VertexHandle, Vec<usize>)> {
    let vertex_edges = vertex_incidences(analyses.iter().map(|a| (a.vertex_start, a.vertex_end)));

    // Keep only vertices with 2+ edges.
    let mut groups: Vec<(VertexHandle, Vec<usize>)> = vertex_edges
        .into_iter()
        .filter(|(_, edges)| edges.len() >= 2)
        .collect();

    // Sort for deterministic output.
    groups.sort_by(|a, b| format!("{:?}", a.0).cmp(&format!("{:?}", b.0)));

    groups
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::analyze_edges;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn make_cube() -> BRep {
        BRep::unit_cube(Tolerance::default()).unwrap()
    }

    /// Collect all edge handles from the BRep.
    fn all_edges(brep: &BRep) -> Vec<popkern_core::topology::EdgeHandle> {
        brep.edges.keys().collect()
    }

    #[test]
    fn detect_chains_cube() {
        let cube = make_cube();
        let edges = all_edges(&cube);

        let analyses = analyze_edges(&cube, &edges, false).unwrap();
        // A cube has 12 edges, all meeting at right angles. No tangent chains.
        let chains = detect_chains(&cube, &analyses, 0.1);

        // Each edge should be its own chain (90-degree corners are not tangent-continuous).
        for chain in &chains {
            assert_eq!(chain.len(), 1, "cube edges should not form tangent chains");
        }
        assert_eq!(chains.len(), analyses.len());
    }

    #[test]
    fn group_vertices_cube_corner() {
        let cube = make_cube();
        let edges = all_edges(&cube);

        let analyses = analyze_edges(&cube, &edges, false).unwrap();
        let groups = group_vertices(&analyses);

        // A cube has 8 vertices, each with 3 edges meeting.
        assert_eq!(groups.len(), 8, "cube should have 8 vertex groups");
        for (_, edge_indices) in &groups {
            assert_eq!(
                edge_indices.len(),
                3,
                "each cube vertex should have 3 edges"
            );
        }
    }

    #[test]
    fn group_vertices_empty() {
        let groups = group_vertices(&[]);
        assert!(groups.is_empty());
    }

    /// PK-296: a closed edge (start == end) contributes ONE incidence
    /// at its closure vertex — the closure vertex is a parameterization
    /// seam, not a junction, and must not read as a phantom 2-edge
    /// corner in the vertex-blend dispatch.
    #[test]
    fn vertex_incidences_counts_closed_edge_once() {
        let cube = make_cube();
        let v: Vec<VertexHandle> = cube.vertices.keys().collect();
        // One closed edge (v0→v0) plus one open edge (v0→v1).
        let incidences = vertex_incidences([(v[0], v[0]), (v[0], v[1])].into_iter());
        assert_eq!(
            incidences[&v[0]],
            vec![0, 1],
            "closure vertex: one incidence from the closed edge, one from the open edge"
        );
        assert_eq!(incidences[&v[1]], vec![1]);
        // A lone closed edge yields valence 1 at its closure vertex —
        // below the >= 2 corner threshold everywhere it is consumed.
        let lone = vertex_incidences([(v[0], v[0])].into_iter());
        assert_eq!(lone[&v[0]], vec![0]);
    }

    #[test]
    fn detect_chains_empty() {
        let brep = BRep::new(Tolerance::default());
        let chains = detect_chains(&brep, &[], 0.1);
        assert!(chains.is_empty());
    }
}
