//! Edge-level continuity classification on a BRep.
//!
//! Classifies each edge by the achieved cross-boundary continuity level
//! between its two adjacent faces, and detects maximal chains of smooth
//! edges sharing a vertex.

use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceEval;
use crate::healing::continuity::measurement::{BoundarySample, cross_boundary_continuity};
use crate::healing::quality::ContinuityLevel;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle};

/// Number of evenly-spaced samples along an edge for continuity measurement.
const EDGE_SAMPLE_COUNT: usize = 8;

/// Continuity classification for a single edge.
#[derive(Debug, Clone, Copy)]
pub struct EdgeClassification {
    /// The classified edge.
    pub edge: EdgeHandle,
    /// Achieved continuity level across this edge.
    pub level: ContinuityLevel,
}

/// A maximal chain of connected edges that all meet a minimum continuity level.
#[derive(Debug, Clone)]
pub struct SmoothChain {
    /// Edges in the chain.
    pub edges: Vec<EdgeHandle>,
    /// Minimum continuity level achieved by every edge in the chain.
    pub min_level: ContinuityLevel,
}

/// Classify the cross-boundary continuity of a single edge.
///
/// Determines the two faces adjacent to the edge via its half-edges,
/// samples the shared boundary curve, projects each 3D sample onto both
/// face surfaces, and measures the achieved continuity level.
///
/// Returns `None` if the edge is non-manifold (missing twin, missing face,
/// or either face lacks surface geometry).
pub fn classify_edge_smoothness(
    brep: &BRep,
    edge: EdgeHandle,
    tolerance: &Tolerance,
) -> Option<EdgeClassification> {
    let edge_data = brep.edges.get(edge)?;
    let he_a = edge_data.half_edge;
    let he_a_data = brep.half_edges.get(he_a)?;

    // Get the twin half-edge and its face.
    let he_b = he_a_data.twin;
    let he_b_data = brep.half_edges.get(he_b)?;

    // Retrieve the two faces.
    let loop_a = brep.loops.get(he_a_data.loop_)?;
    let loop_b = brep.loops.get(he_b_data.loop_)?;
    let face_a = brep.faces.get(loop_a.face)?;
    let face_b = brep.faces.get(loop_b.face)?;

    // Both faces must have surface geometry.
    let surf_a = face_a.surface.as_ref()?;
    let surf_b = face_b.surface.as_ref()?;

    // Sample 3D points along the shared edge boundary.
    let boundary_points = sample_edge_points(brep, edge, EDGE_SAMPLE_COUNT);
    if boundary_points.is_empty() {
        return None;
    }

    // Project each 3D point onto both surfaces to get UV pairs.
    let mut samples = Vec::with_capacity(boundary_points.len());
    for pt in &boundary_points {
        let uv_a = surf_a.geometry.closest_point(pt)?;
        let uv_b = surf_b.geometry.closest_point(pt)?;
        samples.push(BoundarySample {
            u_a: uv_a.0,
            v_a: uv_a.1,
            u_b: uv_b.0,
            v_b: uv_b.1,
        });
    }

    let measurement =
        cross_boundary_continuity(&surf_a.geometry, &surf_b.geometry, &samples, tolerance);

    Some(EdgeClassification {
        edge,
        level: measurement.achieved_level,
    })
}

/// Sample 3D points along an edge's curve.
///
/// If the edge has curve geometry, samples at `count` evenly-spaced parameter
/// values. Otherwise falls back to the endpoint vertex positions.
fn sample_edge_points(brep: &BRep, edge: EdgeHandle, count: usize) -> Vec<crate::geometry::Point3> {
    let Some(edge_data) = brep.edges.get(edge) else {
        return Vec::new();
    };

    if let Some(ref curve) = edge_data.curve {
        let n = count.max(2);
        let (t0, t1) = (curve.t_start, curve.t_end);
        (0..n)
            .map(|i| {
                let t = t0 + (t1 - t0) * (i as f64) / ((n - 1) as f64);
                curve.geometry.point_at(t)
            })
            .collect()
    } else {
        // No curve geometry — use the two endpoint vertex positions.
        let Some(he) = brep.half_edges.get(edge_data.half_edge) else {
            return Vec::new();
        };
        let mut pts = Vec::new();
        if let Some(v) = brep.vertices.get(he.origin) {
            pts.push(v.position);
        }
        let twin = he.twin;
        if let Some(twin_he) = brep.half_edges.get(twin)
            && let Some(v) = brep.vertices.get(twin_he.origin)
        {
            pts.push(v.position);
        }
        pts
    }
}

/// Classify the cross-boundary continuity of every edge in a BRep.
///
/// Iterates all edges and calls [`classify_edge_smoothness`] for each.
/// Edges that cannot be classified (non-manifold, missing geometry) are
/// silently skipped.
pub fn classify_all_edges(brep: &BRep, tolerance: &Tolerance) -> Vec<EdgeClassification> {
    brep.edges
        .keys()
        .filter_map(|eh| classify_edge_smoothness(brep, eh, tolerance))
        .collect()
}

/// Detect maximal chains of connected edges that all achieve at least `min_level`.
///
/// Two edges are considered adjacent if they share a vertex. The function
/// uses a BFS flood-fill to find connected components among qualifying edges,
/// then returns chains sorted longest-first.
pub fn detect_smooth_chains(
    brep: &BRep,
    tolerance: &Tolerance,
    min_level: ContinuityLevel,
) -> Vec<SmoothChain> {
    let classifications = classify_all_edges(brep, tolerance);

    // Filter to edges that meet the minimum level.
    let qualifying: Vec<EdgeHandle> = classifications
        .iter()
        .filter(|c| c.level >= min_level)
        .map(|c| c.edge)
        .collect();

    if qualifying.is_empty() {
        return Vec::new();
    }

    // Build adjacency: collect the vertices for each qualifying edge,
    // then two edges are adjacent if they share a vertex.
    let edge_vertices: Vec<(EdgeHandle, Vec<crate::topology::VertexHandle>)> = qualifying
        .iter()
        .filter_map(|&eh| {
            let edge_data = brep.edges.get(eh)?;
            let he = brep.half_edges.get(edge_data.half_edge)?;
            let twin = brep.half_edges.get(he.twin)?;
            Some((eh, vec![he.origin, twin.origin]))
        })
        .collect();

    // Build vertex -> edge indices map.
    let mut vertex_to_edges: std::collections::HashMap<crate::topology::VertexHandle, Vec<usize>> =
        std::collections::HashMap::new();
    for (i, (_, verts)) in edge_vertices.iter().enumerate() {
        for &vh in verts {
            vertex_to_edges.entry(vh).or_default().push(i);
        }
    }

    // Union-find (disjoint set) for connected components.
    let n = edge_vertices.len();
    let mut parent: Vec<usize> = (0..n).collect();
    let mut rank: Vec<usize> = vec![0; n];

    let find = |parent: &mut Vec<usize>, mut x: usize| -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]]; // path compression
            x = parent[x];
        }
        x
    };

    let union = |parent: &mut Vec<usize>, rank: &mut Vec<usize>, a: usize, b: usize| {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra == rb {
            return;
        }
        if rank[ra] < rank[rb] {
            parent[ra] = rb;
        } else if rank[ra] > rank[rb] {
            parent[rb] = ra;
        } else {
            parent[rb] = ra;
            rank[ra] += 1;
        }
    };

    // Union edges that share a vertex.
    for edges_at_vertex in vertex_to_edges.values() {
        for window in edges_at_vertex.windows(2) {
            union(&mut parent, &mut rank, window[0], window[1]);
        }
    }

    // Collect components.
    let mut components: std::collections::HashMap<usize, Vec<EdgeHandle>> =
        std::collections::HashMap::new();
    for (i, (eh, _)) in edge_vertices.iter().enumerate() {
        let root = find(&mut parent, i);
        components.entry(root).or_default().push(*eh);
    }

    // Build SmoothChain for each component.
    let mut chains: Vec<SmoothChain> = components
        .into_values()
        .map(|edges| SmoothChain { edges, min_level })
        .collect();

    // Sort longest-first.
    chains.sort_by_key(|c| std::cmp::Reverse(c.edges.len()));
    chains
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    /// On a unit cube, all edges should be G0 but not G1 (sharp 90-degree edges).
    #[test]
    fn cube_edges_are_g0_not_g1() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let classifications = classify_all_edges(&brep, &tol);

        // A cube has 12 edges; all should be classifiable.
        assert!(
            !classifications.is_empty(),
            "should classify at least some edges"
        );

        for c in &classifications {
            // Each edge borders two perpendicular faces (90-degree angle),
            // which is well above the default angular tolerance (~1e-4 rad).
            // So all edges should be G0 (they share boundary positionally)
            // but NOT G1.
            assert_eq!(
                c.level,
                ContinuityLevel::G0,
                "cube edge should be G0 (sharp), got {:?}",
                c.level
            );
        }
    }

    /// No G1 smooth chains should exist on a cube.
    #[test]
    fn cube_has_no_g1_chains() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let chains = detect_smooth_chains(&brep, &tol, ContinuityLevel::G1);
        assert!(
            chains.is_empty(),
            "cube should have no G1 smooth chains, found {}",
            chains.len()
        );
    }

    /// All cube edges form G0 smooth chains.
    #[test]
    fn cube_g0_chains_cover_all_edges() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let chains = detect_smooth_chains(&brep, &tol, ContinuityLevel::G0);

        let total_edges: usize = chains.iter().map(|c| c.edges.len()).sum();
        // All 12 cube edges should appear in G0 chains.
        assert_eq!(
            total_edges, 12,
            "expected 12 edges in G0 chains, got {total_edges}"
        );
    }

    #[test]
    fn empty_brep_produces_no_classifications() {
        let brep = BRep::new(Tolerance::default());
        let tol = Tolerance::default();
        let classifications = classify_all_edges(&brep, &tol);
        assert!(classifications.is_empty());
    }

    #[test]
    fn empty_brep_produces_no_chains() {
        let brep = BRep::new(Tolerance::default());
        let tol = Tolerance::default();
        let chains = detect_smooth_chains(&brep, &tol, ContinuityLevel::G0);
        assert!(chains.is_empty());
    }
}
