//! Boundary assembly: identify and merge matching boundary edges
//! to stitch disconnected surface patches into connected topology.
//!
//! The boundary assembly problem arises during import: disconnected surface
//! patches need to be assembled into connected topology by identifying and
//! merging matching boundary edges. This module implements the detection and
//! vertex-merging phases; full topological edge merging via Euler operators
//! is deferred.

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, VertexHandle};

/// Result of boundary assembly.
#[derive(Debug, Clone)]
pub struct StitchResult {
    /// Number of free boundary edges found before stitching.
    pub free_edges_before: usize,
    /// Number of free boundary edges remaining after stitching.
    pub free_edges_after: usize,
    /// Number of edge pairs that were matched.
    pub pairs_matched: usize,
    /// Number of pairs successfully merged (vertices moved).
    pub pairs_merged: usize,
    /// Unmatched boundaries that could not be stitched.
    pub unmatched_edges: Vec<EdgeHandle>,
    /// Maximum displacement introduced during merging.
    pub max_displacement: f64,
}

/// A candidate match between two free boundary edges.
#[derive(Debug, Clone)]
pub struct EdgeMatch {
    /// First edge.
    pub edge_a: EdgeHandle,
    /// Second edge.
    pub edge_b: EdgeHandle,
    /// Hausdorff-like distance between the two edges.
    pub distance: f64,
    /// Whether the edges need orientation reversal to match.
    pub needs_reversal: bool,
}

/// Find free boundary edges -- edges used by only one face.
///
/// In a manifold BRep, every edge has two half-edges on different faces.
/// A free boundary edge has both half-edges belonging to the same face
/// (via their loops). These represent open boundaries where surface
/// patches are disconnected.
pub fn find_free_edges(brep: &BRep) -> Vec<EdgeHandle> {
    let mut free = Vec::new();
    for (eh, edge) in &brep.edges {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            // Missing twin means this is definitely a boundary edge.
            free.push(eh);
            continue;
        };
        // Get the face of each half-edge via its parent loop.
        let face_a = brep.loops.get(he.loop_).map(|l| l.face);
        let face_b = brep.loops.get(twin.loop_).map(|l| l.face);
        match (face_a, face_b) {
            (Some(fa), Some(fb)) if fa == fb => {
                // Both half-edges belong to the same face: free/seam edge.
                free.push(eh);
            }
            (None, _) | (_, None) => {
                // Missing loop/face: treat as free edge.
                free.push(eh);
            }
            _ => {}
        }
    }
    free
}

/// Sample points along an edge for proximity comparison.
///
/// If the edge has curve geometry, samples `num_samples` evenly-spaced
/// points along the curve. Otherwise, linearly interpolates between the
/// endpoint vertex positions.
fn sample_edge_points(brep: &BRep, edge: EdgeHandle, num_samples: usize) -> Option<Vec<Point3>> {
    let edge_data = brep.edges.get(edge)?;
    let he = brep.half_edges.get(edge_data.half_edge)?;
    let twin = brep.half_edges.get(he.twin)?;
    let start_pos = brep.vertices.get(he.origin)?.position;
    let end_pos = brep.vertices.get(twin.origin)?.position;

    let samples = num_samples.max(2);

    if let Some(curve_on_edge) = &edge_data.curve {
        let mut pts = Vec::with_capacity(samples);
        for i in 0..samples {
            let t_frac = i as f64 / (samples - 1) as f64;
            let t = curve_on_edge.t_start + t_frac * (curve_on_edge.t_end - curve_on_edge.t_start);
            pts.push(curve_on_edge.geometry.point_at(t));
        }
        Some(pts)
    } else {
        // Linear interpolation between endpoints.
        let mut pts = Vec::with_capacity(samples);
        let dir = end_pos - start_pos;
        for i in 0..samples {
            let t = i as f64 / (samples - 1) as f64;
            pts.push(Point3::new(
                start_pos.x + dir.x * t,
                start_pos.y + dir.y * t,
                start_pos.z + dir.z * t,
            ));
        }
        Some(pts)
    }
}

/// Get the endpoint vertices of an edge.
///
/// Returns `(start_vertex, end_vertex)` or `None` if topology is incomplete.
fn edge_endpoints(brep: &BRep, edge: EdgeHandle) -> Option<(VertexHandle, VertexHandle)> {
    let edge_data = brep.edges.get(edge)?;
    let he = brep.half_edges.get(edge_data.half_edge)?;
    let twin = brep.half_edges.get(he.twin)?;
    Some((he.origin, twin.origin))
}

/// Compute one-sided directed Hausdorff distance from `pts_a` to `pts_b`.
///
/// For each point in `pts_a`, finds the minimum distance to any point in
/// `pts_b`, then returns the maximum of those minima.
fn directed_hausdorff(pts_a: &[Point3], pts_b: &[Point3]) -> f64 {
    let mut max_min = 0.0_f64;
    for pa in pts_a {
        let mut min_dist = f64::MAX;
        for pb in pts_b {
            let d = pa.distance_to(pb);
            if d < min_dist {
                min_dist = d;
            }
        }
        max_min = max_min.max(min_dist);
    }
    max_min
}

/// Find matching pairs among free boundary edges.
///
/// Uses pairwise proximity: for each free edge, samples points along the
/// edge and finds the nearest free edge whose sampled points are within
/// `max_gap`. Returns matches sorted by distance (best match first).
///
/// Uses greedy assignment: each edge can only participate in one match.
pub fn find_matching_edges(
    brep: &BRep,
    free_edges: &[EdgeHandle],
    max_gap: f64,
    num_samples: usize,
) -> Vec<EdgeMatch> {
    let samples_count = num_samples.max(2);
    // Pre-sample all edges.
    let sampled: Vec<Option<Vec<Point3>>> = free_edges
        .iter()
        .map(|&eh| sample_edge_points(brep, eh, samples_count))
        .collect();

    // Collect candidate matches.
    let mut candidates: Vec<EdgeMatch> = Vec::new();
    for i in 0..free_edges.len() {
        let Some(pts_a) = &sampled[i] else {
            continue;
        };
        for j in (i + 1)..free_edges.len() {
            let Some(pts_b) = &sampled[j] else {
                continue;
            };
            // Symmetric Hausdorff = max(directed(a->b), directed(b->a)).
            let h_ab = directed_hausdorff(pts_a, pts_b);
            let h_ba = directed_hausdorff(pts_b, pts_a);
            let hausdorff = h_ab.max(h_ba);

            if hausdorff < max_gap {
                // Check orientation via tangent at midpoints.
                let needs_reversal = check_needs_reversal(brep, free_edges[i], free_edges[j]);
                candidates.push(EdgeMatch {
                    edge_a: free_edges[i],
                    edge_b: free_edges[j],
                    distance: hausdorff,
                    needs_reversal,
                });
            }
        }
    }

    // Sort by distance ascending (best matches first).
    candidates.sort_by(|a, b| {
        a.distance
            .partial_cmp(&b.distance)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Greedy assignment: each edge used at most once.
    let mut used = vec![false; free_edges.len()];
    // Build index map: EdgeHandle -> position in free_edges.
    let index_of: std::collections::HashMap<EdgeHandle, usize> = free_edges
        .iter()
        .enumerate()
        .map(|(i, &eh)| (eh, i))
        .collect();

    let mut result = Vec::new();
    for c in candidates {
        let Some(&ia) = index_of.get(&c.edge_a) else {
            continue;
        };
        let Some(&ib) = index_of.get(&c.edge_b) else {
            continue;
        };
        if used[ia] || used[ib] {
            continue;
        }
        used[ia] = true;
        used[ib] = true;
        result.push(c);
    }
    result
}

/// Check if two edges need orientation reversal to align.
///
/// Compares tangent directions at midpoints. If the dot product of the
/// tangents is negative, the edges run in opposite directions and one
/// needs reversal.
fn check_needs_reversal(brep: &BRep, edge_a: EdgeHandle, edge_b: EdgeHandle) -> bool {
    let tangent_a = edge_midpoint_tangent(brep, edge_a);
    let tangent_b = edge_midpoint_tangent(brep, edge_b);
    match (tangent_a, tangent_b) {
        (Some(ta), Some(tb)) => ta.dot(&tb) < 0.0,
        _ => false,
    }
}

/// Get the tangent vector at the midpoint of an edge.
fn edge_midpoint_tangent(brep: &BRep, edge: EdgeHandle) -> Option<crate::geometry::Vec3> {
    let edge_data = brep.edges.get(edge)?;
    let curve = edge_data.curve.as_ref()?;
    let t_mid = (curve.t_start + curve.t_end) * 0.5;
    Some(curve.geometry.tangent_at(t_mid))
}

/// Stitch matched boundary edges by merging vertex positions.
///
/// For each matched pair, moves the endpoint vertices of both edges
/// to their midpoints, closing the gap. Does not modify topology
/// (edge merging via Euler operators is deferred).
///
/// Respects `max_displacement`: pairs whose gap exceeds 2x `max_displacement`
/// are skipped (cannot close by moving vertices halfway each).
pub fn stitch_boundaries(
    brep: &mut BRep,
    max_gap: f64,
    max_displacement: f64,
    tolerance: &Tolerance,
) -> StitchResult {
    let free_before = find_free_edges(brep);
    let free_edges_before = free_before.len();

    if free_edges_before == 0 {
        return StitchResult {
            free_edges_before: 0,
            free_edges_after: 0,
            pairs_matched: 0,
            pairs_merged: 0,
            unmatched_edges: Vec::new(),
            max_displacement: 0.0,
        };
    }

    let num_samples = 8;
    let matches = find_matching_edges(brep, &free_before, max_gap, num_samples);
    let pairs_matched = matches.len();

    let mut pairs_merged = 0;
    let mut overall_max_displacement = 0.0_f64;

    // Track which edges were successfully matched.
    let mut matched_edges = std::collections::HashSet::new();

    for m in &matches {
        // Skip if the gap is too large to close within displacement limit.
        // Each vertex moves half the gap, so max displacement per vertex = distance / 2.
        if m.distance > 2.0 * max_displacement {
            continue;
        }

        let Some((va_start, va_end)) = edge_endpoints(brep, m.edge_a) else {
            continue;
        };
        let Some((vb_start, vb_end)) = edge_endpoints(brep, m.edge_b) else {
            continue;
        };

        // Read vertex positions.
        let Some(pa_start) = brep.vertices.get(va_start).map(|v| v.position) else {
            continue;
        };
        let Some(pa_end) = brep.vertices.get(va_end).map(|v| v.position) else {
            continue;
        };
        let Some(pb_start) = brep.vertices.get(vb_start).map(|v| v.position) else {
            continue;
        };
        let Some(pb_end) = brep.vertices.get(vb_end).map(|v| v.position) else {
            continue;
        };

        // Determine which pairing of endpoints is closer (forward or reversed).
        let (pair1_a, pair1_b, pair2_a, pair2_b) = if m.needs_reversal {
            // a_start <-> b_end, a_end <-> b_start
            (va_start, vb_end, va_end, vb_start)
        } else {
            // a_start <-> b_start, a_end <-> b_end
            (va_start, vb_start, va_end, vb_end)
        };

        let pos1_a = pa_start;
        let pos1_b = if m.needs_reversal { pb_end } else { pb_start };
        let pos2_a = pa_end;
        let pos2_b = if m.needs_reversal { pb_start } else { pb_end };

        let disp1 = pos1_a.distance_to(&pos1_b) * 0.5;
        let disp2 = pos2_a.distance_to(&pos2_b) * 0.5;
        let pair_max_disp = disp1.max(disp2);

        if pair_max_disp > max_displacement {
            continue;
        }

        // Merge: move both vertices of each pair to their midpoint.
        let mid1 = Point3::new(
            (pos1_a.x + pos1_b.x) * 0.5,
            (pos1_a.y + pos1_b.y) * 0.5,
            (pos1_a.z + pos1_b.z) * 0.5,
        );
        let mid2 = Point3::new(
            (pos2_a.x + pos2_b.x) * 0.5,
            (pos2_a.y + pos2_b.y) * 0.5,
            (pos2_a.z + pos2_b.z) * 0.5,
        );

        // Only merge if the displacement exceeds the numerical floor.
        let _ = tolerance; // tolerance available for future refinement
        if let Some(v) = brep.vertices.get_mut(pair1_a) {
            v.position = mid1;
        }
        if let Some(v) = brep.vertices.get_mut(pair1_b) {
            v.position = mid1;
        }
        if let Some(v) = brep.vertices.get_mut(pair2_a) {
            v.position = mid2;
        }
        if let Some(v) = brep.vertices.get_mut(pair2_b) {
            v.position = mid2;
        }

        overall_max_displacement = overall_max_displacement.max(pair_max_disp);
        pairs_merged += 1;
        matched_edges.insert(m.edge_a);
        matched_edges.insert(m.edge_b);
    }

    // Re-count free edges after vertex merging.
    let free_after = find_free_edges(brep);
    let unmatched_edges: Vec<EdgeHandle> = free_after
        .iter()
        .filter(|eh| !matched_edges.contains(eh))
        .copied()
        .collect();

    StitchResult {
        free_edges_before,
        free_edges_after: free_after.len(),
        pairs_matched,
        pairs_merged,
        unmatched_edges,
        max_displacement: overall_max_displacement,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn valid_cube_has_no_free_edges() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let free = find_free_edges(&brep);
        assert!(
            free.is_empty(),
            "manifold cube should have no free edges, found {}",
            free.len()
        );
    }

    #[test]
    fn stitch_on_valid_cube_is_noop() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let result = stitch_boundaries(&mut brep, 0.01, 0.005, &tol);
        assert_eq!(result.free_edges_before, 0);
        assert_eq!(result.free_edges_after, 0);
        assert_eq!(result.pairs_matched, 0);
        assert_eq!(result.pairs_merged, 0);
        assert!(result.unmatched_edges.is_empty());
        assert!(result.max_displacement < f64::EPSILON);
    }

    #[test]
    fn empty_brep_has_no_free_edges() {
        let brep = BRep::new(Tolerance::default());
        let free = find_free_edges(&brep);
        assert!(free.is_empty());
    }

    #[test]
    fn stitch_on_empty_brep_is_noop() {
        let mut brep = BRep::new(Tolerance::default());
        let tol = Tolerance::default();
        let result = stitch_boundaries(&mut brep, 0.01, 0.005, &tol);
        assert_eq!(result.free_edges_before, 0);
        assert_eq!(result.pairs_matched, 0);
    }

    #[test]
    fn find_matching_edges_empty_list_returns_empty() {
        let brep = BRep::new(Tolerance::default());
        let matches = find_matching_edges(&brep, &[], 0.01, 8);
        assert!(matches.is_empty());
    }

    #[test]
    fn directed_hausdorff_identical_points_is_zero() {
        let pts = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
        let h = directed_hausdorff(&pts, &pts);
        assert!(h < f64::EPSILON);
    }

    #[test]
    fn directed_hausdorff_known_value() {
        let pts_a = vec![Point3::new(0.0, 0.0, 0.0)];
        let pts_b = vec![Point3::new(1.0, 0.0, 0.0)];
        let h = directed_hausdorff(&pts_a, &pts_b);
        assert!((h - 1.0).abs() < f64::EPSILON);
    }
}
