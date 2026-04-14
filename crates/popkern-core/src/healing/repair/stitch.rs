//! Boundary assembly: identify and merge matching boundary edges
//! to stitch disconnected surface patches into connected topology.
//! Implements detection and vertex-merging; full topological edge merging
//! via Euler operators is deferred.

use std::collections::{HashMap, HashSet};

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, VertexHandle};

use super::spatial_grid::SpatialGrid;

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
    /// Edges that had multiple matching candidates (ambiguous).
    pub ambiguous_edges: Vec<EdgeHandle>,
    /// Maximum displacement introduced during merging.
    pub max_displacement: f64,
}

/// Result of [`find_matching_edges`] including ambiguity information.
#[derive(Debug, Clone)]
pub struct MatchResult {
    /// Matched edge pairs after greedy assignment.
    pub result: Vec<EdgeMatch>,
    /// Edges that had more than one candidate match within max_gap.
    pub ambiguous: Vec<EdgeHandle>,
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
/// A free boundary edge has both half-edges belonging to the same face
/// (via their loops), representing open boundaries.
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
pub(super) fn sample_edge_points(
    brep: &BRep,
    edge: EdgeHandle,
    num_samples: usize,
) -> Option<Vec<Point3>> {
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
/// Uses a spatial grid to accelerate proximity queries: for each free
/// edge, computes its midpoint, inserts it into the grid, then queries
/// nearby edges instead of checking all pairs.
///
/// Uses greedy assignment: each edge can only participate in one match.
/// Also reports edges that had multiple candidate matches (ambiguous).
pub fn find_matching_edges(
    brep: &BRep,
    free_edges: &[EdgeHandle],
    max_gap: f64,
    num_samples: usize,
) -> MatchResult {
    let samples_count = num_samples.max(2);
    let sampled: Vec<Option<Vec<Point3>>> = free_edges
        .iter()
        .map(|&eh| sample_edge_points(brep, eh, samples_count))
        .collect();
    let mut grid = SpatialGrid::new(if max_gap > 0.0 { max_gap } else { 1.0 });
    let midpoints: Vec<Option<Point3>> = sampled
        .iter()
        .enumerate()
        .map(|(i, pts_opt)| {
            let pts = pts_opt.as_ref()?;
            if pts.is_empty() {
                return None;
            }
            let n = pts.len() as f64;
            let mid = Point3::new(
                pts.iter().map(|p| p.x).sum::<f64>() / n,
                pts.iter().map(|p| p.y).sum::<f64>() / n,
                pts.iter().map(|p| p.z).sum::<f64>() / n,
            );
            grid.insert(i, &mid);
            Some(mid)
        })
        .collect();
    let mut match_counts = vec![0_usize; free_edges.len()];
    let mut candidates: Vec<EdgeMatch> = Vec::new();
    let mut checked = HashSet::new();
    for i in 0..free_edges.len() {
        let Some(mid) = &midpoints[i] else {
            continue;
        };
        let Some(pts_a) = &sampled[i] else {
            continue;
        };
        for j in grid.query_nearby(mid) {
            if j <= i {
                continue;
            }
            if !checked.insert((i, j)) {
                continue;
            }
            let Some(pts_b) = &sampled[j] else {
                continue;
            };
            let hausdorff = directed_hausdorff(pts_a, pts_b).max(directed_hausdorff(pts_b, pts_a));
            if hausdorff < max_gap {
                match_counts[i] += 1;
                match_counts[j] += 1;
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
    let ambiguous: Vec<EdgeHandle> = free_edges
        .iter()
        .enumerate()
        .filter(|&(i, _)| match_counts[i] > 1)
        .map(|(_, &eh)| eh)
        .collect();
    candidates.sort_by(|a, b| {
        a.distance
            .partial_cmp(&b.distance)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut used = vec![false; free_edges.len()];
    let index_of: HashMap<EdgeHandle, usize> = free_edges
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
    MatchResult { result, ambiguous }
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
            ambiguous_edges: Vec::new(),
            max_displacement: 0.0,
        };
    }

    let num_samples = 8;
    let match_result = find_matching_edges(brep, &free_before, max_gap, num_samples);
    let matches = &match_result.result;
    let ambiguous_edges = match_result.ambiguous.clone();
    let pairs_matched = matches.len();

    let mut pairs_merged = 0;
    let mut overall_max_displacement = 0.0_f64;

    // Track which edges were successfully matched.
    let mut matched_edges = HashSet::new();

    for m in matches {
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

        let (pair1_a, pair1_b, pair2_a, pair2_b) = if m.needs_reversal {
            (va_start, vb_end, va_end, vb_start)
        } else {
            (va_start, vb_start, va_end, vb_end)
        };
        let (p1a, p1b) = (pa_start, if m.needs_reversal { pb_end } else { pb_start });
        let (p2a, p2b) = (pa_end, if m.needs_reversal { pb_start } else { pb_end });
        let pair_max_disp = (p1a.distance_to(&p1b) * 0.5).max(p2a.distance_to(&p2b) * 0.5);
        if pair_max_disp > max_displacement {
            continue;
        }
        let mid1 = Point3::new(
            (p1a.x + p1b.x) * 0.5,
            (p1a.y + p1b.y) * 0.5,
            (p1a.z + p1b.z) * 0.5,
        );
        let mid2 = Point3::new(
            (p2a.x + p2b.x) * 0.5,
            (p2a.y + p2b.y) * 0.5,
            (p2a.z + p2b.z) * 0.5,
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
        ambiguous_edges,
        max_displacement: overall_max_displacement,
    }
}

#[path = "stitch_overlap.rs"]
mod overlap;

pub use overlap::detect_partial_overlap;

#[cfg(test)]
#[path = "stitch_tests.rs"]
mod tests;
