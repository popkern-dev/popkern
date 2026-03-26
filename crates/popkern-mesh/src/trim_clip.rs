//! Boundary triangle clipping and T-junction prevention.
//!
//! Clips BOUNDARY triangles against trim polyline segments with vertex snapping.
//! Propagates split vertices into neighboring INSIDE triangles to prevent T-junctions.

use std::collections::HashMap;

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use crate::mesh::MeshVertex;
use crate::trim::{TriClass, TrimLoop, classify_triangle, point_inside_trim_loops};

/// Snap tolerance for merging intersection points with existing triangle vertices.
const SNAP_TOL: f64 = 1e-10;

/// Apply trim clipping to a mesh, returning new vertices and triangles.
///
/// Phase 1: classifies each triangle as INSIDE, OUTSIDE, or BOUNDARY.
/// INSIDE triangles are kept, OUTSIDE are discarded, and BOUNDARY triangles
/// are clipped against the trim polyline segments. Phase 2 propagates edge
/// splits into neighboring INSIDE triangles to prevent T-junctions.
pub(crate) fn apply_trim(
    vertices: &[MeshVertex],
    triangles: &[[usize; 3]],
    trim_loops: &[TrimLoop],
    surface: &SurfaceGeometry,
    _orientation: bool,
) -> (Vec<MeshVertex>, Vec<[usize; 3]>) {
    let mut new_verts: Vec<MeshVertex> = vertices.to_vec();
    let mut new_tris: Vec<[usize; 3]> = Vec::new();
    // Map from canonical edge (min, max) to inserted split vertex indices.
    let mut edge_splits: HashMap<(usize, usize), Vec<usize>> = HashMap::new();

    // Phase 1: classify and process each triangle.
    for tri in triangles {
        let tri_uv = [
            vertices[tri[0]].uv,
            vertices[tri[1]].uv,
            vertices[tri[2]].uv,
        ];
        let class = classify_triangle(&tri_uv, trim_loops);

        match class {
            TriClass::Inside => {
                new_tris.push(*tri);
            }
            TriClass::Outside => {
                // Discard.
            }
            TriClass::Boundary => {
                clip_boundary_triangle(
                    tri,
                    vertices,
                    trim_loops,
                    surface,
                    &mut new_verts,
                    &mut new_tris,
                    &mut edge_splits,
                );
            }
        }
    }

    // Phase 2: T-junction prevention.
    // Re-triangulate INSIDE triangles that share split edges.
    let inside_tris: Vec<[usize; 3]> = new_tris
        .iter()
        .copied()
        .filter(|t| {
            let uv = [new_verts[t[0]].uv, new_verts[t[1]].uv, new_verts[t[2]].uv];
            classify_triangle(&uv, trim_loops) == TriClass::Inside
        })
        .collect();

    let mut fixed_tris: Vec<[usize; 3]> = Vec::new();
    // Collect non-inside tris (boundary clipped output) as-is.
    for t in &new_tris {
        let uv = [new_verts[t[0]].uv, new_verts[t[1]].uv, new_verts[t[2]].uv];
        if classify_triangle(&uv, trim_loops) != TriClass::Inside {
            fixed_tris.push(*t);
        }
    }

    for tri in &inside_tris {
        let splits = collect_edge_splits(tri, &edge_splits);
        if splits.iter().all(|(_, s)| s.is_empty()) {
            fixed_tris.push(*tri);
        } else {
            split_triangle_at_edge_vertices(tri, &splits, &mut fixed_tris);
        }
    }

    (new_verts, fixed_tris)
}

/// Clip a single BOUNDARY triangle against all trim polyline segments.
///
/// Finds intersection points between triangle edges and trim segments,
/// snaps near-vertex intersections, creates new mesh vertices, and
/// fan-triangulates the resulting polygon, keeping only sub-triangles
/// whose centroids are inside the trim region.
fn clip_boundary_triangle(
    tri: &[usize; 3],
    vertices: &[MeshVertex],
    trim_loops: &[TrimLoop],
    surface: &SurfaceGeometry,
    new_verts: &mut Vec<MeshVertex>,
    new_tris: &mut Vec<[usize; 3]>,
    edge_splits: &mut HashMap<(usize, usize), Vec<usize>>,
) {
    let tri_uv = [
        vertices[tri[0]].uv,
        vertices[tri[1]].uv,
        vertices[tri[2]].uv,
    ];

    // Collect all intersection points per edge.
    // Per edge: Vec<(t, uv_point)>
    let mut edge_intersections: [Vec<(f64, (f64, f64))>; 3] = [Vec::new(), Vec::new(), Vec::new()];

    for tl in trim_loops {
        let n = tl.points.len();
        if n < 2 {
            continue;
        }
        for i in 0..n {
            let j = (i + 1) % n;
            let hits = clip_triangle_against_segment(&tri_uv, tl.points[i], tl.points[j]);
            for (edge_idx, t, uv_pt) in hits {
                edge_intersections[edge_idx].push((t, uv_pt));
            }
        }
    }

    // Sort intersections per edge by parameter t.
    for edge in &mut edge_intersections {
        edge.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }

    // Deduplicate very close intersections.
    for edge in &mut edge_intersections {
        edge.dedup_by(|a, b| (a.0 - b.0).abs() < SNAP_TOL);
    }

    // Map each intersection to a vertex index (snap or create).
    let edge_indices = [(0usize, 1usize), (1, 2), (2, 0)];
    let mut edge_new_verts: [Vec<usize>; 3] = [Vec::new(), Vec::new(), Vec::new()];

    for (ei, intersections) in edge_intersections.iter().enumerate() {
        let (vi, vj) = edge_indices[ei];
        let v_start = tri[vi];
        let v_end = tri[vj];
        let uv_start = vertices[v_start].uv;
        let uv_end = vertices[v_end].uv;

        for &(t, uv_pt) in intersections {
            // Check snap to start vertex.
            let dist_start =
                ((uv_pt.0 - uv_start.0).powi(2) + (uv_pt.1 - uv_start.1).powi(2)).sqrt();
            if dist_start < SNAP_TOL {
                edge_new_verts[ei].push(v_start);
                continue;
            }

            // Check snap to end vertex.
            let dist_end = ((uv_pt.0 - uv_end.0).powi(2) + (uv_pt.1 - uv_end.1).powi(2)).sqrt();
            if dist_end < SNAP_TOL {
                edge_new_verts[ei].push(v_end);
                continue;
            }

            // Check snap to parameter extremes.
            if t < SNAP_TOL {
                edge_new_verts[ei].push(v_start);
                continue;
            }
            if t > 1.0 - SNAP_TOL {
                edge_new_verts[ei].push(v_end);
                continue;
            }

            // Create new vertex by evaluating the surface at the UV point.
            let pos = surface.point_at(uv_pt.0, uv_pt.1);
            let normal = surface.normal_at(uv_pt.0, uv_pt.1);
            let new_idx = new_verts.len();
            new_verts.push(MeshVertex {
                position: pos,
                normal,
                uv: uv_pt,
            });
            edge_new_verts[ei].push(new_idx);

            // Record in edge_splits using canonical (min, max) key.
            let key = (v_start.min(v_end), v_start.max(v_end));
            edge_splits.entry(key).or_default().push(new_idx);
        }
    }

    // Build ordered vertex ring: v0 -> [edge0 splits] -> v1 -> [edge1 splits] -> v2 -> [edge2 splits]
    let mut ring: Vec<usize> = Vec::new();
    for ei in 0..3 {
        ring.push(tri[edge_indices[ei].0]);
        ring.extend_from_slice(&edge_new_verts[ei]);
    }

    // If ring has fewer than 3 unique vertices, skip.
    if ring.len() < 3 {
        return;
    }

    // Fan-triangulate the ring from the first vertex, keep sub-triangles
    // whose centroids are inside the trim region.
    for i in 1..ring.len() - 1 {
        let a = ring[0];
        let b = ring[i];
        let c = ring[i + 1];

        // Skip degenerate triangles.
        if a == b || b == c || a == c {
            continue;
        }

        let uv_a = new_verts[a].uv;
        let uv_b = new_verts[b].uv;
        let uv_c = new_verts[c].uv;
        let cx = (uv_a.0 + uv_b.0 + uv_c.0) / 3.0;
        let cy = (uv_a.1 + uv_b.1 + uv_c.1) / 3.0;

        if point_inside_trim_loops(cx, cy, trim_loops) {
            new_tris.push([a, b, c]);
        }
    }
}

/// Find intersection points between a triangle's edges and a polyline segment.
///
/// For each of the 3 triangle edges, computes the intersection with segment
/// `(seg_a, seg_b)`. Returns a list of `(edge_index, t_along_edge, uv_point)`.
pub(crate) fn clip_triangle_against_segment(
    tri_uv: &[(f64, f64); 3],
    seg_a: (f64, f64),
    seg_b: (f64, f64),
) -> Vec<(usize, f64, (f64, f64))> {
    let edges = [(0usize, 1usize), (1, 2), (2, 0)];
    let mut result = Vec::new();

    for (ei, &(vi, vj)) in edges.iter().enumerate() {
        if let Some((t, uv_pt)) = segment_segment_intersection(tri_uv[vi], tri_uv[vj], seg_a, seg_b)
        {
            result.push((ei, t, uv_pt));
        }
    }

    result
}

/// Standard 2D segment-segment intersection via cross products.
///
/// Segments are `(a1, a2)` and `(b1, b2)`. Returns `Some((t, point))` where
/// `t` is the parameter along the first segment and `point` is the intersection
/// in UV space. Returns `None` for parallel or non-intersecting segments.
pub(crate) fn segment_segment_intersection(
    a1: (f64, f64),
    a2: (f64, f64),
    b1: (f64, f64),
    b2: (f64, f64),
) -> Option<(f64, (f64, f64))> {
    let d1 = (a2.0 - a1.0, a2.1 - a1.1);
    let d2 = (b2.0 - b1.0, b2.1 - b1.1);

    let denom = d1.0 * d2.1 - d1.1 * d2.0;
    if denom.abs() < 1e-15 {
        return None; // Parallel or coincident.
    }

    let diff = (b1.0 - a1.0, b1.1 - a1.1);
    let t = (diff.0 * d2.1 - diff.1 * d2.0) / denom;
    let u = (diff.0 * d1.1 - diff.1 * d1.0) / denom;

    if !(0.0..=1.0).contains(&t) || !(0.0..=1.0).contains(&u) {
        return None; // Intersection outside both segments.
    }

    let pt = (a1.0 + t * d1.0, a1.1 + t * d1.1);
    Some((t, pt))
}

/// Collect edge splits for a triangle from the global edge-splits map.
///
/// For each of the 3 edges of `tri`, looks up the canonical `(min, max)` key
/// in `edge_splits` and returns the list of split vertex indices. The result
/// is `[(edge_index, split_vertex_indices), ...]`.
pub(crate) fn collect_edge_splits(
    tri: &[usize; 3],
    edge_splits: &HashMap<(usize, usize), Vec<usize>>,
) -> Vec<(usize, Vec<usize>)> {
    let edges = [(0usize, 1usize), (1, 2), (2, 0)];
    let mut result = Vec::new();

    for (ei, &(vi, vj)) in edges.iter().enumerate() {
        let key = (tri[vi].min(tri[vj]), tri[vi].max(tri[vj]));
        let splits = edge_splits.get(&key).cloned().unwrap_or_default();
        result.push((ei, splits));
    }

    result
}

/// Re-triangulate a triangle by inserting split vertices along its edges.
///
/// Builds an ordered vertex ring (vertex -> split points -> next vertex)
/// for each edge, then fan-triangulates from the first vertex. All resulting
/// sub-triangles are appended to `output`.
pub(crate) fn split_triangle_at_edge_vertices(
    tri: &[usize; 3],
    splits: &[(usize, Vec<usize>)],
    output: &mut Vec<[usize; 3]>,
) {
    let edge_indices = [(0usize, 1usize), (1, 2), (2, 0)];

    // Build ordered vertex ring.
    let mut ring: Vec<usize> = Vec::new();
    for (ei, (vi, _vj)) in edge_indices.iter().enumerate() {
        ring.push(tri[*vi]);
        // Find the split list for this edge.
        for (split_ei, split_verts) in splits {
            if *split_ei == ei {
                ring.extend_from_slice(split_verts);
            }
        }
    }

    if ring.len() < 3 {
        output.push(*tri);
        return;
    }

    // Fan-triangulate from the first vertex.
    for i in 1..ring.len() - 1 {
        let a = ring[0];
        let b = ring[i];
        let c = ring[i + 1];

        // Skip degenerate triangles.
        if a == b || b == c || a == c {
            continue;
        }

        output.push([a, b, c]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_core::geometry::{Point3, Vec3};

    /// Create a simple MeshVertex at given position with UV.
    fn make_vertex(x: f64, y: f64, z: f64, u: f64, v: f64) -> MeshVertex {
        MeshVertex {
            position: Point3::new(x, y, z),
            normal: Vec3::new(0.0, 0.0, 1.0),
            uv: (u, v),
        }
    }

    /// Create a Z-up plane surface for testing.
    fn test_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    /// Helper: create a unit-square outer trim loop.
    fn unit_square_loop() -> TrimLoop {
        TrimLoop {
            points: vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
            is_outer: true,
        }
    }

    #[test]
    fn clip_triangle_against_vertical_line() {
        let tri_uv = [(0.0, 0.0), (2.0, 0.0), (1.0, 2.0)];
        let hits = clip_triangle_against_segment(&tri_uv, (1.0, -1.0), (1.0, 3.0));
        assert!(!hits.is_empty(), "should find at least one intersection");
        let edge0_hits: Vec<_> = hits.iter().filter(|h| h.0 == 0).collect();
        assert_eq!(edge0_hits.len(), 1, "edge 0 should have 1 intersection");
        assert!(
            (edge0_hits[0].1 - 0.5).abs() < 1e-10,
            "edge 0 t should be 0.5"
        );
        assert!(
            (edge0_hits[0].2.0 - 1.0).abs() < 1e-10,
            "edge 0 x should be 1.0"
        );
    }

    #[test]
    fn apply_trim_preserves_inside_triangles() {
        let loops = vec![unit_square_loop()];
        let surface = test_plane();

        // Triangle fully inside the unit square.
        let vertices = vec![
            make_vertex(0.3, 0.3, 0.0, 0.3, 0.3),
            make_vertex(0.5, 0.3, 0.0, 0.5, 0.3),
            make_vertex(0.4, 0.5, 0.0, 0.4, 0.5),
        ];
        let triangles = vec![[0, 1, 2]];

        let (out_verts, out_tris) = apply_trim(&vertices, &triangles, &loops, &surface, true);

        assert_eq!(out_tris.len(), 1, "inside triangle should be preserved");
        assert_eq!(out_tris[0], [0, 1, 2]);
        assert_eq!(out_verts.len(), 3);
    }

    #[test]
    fn apply_trim_removes_outside_triangles() {
        let loops = vec![unit_square_loop()];
        let surface = test_plane();

        // Triangle fully outside the unit square.
        let vertices = vec![
            make_vertex(1.5, 1.5, 0.0, 1.5, 1.5),
            make_vertex(2.0, 1.5, 0.0, 2.0, 1.5),
            make_vertex(1.7, 2.0, 0.0, 1.7, 2.0),
        ];
        let triangles = vec![[0, 1, 2]];

        let (_out_verts, out_tris) = apply_trim(&vertices, &triangles, &loops, &surface, true);

        assert!(out_tris.is_empty(), "outside triangle should be removed");
    }

    #[test]
    fn segment_segment_intersection_basic() {
        // Crossing segments.
        let result = segment_segment_intersection((0.0, 0.0), (2.0, 2.0), (0.0, 2.0), (2.0, 0.0));
        assert!(result.is_some(), "crossing segments should intersect");
        let (t, pt) = result.unwrap();
        assert!((t - 0.5).abs() < 1e-10);
        assert!((pt.0 - 1.0).abs() < 1e-10);
        assert!((pt.1 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn segment_segment_intersection_parallel() {
        // Parallel segments.
        let result = segment_segment_intersection((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0));
        assert!(result.is_none(), "parallel segments should not intersect");
    }

    #[test]
    fn segment_segment_intersection_non_overlapping() {
        // Non-overlapping segments.
        let result = segment_segment_intersection((0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 1.0));
        assert!(
            result.is_none(),
            "non-overlapping segments should not intersect"
        );
    }

    #[test]
    fn split_triangle_at_edge_vertices_no_splits() {
        let tri = [0, 1, 2];
        let splits = vec![(0, vec![]), (1, vec![]), (2, vec![])];
        let mut output = Vec::new();
        split_triangle_at_edge_vertices(&tri, &splits, &mut output);
        assert_eq!(output.len(), 1);
        assert_eq!(output[0], [0, 1, 2]);
    }

    #[test]
    fn split_triangle_at_edge_vertices_with_split() {
        // Triangle [0,1,2] with vertex 3 inserted on edge 0 (between 0 and 1).
        let tri = [0, 1, 2];
        let splits = vec![(0, vec![3]), (1, vec![]), (2, vec![])];
        let mut output = Vec::new();
        split_triangle_at_edge_vertices(&tri, &splits, &mut output);
        // Ring: [0, 3, 1, 2] -> fan from 0: [0,3,1], [0,1,2]
        assert_eq!(output.len(), 2);
        assert_eq!(output[0], [0, 3, 1]);
        assert_eq!(output[1], [0, 1, 2]);
    }
}
