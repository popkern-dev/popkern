//! Boundary triangle clipping and T-junction prevention.
//!
//! Clips BOUNDARY triangles against trim polyline segments with vertex snapping.
//! Propagates split vertices into neighboring INSIDE triangles to prevent T-junctions.

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use crate::mesh::MeshVertex;
use crate::trim::{TriClass, TrimLoop, classify_triangle, point_inside_trim_loops};

/// Snap tolerance for merging intersection points with existing triangle vertices.
const SNAP_TOL: f64 = 1e-10;

/// Maximum 3D distance between surface-evaluated and trim-interpolated positions
/// before the interpolated 3D position is used instead.
const UV_3D_ERROR_TOL: f64 = 0.1;

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
    orientation: bool,
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
                    orientation,
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
#[allow(clippy::too_many_arguments)]
fn clip_boundary_triangle(
    tri: &[usize; 3],
    vertices: &[MeshVertex],
    trim_loops: &[TrimLoop],
    surface: &SurfaceGeometry,
    orientation: bool,
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

            // Create new vertex by evaluating the surface at the UV point,
            // with 3D validation: if the surface-evaluated position diverges
            // from the trim polyline's interpolated 3D position, prefer the
            // trim 3D position (projection-based UVs may be imprecise).
            let pos_from_uv = surface.point_at(uv_pt.0, uv_pt.1);
            let (pos, normal) = if let Some(expected_3d) = interpolate_trim_3d(uv_pt, trim_loops) {
                if pos_from_uv.distance_to(&expected_3d) > UV_3D_ERROR_TOL {
                    // UV is unreliable — use the interpolated 3D position
                    // and derive the normal from the nearest valid UV.
                    let fallback_normal = surface.normal_at(uv_pt.0, uv_pt.1);
                    (expected_3d, fallback_normal)
                } else {
                    (pos_from_uv, surface.normal_at(uv_pt.0, uv_pt.1))
                }
            } else {
                (pos_from_uv, surface.normal_at(uv_pt.0, uv_pt.1))
            };
            // Flip normal for reversed-orientation faces (e.g., inner shell walls).
            let normal = if orientation { normal } else { -normal };
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
    // Track which ring entries are boundary intersections (always part
    // of the clipped region) vs original triangle corners (kept only
    // when inside the trim region).
    let mut ring: Vec<(usize, bool)> = Vec::new(); // (vertex, is_intersection)
    for ei in 0..3 {
        ring.push((tri[edge_indices[ei].0], false));
        for &v in &edge_new_verts[ei] {
            ring.push((v, true));
        }
    }
    if ring.len() < 3 {
        return;
    }

    // Clipped polygon = intersection points + inside grid corners +
    // TRIM POLYGON VERTICES lying inside this triangle (the previous
    // construction dropped the corner wedge of every cell containing
    // a trim vertex — measured 9% area loss on a clean synthetic
    // triangle trim). The clipped piece is convex for convex trim
    // corners, so angular ordering around the centroid triangulates
    // it; the centroid-inside test stays as the non-convex backstop.
    let mut clipped: Vec<usize> = Vec::new();
    for &(v, is_intersection) in &ring {
        if is_intersection {
            clipped.push(v);
            continue;
        }
        let uv = new_verts[v].uv;
        if point_inside_trim_loops(uv.0, uv.1, trim_loops) {
            clipped.push(v);
        }
    }
    // Interior trim vertices (trim corners poking into this cell)
    // join the polygon by cheapest insertion, which preserves the
    // ring's walk order — an angular sort around the centroid is NOT
    // order-safe here and produced bowtie slivers spanning the face.
    for tl in trim_loops {
        for &(pu, pv) in &tl.points {
            if !crate::trim::point_in_triangle_uv((pu, pv), &tri_uv)
                || clipped.iter().any(|&v| {
                    let uv = new_verts[v].uv;
                    (uv.0 - pu).abs() < SNAP_TOL && (uv.1 - pv).abs() < SNAP_TOL
                })
            {
                continue;
            }
            let pos = surface.point_at(pu, pv);
            let mut normal = surface.normal_at(pu, pv);
            if !orientation {
                normal = -normal;
            }
            let idx = new_verts.len();
            new_verts.push(MeshVertex {
                position: pos,
                normal,
                uv: (pu, pv),
            });
            if clipped.len() < 2 {
                clipped.push(idx);
                continue;
            }
            let dist =
                |a: (f64, f64), b: (f64, f64)| ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt();
            let mut best_slot = 0;
            let mut best_cost = f64::INFINITY;
            for slot in 0..clipped.len() {
                let a = new_verts[clipped[slot]].uv;
                let b = new_verts[clipped[(slot + 1) % clipped.len()]].uv;
                let cost = dist(a, (pu, pv)) + dist((pu, pv), b) - dist(a, b);
                if cost < best_cost {
                    best_cost = cost;
                    best_slot = slot;
                }
            }
            clipped.insert(best_slot + 1, idx);
        }
    }
    clipped.dedup();
    if clipped.len() >= 2 && clipped[0] == clipped[clipped.len() - 1] {
        clipped.pop();
    }
    if clipped.len() < 3 {
        // No usable clipped polygon (e.g. the trim corner pokes into
        // the cell without edge crossings, or near-line corners
        // misclassify): fall back to fanning the full ring with the
        // centroid filter rather than dropping the cell.
        for i in 1..ring.len() - 1 {
            let a = ring[0].0;
            let b = ring[i].0;
            let c = ring[i + 1].0;
            if a == b || b == c || a == c {
                continue;
            }
            let uv_a = new_verts[a].uv;
            let uv_b = new_verts[b].uv;
            let uv_c = new_verts[c].uv;
            let cx = (uv_a.0 + uv_b.0 + uv_c.0) / 3.0;
            let cy = (uv_a.1 + uv_b.1 + uv_c.1) / 3.0;
            if point_inside_trim_loops(cx, cy, trim_loops) {
                let signed =
                    (uv_b.0 - uv_a.0) * (uv_c.1 - uv_a.1) - (uv_c.0 - uv_a.0) * (uv_b.1 - uv_a.1);
                if signed >= 0.0 {
                    new_tris.push([a, b, c]);
                } else {
                    new_tris.push([a, c, b]);
                }
            }
        }
        return;
    }

    // Fan-triangulate the clipped polygon from an intersection vertex
    // when one exists (guaranteed on the region boundary); keep the
    // centroid test as a backstop for non-convex corner cases.
    let apex = clipped
        .iter()
        .position(|v| ring.iter().any(|&(rv, ri)| ri && rv == *v))
        .unwrap_or(0);
    let m = clipped.len();
    for k in 1..m - 1 {
        let a = clipped[apex];
        let b = clipped[(apex + k) % m];
        let c = clipped[(apex + k + 1) % m];
        if a == b || b == c || a == c {
            continue;
        }
        let uv_a = new_verts[a].uv;
        let uv_b = new_verts[b].uv;
        let uv_c = new_verts[c].uv;
        let cx = (uv_a.0 + uv_b.0 + uv_c.0) / 3.0;
        let cy = (uv_a.1 + uv_b.1 + uv_c.1) / 3.0;
        if point_inside_trim_loops(cx, cy, trim_loops) {
            // Winding: the angular sort's direction is arbitrary;
            // grid triangles are CCW in UV, so clipped output must be
            // too (mixed winding rendered whole boundary bands as
            // backfaces).
            let signed =
                (uv_b.0 - uv_a.0) * (uv_c.1 - uv_a.1) - (uv_c.0 - uv_a.0) * (uv_b.1 - uv_a.1);
            if signed >= 0.0 {
                new_tris.push([a, b, c]);
            } else {
                new_tris.push([a, c, b]);
            }
        }
    }
}

/// Interpolate the 3D position on the nearest trim polyline segment for a UV point.
///
/// Finds the trim segment whose UV projection is closest to `uv_pt` and
/// interpolates the corresponding 3D position along that segment. Returns
/// `None` if no trim loops have 3D data.
fn interpolate_trim_3d(uv_pt: (f64, f64), trim_loops: &[TrimLoop]) -> Option<Point3> {
    let mut best_dist_sq = f64::INFINITY;
    let mut best_pos: Option<Point3> = None;

    for tl in trim_loops {
        if tl.points_3d.len() != tl.points.len() || tl.points.len() < 2 {
            continue;
        }
        let n = tl.points.len();
        for i in 0..n {
            let j = (i + 1) % n;
            let (ax, ay) = tl.points[i];
            let (bx, by) = tl.points[j];
            // Project uv_pt onto segment [a, b] in UV space.
            let dx = bx - ax;
            let dy = by - ay;
            let seg_len_sq = dx * dx + dy * dy;
            if seg_len_sq < 1e-30 {
                continue;
            }
            let t = ((uv_pt.0 - ax) * dx + (uv_pt.1 - ay) * dy) / seg_len_sq;
            let t_clamped = t.clamp(0.0, 1.0);
            let proj_u = ax + t_clamped * dx;
            let proj_v = ay + t_clamped * dy;
            let dist_sq = (uv_pt.0 - proj_u).powi(2) + (uv_pt.1 - proj_v).powi(2);
            if dist_sq < best_dist_sq {
                best_dist_sq = dist_sq;
                let p0 = tl.points_3d[i];
                let p1 = tl.points_3d[j];
                best_pos = Some(Point3::new(
                    p0.x + t_clamped * (p1.x - p0.x),
                    p0.y + t_clamped * (p1.y - p0.y),
                    p0.z + t_clamped * (p1.z - p0.z),
                ));
            }
        }
    }

    best_pos
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
#[path = "trim_clip_tests.rs"]
mod tests;
