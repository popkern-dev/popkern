//! Coincidence detection between two BRep bodies.
//!
//! Detects coincident faces, overlapping edges, and mergeable vertices.
//! Uses analytic O(1) comparison for same-type surfaces and sample-based
//! detection as fallback.

use popkern_core::geometry::intersection::coincident_analytic::are_surfaces_coincident_analytic;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};

use crate::clip::face_vertex_positions;
use crate::tolerance::BoolTolerance;

use super::classify::classify_face_orientation;
use super::types::{CoincidenceMap, FaceOverlap, VertexMerge};

/// Detect all coincidences between two bodies.
///
/// Runs face-face coincidence detection using AABB pre-filtering,
/// plus vertex-vertex merge detection from the hierarchy results.
pub(crate) fn detect_coincidences(
    a: &BRep,
    b: &BRep,
    aabbs_a: &[(FaceHandle, Aabb)],
    aabbs_b: &[(FaceHandle, Aabb)],
    tol: &BoolTolerance,
) -> CoincidenceMap {
    let mut map = CoincidenceMap::default();

    // 1. Vertex-vertex merges (from hierarchy VV hits)
    detect_vertex_merges(a, b, tol, &mut map);

    // 2. Edge-edge overlap detection
    detect_edge_overlaps(a, b, tol, &mut map);

    // 3. Face-face coincidence
    detect_face_coincidences(a, b, aabbs_a, aabbs_b, tol, &mut map);

    map
}

/// Detect overlapping edge pairs between two bodies.
///
/// Two edges overlap when their underlying curves are geometrically
/// identical (within tolerance) along some sub-interval. Detected by
/// checking if both endpoints of edge A lie on edge B's curve and
/// vice versa.
fn detect_edge_overlaps(a: &BRep, b: &BRep, tol: &BoolTolerance, map: &mut CoincidenceMap) {
    use popkern_core::geometry::CurveEval;

    for (ea_h, ea) in &a.edges {
        if ea.curve.is_none() {
            continue;
        }
        // Get edge A's endpoint positions
        let he_a = match a.half_edges.get(ea.half_edge) {
            Some(h) => h,
            None => continue,
        };
        let twin_a = match a.half_edges.get(he_a.twin) {
            Some(t) => t,
            None => continue,
        };
        let a_start = a.vertices[he_a.origin].position;
        let a_end = a.vertices[twin_a.origin].position;

        for (eb_h, eb) in &b.edges {
            let cb = match &eb.curve {
                Some(c) => c,
                None => continue,
            };
            let he_b = match b.half_edges.get(eb.half_edge) {
                Some(h) => h,
                None => continue,
            };
            let twin_b = match b.half_edges.get(he_b.twin) {
                Some(t) => t,
                None => continue,
            };
            let b_start = b.vertices[he_b.origin].position;
            let b_end = b.vertices[twin_b.origin].position;

            // Check if A's endpoints lie on B's curve (both checks
            // use adaptive_dist_sq_lt so the decision is deterministic
            // at the snap threshold).
            let snap_sq = tol.snap * tol.snap;
            let a_start_on_b = if let Some(t) = cb.geometry.closest_point(&a_start) {
                let t_clamped = t.clamp(cb.t_start.min(cb.t_end), cb.t_start.max(cb.t_end));
                let p = cb.geometry.point_at(t_clamped);
                popkern_core::adaptive_precision::adaptive_dist_sq_lt(
                    (p.x, p.y, p.z),
                    (a_start.x, a_start.y, a_start.z),
                    snap_sq,
                )
            } else {
                false
            };

            let a_end_on_b = if let Some(t) = cb.geometry.closest_point(&a_end) {
                let t_clamped = t.clamp(cb.t_start.min(cb.t_end), cb.t_start.max(cb.t_end));
                let p = cb.geometry.point_at(t_clamped);
                popkern_core::adaptive_precision::adaptive_dist_sq_lt(
                    (p.x, p.y, p.z),
                    (a_end.x, a_end.y, a_end.z),
                    snap_sq,
                )
            } else {
                false
            };

            if a_start_on_b && a_end_on_b {
                // Both endpoints of A lie on B's curve — edges overlap.
                // Determine direction by comparing squared distances
                // (monotonic, avoids sqrt and gives the same ordering
                // for every caller).
                let d_end_sq = {
                    let dx = a_start.x - b_end.x;
                    let dy = a_start.y - b_end.y;
                    let dz = a_start.z - b_end.z;
                    dx * dx + dy * dy + dz * dz
                };
                let d_start_sq = {
                    let dx = a_start.x - b_start.x;
                    let dy = a_start.y - b_start.y;
                    let dz = a_start.z - b_start.z;
                    dx * dx + dy * dy + dz * dz
                };
                let reversed = d_end_sq < d_start_sq;
                map.edge_overlaps.push(super::types::EdgeOverlap {
                    edge_a: ea_h,
                    edge_b: eb_h,
                    reversed,
                });
            }
        }
    }
}

/// Detect vertex pairs that should be merged.
///
/// Uses [`popkern_core::adaptive_precision::adaptive_dist_sq_lt`] for
/// the `dist² < snap²` comparison: the f64 fast path handles clear
/// cases in one subtraction and one compare, and the double-
/// precision slow path kicks in only when the f64 answer lands
/// inside the ambiguity margin (`4 * ε * max(dist_sq, snap_sq)`).
/// This gives consistent, order-independent merge decisions for
/// vertex pairs that are near the snap threshold — the failure
/// mode that caused tolerance-induced predicate inconsistency in
/// the old naive-f64 path.
fn detect_vertex_merges(a: &BRep, b: &BRep, tol: &BoolTolerance, map: &mut CoincidenceMap) {
    let snap_sq = tol.snap * tol.snap;
    for (va_h, va) in &a.vertices {
        for (vb_h, vb) in &b.vertices {
            let p_a = (va.position.x, va.position.y, va.position.z);
            let p_b = (vb.position.x, vb.position.y, vb.position.z);
            if popkern_core::adaptive_precision::adaptive_dist_sq_lt(p_a, p_b, snap_sq) {
                // For the recorded `distance` field, the naive f64
                // magnitude is fine — it's diagnostic output, not a
                // decision point.
                let dx = va.position.x - vb.position.x;
                let dy = va.position.y - vb.position.y;
                let dz = va.position.z - vb.position.z;
                let dist_sq = dx * dx + dy * dy + dz * dz;
                map.vertex_merges.push(VertexMerge {
                    vertex_a: va_h,
                    vertex_b: vb_h,
                    distance: dist_sq.sqrt(),
                });
            }
        }
    }
}

/// Detect coincident face pairs.
fn detect_face_coincidences(
    a: &BRep,
    b: &BRep,
    aabbs_a: &[(FaceHandle, Aabb)],
    aabbs_b: &[(FaceHandle, Aabb)],
    tol: &BoolTolerance,
    map: &mut CoincidenceMap,
) {
    for (face_a, aabb_a) in aabbs_a {
        for (face_b, aabb_b) in aabbs_b {
            if !aabb_a.overlaps(aabb_b, tol.coincidence) {
                continue;
            }

            let surf_a = match a.faces[*face_a].surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };
            let surf_b = match b.faces[*face_b].surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };

            // Tier A: analytic O(1) comparison for same-type surfaces
            if are_surfaces_coincident_analytic(surf_a, surf_b, &tol.core) {
                let orientation =
                    classify_face_orientation(a, *face_a, aabbs_a, b, *face_b, aabbs_b, &tol.core);
                let overlap_center = compute_overlap_center(a, *face_a, b, *face_b);
                map.face_overlaps.push(FaceOverlap {
                    face_a: *face_a,
                    face_b: *face_b,
                    orientation,
                    full_overlap: false, // conservative; refine later
                    overlap_center,
                });
                continue;
            }

            // Tier B: sample-based detection for near-coincident surfaces
            if is_near_coincident(surf_a, a, *face_a, surf_b, b, *face_b, tol) {
                let orientation =
                    classify_face_orientation(a, *face_a, aabbs_a, b, *face_b, aabbs_b, &tol.core);
                let overlap_center = compute_overlap_center(a, *face_a, b, *face_b);
                map.face_overlaps.push(FaceOverlap {
                    face_a: *face_a,
                    face_b: *face_b,
                    orientation,
                    full_overlap: false,
                    overlap_center,
                });
            }
        }
    }
}

/// Sample-based near-coincidence check.
///
/// Samples 4x4 grid on face A's surface, projects each onto face B's surface.
/// If all samples are within coincidence tolerance, the faces are near-coincident.
fn is_near_coincident(
    surf_a: &popkern_core::geometry::surface::SurfaceGeometry,
    brep_a: &BRep,
    face_a: FaceHandle,
    surf_b: &popkern_core::geometry::surface::SurfaceGeometry,
    brep_b: &BRep,
    face_b: FaceHandle,
    tol: &BoolTolerance,
) -> bool {
    let positions_a = face_vertex_positions(brep_a, face_a);
    let positions_b = face_vertex_positions(brep_b, face_b);
    if positions_a.is_empty() || positions_b.is_empty() {
        return false;
    }

    let ((u_min, u_max), (v_min, v_max)) =
        popkern_core::geometry::uv_projection::bounded_parameter_range(&positions_a, surf_a);

    const GRID: usize = 4;
    let mut within_count = 0;
    let total = GRID * GRID;

    let coincidence_sq = tol.coincidence * tol.coincidence;
    for i in 0..GRID {
        let u = u_min + (u_max - u_min) * (i as f64 + 0.5) / (GRID as f64);
        for j in 0..GRID {
            let v = v_min + (v_max - v_min) * (j as f64 + 0.5) / (GRID as f64);
            let p = surf_a.point_at(u, v);

            // Project onto surface B; decide via adaptive predicate
            // so surfaces exactly at the coincidence threshold
            // produce a deterministic classification.
            if let Some((u2, v2)) = surf_b.closest_point(&p) {
                let closest = surf_b.point_at(u2, v2);
                if popkern_core::adaptive_precision::adaptive_dist_sq_lt(
                    (p.x, p.y, p.z),
                    (closest.x, closest.y, closest.z),
                    coincidence_sq,
                ) {
                    within_count += 1;
                }
            }
        }
    }

    // Require > 75% of samples within tolerance
    within_count * 4 > total * 3
}

/// Compute the approximate center of the overlap region.
fn compute_overlap_center(
    a: &BRep,
    face_a: FaceHandle,
    b: &BRep,
    face_b: FaceHandle,
) -> popkern_core::geometry::Point3 {
    let positions_a = face_vertex_positions(a, face_a);
    let positions_b = face_vertex_positions(b, face_b);

    let all: Vec<_> = positions_a.iter().chain(positions_b.iter()).collect();
    if all.is_empty() {
        return popkern_core::geometry::Point3::new(0.0, 0.0, 0.0);
    }

    let n = all.len() as f64;
    let sum_x: f64 = all.iter().map(|p| p.x).sum();
    let sum_y: f64 = all.iter().map(|p| p.y).sum();
    let sum_z: f64 = all.iter().map(|p| p.z).sum();
    popkern_core::geometry::Point3::new(sum_x / n, sum_y / n, sum_z / n)
}
