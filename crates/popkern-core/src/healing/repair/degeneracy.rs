//! Degeneracy repair: remove zero-length edges and collapse small faces.
//!
//! Handles degenerate edges at surface singularities (poles, collapsed
//! parameter lines) by detecting edges whose endpoints are coincident
//! and merging the vertex positions.

use crate::geometry::surface::SurfaceEval;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, FaceHandle};

/// Result of a degeneracy repair operation.
#[derive(Debug, Clone)]
pub struct DegeneracyRepairResult {
    /// Number of zero-length edges whose vertex positions were merged.
    pub zero_length_edges_removed: usize,
    /// Number of small faces collapsed (currently always 0; reserved for
    /// future iterations).
    pub small_faces_collapsed: usize,
    /// Number of zero-area faces detected (flagged for removal).
    pub zero_area_faces_detected: usize,
    /// Number of cusps detected (surface points where curvature is infinite).
    pub cusps_detected: usize,
    /// Number of self-tangent face locations detected.
    pub self_tangent_faces_detected: usize,
}

/// Remove zero-length edges from a BRep.
///
/// For each edge whose two endpoint vertices are coincident (distance <
/// tolerance), removes the edge geometrically by merging the two vertex
/// positions into one.  This handles degenerate edges at surface
/// singularities (poles, collapsed parameter lines).
///
/// Edges where both half-edges originate from the same vertex (seam edges
/// with identical origin) are skipped -- they are topological artefacts,
/// not degeneracies.
///
/// **Note:** this merges vertex *positions* only.  Full topological removal
/// of the edge (via Euler operators) is deferred to a future iteration
/// because it requires careful half-edge relinking.
///
/// Returns the number of edges whose vertex positions were merged.
pub fn remove_zero_length_edges(brep: &mut BRep, tolerance: &Tolerance) -> usize {
    // Collect candidate edges first to avoid borrowing conflicts.
    let candidates: Vec<EdgeHandle> = brep.edges.keys().collect();
    let linear_tol = tolerance.linear();
    let mut count = 0usize;

    for edge_handle in candidates {
        let Some(edge_data) = brep.edges.get(edge_handle) else {
            continue;
        };
        let he_handle = edge_data.half_edge;

        let Some(he) = brep.half_edges.get(he_handle) else {
            continue;
        };
        let start_vertex = he.origin;
        let twin_handle = he.twin;

        let Some(twin) = brep.half_edges.get(twin_handle) else {
            continue;
        };
        let end_vertex = twin.origin;

        // Skip seam edges where both half-edges originate from the same vertex.
        if start_vertex == end_vertex {
            continue;
        }

        let Some(start_data) = brep.vertices.get(start_vertex) else {
            continue;
        };
        let start_pos = start_data.position;

        let Some(end_data) = brep.vertices.get(end_vertex) else {
            continue;
        };
        let end_pos = end_data.position;

        // Check if the two endpoints are within linear tolerance.
        if start_pos.distance_to(&end_pos) >= linear_tol {
            continue;
        }

        // Merge: move the end vertex to the start vertex's position.
        if let Some(v) = brep.vertices.get_mut(end_vertex) {
            v.position = start_pos;
        }
        count += 1;
    }

    count
}

/// Repair all degeneracy defects in a BRep.
///
/// Removes zero-length edges (by merging vertex positions) and collapses
/// small faces.  Small-face collapse is reserved for a future iteration;
/// currently only zero-length edge repair is performed.
pub fn repair_degeneracies(brep: &mut BRep, tolerance: &Tolerance) -> DegeneracyRepairResult {
    let zero_length_edges_removed = remove_zero_length_edges(brep, tolerance);

    DegeneracyRepairResult {
        zero_length_edges_removed,
        small_faces_collapsed: 0,
        zero_area_faces_detected: 0,
        cusps_detected: 0,
        self_tangent_faces_detected: 0,
    }
}

/// Detect and repair zero-area face regions.
///
/// For each face, computes the area element from the first fundamental
/// form. If the total integrated area falls below a threshold, the face
/// is considered degenerate. Degenerate faces are flagged for removal.
///
/// Returns the list of zero-area face handles.
pub fn detect_zero_area_faces(brep: &BRep, tolerance: &Tolerance) -> Vec<FaceHandle> {
    let linear_tol = tolerance.linear();
    let area_threshold = linear_tol * linear_tol;
    let mut zero_area = Vec::new();

    for (face_handle, face_data) in &brep.faces {
        let surface_on_face = match face_data.surface.as_ref() {
            Some(s) => s,
            None => continue,
        };
        let surface = &surface_on_face.geometry;

        // Get bounded parameter range from face vertices.
        let positions = brep.face_outer_loop_positions(face_handle);
        if positions.is_empty() {
            continue;
        }
        let ((u_min, u_max), (v_min, v_max)) =
            crate::geometry::uv_projection::tight_bounded_parameter_range(&positions, surface);

        // Sample first fundamental form on a 4x4 grid and integrate area.
        let n = 4usize;
        let du = (u_max - u_min) / (n as f64);
        let dv = (v_max - v_min) / (n as f64);
        let cell_area = du * dv;
        let mut total_area = 0.0_f64;

        for i in 0..n {
            let u = u_min + du * (i as f64 + 0.5);
            for j in 0..n {
                let v = v_min + dv * (j as f64 + 0.5);
                let su = surface.partial_u(u, v);
                let sv = surface.partial_v(u, v);
                // Area element = |Su x Sv|
                let cross = su.cross(&sv);
                let area_element = cross.length();
                total_area += area_element * cell_area;
            }
        }

        if total_area < area_threshold {
            zero_area.push(face_handle);
        }
    }

    zero_area
}

/// Detect cusps -- surface points where curvature becomes infinite.
///
/// A cusp occurs where the surface folds to a point (the derivative
/// vanishes in one direction). Detected by finding locations where
/// the first fundamental form determinant approaches zero while the
/// second fundamental form is non-zero.
///
/// Returns a list of `(face, u, v)` triples locating each detected cusp.
pub fn detect_cusps(brep: &BRep, tolerance: &Tolerance) -> Vec<(FaceHandle, f64, f64)> {
    let linear_tol = tolerance.linear();
    let area_threshold = linear_tol * 0.01;
    let mut cusps = Vec::new();

    for (face_handle, face_data) in &brep.faces {
        let surface_on_face = match face_data.surface.as_ref() {
            Some(s) => s,
            None => continue,
        };
        let surface = &surface_on_face.geometry;

        // Get bounded parameter range from face vertices.
        let positions = brep.face_outer_loop_positions(face_handle);
        if positions.is_empty() {
            continue;
        }
        let ((u_min, u_max), (v_min, v_max)) =
            crate::geometry::uv_projection::tight_bounded_parameter_range(&positions, surface);

        // Sample on a 6x6 grid.
        let n = 6usize;
        for i in 0..n {
            let u = u_min + (u_max - u_min) * (i as f64 + 0.5) / (n as f64);
            for j in 0..n {
                let v = v_min + (v_max - v_min) * (j as f64 + 0.5) / (n as f64);

                // First fundamental form: E = Su·Su, F = Su·Sv, G = Sv·Sv
                let su = surface.partial_u(u, v);
                let sv = surface.partial_v(u, v);
                let cross = su.cross(&sv);
                let area_element = cross.length();

                if area_element >= area_threshold {
                    continue; // Not degenerate at this point.
                }

                // Area element is near zero -- derivative vanishes.
                // Check second fundamental form: if curvature is non-zero
                // at a degenerate point, it's a cusp.
                let suu = surface.second_partial_uu(u, v);
                let suv = surface.second_partial_uv(u, v);
                let svv = surface.second_partial_vv(u, v);

                // Second fundamental form coefficients (projected onto
                // normal direction). Use un-normalized normal to detect
                // magnitude.
                let second_mag = suu.length() + suv.length() + svv.length();
                if second_mag > linear_tol {
                    // Non-zero curvature at a degenerate point => cusp.
                    cusps.push((face_handle, u, v));
                }
            }
        }
    }

    cusps
}

/// Repair degeneracies: remove zero-length edges, flag zero-area faces
/// and cusps.
///
/// Extends the existing [`repair_degeneracies`] with zero-area and cusp
/// detection. This is the comprehensive variant.
pub fn repair_degeneracies_full(brep: &mut BRep, tolerance: &Tolerance) -> DegeneracyRepairResult {
    let zero_length_edges_removed = remove_zero_length_edges(brep, tolerance);
    let zero_area = detect_zero_area_faces(brep, tolerance);
    let cusps = detect_cusps(brep, tolerance);
    let self_tangent = detect_self_tangent_faces(brep, tolerance);

    DegeneracyRepairResult {
        zero_length_edges_removed,
        small_faces_collapsed: 0,
        zero_area_faces_detected: zero_area.len(),
        cusps_detected: cusps.len(),
        self_tangent_faces_detected: self_tangent.len(),
    }
}

/// Detect self-tangent surfaces -- surfaces that touch themselves without crossing.
///
/// Returns a list of `(face, u1, v1, u2, v2)` tuples identifying each
/// detected self-tangent location in parameter space.
pub fn detect_self_tangent_faces(
    brep: &BRep,
    tolerance: &Tolerance,
) -> Vec<(FaceHandle, f64, f64, f64, f64)> {
    let linear_tol = tolerance.linear();
    let mut results = Vec::new();
    for (face_handle, face_data) in &brep.faces {
        let sof = match face_data.surface.as_ref() {
            Some(s) => s,
            None => continue,
        };
        if !matches!(
            sof.geometry,
            crate::geometry::surface::SurfaceGeometry::Nurbs(_)
        ) {
            continue;
        }
        let positions = brep.face_outer_loop_positions(face_handle);
        if positions.is_empty() {
            continue;
        }
        let ((u_min, u_max), (v_min, v_max)) =
            crate::geometry::uv_projection::tight_bounded_parameter_range(
                &positions,
                &sof.geometry,
            );
        let u_span = u_max - u_min;
        let v_span = v_max - v_min;
        if u_span < 1e-15 || v_span < 1e-15 {
            continue;
        }
        let grid_n = 8usize;
        let mut grid_pts = Vec::with_capacity(grid_n * grid_n);
        let mut grid_uvs: Vec<(usize, usize, f64, f64)> = Vec::with_capacity(grid_n * grid_n);
        for i in 0..grid_n {
            let u = u_min + u_span * (i as f64 + 0.5) / (grid_n as f64);
            for j in 0..grid_n {
                let v = v_min + v_span * (j as f64 + 0.5) / (grid_n as f64);
                grid_pts.push(sof.geometry.point_at(u, v));
                grid_uvs.push((i, j, u, v));
            }
        }
        let uv_min_dist = 0.1 * (u_span * u_span + v_span * v_span).sqrt();
        for a in 0..grid_pts.len() {
            let (ai, aj, au, av) = grid_uvs[a];
            for b in (a + 1)..grid_pts.len() {
                let (bi, bj, bu, bv) = grid_uvs[b];
                if ai.abs_diff(bi) <= 1 && aj.abs_diff(bj) <= 1 {
                    continue;
                }
                let du = au - bu;
                let dv = av - bv;
                if (du * du + dv * dv).sqrt() < uv_min_dist {
                    continue;
                }
                let pa = grid_pts[a];
                let pb = grid_pts[b];
                if crate::adaptive_precision::adaptive_dist_sq_lt(
                    (pa.x, pa.y, pa.z),
                    (pb.x, pb.y, pb.z),
                    linear_tol * linear_tol,
                ) {
                    results.push((face_handle, au, av, bu, bv));
                }
            }
        }
    }
    results
}

#[cfg(test)]
#[path = "degeneracy_tests.rs"]
mod tests;
