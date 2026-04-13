//! Helper functions for topology manipulation during blend construction.
//!
//! Provides perpendicular edge splitting, consumed-face identification,
//! and fillet-face merging used by the main `apply_blend_geometry` orchestrator.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, LoopHandle, VertexHandle};

pub use super::merge::{
    identify_consumed_face, identify_consumed_face_contains, merge_consumed_faces,
};

/// Split perpendicular edges at the fillet offset distance.
///
/// For each face adjacent to the fillet edge, finds the boundary edges that
/// share a vertex with the fillet edge, and splits them at the offset distance
/// (fillet radius) from the shared vertex. This creates the anchor vertices
/// where the trim curves will start and end.
///
/// The 3D trim curves from the geometry pipeline may extend beyond the face
/// boundary, so instead we compute split positions geometrically from the
/// fillet edge endpoints and radius.
pub fn split_perpendicular_edges_by_radius(
    brep: &mut BRep,
    face_1: FaceHandle,
    face_2: FaceHandle,
    fillet_edge: EdgeHandle,
    radius: f64,
    _tolerance: &Tolerance,
) -> Result<SplitInfo, String> {
    let edge_data = brep
        .edges
        .get(fillet_edge)
        .ok_or_else(|| format!("fillet edge {fillet_edge:?} not found"))?;
    let he_fwd = edge_data.half_edge;
    let he_fwd_data = brep
        .half_edges
        .get(he_fwd)
        .ok_or_else(|| "half-edge not found".to_string())?;
    let he_rev = he_fwd_data.twin;
    let v_start = he_fwd_data.origin;
    let v_end = brep
        .half_edges
        .get(he_rev)
        .ok_or_else(|| "twin not found".to_string())?
        .origin;

    // Closed fillet edges (v_start == v_end) — e.g. the rim of a
    // cylindrical hole — are not yet supported. The construct phase
    // builds anchor vertices on perpendicular boundary edges and
    // there are no such edges when the fillet edge is a closed loop.
    // Return a clean diagnostic so the caller sees an actionable
    // error instead of "no perpendicular edge found".
    if v_start == v_end {
        return Err(format!(
            "closed fillet edges (v_start == v_end on edge {fillet_edge:?}) are not yet \
             supported by the construct phase — the cylindrical-hole rim case requires a \
             separate non-anchored construction path"
        ));
    }

    // Split perpendicular edges on face_1.
    let result_1 = split_face_perp_edges(brep, face_1, fillet_edge, v_start, v_end, radius)?;

    // Split perpendicular edges on face_2.
    let result_2 = split_face_perp_edges(brep, face_2, fillet_edge, v_start, v_end, radius)?;

    // Collect all edge splits for handle remapping.
    let mut edge_splits = result_1.edge_splits;
    edge_splits.extend(result_2.edge_splits);

    // Build trim 3D curves from the new vertex positions (straight lines for now).
    let p1_a = brep.vertices[result_1.anchor_start].position;
    let p1_b = brep.vertices[result_1.anchor_end].position;
    let p2_a = brep.vertices[result_2.anchor_start].position;
    let p2_b = brep.vertices[result_2.anchor_end].position;

    let trim_3d_1 = line_nurbs(&p1_a, &p1_b);
    let trim_3d_2 = line_nurbs(&p2_a, &p2_b);

    Ok(SplitInfo {
        trim_3d_1,
        trim_3d_2,
        face_1_anchors: (result_1.anchor_start, result_1.anchor_end),
        face_2_anchors: (result_2.anchor_start, result_2.anchor_end),
        edge_splits,
    })
}

/// Information returned by perpendicular edge splitting.
#[derive(Debug)]
pub struct SplitInfo {
    /// Corrected 3D trim curve for face 1 (kept for compatibility).
    #[allow(dead_code)]
    pub trim_3d_1: NurbsCurve,
    /// Corrected 3D trim curve for face 2 (kept for compatibility).
    #[allow(dead_code)]
    pub trim_3d_2: NurbsCurve,
    /// Anchor vertices on face 1 (start, end matching fillet edge direction).
    pub face_1_anchors: (VertexHandle, VertexHandle),
    /// Anchor vertices on face 2 (start, end matching fillet edge direction).
    pub face_2_anchors: (VertexHandle, VertexHandle),
    /// Edge splits performed: `(old_edge, edge_a, edge_b)`.
    #[allow(dead_code)]
    pub edge_splits: Vec<(EdgeHandle, EdgeHandle, EdgeHandle)>,
}

/// Build a simple linear NURBS curve between two points.
pub(crate) fn line_nurbs(a: &Point3, b: &Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![*a, *b],
        weights: vec![1.0, 1.0],
    }
}

/// Build a pcurve (UV-space line) between two points projected onto a surface.
fn compute_pcurve_for_trim(surface: &SurfaceGeometry, pt_a: &Point3, pt_b: &Point3) -> NurbsCurve {
    let (u_a, v_a) = project_to_uv(surface, pt_a);
    let (u_b, v_b) = project_to_uv(surface, pt_b);
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(u_a, v_a, 0.0), Point3::new(u_b, v_b, 0.0)],
        weights: vec![1.0, 1.0],
    }
}

/// Project a 3D point onto a surface's UV parameter space.
pub(crate) fn project_to_uv(surface: &SurfaceGeometry, point: &Point3) -> (f64, f64) {
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    match project_point_to_uv(surface, point) {
        Some(uv) => uv,
        None => {
            // Fallback: grid search.
            use popkern_core::geometry::surface::SurfaceEval;
            let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
            let u_lo = if u_min.is_finite() { u_min } else { -10.0 };
            let u_hi = if u_max.is_finite() { u_max } else { 10.0 };
            let v_lo = if v_min.is_finite() { v_min } else { -10.0 };
            let v_hi = if v_max.is_finite() { v_max } else { 10.0 };
            let n = 20;
            let mut best = ((u_lo + u_hi) * 0.5, (v_lo + v_hi) * 0.5);
            let mut best_d = f64::INFINITY;
            for iu in 0..n {
                let u = u_lo + (u_hi - u_lo) * (iu as f64 / (n - 1) as f64);
                for iv in 0..n {
                    let v = v_lo + (v_hi - v_lo) * (iv as f64 / (n - 1) as f64);
                    let d = surface.point_at(u, v).distance_to(point);
                    if d < best_d {
                        best_d = d;
                        best = (u, v);
                    }
                }
            }
            best
        }
    }
}

/// Compute a pcurve for a trim curve on the given face surface.
pub fn pcurve_for_face(
    brep: &BRep,
    face: FaceHandle,
    trim_start: &Point3,
    trim_end: &Point3,
) -> Result<NurbsCurve, String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} not found"))?;
    let surface = face_data
        .surface
        .as_ref()
        .ok_or_else(|| format!("face {face:?} has no surface"))?;
    Ok(compute_pcurve_for_trim(
        &surface.geometry,
        trim_start,
        trim_end,
    ))
}

/// Result of splitting perpendicular edges on a single face.
struct FacePerpSplitResult {
    /// Anchor vertex near the fillet edge start vertex.
    anchor_start: VertexHandle,
    /// Anchor vertex near the fillet edge end vertex.
    anchor_end: VertexHandle,
    /// Edge splits: `(old_edge, edge_a, edge_b)`.
    edge_splits: Vec<(EdgeHandle, EdgeHandle, EdgeHandle)>,
}

/// For a given face, find the two perpendicular boundary edges that share
/// vertices with the fillet edge, and split them at `radius` distance from
/// the shared vertex. Returns the two new anchor vertices and the list of
/// edge splits performed as `(old_edge, edge_a, edge_b)`.
///
/// Walks every loop on the face (outer + inner). The fillet edge can sit
/// on an inner loop when the face has a hole — for example, when the
/// fillet edge is the rim of a cylindrical hole on a planar face. Walking
/// only the outer loop in that case finds no incident edges and produces
/// the historical "no perpendicular edge found" failure.
fn split_face_perp_edges(
    brep: &mut BRep,
    face: FaceHandle,
    fillet_edge: EdgeHandle,
    v_start: VertexHandle,
    v_end: VertexHandle,
    radius: f64,
) -> Result<FacePerpSplitResult, String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} not found"))?;

    // Walk every loop on the face — outer + inner. The fillet edge may
    // be on an inner loop (hole boundary) when the support face has a
    // hole, e.g. a cube top face with a cylindrical drilling.
    let mut loop_handles: Vec<LoopHandle> = Vec::with_capacity(1 + face_data.inner_loops.len());
    loop_handles.push(face_data.outer_loop);
    loop_handles.extend(face_data.inner_loops.iter().copied());

    let mut perp_start: Option<(EdgeHandle, f64)> = None;
    let mut perp_end: Option<(EdgeHandle, f64)> = None;

    let max_iters = brep.half_edges.len() + 1;
    for loop_handle in loop_handles {
        let start_he = brep
            .loops
            .get(loop_handle)
            .ok_or_else(|| format!("loop {loop_handle:?} not found"))?
            .half_edge;

        let mut cursor = start_he;
        let mut iters = 0;
        loop {
            let he = brep
                .half_edges
                .get(cursor)
                .ok_or_else(|| format!("half-edge {cursor:?} not found"))?;
            let next_he = brep
                .half_edges
                .get(he.next)
                .ok_or_else(|| "next half-edge not found".to_string())?;

            if he.edge != fillet_edge {
                let origin = he.origin;
                let dest = next_he.origin;
                let p0 = brep.vertices[origin].position;
                let p1 = brep.vertices[dest].position;
                let edge_len = p0.distance_to(&p1);

                // `split_edge(edge, t)` interprets `t` against the
                // edge's canonical primary half-edge — *not* the
                // half-edge we're currently walking. When the walked
                // half-edge is the twin, the parameter must be
                // inverted so the new vertex still lands at
                // `radius` from the matched fillet endpoint.
                let canonical_he_h = brep.edges[he.edge].half_edge;
                let canonical_origin = brep.half_edges[canonical_he_h].origin;
                let walked_matches_canonical = canonical_origin == origin;
                let invert = !walked_matches_canonical;

                let normalise = |t_from_origin: f64| -> f64 {
                    let raw = if invert {
                        1.0 - t_from_origin
                    } else {
                        t_from_origin
                    };
                    raw.clamp(0.001, 0.999)
                };

                if origin == v_start && edge_len > 1e-15 {
                    let t_from_origin = (radius / edge_len).min(0.999);
                    if t_from_origin > 0.001 && perp_start.is_none() {
                        perp_start = Some((he.edge, normalise(t_from_origin)));
                    }
                } else if dest == v_start && edge_len > 1e-15 {
                    // The matched endpoint is at the END of the edge
                    // (in walked-direction terms), so t-from-walked-
                    // origin is `1 - r/L`. The normalise() call
                    // re-orients it for the canonical direction.
                    let t_from_origin = 1.0 - (radius / edge_len).min(0.999);
                    if t_from_origin > 0.001 && t_from_origin < 0.999 && perp_start.is_none() {
                        perp_start = Some((he.edge, normalise(t_from_origin)));
                    }
                }

                if origin == v_end && edge_len > 1e-15 {
                    let t_from_origin = (radius / edge_len).min(0.999);
                    if t_from_origin > 0.001 && perp_end.is_none() {
                        perp_end = Some((he.edge, normalise(t_from_origin)));
                    }
                } else if dest == v_end && edge_len > 1e-15 {
                    let t_from_origin = 1.0 - (radius / edge_len).min(0.999);
                    if t_from_origin > 0.001 && t_from_origin < 0.999 && perp_end.is_none() {
                        perp_end = Some((he.edge, normalise(t_from_origin)));
                    }
                }
            }

            cursor = he.next;
            iters += 1;
            if cursor == start_he || iters >= max_iters {
                break;
            }
        }
    }

    let (edge_a, t_a) = perp_start
        .ok_or_else(|| format!("no perpendicular edge found at v_start on face {face:?}"))?;
    let (edge_b, t_b) =
        perp_end.ok_or_else(|| format!("no perpendicular edge found at v_end on face {face:?}"))?;

    // Split the edges to create anchor vertices.
    let split_a = brep
        .split_edge(edge_a, t_a)
        .map_err(|e| format!("split_edge at v_start perp failed: {e}"))?;
    let split_b = brep
        .split_edge(edge_b, t_b)
        .map_err(|e| format!("split_edge at v_end perp failed: {e}"))?;

    let edge_splits = vec![
        (edge_a, split_a.edge_a, split_a.edge_b),
        (edge_b, split_b.edge_a, split_b.edge_b),
    ];

    Ok(FacePerpSplitResult {
        anchor_start: split_a.vertex,
        anchor_end: split_b.vertex,
        edge_splits,
    })
}

/// Extract the fillet radius from the edge geometry analysis.
///
/// Uses the distance between a trim curve endpoint and the nearest fillet
/// edge vertex as a fallback radius estimate.
#[allow(dead_code)]
pub fn estimate_radius_from_geometry(
    brep: &BRep,
    fillet_edge: EdgeHandle,
    trim_3d: &NurbsCurve,
) -> f64 {
    let edge_data = match brep.edges.get(fillet_edge) {
        Some(e) => e,
        None => return 0.0,
    };
    let he = match brep.half_edges.get(edge_data.half_edge) {
        Some(h) => h,
        None => return 0.0,
    };
    let v_start_pos = brep.vertices[he.origin].position;

    let (t_start, _) = trim_3d.parameter_range();
    let trim_start = trim_3d.point_at(t_start);

    v_start_pos.distance_to(&trim_start)
}

#[cfg(test)]
#[path = "helpers_tests.rs"]
mod tests;
