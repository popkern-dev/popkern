//! Helper functions for topology manipulation during blend construction.
//!
//! Provides perpendicular edge splitting, consumed-face identification,
//! and fillet-face merging used by the main `apply_blend_geometry` orchestrator.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

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

    // Split perpendicular edges on face_1.
    let (v1_a, v1_b) = split_face_perp_edges(brep, face_1, fillet_edge, v_start, v_end, radius)?;

    // Split perpendicular edges on face_2.
    let (v2_a, v2_b) = split_face_perp_edges(brep, face_2, fillet_edge, v_start, v_end, radius)?;

    // Build trim 3D curves from the new vertex positions (straight lines for now).
    let p1_a = brep.vertices[v1_a].position;
    let p1_b = brep.vertices[v1_b].position;
    let p2_a = brep.vertices[v2_a].position;
    let p2_b = brep.vertices[v2_b].position;

    let trim_3d_1 = line_nurbs(&p1_a, &p1_b);
    let trim_3d_2 = line_nurbs(&p2_a, &p2_b);

    Ok(SplitInfo {
        trim_3d_1,
        trim_3d_2,
        face_1_anchors: (v1_a, v1_b),
        face_2_anchors: (v2_a, v2_b),
    })
}

/// Information returned by perpendicular edge splitting.
#[derive(Debug)]
pub struct SplitInfo {
    /// Corrected 3D trim curve for face 1.
    pub trim_3d_1: NurbsCurve,
    /// Corrected 3D trim curve for face 2.
    pub trim_3d_2: NurbsCurve,
    /// Anchor vertices on face 1 (start, end matching fillet edge direction).
    pub face_1_anchors: (VertexHandle, VertexHandle),
    /// Anchor vertices on face 2 (start, end matching fillet edge direction).
    pub face_2_anchors: (VertexHandle, VertexHandle),
}

/// Build a simple linear NURBS curve between two points.
fn line_nurbs(a: &Point3, b: &Point3) -> NurbsCurve {
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
fn project_to_uv(surface: &SurfaceGeometry, point: &Point3) -> (f64, f64) {
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

/// For a given face, find the two perpendicular boundary edges that share
/// vertices with the fillet edge, and split them at `radius` distance from
/// the shared vertex. Returns the two new anchor vertices.
fn split_face_perp_edges(
    brep: &mut BRep,
    face: FaceHandle,
    fillet_edge: EdgeHandle,
    v_start: VertexHandle,
    v_end: VertexHandle,
    radius: f64,
) -> Result<(VertexHandle, VertexHandle), String> {
    // Walk the face loop to find edges incident to v_start and v_end
    // that are NOT the fillet edge.
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} not found"))?;
    let outer_loop = face_data.outer_loop;
    let start_he = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| format!("loop {outer_loop:?} not found"))?
        .half_edge;

    let mut perp_start: Option<(EdgeHandle, f64)> = None;
    let mut perp_end: Option<(EdgeHandle, f64)> = None;

    let mut cursor = start_he;
    let max_iters = brep.half_edges.len() + 1;
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

            if origin == v_start && edge_len > 1e-15 {
                // This edge leaves v_start. Split at radius/edge_len from origin.
                let t = (radius / edge_len).min(0.999);
                if t > 0.001 {
                    perp_start = Some((he.edge, t));
                }
            } else if dest == v_start && edge_len > 1e-15 {
                // This edge arrives at v_start. Split at 1 - radius/edge_len.
                let t = 1.0 - (radius / edge_len).min(0.999);
                if t > 0.001 && t < 0.999 {
                    perp_start = Some((he.edge, t));
                }
            }

            if origin == v_end && edge_len > 1e-15 {
                let t = (radius / edge_len).min(0.999);
                if t > 0.001 {
                    perp_end = Some((he.edge, t));
                }
            } else if dest == v_end && edge_len > 1e-15 {
                let t = 1.0 - (radius / edge_len).min(0.999);
                if t > 0.001 && t < 0.999 {
                    perp_end = Some((he.edge, t));
                }
            }
        }

        cursor = he.next;
        iters += 1;
        if cursor == start_he || iters >= max_iters {
            break;
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

    Ok((split_a.vertex, split_b.vertex))
}

/// Determine which of the two split faces is the "consumed" face.
///
/// The consumed face is the one whose outer loop contains the original
/// fillet edge. Returns `(kept_face, consumed_face)`.
pub fn identify_consumed_face(
    brep: &BRep,
    face_a: FaceHandle,
    face_b: FaceHandle,
    fillet_edge: EdgeHandle,
) -> Result<(FaceHandle, FaceHandle), String> {
    if face_loop_contains_edge(brep, face_b, fillet_edge)? {
        Ok((face_a, face_b))
    } else if face_loop_contains_edge(brep, face_a, fillet_edge)? {
        Ok((face_b, face_a))
    } else {
        Err(format!(
            "neither split face contains fillet edge {fillet_edge:?}"
        ))
    }
}

/// Check if a face's outer loop contains a half-edge belonging to the given edge.
fn face_loop_contains_edge(
    brep: &BRep,
    face: FaceHandle,
    edge: EdgeHandle,
) -> Result<bool, String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} not found"))?;
    let outer_loop = face_data.outer_loop;
    let start_he = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| format!("loop {outer_loop:?} not found"))?
        .half_edge;

    let mut cursor = start_he;
    let max_iters = brep.half_edges.len() + 1;
    let mut iters = 0;
    loop {
        let he = brep
            .half_edges
            .get(cursor)
            .ok_or_else(|| format!("half-edge {cursor:?} not found"))?;
        if he.edge == edge {
            return Ok(true);
        }
        cursor = he.next;
        iters += 1;
        if cursor == start_he || iters >= max_iters {
            break;
        }
    }
    Ok(false)
}

/// Find the half-edge of the given edge that lies in the given face's loop.
pub fn find_edge_halfedge_on_face(
    brep: &BRep,
    edge: EdgeHandle,
    face: FaceHandle,
) -> Result<HalfEdgeHandle, String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} not found"))?;
    let outer_loop = face_data.outer_loop;
    let start_he = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| format!("loop {outer_loop:?} not found"))?
        .half_edge;

    let mut cursor = start_he;
    let max_iters = brep.half_edges.len() + 1;
    let mut iters = 0;
    loop {
        let he = brep
            .half_edges
            .get(cursor)
            .ok_or_else(|| format!("half-edge {cursor:?} not found"))?;
        if he.edge == edge {
            return Ok(cursor);
        }
        cursor = he.next;
        iters += 1;
        if cursor == start_he || iters >= max_iters {
            break;
        }
    }
    Err(format!(
        "edge {edge:?} not found in face {face:?}'s outer loop"
    ))
}

/// Merge two consumed faces by removing the original fillet edge via `kef`.
///
/// After merging, the surviving face is assigned the fillet surface geometry.
/// Returns the handle of the merged (fillet) face.
pub fn merge_consumed_faces(
    brep: &mut BRep,
    consumed_1: FaceHandle,
    consumed_2: FaceHandle,
    fillet_edge: EdgeHandle,
    surface: &SurfaceGeometry,
) -> Result<FaceHandle, String> {
    // Find the half-edge of the fillet edge that lies on consumed_1's loop.
    let he_on_consumed_1 = find_edge_halfedge_on_face(brep, fillet_edge, consumed_1)?;

    // Verify the twin is on consumed_2 (sanity check before kef).
    let twin_handle = brep.half_edges[he_on_consumed_1].twin;
    let twin_face = brep.loops[brep.half_edges[twin_handle].loop_].face;
    if twin_face != consumed_2 {
        return Err(format!(
            "fillet edge twin is on face {twin_face:?}, expected {consumed_2:?}"
        ));
    }

    // kef removes the edge and merges the two faces.
    // It kills the twin's loop/face (consumed_2) and keeps consumed_1.
    brep.kef(he_on_consumed_1)
        .map_err(|e| format!("kef failed: {e}"))?;

    // After kef, consumed_2 is removed and consumed_1 survives with the merged loop.
    brep.faces[consumed_1].surface = Some(SurfaceOnFace {
        geometry: surface.clone(),
        orientation: true,
    });

    Ok(consumed_1)
}

/// Extract the fillet radius from the edge geometry analysis.
///
/// Uses the distance between a trim curve endpoint and the nearest fillet
/// edge vertex as a fallback radius estimate.
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
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    #[test]
    fn identify_consumed_face_finds_correct_side() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let edge = cube.edges.keys().next().unwrap();

        let edge_data = &cube.edges[edge];
        let he = &cube.half_edges[edge_data.half_edge];
        let face_of_he = cube.loops[he.loop_].face;

        assert!(face_loop_contains_edge(&cube, face_of_he, edge).unwrap());
    }

    #[test]
    fn find_edge_halfedge_on_face_works() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let edge = cube.edges.keys().next().unwrap();

        let edge_data = &cube.edges[edge];
        let he = &cube.half_edges[edge_data.half_edge];
        let face_of_he = cube.loops[he.loop_].face;

        let found = find_edge_halfedge_on_face(&cube, edge, face_of_he).unwrap();
        assert_eq!(cube.half_edges[found].edge, edge);
    }

    #[test]
    fn split_perp_edges_on_cube() {
        let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let edge = cube.edges.keys().next().unwrap();

        // Find the two adjacent faces.
        let edge_data = &cube.edges[edge];
        let he = &cube.half_edges[edge_data.half_edge];
        let face_1 = cube.loops[he.loop_].face;
        let twin = &cube.half_edges[he.twin];
        let face_2 = cube.loops[twin.loop_].face;

        let result = split_perpendicular_edges_by_radius(
            &mut cube,
            face_1,
            face_2,
            edge,
            0.1,
            &Tolerance::default(),
        );
        assert!(
            result.is_ok(),
            "split perp edges failed: {:?}",
            result.err()
        );
    }

    #[test]
    fn line_nurbs_produces_valid_curve() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let curve = line_nurbs(&a, &b);
        let (t0, t1) = curve.parameter_range();
        let p0 = curve.point_at(t0);
        let p1 = curve.point_at(t1);
        assert!(p0.distance_to(&a) < 1e-10);
        assert!(p1.distance_to(&b) < 1e-10);
    }
}
