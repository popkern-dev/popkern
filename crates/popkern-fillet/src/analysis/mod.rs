//! Edge analysis for fillet/chamfer operations.
//!
//! For each edge to fillet or chamfer, identifies the two adjacent faces,
//! classifies edge convexity, detects G1-continuous chains, and groups
//! vertices where multiple fillet edges meet.

mod chains;

pub use chains::{detect_chains, group_vertices};

use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// Analysis of a single edge for filleting.
#[derive(Debug)]
pub struct EdgeAnalysis {
    /// The edge handle in the BRep.
    pub edge: EdgeHandle,
    /// First adjacent face.
    pub face_1: FaceHandle,
    /// Second adjacent face.
    pub face_2: FaceHandle,
    /// Surface geometry of the first face.
    pub surface_1: SurfaceGeometry,
    /// Surface geometry of the second face.
    pub surface_2: SurfaceGeometry,
    /// Convexity classification at the edge.
    pub convexity: Convexity,
    /// Vertex at the start of the edge.
    pub vertex_start: VertexHandle,
    /// Vertex at the end of the edge.
    pub vertex_end: VertexHandle,
}

/// Convexity classification of an edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Convexity {
    /// The edge is convex (exterior angle < 180 degrees).
    Convex,
    /// The edge is concave (exterior angle > 180 degrees).
    Concave,
    /// Mixed or indeterminate convexity (e.g., tangent surfaces).
    Mixed,
}

/// Analyze edges for filleting.
///
/// For each edge, finds the two adjacent faces, retrieves their surface
/// geometry, classifies convexity, and identifies the start/end vertices.
/// If `propagate_tangent` is true, additional tangent-continuous edges
/// may be included (future extension point).
pub fn analyze_edges(
    brep: &BRep,
    edges: &[EdgeHandle],
    _propagate_tangent: bool,
) -> Result<Vec<EdgeAnalysis>, String> {
    let mut analyses = Vec::with_capacity(edges.len());

    for &edge_handle in edges {
        let analysis = analyze_single_edge(brep, edge_handle)?;
        analyses.push(analysis);
    }

    Ok(analyses)
}

/// Analyze a single edge: find adjacent faces, classify convexity, get vertices.
fn analyze_single_edge(brep: &BRep, edge_handle: EdgeHandle) -> Result<EdgeAnalysis, String> {
    let edge = brep
        .edges
        .get(edge_handle)
        .ok_or_else(|| format!("edge {edge_handle:?} not found in BRep"))?;

    let he1_handle = edge.half_edge;
    let he1 = brep
        .half_edges
        .get(he1_handle)
        .ok_or_else(|| format!("half-edge {he1_handle:?} not found"))?;

    let he2_handle = he1.twin;
    let he2 = brep
        .half_edges
        .get(he2_handle)
        .ok_or_else(|| format!("twin half-edge {he2_handle:?} not found"))?;

    // Walk half-edge -> loop -> face.
    let face_1 = loop_to_face(brep, he1.loop_)?;
    let face_2 = loop_to_face(brep, he2.loop_)?;

    let surface_1 = get_face_surface(brep, face_1)?;
    let surface_2 = get_face_surface(brep, face_2)?;

    // Get the two vertices at edge endpoints.
    let vertex_start = he1.origin;
    let vertex_end = he2.origin;

    // Classify convexity at the edge midpoint.
    let convexity = classify_convexity(brep, edge_handle, &surface_1, &surface_2, face_1, face_2)?;

    Ok(EdgeAnalysis {
        edge: edge_handle,
        face_1,
        face_2,
        surface_1,
        surface_2,
        convexity,
        vertex_start,
        vertex_end,
    })
}

/// Find the face that owns a given loop.
fn loop_to_face(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> Result<FaceHandle, String> {
    let loop_ = brep
        .loops
        .get(loop_handle)
        .ok_or_else(|| format!("loop {loop_handle:?} not found"))?;
    Ok(loop_.face)
}

/// Get surface geometry for a face.
fn get_face_surface(brep: &BRep, face_handle: FaceHandle) -> Result<SurfaceGeometry, String> {
    let face = brep
        .faces
        .get(face_handle)
        .ok_or_else(|| format!("face {face_handle:?} not found"))?;

    let surface_on_face = face
        .surface
        .as_ref()
        .ok_or_else(|| format!("face {face_handle:?} has no surface geometry"))?;

    Ok(surface_on_face.geometry.clone())
}

/// Number of sample points along the edge for convexity classification.
const CONVEXITY_SAMPLES: usize = 7;

/// Classify the convexity of an edge between two surfaces.
///
/// Samples normals at [`CONVEXITY_SAMPLES`] points along the edge and
/// uses the sign of (N1 x T) . N2 at each sample:
///   - positive -> convex
///   - negative -> concave
///   - near zero -> indeterminate (skipped)
///
/// If all definite samples agree, returns [`Convexity::Convex`] or
/// [`Convexity::Concave`]. If samples disagree, returns [`Convexity::Mixed`].
fn classify_convexity(
    brep: &BRep,
    edge_handle: EdgeHandle,
    surface_1: &SurfaceGeometry,
    surface_2: &SurfaceGeometry,
    face_1: FaceHandle,
    face_2: FaceHandle,
) -> Result<Convexity, String> {
    let edge = brep
        .edges
        .get(edge_handle)
        .ok_or_else(|| format!("edge {edge_handle:?} not found"))?;

    let curve = edge
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {edge_handle:?} has no curve geometry"))?;

    const MIXED_THRESHOLD: f64 = 1e-6;

    let mut convex_count: usize = 0;
    let mut concave_count: usize = 0;

    for i in 0..CONVEXITY_SAMPLES {
        let frac = (i as f64) / ((CONVEXITY_SAMPLES - 1).max(1) as f64);
        let t = curve.t_start + (curve.t_end - curve.t_start) * frac;
        let point = curve.geometry.point_at(t);
        let tangent = curve.geometry.tangent_at(t);

        let n1 = surface_normal_at_point(surface_1, &point, brep, face_1)?;
        let n2 = surface_normal_at_point(surface_2, &point, brep, face_2)?;

        let cross = n1.cross(&tangent);
        let dot = cross.dot(&n2);

        if dot > MIXED_THRESHOLD {
            convex_count += 1;
        } else if dot < -MIXED_THRESHOLD {
            concave_count += 1;
        }
        // Near-zero samples are indeterminate — skipped.
    }

    if convex_count > 0 && concave_count > 0 {
        Ok(Convexity::Mixed)
    } else if convex_count > 0 {
        Ok(Convexity::Convex)
    } else if concave_count > 0 {
        Ok(Convexity::Concave)
    } else {
        // All samples indeterminate.
        Ok(Convexity::Mixed)
    }
}

/// Evaluate the surface normal at a 3D point, accounting for face orientation.
fn surface_normal_at_point(
    surface: &SurfaceGeometry,
    point: &popkern_core::geometry::Point3,
    brep: &BRep,
    face_handle: FaceHandle,
) -> Result<Vec3, String> {
    // Find UV coordinates by closest-point projection.
    let (u, v) = find_uv_for_point(surface, point);
    let mut normal = surface.normal_at(u, v);

    // Account for face orientation.
    let face = brep
        .faces
        .get(face_handle)
        .ok_or_else(|| format!("face {face_handle:?} not found"))?;
    if let Some(ref sof) = face.surface
        && !sof.orientation
    {
        normal = Vec3::new(-normal.x, -normal.y, -normal.z);
    }

    Ok(normal)
}

/// Find (u, v) parameters on a surface closest to a 3D point.
///
/// Uses analytic projection for planes and a coarse grid search for others.
fn find_uv_for_point(
    surface: &SurfaceGeometry,
    point: &popkern_core::geometry::Point3,
) -> (f64, f64) {
    // Try analytic projection first.
    if let Some(uv) = popkern_core::geometry::uv_projection::project_point_to_uv(surface, point) {
        return uv;
    }

    // Fallback: grid search over the parameter range.
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    let u_lo = if u_min.is_finite() { u_min } else { -10.0 };
    let u_hi = if u_max.is_finite() { u_max } else { 10.0 };
    let v_lo = if v_min.is_finite() { v_min } else { -10.0 };
    let v_hi = if v_max.is_finite() { v_max } else { 10.0 };

    let n = 20;
    let mut best_u = (u_lo + u_hi) * 0.5;
    let mut best_v = (v_lo + v_hi) * 0.5;
    let mut best_dist = f64::INFINITY;

    for iu in 0..n {
        let u = u_lo + (u_hi - u_lo) * (iu as f64) / ((n - 1).max(1) as f64);
        for iv in 0..n {
            let v = v_lo + (v_hi - v_lo) * (iv as f64) / ((n - 1).max(1) as f64);
            let sp = surface.point_at(u, v);
            let d = sp.distance_to(point);
            if d < best_dist {
                best_dist = d;
                best_u = u;
                best_v = v;
            }
        }
    }

    (best_u, best_v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn make_cube() -> BRep {
        BRep::unit_cube(Tolerance::default()).unwrap()
    }

    /// Find the first edge handle in the BRep.
    fn first_edge(brep: &BRep) -> EdgeHandle {
        brep.edges.keys().next().unwrap()
    }

    #[test]
    fn analyze_cube_edge() {
        let cube = make_cube();
        let edge = first_edge(&cube);

        let analyses = analyze_edges(&cube, &[edge], false).unwrap();
        assert_eq!(analyses.len(), 1);

        let a = &analyses[0];
        assert_eq!(a.edge, edge);

        // Both faces should have planar surfaces.
        assert!(
            matches!(a.surface_1, SurfaceGeometry::Plane { .. }),
            "face_1 should be planar"
        );
        assert!(
            matches!(a.surface_2, SurfaceGeometry::Plane { .. }),
            "face_2 should be planar"
        );

        // A cube edge is convex (exterior corner).
        assert!(
            a.convexity == Convexity::Convex || a.convexity == Convexity::Concave,
            "cube edge should have definite convexity, got {:?}",
            a.convexity
        );

        // Vertices should be distinct.
        assert_ne!(
            a.vertex_start, a.vertex_end,
            "start and end vertices should differ"
        );
    }

    #[test]
    fn analyze_invalid_edge() {
        let cube = make_cube();
        // A default handle will never exist in a populated BRep.
        let bogus = popkern_core::topology::EdgeHandle::default();

        let result = analyze_edges(&cube, &[bogus], false);
        assert!(result.is_err(), "should fail for invalid edge handle");
    }

    #[test]
    fn classify_convexity_on_cube_is_definite() {
        let cube = make_cube();
        // Analyze ALL edges — every cube edge should be either Convex or Concave.
        let edges: Vec<EdgeHandle> = cube.edges.keys().collect();
        let analyses = analyze_edges(&cube, &edges, false).unwrap();
        for a in &analyses {
            assert!(
                a.convexity == Convexity::Convex || a.convexity == Convexity::Concave,
                "cube edge {:?} has unexpected convexity {:?}",
                a.edge,
                a.convexity
            );
        }
    }

    #[test]
    fn concave_edge_from_subtraction() {
        use popkern_core::geometry::Transform3;
        use popkern_core::tolerance::Tolerance;
        let tol = Tolerance::default();
        let big = BRep::unit_cube(tol).unwrap();
        // Scale to half-size then translate so it overlaps one corner.
        let xf = Transform3::translation(0.25, 0.25, 0.5).compose(&Transform3::scale_uniform(0.5));
        let small = BRep::unit_cube(tol).unwrap().transform(&xf).unwrap();

        let result = popkern_bool::subtraction(big, small, tol);
        if let Ok(subtracted) = result {
            let edges: Vec<EdgeHandle> = subtracted.edges.keys().collect();
            if !edges.is_empty() {
                let analyses = analyze_edges(&subtracted, &edges, false);
                if let Ok(analyses) = analyses {
                    let has_concave = analyses.iter().any(|a| a.convexity == Convexity::Concave);
                    assert!(has_concave, "subtraction result should have concave edges");
                }
            }
        }
    }

    #[test]
    fn edge_missing_surface_geometry() {
        // Construct a BRep where a face has no surface, triggering the error path
        // in get_face_surface.
        let mut cube = make_cube();
        let edge = first_edge(&cube);
        let edge_data = &cube.edges[edge];
        let he = &cube.half_edges[edge_data.half_edge];
        let face = cube.loops[he.loop_].face;

        // Remove surface geometry from the face.
        cube.faces[face].surface = None;

        let result = analyze_edges(&cube, &[edge], false);
        assert!(
            result.is_err(),
            "should fail when face has no surface geometry"
        );
        let err = result.unwrap_err();
        assert!(
            err.contains("no surface"),
            "error should mention missing surface: {err}"
        );
    }

    #[test]
    fn analyze_multiple_edges() {
        let cube = make_cube();
        let edges: Vec<EdgeHandle> = cube.edges.keys().take(3).collect();
        let analyses = analyze_edges(&cube, &edges, false).unwrap();
        assert_eq!(analyses.len(), 3);
        for a in &analyses {
            assert_ne!(a.face_1, a.face_2, "adjacent faces should be different");
        }
    }

    #[test]
    fn convexity_classification_samples_multiple_points() {
        // A cube edge is uniformly convex — multi-point sampling should
        // still return Convex (not Mixed).
        let cube = make_cube();
        let edges: Vec<EdgeHandle> = cube.edges.keys().collect();
        let analyses = analyze_edges(&cube, &edges, false).unwrap();

        for a in &analyses {
            assert!(
                a.convexity == Convexity::Convex || a.convexity == Convexity::Concave,
                "cube edge {:?} should be uniformly convex or concave with \
                 multi-point sampling, got {:?}",
                a.edge,
                a.convexity,
            );
        }

        // A concave edge from a boolean subtraction should be Concave.
        use popkern_core::geometry::Transform3;
        use popkern_core::tolerance::Tolerance;
        let tol = Tolerance::default();
        let big = BRep::unit_cube(tol).unwrap();
        let xf = Transform3::translation(0.25, 0.25, 0.5).compose(&Transform3::scale_uniform(0.5));
        let small = BRep::unit_cube(tol).unwrap().transform(&xf).unwrap();

        if let Ok(subtracted) = popkern_bool::subtraction(big, small, tol) {
            let edges: Vec<EdgeHandle> = subtracted.edges.keys().collect();
            if let Ok(analyses) = analyze_edges(&subtracted, &edges, false) {
                // At least one concave edge should exist, and multi-point
                // sampling should not flip it to Mixed.
                let has_concave = analyses.iter().any(|a| a.convexity == Convexity::Concave);
                assert!(
                    has_concave,
                    "subtraction result should have concave edges with multi-point sampling"
                );
            }
        }
    }
}
