//! Edge analysis for fillet/chamfer operations.
//!
//! For each edge to fillet or chamfer, identifies the two adjacent faces,
//! classifies edge convexity, detects G1-continuous chains, and groups
//! vertices where multiple fillet edges meet.

mod chains;
mod convexity;

pub use chains::{detect_chains, group_vertices};

use popkern_core::geometry::Vec3;
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
    /// Face orientation for surface 1. When `false`, the effective face
    /// normal is opposite to the surface geometry's normal.
    pub orientation_1: bool,
    /// Face orientation for surface 2.
    pub orientation_2: bool,
    /// Whether the shell containing these faces is globally
    /// OUTWARD-oriented (measured via signed volume). Primitives
    /// build inward shells; boolean/sweep output builds outward
    /// ones. True face outwardness = `orientation_i == shell_outward`.
    pub shell_outward: bool,
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

    // Seam edges have the same face on both sides — cannot be filleted.
    if face_1 == face_2 {
        return Err(format!(
            "edge {edge_handle:?} is a seam edge (same face on both sides) and cannot be filleted"
        ));
    }

    let (surface_1, stored_orientation_1) = get_face_surface_with_orientation(brep, face_1)?;
    let (surface_2, stored_orientation_2) = get_face_surface_with_orientation(brep, face_2)?;
    // EFFECTIVE orientation: does the surface's natural normal point
    // OUT of the material? Derived from the outer loop's UV winding
    // (outer loops wind CCW in the chart iff the natural normal is
    // outward — the manifold convention), not from the stored flag:
    // face producers disagree on the flag (fillet strips carry
    // radial-out cylinder normals into concave grooves; sweep caps
    // and boolean output vary), and the rolling-ball side depends on
    // the PHYSICAL outward direction. Falls back to the stored flag
    // when pcurves are unavailable.
    // orientation_1/2 keep legacy STORED-flag semantics (a large
    // consumer surface depends on them); shell_outward carries the
    // measured global shell orientation so ball placement can
    // compute true outwardness (stored == shell_outward).
    let shell_outward = shell_outward_sign(brep) > 0.0;
    let orientation_1 = stored_orientation_1;
    let orientation_2 = stored_orientation_2;

    // Get the two vertices at edge endpoints.
    let vertex_start = he1.origin;
    let vertex_end = he2.origin;

    // Classify convexity at the edge midpoint.
    let convexity =
        convexity::classify_convexity(brep, edge_handle, &surface_1, &surface_2, face_1, face_2)?;

    Ok(EdgeAnalysis {
        edge: edge_handle,
        face_1,
        face_2,
        surface_1,
        surface_2,
        orientation_1,
        orientation_2,
        shell_outward,
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
#[allow(dead_code)]
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

/// Get the surface geometry AND orientation flag from a face.
fn get_face_surface_with_orientation(
    brep: &BRep,
    face_handle: FaceHandle,
) -> Result<(SurfaceGeometry, bool), String> {
    let face = brep
        .faces
        .get(face_handle)
        .ok_or_else(|| format!("face {face_handle:?} not found"))?;
    let surface_on_face = face
        .surface
        .as_ref()
        .ok_or_else(|| format!("face {face_handle:?} has no surface geometry"))?;
    Ok((
        surface_on_face.geometry.clone(),
        surface_on_face.orientation,
    ))
}

/// Number of sample points along the edge for convexity classification.
/// Evaluate the surface normal at a 3D point, accounting for face orientation.
pub(super) fn surface_normal_at_point_natural(
    surface: &SurfaceGeometry,
    point: &popkern_core::geometry::Point3,
) -> Vec3 {
    let (u, v) = find_uv_for_point(surface, point);
    surface.normal_at(u, v)
}

#[allow(dead_code)]
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

/// +1.0 when the shell's loop windings + orientation flags describe
/// an OUTWARD-oriented boundary, −1.0 for inward, via the divergence
/// theorem over outer-loop polygon fans (sign-exact for sane closed
/// shells). Primitive constructors currently build inward shells
/// (extrude keeps the base winding); boolean/sweep output builds
/// outward ones — consumers needing true outwardness must fold this
/// in until constructors are normalized (PK-100).
pub(super) fn shell_outward_sign(brep: &BRep) -> f64 {
    use popkern_core::geometry::surface::SurfaceEval;
    // Divergence theorem per face using the FLAG-CORRECTED NATURAL
    // normal — deliberately NOT the loop winding: producers disagree
    // on winding-vs-flag consistency (the cube winds inward-wound
    // loops, the square pyramid outward ones, both with flag=true),
    // so the loop fan supplies only area and centroid, while the
    // direction comes from what the flags claim. V > 0 ⇔ the
    // flag-corrected normals describe an outward shell.
    let mut vol = 0.0;
    for (fh, face) in brep.faces.iter() {
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        let flag_sign = if sof.orientation { 1.0 } else { -1.0 };
        let pts = brep.face_outer_loop_positions(fh);
        if pts.len() < 3 {
            continue;
        }
        let m = pts.len() as f64;
        let centroid = popkern_core::geometry::Point3::new(
            pts.iter().map(|p| p.x).sum::<f64>() / m,
            pts.iter().map(|p| p.y).sum::<f64>() / m,
            pts.iter().map(|p| p.z).sum::<f64>() / m,
        );
        let mut area_vec = Vec3::new(0.0, 0.0, 0.0);
        for i in 1..pts.len() - 1 {
            let v1 = pts[i] - pts[0];
            let v2 = pts[i + 1] - pts[0];
            area_vec = area_vec + v1.cross(&v2) * 0.5;
        }
        let area = area_vec.length();
        if area < 1e-16 {
            continue;
        }
        let (u, v) = find_uv_for_point(&sof.geometry, &centroid);
        let n = sof.geometry.normal_at(u, v) * flag_sign;
        let x = Vec3::new(centroid.x, centroid.y, centroid.z);
        vol += n.dot(&x) * area;
    }
    if std::env::var_os("PK_CONVEXITY_TRACE").is_some() {
        eprintln!("    [shell] signed 6V = {vol:.3e}");
    }
    if vol >= 0.0 { 1.0 } else { -1.0 }
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
