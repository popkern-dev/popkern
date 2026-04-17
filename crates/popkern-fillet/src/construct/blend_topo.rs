//! Topology walking helpers for vertex blend insertion.
//!
//! Provides functions to walk around a vertex in the BRep topology,
//! collect adjacent faces, anchor vertices, and incident edges.
//! These are used by [`super::blend_insert`] during blend face insertion.

use std::collections::HashSet;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::construct::EdgeGeometry;
use crate::fillet_surface::FilletSurfaceResult;

/// Repair a vertex's outgoing half-edge if it has become stale.
///
/// After edge splitting and face splitting, a vertex's `outgoing` field
/// may point to a half-edge that was removed. This scans all half-edges
/// to find one originating from the vertex and updates the field.
pub(crate) fn repair_vertex_outgoing(brep: &mut BRep, vertex: VertexHandle) -> Result<(), String> {
    let v = brep
        .vertices
        .get(vertex)
        .ok_or_else(|| format!("vertex {vertex:?} not found"))?;

    if brep.half_edges.contains_key(v.outgoing) {
        return Ok(()); // Already valid.
    }

    // Scan for any half-edge originating from this vertex.
    let mut found = None;
    for (he_handle, he) in &brep.half_edges {
        if he.origin == vertex {
            found = Some(he_handle);
            break;
        }
    }

    let new_outgoing = found.ok_or_else(|| {
        let count = brep
            .half_edges
            .values()
            .filter(|he| he.origin == vertex)
            .count();
        format!("vertex {vertex:?} has no half-edges originating from it (scan found {count})")
    })?;
    brep.vertices[vertex].outgoing = new_outgoing;
    Ok(())
}

/// Collect faces around a vertex with the half-edge arriving at V.
///
/// Returns `(face, he_arriving_at_v)` for each face incident to V.
pub(crate) fn collect_faces_around_vertex(
    brep: &BRep,
    vertex: VertexHandle,
) -> Result<Vec<(FaceHandle, HalfEdgeHandle)>, String> {
    let outgoing = brep
        .vertices
        .get(vertex)
        .ok_or_else(|| format!("vertex {vertex:?} not found"))?
        .outgoing;

    if !brep.half_edges.contains_key(outgoing) {
        return Err(format!("vertex {vertex:?} has no valid outgoing half-edge"));
    }

    let mut faces = Vec::new();
    let mut seen_faces = HashSet::new();
    let mut cursor = outgoing;
    let max_iters = brep.half_edges.len() + 1;
    let mut iters = 0;

    loop {
        let he = brep
            .half_edges
            .get(cursor)
            .ok_or_else(|| format!("half-edge {cursor:?} not found"))?;

        let face = brep.loops[he.loop_].face;
        if seen_faces.insert(face) {
            faces.push((face, he.prev));
        }

        cursor = brep.half_edges[he.twin].next;
        iters += 1;
        if cursor == outgoing || iters >= max_iters {
            break;
        }
    }

    Ok(faces)
}

/// Compute the V-side trim-endpoint positions that lie on `face`'s
/// supporting plane. Each incident fillet contributes exactly one
/// such point (its trim curve that's co-planar with `face`, closest
/// point to V). For a trihedral corner's planar adjacent face this
/// yields 2 target positions — the 2 fillets bounding that face.
///
/// **Why not distance-based matching?** A vertex on a cube edge
/// at distance `r` from V matches *some* fillet's radius but is
/// not necessarily a trim endpoint of any fillet incident to the
/// current face. Example: on the mixed-radius cube's `y=0` corner
/// region, `V14 = (0.10, 0, 0)` sits at distance 0.10 (matching
/// fillet-1's radius) but fillet-1 bounds `z=0` and `x=0`, not
/// `y=0`. Distance-based matching picked V14 as a y=0 anchor and
/// cut the wrong chord. The plane-containment test filters on the
/// actual face the trim lies on.
/// Radii of fillets whose trim curves lie on `face`'s plane — i.e.
/// the fillets actually bordering this face at `vertex`.
///
/// A trim curve is *on* `face` iff all its sample points lie on
/// `face`'s plane (mean plane-distance < `tol_plane`). This is far
/// more robust than checking a single closest-to-V point whose
/// sampling error bleeds into the face normal direction. Mean-of-8
/// averages out sampling noise, while a trim curve that truly does
/// NOT lie on this plane shows a mean distance ≥ the dihedral
/// angle's offset (much larger than sampling noise).
pub(crate) fn collect_radii_for_face(
    brep: &BRep,
    face: FaceHandle,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Vec<f64> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let Some(face_data) = brep.faces.get(face) else {
        return Vec::new();
    };
    let Some(sof) = face_data.surface.as_ref() else {
        return Vec::new();
    };
    let (plane_origin, plane_normal) = match &sof.geometry {
        SurfaceGeometry::Plane { origin, normal, .. } => (*origin, *normal),
        _ => return Vec::new(),
    };

    let tol_plane = (tolerance.linear() * 100.0).max(1e-4);

    let mut radii = Vec::with_capacity(edge_geometries.len());
    for eg in edge_geometries {
        for trim in [&eg.trim_3d_1, &eg.trim_3d_2] {
            let (t0, t1) = trim.parameter_range();
            let n_samples = 8usize;
            let mut mean_dist = 0.0_f64;
            for i in 0..n_samples {
                let t = t0 + (t1 - t0) * (i as f64 + 0.5) / (n_samples as f64);
                let p = trim.point_at(t);
                let delta = p - plane_origin;
                mean_dist += delta.dot(&plane_normal).abs();
            }
            mean_dist /= n_samples as f64;

            if mean_dist < tol_plane {
                radii.push(eg.effective_radius);
                break; // At most one of eg's two trims lies on `face`.
            }
        }
    }
    radii
}

/// Find the per-face fillet-trim-endpoint anchors for a face
/// incident to `vertex`. For every fillet whose trim curves lie on
/// `face`'s plane (identified by [`collect_radii_for_face`]), picks
/// the face-loop vertex at distance `r` from V. Matching by
/// (plane-filtered radius, face-loop membership) together is what
/// makes this robust: a vertex at distance `r` from V that's *on*
/// `face` but NOT actually the trim endpoint of a fillet on
/// `face` (e.g. a cube-sub perp-split from a different fillet's
/// radius match) is excluded because the matching radius only
/// lists fillets bordering `face`.
///
/// For a trihedral corner's planar adjacent face the return is
/// exactly two vertices — the two trim endpoints of the two
/// fillets that bound this face.
pub(crate) fn find_face_trim_anchors(
    brep: &BRep,
    face: FaceHandle,
    vertex: VertexHandle,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Result<Vec<VertexHandle>, String> {
    let radii = collect_radii_for_face(brep, face, edge_geometries, tolerance);
    if radii.is_empty() {
        return Ok(Vec::new());
    }

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
    let v_pos = brep
        .vertices
        .get(vertex)
        .ok_or_else(|| format!("vertex {vertex:?} not found"))?
        .position;

    // Match band: absolute BRep tolerance + 0.1%-of-radius
    // padding. Tight enough to exclude chord-crossing points at
    // `√(r_a² + r_b²) ≥ r·√2` — the gap between a true perp-split
    // vertex at `r` and the nearest crossing point is O(r),
    // orders of magnitude larger than this band.
    let min_r = radii
        .iter()
        .copied()
        .filter(|r| r.is_finite() && *r > 0.0)
        .fold(f64::INFINITY, f64::min);
    let band = (tolerance.linear() * 100.0).max(min_r * 1e-3);

    let mut anchors: Vec<VertexHandle> = Vec::new();
    let mut cursor = start_he;
    let max_iters = brep.half_edges.len() + 1;
    for _ in 0..max_iters {
        let vh = brep.half_edges[cursor].origin;
        if vh != vertex {
            let vpos = brep.vertices[vh].position;
            let d = v_pos.distance_to(&vpos);
            if radii.iter().any(|&r| (d - r).abs() < band) && !anchors.contains(&vh) {
                anchors.push(vh);
            }
        }
        cursor = brep.half_edges[cursor].next;
        if cursor == start_he {
            break;
        }
    }

    Ok(anchors)
}

/// Find a half-edge originating from the given vertex in a face's loop.
pub(crate) fn find_halfedge_from_vertex_in_face(
    brep: &BRep,
    face: FaceHandle,
    vertex: VertexHandle,
) -> Result<HalfEdgeHandle, String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} not found"))?;
    let outer_loop = face_data.outer_loop;
    let start_he = brep.loops[outer_loop].half_edge;

    let mut cursor = start_he;
    let max_iters = brep.half_edges.len() + 1;
    let mut iters = 0;
    loop {
        if brep.half_edges[cursor].origin == vertex {
            return Ok(cursor);
        }
        cursor = brep.half_edges[cursor].next;
        iters += 1;
        if cursor == start_he || iters >= max_iters {
            break;
        }
    }

    Err(format!(
        "no half-edge from vertex {vertex:?} found in face {face:?}"
    ))
}

/// Determine which of two faces is the corner triangle containing the vertex.
pub(crate) fn identify_corner_face(
    brep: &BRep,
    face_a: FaceHandle,
    face_b: FaceHandle,
    vertex: VertexHandle,
) -> Result<FaceHandle, String> {
    if face_loop_contains_vertex(brep, face_b, vertex)? {
        Ok(face_b)
    } else if face_loop_contains_vertex(brep, face_a, vertex)? {
        Ok(face_a)
    } else {
        Err(format!(
            "neither face {face_a:?} nor {face_b:?} contains vertex {vertex:?}"
        ))
    }
}

/// Check if a face's outer loop contains a vertex.
pub(crate) fn face_loop_contains_vertex(
    brep: &BRep,
    face: FaceHandle,
    vertex: VertexHandle,
) -> Result<bool, String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} not found"))?;
    let start_he = brep.loops[face_data.outer_loop].half_edge;

    let mut cursor = start_he;
    let max_iters = brep.half_edges.len() + 1;
    let mut iters = 0;
    loop {
        if brep.half_edges[cursor].origin == vertex {
            return Ok(true);
        }
        cursor = brep.half_edges[cursor].next;
        iters += 1;
        if cursor == start_he || iters >= max_iters {
            break;
        }
    }
    Ok(false)
}

/// Collect all edges incident to a vertex.
pub(crate) fn collect_edges_through_vertex(
    brep: &BRep,
    vertex: VertexHandle,
) -> Result<Vec<EdgeHandle>, String> {
    let outgoing = brep
        .vertices
        .get(vertex)
        .ok_or_else(|| format!("vertex {vertex:?} not found"))?
        .outgoing;

    if !brep.half_edges.contains_key(outgoing) {
        return Err(format!("vertex {vertex:?} has no valid outgoing half-edge"));
    }

    let mut edges = Vec::new();
    let mut seen = HashSet::new();
    let mut cursor = outgoing;
    let max_iters = brep.half_edges.len() + 1;
    let mut iters = 0;

    loop {
        let he = brep
            .half_edges
            .get(cursor)
            .ok_or_else(|| format!("half-edge {cursor:?} not found"))?;

        if seen.insert(he.edge) {
            edges.push(he.edge);
        }

        cursor = brep.half_edges[he.twin].next;
        iters += 1;
        if cursor == outgoing || iters >= max_iters {
            break;
        }
    }

    Ok(edges)
}

/// Build `FilletSurfaceResult` list from edge geometries for blend boundary
/// computation.
pub(crate) fn build_fillet_surface_list(
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
) -> Vec<FilletSurfaceResult> {
    edge_indices
        .iter()
        .filter_map(|&idx| {
            let eg = edge_geometries.get(idx)?;
            Some(FilletSurfaceResult {
                surface: eg.surface.clone(),
                trim_curve_1: eg.trim_3d_1.clone(),
                trim_curve_2: eg.trim_3d_2.clone(),
                continuity_downgrade: None,
            })
        })
        .collect()
}

#[cfg(test)]
#[path = "blend_topo_tests.rs"]
mod tests;
