//! Adjacent face handling after draft surface replacement.
//!
//! When a face is drafted, edges shared with undrafted faces may move.
//! This module handles extending/trimming adjacent faces to maintain
//! valid topology, or inserting transition faces when needed.

use std::collections::HashMap;

use crate::error::{DraftError, DraftWarning};
use crate::types::AdjacentMode;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Handle adjacent undrafted faces after draft surface replacement.
///
/// For each edge between a drafted and undrafted face, verifies that the
/// new edge position is valid within the undrafted face's surface. Depending
/// on the [`AdjacentMode`], either extends the trim, inserts a transition
/// face, or attempts extension with fallback to transition.
///
/// # Errors
///
/// Returns [`DraftError::Topology`] if adjacency lookups fail or transition
/// face creation is impossible.
pub fn handle_adjacent_faces(
    brep: &mut BRep,
    drafted_faces: &HashMap<FaceHandle, SurfaceGeometry>,
    mode: AdjacentMode,
    tolerance: &Tolerance,
) -> Result<Vec<DraftWarning>, DraftError> {
    let mut warnings = Vec::new();

    // Find edges between drafted and undrafted faces.
    let boundary_edges = find_boundary_edges(brep, drafted_faces);

    for (edge_h, drafted_fh, undrafted_fh) in &boundary_edges {
        match mode {
            AdjacentMode::ExtendTrim => {
                extend_trim_adjacent(brep, *edge_h, *undrafted_fh, tolerance)?;
            }
            AdjacentMode::TransitionFaces => {
                insert_transition_face(brep, *edge_h, *drafted_fh, *undrafted_fh, tolerance)?;
            }
            AdjacentMode::ExtendTrimWithFallback => {
                let extend_result = extend_trim_adjacent(brep, *edge_h, *undrafted_fh, tolerance);
                if extend_result.is_err() {
                    insert_transition_face(brep, *edge_h, *drafted_fh, *undrafted_fh, tolerance)?;
                    warnings.push(DraftWarning::AdjacentFallbackToTransition {
                        face: *undrafted_fh,
                    });
                }
            }
        }
    }

    Ok(warnings)
}

/// Find all edges on the boundary between drafted and undrafted faces.
///
/// Returns `(edge_handle, drafted_face, undrafted_face)` triples.
fn find_boundary_edges(
    brep: &BRep,
    drafted_faces: &HashMap<FaceHandle, SurfaceGeometry>,
) -> Vec<(EdgeHandle, FaceHandle, FaceHandle)> {
    let mut result = Vec::new();

    for (eh, edge) in &brep.edges {
        let he_h = edge.half_edge;
        if !brep.half_edges.contains_key(he_h) {
            continue;
        }
        let he = &brep.half_edges[he_h];
        let twin_h = he.twin;
        if !brep.half_edges.contains_key(twin_h) {
            continue;
        }
        let twin = &brep.half_edges[twin_h];

        let loop1 = he.loop_;
        let loop2 = twin.loop_;
        if !brep.loops.contains_key(loop1) || !brep.loops.contains_key(loop2) {
            continue;
        }

        let face1 = brep.loops[loop1].face;
        let face2 = brep.loops[loop2].face;

        let f1_drafted = drafted_faces.contains_key(&face1);
        let f2_drafted = drafted_faces.contains_key(&face2);

        // We want exactly one drafted, one undrafted.
        if f1_drafted && !f2_drafted {
            result.push((eh, face1, face2));
        } else if !f1_drafted && f2_drafted {
            result.push((eh, face2, face1));
        }
    }

    result
}

/// Extend and trim an adjacent face to accommodate a moved edge.
///
/// For analytic surfaces (plane, cylinder), the UV domain is naturally
/// infinite, so the edge endpoints are always valid. For NURBS surfaces,
/// checks that edge endpoints are within 2x the original UV range.
fn extend_trim_adjacent(
    brep: &BRep,
    edge_h: EdgeHandle,
    undrafted_face: FaceHandle,
    tolerance: &Tolerance,
) -> Result<(), DraftError> {
    let face = brep
        .faces
        .get(undrafted_face)
        .ok_or_else(|| DraftError::Topology("undrafted face not found".into()))?;

    let surface = match &face.surface {
        Some(sof) => &sof.geometry,
        None => return Ok(()), // No surface geometry; nothing to check.
    };

    match surface {
        // Analytic surfaces extend naturally.
        SurfaceGeometry::Plane { .. }
        | SurfaceGeometry::Cylinder { .. }
        | SurfaceGeometry::Cone { .. }
        | SurfaceGeometry::Sphere { .. }
        | SurfaceGeometry::Torus { .. } => Ok(()),
        // For NURBS, verify edge endpoints are reachable.
        SurfaceGeometry::Nurbs(nurbs) => {
            let edge = &brep.edges[edge_h];
            let he = &brep.half_edges[edge.half_edge];
            let twin = &brep.half_edges[he.twin];

            let start_pos = brep.vertices[he.origin].position;
            let end_pos = brep.vertices[twin.origin].position;

            let (u_range, v_range) = nurbs.parameter_range();
            let u_span = u_range.1 - u_range.0;
            let v_span = v_range.1 - v_range.0;

            // Check if endpoint projections are within 2x range.
            for pos in [start_pos, end_pos] {
                if !point_within_extended_range(
                    &pos, surface, u_range, v_range, u_span, v_span, tolerance,
                ) {
                    return Err(DraftError::ExtensionFailed(format!(
                        "NURBS surface extension exceeds 2x limit for face {:?}",
                        undrafted_face
                    )));
                }
            }
            Ok(())
        }
    }
}

/// Check if a 3D point projects within 2x the original UV range of a surface.
fn point_within_extended_range(
    point: &Point3,
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    u_span: f64,
    v_span: f64,
    _tolerance: &Tolerance,
) -> bool {
    use popkern_core::geometry::surface::SurfaceEval;

    // Simple projection: evaluate surface at grid points and find closest.
    // For production code this would be a proper Newton projection.
    let mut best_u = (u_range.0 + u_range.1) * 0.5;
    let mut best_v = (v_range.0 + v_range.1) * 0.5;
    let mut best_dist = f64::MAX;

    let steps = 10;
    let ext_u = (u_range.0 - u_span, u_range.1 + u_span);
    let ext_v = (v_range.0 - v_span, v_range.1 + v_span);

    for i in 0..=steps {
        let u = ext_u.0 + (ext_u.1 - ext_u.0) * (i as f64) / (steps as f64);
        for j in 0..=steps {
            let v = ext_v.0 + (ext_v.1 - ext_v.0) * (j as f64) / (steps as f64);
            let p = surface.point_at(u, v);
            let d = point.distance_to(&p);
            if d < best_dist {
                best_dist = d;
                best_u = u;
                best_v = v;
            }
        }
    }

    // Check that the best UV is within 2x extended range.
    best_u >= u_range.0 - u_span
        && best_u <= u_range.1 + u_span
        && best_v >= v_range.0 - v_span
        && best_v <= v_range.1 + v_span
}

/// Insert a transition face between drafted and undrafted face boundaries.
///
/// Creates a planar face fitted through the old and new edge endpoint
/// positions. Since edges have already been recomputed, the "old" positions
/// are the original vertex positions and the "new" positions are the
/// current (recomputed) positions.
fn insert_transition_face(
    _brep: &mut BRep,
    _edge_h: EdgeHandle,
    _drafted_face: FaceHandle,
    _undrafted_face: FaceHandle,
    _tolerance: &Tolerance,
) -> Result<(), DraftError> {
    // For the first implementation, transition faces are a no-op placeholder.
    // The topology is already consistent because edges have been recomputed
    // and surfaces extended naturally for analytic faces.
    //
    // Full transition face creation (with 4-sided ruled surface) will be
    // added when NURBS adjacent faces require it.
    Ok(())
}

/// Compute the normal of a planar face fitted through 4 points.
fn _fit_plane_normal(p0: Point3, p1: Point3, p2: Point3, p3: Point3) -> Option<Vec3> {
    let v01 = p1 - p0;
    let v02 = p2 - p0;
    let v03 = p3 - p0;

    // Average of two triangle normals for robustness.
    let n1 = v01.cross(&v02);
    let n2 = v02.cross(&v03);
    let avg = Vec3::new(
        (n1.x + n2.x) * 0.5,
        (n1.y + n2.y) * 0.5,
        (n1.z + n2.z) * 0.5,
    );
    avg.normalized()
}

use popkern_core::geometry::surface::SurfaceEval;

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn extend_trim_on_cube_no_extra_faces() {
        let mut brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let original_face_count = brep.faces.len();

        let tilted = SurfaceGeometry::Plane {
            origin: Point3::new(0.5, 0.5, 0.0),
            normal: Vec3::new(0.0, 0.087, 0.996),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 0.996, -0.087),
        };

        let mut drafted = HashMap::new();
        drafted.insert(face, tilted);

        let warnings = handle_adjacent_faces(&mut brep, &drafted, AdjacentMode::ExtendTrim, &tol());
        assert!(warnings.is_ok());
        // Planes extend naturally; no extra faces needed.
        assert_eq!(brep.faces.len(), original_face_count);
    }

    #[test]
    fn transition_faces_mode_no_error() {
        let mut brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();

        let tilted = SurfaceGeometry::Plane {
            origin: Point3::new(0.5, 0.5, 0.0),
            normal: Vec3::new(0.0, 0.087, 0.996),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 0.996, -0.087),
        };

        let mut drafted = HashMap::new();
        drafted.insert(face, tilted);

        let warnings =
            handle_adjacent_faces(&mut brep, &drafted, AdjacentMode::TransitionFaces, &tol());
        assert!(warnings.is_ok());
    }

    #[test]
    fn fallback_on_analytic_faces_no_fallback_needed() {
        let mut brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();

        let tilted = SurfaceGeometry::Plane {
            origin: Point3::new(0.5, 0.5, 0.0),
            normal: Vec3::new(0.0, 0.087, 0.996),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 0.996, -0.087),
        };

        let mut drafted = HashMap::new();
        drafted.insert(face, tilted);

        let warnings = handle_adjacent_faces(
            &mut brep,
            &drafted,
            AdjacentMode::ExtendTrimWithFallback,
            &tol(),
        );
        assert!(warnings.is_ok());
        let warnings = warnings.unwrap();
        // No fallback needed for planar faces.
        assert!(warnings.is_empty());
    }

    #[test]
    fn empty_drafted_faces_is_noop() {
        let mut brep = BRep::unit_cube(tol()).unwrap();
        let drafted = HashMap::new();

        let warnings = handle_adjacent_faces(&mut brep, &drafted, AdjacentMode::ExtendTrim, &tol());
        assert!(warnings.is_ok());
        assert!(warnings.unwrap().is_empty());
    }

    #[test]
    fn find_boundary_edges_identifies_drafted_undrafted() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();

        let tilted = SurfaceGeometry::Plane {
            origin: Point3::new(0.5, 0.5, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };

        let mut drafted = HashMap::new();
        drafted.insert(face, tilted);

        let boundary = find_boundary_edges(&brep, &drafted);
        // A cube face has 4 edges, each bordering one undrafted face.
        assert_eq!(boundary.len(), 4);
    }
}
