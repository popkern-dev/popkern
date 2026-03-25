//! Side wall construction at shell openings.
//!
//! When faces are removed from a shell, the resulting gap between the outer
//! and inner skins must be closed with side wall faces. This module identifies
//! boundary edges (edges where a removed face meets a kept face) and builds
//! quadrilateral side wall faces connecting the outer and inner skins.

use std::collections::{HashMap, HashSet};

use crate::error::ShellError;
use crate::skin_map::SkinMap;
use popkern_core::constructors::BRepBuilder;
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, ShellHandle, VertexHandle,
};

/// Result of side wall construction.
#[derive(Debug, Clone)]
pub struct SideWallResult {
    /// New faces created as side walls.
    pub faces: Vec<FaceHandle>,
}

/// Build side wall faces connecting outer and inner skins at shell openings.
///
/// For each removed face, walks its outer loop edges. An edge whose opposite
/// face is NOT in the removed set is a boundary edge. For each boundary edge,
/// a quad face is created connecting the outer and inner vertices.
///
/// # Errors
///
/// Returns [`ShellError::Topology`] if handle lookups fail.
#[allow(clippy::too_many_arguments)]
pub fn build_side_walls(
    builder: &mut BRepBuilder,
    outer: &BRep,
    removed_faces: &[FaceHandle],
    skin_map: &SkinMap,
    outer_vertex_map: &HashMap<VertexHandle, VertexHandle>,
    inner_vertex_map: &HashMap<VertexHandle, VertexHandle>,
    shell_handle: ShellHandle,
    tolerance: &Tolerance,
) -> Result<SideWallResult, ShellError> {
    let removed_set: HashSet<FaceHandle> = removed_faces.iter().copied().collect();
    let mut faces = Vec::new();
    // Track which outer edges already have a side wall to avoid duplicates.
    let mut processed_edges: HashSet<EdgeHandle> = HashSet::new();

    for &removed_fh in removed_faces {
        let boundary_edges = find_boundary_edges(outer, removed_fh, &removed_set)?;

        for (outer_eh, outer_he) in boundary_edges {
            if processed_edges.contains(&outer_eh) {
                continue;
            }
            processed_edges.insert(outer_eh);

            let face = build_one_side_wall(
                builder,
                outer,
                outer_he,
                skin_map,
                outer_vertex_map,
                inner_vertex_map,
                shell_handle,
                tolerance,
            )?;
            faces.push(face);
        }
    }

    Ok(SideWallResult { faces })
}

/// Find boundary edges of a removed face.
///
/// Returns `(EdgeHandle, HalfEdgeHandle)` pairs where the half-edge belongs
/// to the removed face and the twin's face is NOT removed.
fn find_boundary_edges(
    outer: &BRep,
    removed_fh: FaceHandle,
    removed_set: &HashSet<FaceHandle>,
) -> Result<Vec<(EdgeHandle, HalfEdgeHandle)>, ShellError> {
    let face = &outer.faces[removed_fh];
    let start_he = outer.loops[face.outer_loop].half_edge;
    let mut result = Vec::new();
    let mut current = start_he;

    loop {
        let he = &outer.half_edges[current];
        let twin = &outer.half_edges[he.twin];
        let twin_loop = twin.loop_;

        if outer.loops.contains_key(twin_loop) {
            let twin_face = outer.loops[twin_loop].face;
            if !removed_set.contains(&twin_face) {
                result.push((he.edge, current));
            }
        }

        current = he.next;
        if current == start_he {
            break;
        }
    }

    Ok(result)
}

/// Build a single side wall quad face.
///
/// Given a boundary half-edge on the outer skin (belonging to a removed face),
/// creates a quad connecting outer_v0 -> outer_v1 -> inner_v1 -> inner_v0.
#[allow(clippy::too_many_arguments)]
fn build_one_side_wall(
    builder: &mut BRepBuilder,
    outer: &BRep,
    outer_he: HalfEdgeHandle,
    skin_map: &SkinMap,
    outer_vertex_map: &HashMap<VertexHandle, VertexHandle>,
    inner_vertex_map: &HashMap<VertexHandle, VertexHandle>,
    shell_handle: ShellHandle,
    _tolerance: &Tolerance,
) -> Result<FaceHandle, ShellError> {
    let he = &outer.half_edges[outer_he];
    let twin = &outer.half_edges[he.twin];

    // Outer vertices (in the original outer BRep).
    let outer_v0_orig = he.origin;
    let outer_v1_orig = twin.origin;

    // Map to builder vertices.
    let builder_outer_v0 = *outer_vertex_map
        .get(&outer_v0_orig)
        .ok_or_else(|| ShellError::Topology("outer vertex not in vertex map".into()))?;
    let builder_outer_v1 = *outer_vertex_map
        .get(&outer_v1_orig)
        .ok_or_else(|| ShellError::Topology("outer vertex not in vertex map".into()))?;

    // Inner vertices: map outer -> inner via skin_map, then to builder.
    let inner_v0_orig = *skin_map
        .vertices
        .get(&outer_v0_orig)
        .ok_or_else(|| ShellError::Topology("vertex not in skin map".into()))?;
    let inner_v1_orig = *skin_map
        .vertices
        .get(&outer_v1_orig)
        .ok_or_else(|| ShellError::Topology("vertex not in skin map".into()))?;

    let builder_inner_v0 = *inner_vertex_map
        .get(&inner_v0_orig)
        .ok_or_else(|| ShellError::Topology("inner vertex not in vertex map".into()))?;
    let builder_inner_v1 = *inner_vertex_map
        .get(&inner_v1_orig)
        .ok_or_else(|| ShellError::Topology("inner vertex not in vertex map".into()))?;

    // Get 3D positions for geometry.
    let p0 = builder.brep().vertices[builder_outer_v0].position;
    let p1 = builder.brep().vertices[builder_outer_v1].position;
    let p2 = builder.brep().vertices[builder_inner_v1].position;
    let p3 = builder.brep().vertices[builder_inner_v0].position;

    // Create 4 edges for the quad: outer, thickness_end, inner_rev, thickness_start.
    let (e_top, he_top_fwd, _he_top_rev) = builder.add_edge(builder_outer_v0, builder_outer_v1);
    let (e_right, he_right_fwd, _he_right_rev) =
        builder.add_edge(builder_outer_v1, builder_inner_v1);
    let (e_bottom, he_bottom_fwd, _he_bottom_rev) =
        builder.add_edge(builder_inner_v1, builder_inner_v0);
    let (e_left, he_left_fwd, _he_left_rev) = builder.add_edge(builder_inner_v0, builder_outer_v0);

    // Set edge curves (straight lines).
    set_line_curve(builder, e_top, p0, p1);
    set_line_curve(builder, e_right, p1, p2);
    set_line_curve(builder, e_bottom, p2, p3);
    set_line_curve(builder, e_left, p3, p0);

    // Create the face and its loop.
    let face = builder.add_face(shell_handle);
    builder
        .add_loop(
            &[he_top_fwd, he_right_fwd, he_bottom_fwd, he_left_fwd],
            face,
        )
        .map_err(ShellError::Core)?;

    // Compute and set a planar surface for the quad.
    if let Some(surface) = plane_from_quad(p0, p1, p2, p3) {
        builder.set_face_surface(
            face,
            SurfaceOnFace {
                geometry: surface,
                orientation: true,
            },
        );
    }

    Ok(face)
}

/// Set a line curve on an edge.
fn set_line_curve(builder: &mut BRepBuilder, edge: EdgeHandle, start: Point3, end: Point3) {
    let direction = end - start;
    builder.set_edge_curve(
        edge,
        CurveOnEdge {
            geometry: CurveGeometry::Line {
                origin: start,
                direction,
            },
            t_start: 0.0,
            t_end: 1.0,
        },
    );
}

/// Compute a planar surface from 4 quad vertices.
///
/// Uses the first 3 vertices to define the plane. Returns `None` if the
/// points are collinear.
fn plane_from_quad(p0: Point3, p1: Point3, p2: Point3, _p3: Point3) -> Option<SurfaceGeometry> {
    let raw_u = p1 - p0;
    let raw_v = p2 - p0;
    let normal = raw_u.cross(&raw_v);
    let normal = normal.normalized()?;
    let u_axis = raw_u.normalized()?;
    let v_proj = raw_v - u_axis * raw_v.dot(&u_axis);
    let v_axis = v_proj.normalized()?;
    Some(SurfaceGeometry::Plane {
        origin: p0,
        normal,
        u_axis,
        v_axis,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn find_boundary_edges_one_removed_face() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let first_face = brep.faces.keys().next().unwrap();
        let removed_set: HashSet<FaceHandle> = [first_face].into_iter().collect();
        let edges = find_boundary_edges(&brep, first_face, &removed_set).unwrap();
        // A cube face has 4 edges, all boundary (adjacent faces are kept).
        assert_eq!(edges.len(), 4);
    }

    #[test]
    fn find_boundary_edges_two_adjacent_removed() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
        // Remove first two faces. Some edges between them won't be boundary.
        let removed_set: HashSet<FaceHandle> =
            [face_handles[0], face_handles[1]].into_iter().collect();

        let edges_f0 = find_boundary_edges(&brep, face_handles[0], &removed_set).unwrap();
        let edges_f1 = find_boundary_edges(&brep, face_handles[1], &removed_set).unwrap();

        // Each face has 4 edges, but edges shared between the two removed faces
        // are not boundary edges. Depending on adjacency, we expect < 4 for at
        // least one of them (if adjacent) or 4 each (if non-adjacent).
        let total = edges_f0.len() + edges_f1.len();
        assert!(total >= 6 && total <= 8, "total boundary edges: {total}");
    }

    #[test]
    fn plane_from_quad_basic() {
        let p0 = Point3::new(0.0, 0.0, 0.0);
        let p1 = Point3::new(1.0, 0.0, 0.0);
        let p2 = Point3::new(1.0, 0.0, 1.0);
        let p3 = Point3::new(0.0, 0.0, 1.0);
        let surface = plane_from_quad(p0, p1, p2, p3).unwrap();
        match surface {
            SurfaceGeometry::Plane { normal, .. } => {
                // Normal should be along y-axis.
                assert!(normal.y.abs() > 0.9);
            }
            _ => panic!("expected Plane"),
        }
    }

    #[test]
    fn plane_from_quad_collinear_returns_none() {
        let p0 = Point3::new(0.0, 0.0, 0.0);
        let p1 = Point3::new(1.0, 0.0, 0.0);
        let p2 = Point3::new(2.0, 0.0, 0.0);
        let p3 = Point3::new(3.0, 0.0, 0.0);
        assert!(plane_from_quad(p0, p1, p2, p3).is_none());
    }
}
