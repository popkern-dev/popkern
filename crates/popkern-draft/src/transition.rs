//! Transition face creation for draft boundary gaps.
//!
//! When a draft moves vertices but the adjacent undrafted face cannot
//! extend to match, a small ruled-surface quad face is inserted to bridge
//! the gap between old and new vertex positions.

use std::collections::HashMap;

use crate::error::DraftError;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, Edge, EdgeHandle, Face, FaceHandle, HalfEdge, Loop, Vertex, VertexHandle,
};

/// Insert a transition face between drafted and undrafted face boundaries.
///
/// When a draft moves vertices, this function creates a small ruled-surface
/// quad face bridging the old vertex positions to the new (recomputed) ones.
/// If neither vertex moved beyond tolerance, no face is created.
///
/// The transition face has 4 vertices (old_start, old_end, new_end,
/// new_start), 4 edges, and a planar surface fitted through the quad.
pub(crate) fn insert_transition_face(
    brep: &mut BRep,
    edge_h: EdgeHandle,
    _drafted_face: FaceHandle,
    _undrafted_face: FaceHandle,
    original_positions: &HashMap<VertexHandle, Point3>,
    tolerance: &Tolerance,
) -> Result<(), DraftError> {
    let edge = brep
        .edges
        .get(edge_h)
        .ok_or_else(|| DraftError::Topology("edge not found".into()))?;
    let he_h = edge.half_edge;
    let he = brep
        .half_edges
        .get(he_h)
        .ok_or_else(|| DraftError::Topology("half-edge not found".into()))?;
    let twin_h = he.twin;
    let twin = brep
        .half_edges
        .get(twin_h)
        .ok_or_else(|| DraftError::Topology("twin half-edge not found".into()))?;

    let v_start = he.origin;
    let v_end = twin.origin;

    let new_start = brep.vertices[v_start].position;
    let new_end = brep.vertices[v_end].position;

    let old_start = original_positions
        .get(&v_start)
        .copied()
        .unwrap_or(new_start);
    let old_end = original_positions.get(&v_end).copied().unwrap_or(new_end);

    // Check if either vertex actually moved.
    let start_moved = old_start.distance_to(&new_start) > tolerance.linear();
    let end_moved = old_end.distance_to(&new_end) > tolerance.linear();

    if !start_moved && !end_moved {
        // No gap to bridge — no transition face needed.
        return Ok(());
    }

    // Build the transition quad: old_start -> old_end -> new_end -> new_start.
    // Fit a plane through the 4 points.
    let normal = match fit_plane_normal(old_start, old_end, new_end, new_start) {
        Some(n) => n,
        None => return Ok(()), // Degenerate quad, skip.
    };

    let centroid = Point3::new(
        (old_start.x + old_end.x + new_end.x + new_start.x) * 0.25,
        (old_start.y + old_end.y + new_end.y + new_start.y) * 0.25,
        (old_start.z + old_end.z + new_end.z + new_start.z) * 0.25,
    );

    // Compute u_axis from old_start -> old_end.
    let u_dir = old_end - old_start;
    let u_axis = u_dir.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v_axis = normal.cross(&u_axis);

    let surface = SurfaceGeometry::Plane {
        origin: centroid,
        normal,
        u_axis,
        v_axis,
    };

    build_quad_face(
        brep, old_start, old_end, new_end, new_start, surface, tolerance,
    )
}

/// Build a quad face from 4 corner positions with the given surface geometry.
///
/// Creates 4 vertices, 4 edges (each with a half-edge pair), one face, and
/// one loop wiring the half-edges in order: p0 -> p1 -> p2 -> p3.
fn build_quad_face(
    brep: &mut BRep,
    p0: Point3,
    p1: Point3,
    p2: Point3,
    p3: Point3,
    surface: SurfaceGeometry,
    tolerance: &Tolerance,
) -> Result<(), DraftError> {
    let vh0 = brep.vertices.insert(Vertex {
        position: p0,
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let vh1 = brep.vertices.insert(Vertex {
        position: p1,
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let vh2 = brep.vertices.insert(Vertex {
        position: p2,
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let vh3 = brep.vertices.insert(Vertex {
        position: p3,
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });

    let shell_h = brep
        .faces
        .values()
        .next()
        .map(|f| f.shell)
        .ok_or_else(|| DraftError::Topology("no faces in BRep".into()))?;

    let face_h = brep.faces.insert(Face {
        outer_loop: Default::default(),
        inner_loops: Vec::new(),
        shell: shell_h,
        surface: Some(SurfaceOnFace {
            geometry: surface,
            orientation: true,
        }),
        name: None,
        tolerance: None,
        role: Default::default(),
    });
    let loop_h = brep.loops.insert(Loop {
        half_edge: Default::default(),
        face: face_h,
    });
    brep.faces[face_h].outer_loop = loop_h;

    let verts = [vh0, vh1, vh2, vh3];
    let mut he_handles = Vec::with_capacity(4);
    let mut twin_handles = Vec::with_capacity(4);

    for i in 0..4 {
        let origin = verts[i];
        let dest = verts[(i + 1) % 4];

        let he = brep.half_edges.insert(HalfEdge {
            origin,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: Default::default(),
            loop_: loop_h,
            pcurve: None,
        });
        let tw = brep.half_edges.insert(HalfEdge {
            origin: dest,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: Default::default(),
            loop_: Default::default(),
            pcurve: None,
        });

        brep.half_edges[he].twin = tw;
        brep.half_edges[tw].twin = he;

        let start_pos = brep.vertices[origin].position;
        let end_pos = brep.vertices[dest].position;
        let dir = end_pos - start_pos;
        let dir_len = dir.length();
        let curve = if dir_len > tolerance.linear() {
            let direction = dir * (1.0 / dir_len);
            Some(CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: start_pos,
                    direction,
                },
                t_start: 0.0,
                t_end: dir_len,
            })
        } else {
            None
        };

        let edge = brep.edges.insert(Edge {
            half_edge: he,
            curve,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });
        brep.half_edges[he].edge = edge;
        brep.half_edges[tw].edge = edge;

        he_handles.push(he);
        twin_handles.push(tw);
    }

    // Wire loop: next/prev for interior half-edges.
    for i in 0..4 {
        let next = (i + 1) % 4;
        let prev = (i + 3) % 4;
        brep.half_edges[he_handles[i]].next = he_handles[next];
        brep.half_edges[he_handles[i]].prev = he_handles[prev];
    }

    // Wire twin next/prev into a reverse loop.
    for i in 0..4 {
        let next = (i + 3) % 4;
        let prev = (i + 1) % 4;
        brep.half_edges[twin_handles[i]].next = twin_handles[next];
        brep.half_edges[twin_handles[i]].prev = twin_handles[prev];
    }

    brep.loops[loop_h].half_edge = he_handles[0];
    for i in 0..4 {
        brep.vertices[verts[i]].outgoing = he_handles[i];
    }

    Ok(())
}

/// Compute the normal of a planar face fitted through 4 points.
fn fit_plane_normal(p0: Point3, p1: Point3, p2: Point3, p3: Point3) -> Option<Vec3> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn fit_plane_normal_unit_square() {
        let p0 = Point3::new(0.0, 0.0, 0.0);
        let p1 = Point3::new(1.0, 0.0, 0.0);
        let p2 = Point3::new(1.0, 1.0, 0.0);
        let p3 = Point3::new(0.0, 1.0, 0.0);
        let n = fit_plane_normal(p0, p1, p2, p3);
        assert!(n.is_some());
        let n = n.unwrap();
        assert!(n.z.abs() > 0.99, "expected Z-facing normal, got {:?}", n);
    }

    #[test]
    fn fit_plane_normal_degenerate() {
        let p = Point3::new(0.0, 0.0, 0.0);
        assert!(fit_plane_normal(p, p, p, p).is_none());
    }

    #[test]
    fn build_quad_face_on_empty_brep_fails() {
        let mut brep = BRep::new(tol());
        let surface = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let result = build_quad_face(
            &mut brep,
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            surface,
            &tol(),
        );
        // No faces in BRep to get shell from.
        assert!(result.is_err());
    }

    #[test]
    fn build_quad_face_on_cube() {
        let mut brep = BRep::unit_cube(tol()).unwrap();
        let original_faces = brep.faces.len();
        let surface = SurfaceGeometry::Plane {
            origin: Point3::new(0.5, 0.5, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let result = build_quad_face(
            &mut brep,
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            surface,
            &tol(),
        );
        assert!(result.is_ok());
        assert_eq!(brep.faces.len(), original_faces + 1);
    }
}
