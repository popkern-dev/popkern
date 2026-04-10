//! Edge curve recomputation for drafted faces.
//!
//! For each edge between a drafted and an undrafted face (or between two
//! drafted faces), computes the new edge curve from the intersection of
//! the two adjacent surfaces.

use std::collections::HashMap;

use crate::error::DraftError;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Recompute edge curves affected by drafted face surfaces.
///
/// For each edge adjacent to at least one drafted face, computes the new
/// curve from the intersection of the two adjacent surfaces.
///
/// # Errors
///
/// Returns [`DraftError::Topology`] if adjacency lookups fail.
pub fn recompute_draft_edges(
    brep: &mut BRep,
    drafted: &HashMap<FaceHandle, SurfaceGeometry>,
    tolerance: &Tolerance,
) -> Result<(), DraftError> {
    if drafted.is_empty() {
        return Ok(());
    }

    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();

    for eh in edge_handles {
        let edge = &brep.edges[eh];
        let he_handle = edge.half_edge;

        if !brep.half_edges.contains_key(he_handle) {
            continue;
        }

        let he = &brep.half_edges[he_handle];
        let twin_handle = he.twin;

        if !brep.half_edges.contains_key(twin_handle) {
            continue;
        }

        let twin = &brep.half_edges[twin_handle];

        let loop1 = he.loop_;
        let loop2 = twin.loop_;

        if !brep.loops.contains_key(loop1) || !brep.loops.contains_key(loop2) {
            continue;
        }

        let face1 = brep.loops[loop1].face;
        let face2 = brep.loops[loop2].face;

        // Check if at least one face is drafted.
        let f1_drafted = drafted.contains_key(&face1);
        let f2_drafted = drafted.contains_key(&face2);
        if !f1_drafted && !f2_drafted {
            continue;
        }

        // Get surfaces (drafted if available, otherwise original).
        let surf1 = get_surface(brep, face1, drafted);
        let surf2 = get_surface(brep, face2, drafted);

        let (surf1, surf2) = match (surf1, surf2) {
            (Some(s1), Some(s2)) => (s1, s2),
            _ => continue,
        };

        // Compute new edge curve.
        if let Some(curve_on_edge) = compute_edge_curve(brep, eh, &surf1, &surf2, tolerance)? {
            brep.edges[eh].curve = Some(curve_on_edge);
        }
    }

    Ok(())
}

/// Get the surface for a face, preferring the drafted version.
fn get_surface(
    brep: &BRep,
    face: FaceHandle,
    drafted: &HashMap<FaceHandle, SurfaceGeometry>,
) -> Option<SurfaceGeometry> {
    if let Some(surf) = drafted.get(&face) {
        return Some(surf.clone());
    }
    brep.faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|sof| sof.geometry.clone())
}

/// Compute the edge curve from the intersection of two surfaces.
fn compute_edge_curve(
    brep: &BRep,
    edge_handle: EdgeHandle,
    surf1: &SurfaceGeometry,
    surf2: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> Result<Option<CurveOnEdge>, DraftError> {
    let edge = &brep.edges[edge_handle];
    let he = &brep.half_edges[edge.half_edge];
    let twin = &brep.half_edges[he.twin];
    let start_pos = brep.vertices[he.origin].position;
    let end_pos = brep.vertices[twin.origin].position;

    match (surf1, surf2) {
        (
            SurfaceGeometry::Plane {
                origin: o1,
                normal: n1,
                ..
            },
            SurfaceGeometry::Plane {
                origin: o2,
                normal: n2,
                ..
            },
        ) => {
            let direction = n1.cross(n2);
            let dir_len = direction.length();
            if dir_len < tolerance.solver() {
                return Ok(None);
            }
            let direction = direction * (1.0 / dir_len);

            let d1 = n1.dot(&Vec3::new(o1.x, o1.y, o1.z));
            let d2 = n2.dot(&Vec3::new(o2.x, o2.y, o2.z));
            let point_on_line = solve_point_on_plane_plane(*n1, d1, *n2, d2, direction)?;

            let t_start = (start_pos - point_on_line).dot(&direction);
            let t_end = (end_pos - point_on_line).dot(&direction);

            Ok(Some(CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: point_on_line,
                    direction,
                },
                t_start,
                t_end,
            }))
        }
        _ => {
            // For non-planar surfaces, keep the existing edge curve.
            Ok(None)
        }
    }
}

/// Find a point on the intersection of two planes.
fn solve_point_on_plane_plane(
    n1: Vec3,
    d1: f64,
    n2: Vec3,
    d2: f64,
    direction: Vec3,
) -> Result<popkern_core::geometry::Point3, DraftError> {
    let ax = direction.x.abs();
    let ay = direction.y.abs();
    let az = direction.z.abs();

    let (x, y, z) = if ax >= ay && ax >= az {
        let det = n1.y * n2.z - n1.z * n2.y;
        if det.abs() < 1e-10 {
            return Err(DraftError::Topology(
                "degenerate plane-plane intersection".into(),
            ));
        }
        let y = (d1 * n2.z - n1.z * d2) / det;
        let z = (n1.y * d2 - d1 * n2.y) / det;
        (0.0, y, z)
    } else if ay >= ax && ay >= az {
        let det = n1.x * n2.z - n1.z * n2.x;
        if det.abs() < 1e-10 {
            return Err(DraftError::Topology(
                "degenerate plane-plane intersection".into(),
            ));
        }
        let x = (d1 * n2.z - n1.z * d2) / det;
        let z = (n1.x * d2 - d1 * n2.x) / det;
        (x, 0.0, z)
    } else {
        let det = n1.x * n2.y - n1.y * n2.x;
        if det.abs() < 1e-10 {
            return Err(DraftError::Topology(
                "degenerate plane-plane intersection".into(),
            ));
        }
        let x = (d1 * n2.y - n1.y * d2) / det;
        let y = (n1.x * d2 - d1 * n2.x) / det;
        (x, y, 0.0)
    };

    Ok(popkern_core::geometry::Point3::new(x, y, z))
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn recompute_empty_drafted_is_noop() {
        let mut brep = BRep::unit_cube(tol()).unwrap();
        let drafted = HashMap::new();
        recompute_draft_edges(&mut brep, &drafted, &tol()).unwrap();
    }

    #[test]
    fn recompute_drafted_cube_face_edges() {
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

        let result = recompute_draft_edges(&mut brep, &drafted, &tol());
        assert!(result.is_ok());
    }

    #[test]
    fn solve_point_on_plane_plane_basic() {
        let n1 = Vec3::new(0.0, 0.0, 1.0);
        let n2 = Vec3::new(0.0, 1.0, 0.0);
        let direction = n1.cross(&n2);
        let dir_norm = direction.normalized().unwrap();

        let result = solve_point_on_plane_plane(n1, 0.5, n2, 0.3, dir_norm).unwrap();
        // n1.p = 0.5 => z = 0.5
        // n2.p = 0.3 => y = 0.3
        assert!((result.z - 0.5).abs() < 1e-10);
        assert!((result.y - 0.3).abs() < 1e-10);
    }
}
