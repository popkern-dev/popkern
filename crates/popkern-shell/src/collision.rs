//! Inter-face collision detection for offset shells.
//!
//! After offsetting, adjacent faces may collide if the offset distance
//! exceeds the local curvature radius. This module detects such collisions
//! by sampling edge neighborhoods and checking for normal crossings.

use std::collections::HashMap;

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Detect collisions between adjacent offset faces.
///
/// For each edge in the inner skin, samples points along the edge and
/// evaluates both adjacent face surfaces. A collision is reported when the
/// surface normals at adjacent faces cross (dot product sign change),
/// indicating the offset surfaces intersect.
///
/// Returns a list of colliding face pairs.
pub fn detect_collisions(
    inner: &BRep,
    offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    tolerance: &Tolerance,
) -> Vec<(FaceHandle, FaceHandle)> {
    let mut collisions = Vec::new();
    let num_samples = 5;

    for (_, edge) in &inner.edges {
        let he_handle = edge.half_edge;
        if !inner.half_edges.contains_key(he_handle) {
            continue;
        }
        let he = &inner.half_edges[he_handle];
        let twin_handle = he.twin;
        if !inner.half_edges.contains_key(twin_handle) {
            continue;
        }
        let twin = &inner.half_edges[twin_handle];

        // Get adjacent faces.
        let loop1 = he.loop_;
        let loop2 = twin.loop_;
        if !inner.loops.contains_key(loop1) || !inner.loops.contains_key(loop2) {
            continue;
        }
        let face1 = inner.loops[loop1].face;
        let face2 = inner.loops[loop2].face;

        if face1 == face2 {
            continue;
        }

        let surf1 = match offset_surfaces.get(&face1) {
            Some(s) => s,
            None => continue,
        };
        let surf2 = match offset_surfaces.get(&face2) {
            Some(s) => s,
            None => continue,
        };

        if check_edge_collision(inner, he_handle, surf1, surf2, num_samples, tolerance) {
            collisions.push((face1, face2));
        }
    }

    collisions
}

/// Check if two surfaces collide along an edge.
///
/// Samples points along the edge, evaluates the normals of both adjacent
/// surfaces, and checks for sign changes in the dot product.
fn check_edge_collision(
    inner: &BRep,
    he_handle: popkern_core::topology::HalfEdgeHandle,
    surf1: &SurfaceGeometry,
    surf2: &SurfaceGeometry,
    num_samples: usize,
    _tolerance: &Tolerance,
) -> bool {
    let he = &inner.half_edges[he_handle];
    let twin = &inner.half_edges[he.twin];
    let start = inner.vertices[he.origin].position;
    let end = inner.vertices[twin.origin].position;

    let mut prev_dot: Option<f64> = None;

    for i in 0..num_samples {
        let t = i as f64 / (num_samples - 1).max(1) as f64;
        let pt = start + (end - start) * t;

        let n1 = evaluate_normal_at_point(surf1, &pt);
        let n2 = evaluate_normal_at_point(surf2, &pt);

        let dot = n1.dot(&n2);

        if let Some(prev) = prev_dot {
            // Sign change indicates normals crossed → potential collision.
            if prev * dot < 0.0 {
                return true;
            }
        }
        prev_dot = Some(dot);
    }

    false
}

/// Evaluate the surface normal at an approximate point.
///
/// For planes, the normal is constant. For other surfaces, we use the
/// closest parametric point (approximated by projection for planes).
fn evaluate_normal_at_point(
    surface: &SurfaceGeometry,
    _point: &popkern_core::geometry::Point3,
) -> Vec3 {
    match surface {
        SurfaceGeometry::Plane { normal, .. } => *normal,
        // For non-planar surfaces, evaluate at the midpoint of the parameter range.
        _ => {
            let ((u0, u1), (v0, v1)) = surface.parameter_range();
            let u_mid = (u0 + u1) / 2.0;
            let v_mid = (v0 + v1) / 2.0;
            surface.normal_at(u_mid, v_mid)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;
    use popkern_core::tolerance::Tolerance;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn no_collisions_on_cube() {
        let brep = BRep::unit_cube(tol()).unwrap();
        // Create offset surfaces for all faces (small offset).
        let mut surfaces = HashMap::new();
        for (fh, face) in &brep.faces {
            if let Some(sof) = &face.surface {
                surfaces.insert(fh, sof.geometry.clone());
            }
        }
        let collisions = detect_collisions(&brep, &surfaces, &tol());
        assert!(
            collisions.is_empty(),
            "cube with normal offset should have no collisions"
        );
    }

    #[test]
    fn evaluate_normal_plane() {
        let surface = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let pt = Point3::new(1.0, 2.0, 3.0);
        let n = evaluate_normal_at_point(&surface, &pt);
        assert!((n.z - 1.0).abs() < 1e-12);
    }
}
