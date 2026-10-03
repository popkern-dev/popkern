//! The welded open-edge census: the watertightness measure every gate
//! uses (the part facade, the CLI, popkern live, the gallery).
//!
//! Vertices are welded at 1e-7 and an edge used by exactly one
//! triangle is open. A triangle with a non-finite corner is a hole, not
//! a degenerate sliver: its three edges count as open. (Before this
//! rule a mesh of NaN positions welded every vertex to one key, dropped
//! every triangle as degenerate, and reported 0 open edges — a broken
//! STEP import passed as watertight.)

use std::collections::HashMap;

use popkern_core::geometry::Point3;

use crate::mesh::TriMesh;

/// Open edges of `mesh` after welding at 1e-7 (0 = watertight); every
/// triangle with a non-finite corner adds three.
pub fn welded_open_edges(mesh: &TriMesh) -> usize {
    let mut key: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut id = |p: Point3| -> usize {
        let k = (
            (p.x * 1.0e7).round() as i64,
            (p.y * 1.0e7).round() as i64,
            (p.z * 1.0e7).round() as i64,
        );
        let next = key.len();
        *key.entry(k).or_insert(next)
    };
    let finite = |p: &Point3| p.x.is_finite() && p.y.is_finite() && p.z.is_finite();
    let mut holes = 0;
    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        let corners = [
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        ];
        if !corners.iter().all(finite) {
            holes += 3;
            continue;
        }
        let w = [id(corners[0]), id(corners[1]), id(corners[2])];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    holes + counts.values().filter(|&&c| c == 1).count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::MeshVertex;
    use popkern_core::geometry::Vec3;

    fn vertex(x: f64, y: f64, z: f64) -> MeshVertex {
        MeshVertex {
            position: Point3::new(x, y, z),
            normal: Vec3::new(0.0, 0.0, 1.0),
            uv: (0.0, 0.0),
        }
    }

    #[test]
    fn a_nan_mesh_is_not_watertight() {
        let mesh = TriMesh {
            vertices: vec![
                vertex(f64::NAN, 0.0, 0.0),
                vertex(f64::NAN, 1.0, 0.0),
                vertex(f64::NAN, 0.0, 1.0),
            ],
            triangles: vec![[0, 1, 2]],
            ..TriMesh::default()
        };
        assert_eq!(welded_open_edges(&mesh), 3);
    }
}
