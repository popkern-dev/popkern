//! Topological vertices are ONE point in every face mesh.
//!
//! Edge samples snap their ends to the shared vertex positions
//! (`boundary_first_edge_samples`), but a face can place a vertex's
//! point by its own evaluation instead — a corner patch whose apex sits
//! on the vertex evaluates it from its surface — and the two copies
//! then differ in the last bits. The boss plate's four corner patches
//! put their corners ~5e-16 off the top fillets' (16 open edges in the
//! bit-identical weld). After every face is meshed, a mesh vertex
//! within float noise of a topological vertex takes that vertex's
//! exact position, so the faces meeting there share it by
//! construction.

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::topology::BRep;

use crate::mesh::FaceMesh;

/// Float noise, relative to the model's size: far below any tolerance
/// or chord height, far above f64 rounding of a vertex evaluated twice.
const NOISE: f64 = 1e-12;

/// Snap every face-mesh vertex within `NOISE × size` of a topological
/// vertex onto it. Returns how many moved.
pub(crate) fn snap_to_vertices(brep: &BRep, face_meshes: &mut [FaceMesh]) -> usize {
    let points: Vec<Point3> = brep.vertices.values().map(|v| v.position).collect();
    if points.is_empty() {
        return 0;
    }
    let size = points
        .iter()
        .map(|p| p.x.abs().max(p.y.abs()).max(p.z.abs()))
        .fold(0.0_f64, f64::max)
        .max(1e-300);
    let eps = NOISE * size;
    // A grid of cell 4·eps: a point within eps of a vertex lies in the
    // vertex's cell or a neighbour.
    let cell = 4.0 * eps;
    let key = |p: &Point3| {
        (
            (p.x / cell).floor() as i64,
            (p.y / cell).floor() as i64,
            (p.z / cell).floor() as i64,
        )
    };
    let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    for (i, p) in points.iter().enumerate() {
        grid.entry(key(p)).or_default().push(i);
    }
    let mut moved = 0;
    for fm in face_meshes.iter_mut() {
        for v in fm.vertices.iter_mut() {
            let (kx, ky, kz) = key(&v.position);
            let mut best: Option<(f64, usize)> = None;
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let Some(list) = grid.get(&(kx + dx, ky + dy, kz + dz)) else {
                            continue;
                        };
                        for &i in list {
                            let d = points[i].distance_to(&v.position);
                            if d <= eps && best.is_none_or(|(bd, _)| d < bd) {
                                best = Some((d, i));
                            }
                        }
                    }
                }
            }
            if let Some((d, i)) = best
                && d > 0.0
            {
                v.position = points[i];
                moved += 1;
            }
        }
    }
    moved
}
