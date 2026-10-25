//! Welding a face-by-face tessellation into one indexed mesh per solid.
//!
//! popkern's tessellation does not share vertices across faces. Here
//! positions are quantised at 1e-7 (the same weld as
//! `popkern_mesh::welded_open_edges`), first occurrence wins, in
//! triangle order. Triangles that collapse under the weld are dropped.
//! An edge is open unless its reverse is used exactly as often.

use std::collections::BTreeMap;

use popkern_analysis::FaceId;
use popkern_mesh::TriMesh;

use crate::adapter::AdapterError;

/// Weld quantum, matching popkern's watertight census.
pub const WELD_QUANTUM: f64 = 1.0e-7;

/// One solid's welded mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct WeldedMesh {
    /// Welded vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Triangles over `vertices`.
    pub triangles: Vec<[u32; 3]>,
    /// Face of each triangle.
    pub triangle_face: Vec<FaceId>,
    /// Directed edges without a matching reverse edge.
    pub open_edges: usize,
}

fn key(p: [f64; 3]) -> Option<[i64; 3]> {
    if p.iter().any(|c| !c.is_finite()) {
        return None;
    }
    Some(p.map(|c| (c / WELD_QUANTUM).round() as i64))
}

/// Welds the triangles `tagged` (index into `mesh.triangles`, face) of
/// `mesh`. A triangle with a non-finite corner is a hole: it is dropped
/// and its three edges count as open.
pub fn weld(mesh: &TriMesh, tagged: &[(usize, FaceId)]) -> Result<WeldedMesh, AdapterError> {
    let mut index: BTreeMap<[i64; 3], u32> = BTreeMap::new();
    let mut vertices: Vec<[f64; 3]> = Vec::new();
    let mut triangles = Vec::new();
    let mut triangle_face = Vec::new();
    let mut holes = 0usize;
    for &(t, face) in tagged {
        let Some(corners) = mesh.triangles.get(t) else {
            return Err(AdapterError::Mesh {
                detail: format!("triangle {t} out of range"),
            });
        };
        let mut welded = [0u32; 3];
        let mut finite = true;
        for (slot, &v) in welded.iter_mut().zip(corners) {
            let Some(vertex) = mesh.vertices.get(v) else {
                return Err(AdapterError::Mesh {
                    detail: format!("vertex {v} out of range"),
                });
            };
            let p = [vertex.position.x, vertex.position.y, vertex.position.z];
            let Some(k) = key(p) else {
                finite = false;
                break;
            };
            let next = u32::try_from(vertices.len()).map_err(|_| AdapterError::TooLarge)?;
            *slot = *index.entry(k).or_insert_with(|| {
                vertices.push(p);
                next
            });
        }
        if !finite {
            holes += 1;
            continue;
        }
        let [a, b, c] = welded;
        if a == b || b == c || c == a {
            continue;
        }
        triangles.push(welded);
        triangle_face.push(face);
    }
    Ok(WeldedMesh {
        open_edges: open_edges(&triangles) + 3 * holes,
        vertices,
        triangles,
        triangle_face,
    })
}

/// Directed edges whose use count differs from their reverse's.
pub fn open_edges(triangles: &[[u32; 3]]) -> usize {
    let mut balance: BTreeMap<(u32, u32), i64> = BTreeMap::new();
    for t in triangles {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            if a < b {
                *balance.entry((a, b)).or_insert(0) += 1;
            } else {
                *balance.entry((b, a)).or_insert(0) -= 1;
            }
        }
    }
    balance.values().map(|v| v.unsigned_abs() as usize).sum()
}
