//! Welded mesh gauges for hardening matrices (same quantized weld as
//! popkern-part's watertight gate — raw tessellation leaves per-face
//! islands whose borders read as ~60 false opens on a plain cube
//! fillet).

use popkern_core::topology::BRep;
use std::collections::HashMap;

pub(crate) fn welded_open_edges(brep: &BRep) -> usize {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .expect("tessellate");
    let mut key: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut id = |p: popkern_core::geometry::Point3| -> usize {
        let k = (
            (p.x * 1.0e7).round() as i64,
            (p.y * 1.0e7).round() as i64,
            (p.z * 1.0e7).round() as i64,
        );
        let next = key.len();
        *key.entry(k).or_insert(next)
    };
    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        let w = [
            id(mesh.vertices[t[0]].position),
            id(mesh.vertices[t[1]].position),
            id(mesh.vertices[t[2]].position),
        ];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    counts.values().filter(|&&c| c == 1).count()
}

#[allow(dead_code)]
pub(crate) fn mesh_volume(brep: &BRep) -> f64 {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .expect("tessellate");
    let mut vol = 0.0_f64;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        vol += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    (vol / 6.0).abs()
}
