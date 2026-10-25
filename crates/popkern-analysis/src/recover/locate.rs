//! Locating points in the active cells: display vertices (with their
//! inward normals, which pick the fragment a surface vertex bounds) and
//! interior points.

use crate::grid::Grid;
use crate::model::SolidMesh;
use crate::moments::SolidMoments;

/// Unit inward normal at each mesh vertex (area-weighted over its
/// triangles), to tell which fragment of a split cell a surface vertex
/// bounds.
pub(super) fn vertex_inward(mesh: &SolidMesh) -> Vec<[f64; 3]> {
    let mut n = vec![[0.0; 3]; mesh.vertices.len()];
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|i| mesh.vertices[i as usize]);
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let x = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        for &i in t {
            for k in 0..3 {
                n[i as usize][k] -= x[k];
            }
        }
    }
    n.iter()
        .map(|x| {
            let l = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
            if l > 0.0 { x.map(|c| c / l) } else { [0.0; 3] }
        })
        .collect()
}

/// The active cell (index into the moments) containing `p`, and `p` in
/// its local coordinates; `inward` (the surface's inward normal at a
/// surface point) tells which fragment of a split cell it bounds.
pub(super) fn locate(
    grid: &Grid,
    moments: &SolidMoments,
    p: [f64; 3],
    inward: Option<[f64; 3]>,
) -> Option<(usize, [f64; 3])> {
    let h = grid.cell_size();
    let o = grid.origin();
    let cand = |a: usize| -> Vec<usize> {
        let t = (p[a] - o[a]) / h;
        let mut out = Vec::new();
        for v in [t, t - 1e-9, t + 1e-9] {
            let f = v.floor();
            if f >= 0.0 && f < grid.dims()[a] as f64 && !out.contains(&(f as usize)) {
                out.push(f as usize);
            }
        }
        out
    };
    let (cx, cy, cz) = (cand(0), cand(1), cand(2));
    for &k in &cz {
        for &j in &cy {
            for &i in &cx {
                let idx = grid.cell_index([i, j, k])?;
                if let Some(pos) = moments.entry_at(grid, idx, p, inward) {
                    let lo = grid.cell_min([i, j, k]);
                    return Some((pos, [0, 1, 2].map(|a| 2.0 * (p[a] - lo[a]) / h - 1.0)));
                }
            }
        }
    }
    None
}
