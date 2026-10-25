//! The sub-cell lattice: point containment along grid columns, the
//! midpoint-rule fallback for cell moments, and fragment connectivity.
//!
//! Containment of a point counts signed crossings of a ray along +x
//! with the facet pieces of the point's column beyond it: each piece
//! whose (y, z) projection contains the point and whose plane is
//! crossed past the point adds sign(n_x). The sum is 1 inside a closed
//! outward solid and 0 outside. Projected edges are half-open, so a
//! ray through an edge shared by two pieces counts it once.

use crate::clip::{ClippedSolid, FacetPiece};
use crate::grid::Grid;

use super::moment_index;

/// x^n by repeated multiplication (the same bits on every target, unlike
/// `powi`, which lowers to a runtime helper).
fn ipow(x: f64, n: usize) -> f64 {
    let mut r = 1.0;
    for _ in 0..n {
        r *= x;
    }
    r
}

/// True when the 2-D point (y, z) lies in the (y, z) projection of
/// `piece` under the half-open crossing rule.
fn projects_onto(piece: &FacetPiece, y: f64, z: f64) -> bool {
    let v = piece.polygon.vertices();
    let mut inside = false;
    for i in 0..v.len() {
        let (a, b) = (v[i], v[(i + 1) % v.len()]);
        let (ay, az, by, bz) = (a[1], a[2], b[1], b[2]);
        if (az <= z) != (bz <= z) {
            let t = (z - az) / (bz - az);
            if y < ay + t * (by - ay) {
                inside = !inside;
            }
        }
    }
    inside
}

/// The pieces of the column through cell `ijk`, from that cell to the
/// end of the column (pieces are sorted by cell).
fn column_from<'a>(grid: &Grid, clipped: &'a ClippedSolid, ijk: [usize; 3]) -> &'a [FacetPiece] {
    let first = grid.cell_index(ijk).unwrap_or(usize::MAX);
    let last = grid
        .cell_index([grid.dims()[0] - 1, ijk[1], ijk[2]])
        .unwrap_or(0);
    let pieces = &clipped.pieces;
    let start = pieces.partition_point(|p| p.cell < first);
    let end = pieces.partition_point(|p| p.cell <= last);
    &pieces[start..end.max(start)]
}

/// Whether `p` (inside cell `ijk`) is inside the solid.
pub fn contains(grid: &Grid, clipped: &ClippedSolid, ijk: [usize; 3], p: [f64; 3]) -> bool {
    let mut winding = 0i64;
    for piece in column_from(grid, clipped, ijk) {
        let n = piece.normal;
        if n[0] == 0.0 || !projects_onto(piece, p[1], p[2]) {
            continue;
        }
        let v0 = piece.polygon.vertices()[0];
        let x = v0[0] - (n[1] * (p[1] - v0[1]) + n[2] * (p[2] - v0[2])) / n[0];
        if x > p[0] {
            winding += if n[0] > 0.0 { 1 } else { -1 };
        }
    }
    winding > 0
}

/// Sub-cell centres of cell `ijk` on a 2^depth lattice, x fastest, with
/// their lattice coordinates.
pub(crate) fn lattice(grid: &Grid, ijk: [usize; 3], depth: u8) -> Vec<([usize; 3], [f64; 3])> {
    let n = 1usize << depth;
    let s = grid.cell_size() / n as f64;
    let lo = grid.cell_min(ijk);
    let mut out = Vec::with_capacity(n * n * n);
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let l = [i, j, k];
                out.push((l, [0, 1, 2].map(|a| lo[a] + (l[a] as f64 + 0.5) * s)));
            }
        }
    }
    out
}

/// The 27 moments of cell `ijk` by the midpoint rule on the lattice:
/// the fallback when the exact moments are inconsistent.
pub fn lattice_moments(
    grid: &Grid,
    clipped: &ClippedSolid,
    ijk: [usize; 3],
    depth: u8,
) -> [f64; 27] {
    let h = grid.cell_size();
    let n = (1usize << depth) as f64;
    let s = h / n;
    let weight = s * s * s;
    let lo = grid.cell_min(ijk);
    let mut m = [0.0; 27];
    for (_, p) in lattice(grid, ijk, depth) {
        if !contains(grid, clipped, ijk, p) {
            continue;
        }
        let l = [0, 1, 2].map(|a| 2.0 * (p[a] - lo[a]) / h - 1.0);
        for c in 0..3 {
            for b in 0..3 {
                for a in 0..3 {
                    m[moment_index(a, b, c)] +=
                        weight * ipow(l[0], a) * ipow(l[1], b) * ipow(l[2], c);
                }
            }
        }
    }
    m
}

/// Labels of the inside lattice points of cell `ijk` by 6-connected
/// component (0 outside, k + 1 for component k, components numbered in
/// lattice order), and the number of components.
pub fn component_labels(
    grid: &Grid,
    clipped: &ClippedSolid,
    ijk: [usize; 3],
    depth: u8,
) -> (Vec<u8>, usize) {
    let inside: Vec<bool> = lattice(grid, ijk, depth)
        .into_iter()
        .map(|(_, p)| contains(grid, clipped, ijk, p))
        .collect();
    label_components(&inside, depth)
}

/// 6-connected components of `inside` on the 2^depth lattice.
pub(crate) fn label_components(inside: &[bool], depth: u8) -> (Vec<u8>, usize) {
    let n = 1usize << depth;
    let at = |l: [usize; 3]| l[0] + n * (l[1] + n * l[2]);
    let mut label = vec![0u8; inside.len()];
    let mut components = 0usize;
    for start in 0..inside.len() {
        if !inside[start] || label[start] != 0 {
            continue;
        }
        components += 1;
        let tag = u8::try_from(components).unwrap_or(u8::MAX);
        let mut stack = vec![start];
        label[start] = tag;
        while let Some(i) = stack.pop() {
            let l = [i % n, (i / n) % n, i / (n * n)];
            for axis in 0..3 {
                for step in [-1i64, 1] {
                    let v = l[axis] as i64 + step;
                    if v < 0 || v >= n as i64 {
                        continue;
                    }
                    let mut nl = l;
                    nl[axis] = v as usize;
                    let j = at(nl);
                    if inside[j] && label[j] == 0 {
                        label[j] = tag;
                        stack.push(j);
                    }
                }
            }
        }
    }
    (label, components)
}

/// Number of 6-connected components of the inside lattice points of
/// cell `ijk`: more than one means the cell holds disconnected parts of
/// the solid (a slit narrower than the cell).
pub fn fragment_components(
    grid: &Grid,
    clipped: &ClippedSolid,
    ijk: [usize; 3],
    depth: u8,
) -> usize {
    component_labels(grid, clipped, ijk, depth).1
}
