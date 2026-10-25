//! Which fragments of two face-adjacent cells share material across
//! their common face.
//!
//! - Material crosses a face exactly where clipped pieces have edges on
//!   it, and a triangle crossing the face leaves the same segment in
//!   both cells.
//! - A piece lying in the face plane bounds material on one side. When
//!   the clipper gives it to the cell on the other side (a solid face
//!   on a grid plane), that cell holds the piece with no volume: it
//!   bonds to the material its normal points away from, so loads and
//!   supports on the piece act on the material's dofs.
//! - A piece lying on another wall of both cells (a solid face on a grid
//!   plane, crossing the face) continues as one planar sheet when the
//!   two cells' pieces share an edge on the face: the material behind
//!   the sheet is connected there. Only pieces whose fragment is that
//!   material count: in a face-only cell (bonded to it), or with the
//!   material inside the piece's own cell. A piece held by the cell
//!   beyond its material (a slot's root wall on a grid plane, in the
//!   prong's cell) carries that cell's fragment and would weld across
//!   the slot. Face-only cells along such a face touch this way.
//! - A face without crossing edges is all inside or all outside except
//!   where its material is bounded by walls on grid planes; a 4 × 4
//!   grid of face samples decides.

use std::collections::{BTreeMap, BTreeSet};

use crate::clip::ClippedSolid;
use crate::grid::Grid;

use super::super::lattice::contains;
use super::{CellLabels, label_at};

/// Shared inputs of the touch evaluation.
pub(super) struct Context<'a> {
    /// The grid.
    pub grid: &'a Grid,
    /// The clipped solid.
    pub clipped: &'a ClippedSolid,
    /// Lattice labels of cells that were labelled.
    pub labels: &'a BTreeMap<usize, CellLabels>,
    /// Lattice depth.
    pub depth: u8,
    /// Fragment of each piece.
    pub piece_fragment: &'a [u8],
    /// Whether a cell is split into several entries.
    pub split: &'a dyn Fn(usize) -> bool,
    /// Whether a cell holds pieces but no volume (a solid face on a grid
    /// plane, clipped into the cell beyond it).
    pub face_only: &'a dyn Fn(usize) -> bool,
}

/// A face between `lower` and `upper` (cell index, grid coordinates),
/// normal to `axis`.
pub(super) struct Face {
    /// The cell below the face.
    pub lower: (usize, [usize; 3]),
    /// The cell above the face.
    pub upper: (usize, [usize; 3]),
    /// The face's normal axis.
    pub axis: usize,
}

/// A piece edge on the face plane: its end points and fragment.
type FaceEdge = ([f64; 3], [f64; 3], u8);

/// Pairs (lower fragment, upper fragment) touching across `face`.
pub(super) fn touching(ctx: &Context<'_>, face: Face) -> BTreeSet<(u8, u8)> {
    let Face {
        lower: (c, ijk),
        upper: (d, o),
        axis,
    } = face;
    let grid = ctx.grid;
    let h = grid.cell_size();
    let tol = 1.0e-9 * h;
    let face_lo = grid.cell_min(o);
    let plane = face_lo[axis];
    // A split cell's fragment at `p` is its nearest label.
    let owner = |cell: usize, cijk: [usize; 3], p: [f64; 3]| -> u8 {
        ctx.labels
            .get(&cell)
            .filter(|_| (ctx.split)(cell))
            .map_or(0, |l| {
                label_at(grid, cijk, ctx.depth, &l.lattice, &l.slivers, p).saturating_sub(1)
            })
    };
    let mut bonded: BTreeSet<(u8, u8)> = BTreeSet::new();
    let face_hi = face_lo.map(|x| x + h);
    // A piece on a wall of both cells whose fragment is the material
    // behind it.
    let sheet_piece = |cell: usize, v: &[[f64; 3]], normal: [f64; 3]| {
        (0..3).filter(|&u| u != axis).any(|u| {
            [(face_lo[u], -1.0), (face_hi[u], 1.0)]
                .iter()
                .any(|&(f, out)| {
                    v.iter().all(|x| (x[u] - f).abs() <= tol)
                        && ((ctx.face_only)(cell) || normal[u] * out > 0.0)
                })
        })
    };
    // An edge along the face's own boundary (a wall lying on a grid
    // plane, clipped at the face) is contact along a line, not across
    // the face.
    let along_face_boundary = |a: [f64; 3], b: [f64; 3]| {
        (0..3).filter(|&u| u != axis).any(|u| {
            [face_lo[u], face_hi[u]]
                .iter()
                .any(|&f| (a[u] - f).abs() <= tol && (b[u] - f).abs() <= tol)
        })
    };
    let mut on_face = |cell: usize, upper: bool| -> (Vec<FaceEdge>, Vec<FaceEdge>) {
        let pieces = &ctx.clipped.pieces;
        let lo = pieces.partition_point(|p| p.cell < cell);
        let hi = pieces.partition_point(|p| p.cell <= cell);
        let (mut edges, mut sheet) = (Vec::new(), Vec::new());
        for (i, piece) in pieces[lo..hi].iter().enumerate() {
            let v = piece.polygon.vertices();
            let fragment = ctx.piece_fragment.get(lo + i).copied().unwrap_or(0);
            if v.iter().all(|x| (x[axis] - plane).abs() <= tol) {
                // In the face plane: material lies where the normal
                // points away from. Held by the cell on the other side,
                // the piece bonds to that material.
                let n = piece.normal[axis];
                let nv = v.len().max(1) as f64;
                let mut p = [0, 1, 2].map(|a| v.iter().map(|x| x[a]).sum::<f64>() / nv);
                if upper && n > 0.0 {
                    p[axis] -= 1.0e-6 * h;
                    bonded.insert((owner(c, ijk, p), fragment));
                } else if !upper && n < 0.0 {
                    p[axis] += 1.0e-6 * h;
                    bonded.insert((fragment, owner(d, o, p)));
                } else {
                    // Its material is in this cell: its edges inside the
                    // face bound where that material crosses it (a
                    // pocket's end on the grid plane, the material beside
                    // the pocket crossing below it; PH-067).
                    for k in 0..v.len() {
                        let (a, b) = (v[k], v[(k + 1) % v.len()]);
                        if !along_face_boundary(a, b) {
                            edges.push((a, b, fragment));
                        }
                    }
                }
                continue;
            }
            for k in 0..v.len() {
                let (a, b) = (v[k], v[(k + 1) % v.len()]);
                if (a[axis] - plane).abs() > tol || (b[axis] - plane).abs() > tol {
                    continue;
                }
                if !along_face_boundary(a, b) {
                    edges.push((a, b, fragment));
                } else if sheet_piece(cell, v, piece.normal) {
                    sheet.push((a, b, fragment));
                }
            }
        }
        (edges, sheet)
    };
    let ((ec, sc), (ed, sd)) = (on_face(c, false), on_face(d, true));
    let same = |p: [f64; 3], q: [f64; 3]| (0..3).all(|a| (p[a] - q[a]).abs() <= tol);
    let mut touching: BTreeSet<(u8, u8)> = BTreeSet::new();
    for (lower, upper) in [(&ec, &ed), (&sc, &sd)] {
        for (a1, b1, f1) in lower {
            for (a2, b2, f2) in upper {
                if (same(*a1, *a2) && same(*b1, *b2)) || (same(*a1, *b2) && same(*b1, *a2)) {
                    touching.insert((*f1, *f2));
                }
            }
        }
    }
    if touching.is_empty() {
        // Samples just inside the upper cell, off the column boundary
        // the containment ray's half-open rule decides.
        let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
        for j in 0..4 {
            for i in 0..4 {
                let mut p = face_lo;
                p[axis] += 1.0e-6 * h;
                p[u] += (i as f64 + 0.5) * 0.25 * h;
                p[v] += (j as f64 + 0.5) * 0.25 * h;
                if contains(grid, ctx.clipped, o, p) {
                    touching.insert((owner(c, ijk, p), owner(d, o, p)));
                }
            }
        }
    }
    // Slivers are thinner than the face samples' spacing: one reaching
    // the face touches the material across it at the middle of each of
    // its material boxes' footprints on the face, where there is any
    // (PH-060).
    let lattice_count = |l: &CellLabels| l.lattice.iter().copied().max().unwrap_or(0);
    for (cell, other, oijk, step, upper) in [(c, d, o, 1.0, false), (d, c, ijk, -1.0, true)] {
        let Some(l) = ctx.labels.get(&cell).filter(|_| (ctx.split)(cell)) else {
            continue;
        };
        for (k, sl) in l.slivers.iter().enumerate() {
            for (lo, hi) in sl.boxes_on(axis, plane, tol) {
                let mut p = [0, 1, 2].map(|a| 0.5 * (lo[a] + hi[a]));
                p[axis] = plane + step * 1.0e-6 * h;
                if !contains(grid, ctx.clipped, oijk, p) {
                    continue;
                }
                let mine = lattice_count(l).saturating_add(k as u8);
                let theirs = owner(other, oijk, p);
                touching.insert(if upper {
                    (theirs, mine)
                } else {
                    (mine, theirs)
                });
            }
        }
    }
    touching.extend(bonded);
    touching
}
