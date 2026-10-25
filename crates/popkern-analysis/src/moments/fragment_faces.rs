//! Cells holding only pieces of a solid face that lies on one of their
//! faces: a face on a grid plane, given by the clipper to the empty cell
//! beyond it. Such a cell has no volume; its pieces load and support the
//! material across the face. When that material is split, the cell is
//! split the same way, one zero-volume entry per fragment its pieces
//! bound, so each piece bonds to its own fragment and the empty cell
//! never welds two fragments.

use std::collections::BTreeMap;

use crate::clip::ClippedSolid;
use crate::grid::Grid;

use super::super::CellMoments;
use super::{CellLabels, label_at};

/// Volume below which a cell holds no material, relative to h³.
const EMPTY_VOLUME: f64 = 1.0e-12;

/// For each face-only cell whose pieces bound more than one fragment,
/// the fragment index of each of its pieces (in piece order).
pub(super) fn face_only_splits(
    grid: &Grid,
    clipped: &ClippedSolid,
    cells: &[CellMoments],
    labels: &BTreeMap<usize, CellLabels>,
    depth: u8,
) -> BTreeMap<usize, Vec<u8>> {
    let h = grid.cell_size();
    let tol = 1.0e-9 * h;
    let pieces = &clipped.pieces;
    let mut out = BTreeMap::new();
    for c in cells {
        if c.m[0].abs() > EMPTY_VOLUME * h * h * h {
            continue;
        }
        let Some(ijk) = grid.cell_ijk(c.cell) else {
            continue;
        };
        let lo = grid.cell_min(ijk);
        let first = pieces.partition_point(|p| p.cell < c.cell);
        let end = pieces.partition_point(|p| p.cell <= c.cell);
        // The fragment of the material each piece bounds, or None when a
        // piece is not on a face with the material beyond it.
        let owners: Option<Vec<u8>> = pieces[first..end]
            .iter()
            .map(|piece| {
                let v = piece.polygon.vertices();
                let axis = (0..3).find(|&a| {
                    v.iter().all(|x| (x[a] - lo[a]).abs() <= tol)
                        || v.iter().all(|x| (x[a] - lo[a] - h).abs() <= tol)
                })?;
                let lower = v.iter().all(|x| (x[axis] - lo[axis]).abs() <= tol);
                let n = piece.normal[axis];
                // Material beyond the face: the normal points into this
                // cell.
                let step = match (lower, n > 0.0, n < 0.0) {
                    (true, true, _) => -1i64,
                    (false, _, true) => 1,
                    _ => return None,
                };
                let mut o = ijk;
                o[axis] = (ijk[axis] as i64 + step).try_into().ok()?;
                let neighbour = grid.cell_index(o)?;
                let Some(l) = labels.get(&neighbour) else {
                    return Some(0);
                };
                if l.count() < 2 {
                    return Some(0);
                }
                let nv = v.len().max(1) as f64;
                let mut p = [0, 1, 2].map(|a| v.iter().map(|x| x[a]).sum::<f64>() / nv);
                p[axis] += step as f64 * 1.0e-6 * h;
                Some(label_at(grid, o, depth, &l.lattice, &l.slivers, p).saturating_sub(1))
            })
            .collect();
        let Some(owners) = owners else {
            continue;
        };
        let mut distinct = owners.clone();
        distinct.sort_unstable();
        distinct.dedup();
        if distinct.len() < 2 {
            continue;
        }
        let index: Vec<u8> = owners
            .iter()
            .map(|o| distinct.iter().position(|d| d == o).unwrap_or(0) as u8)
            .collect();
        out.insert(c.cell, index);
    }
    out
}
