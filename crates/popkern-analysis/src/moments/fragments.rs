//! Disconnected parts of one solid closer than a cell (cross-cutting
//! Thread 3, step 2): a slot narrower than a cell must not weld the
//! material on either side through shared dofs.
//!
//! - A cut cell whose inside lattice points form several 6-connected
//!   components is split into one entry per component (a fragment).
//!   Each fragment's moments come from the lattice midpoint rule over
//!   its component, plus a share of the lattice's residual against the
//!   exact cell moments by lattice volume, so the fragments sum to the
//!   exact moments.
//! - Material on two sides of a cell face touches when lattice points
//!   of both entries face each other across it. Touching is evaluated
//!   for every face pair among candidate cells: cut cells holding
//!   pieces that face each other across empty space (opposing normals,
//!   each in front of the other: a slot's walls, not a thin wall's
//!   faces), inside the cell or across a face (a slot centred on a cell
//!   face splits no cell but separates two).
//!   Elsewhere, face neighbours touch.
//! - Each grid node gets one dof copy per group of entries around it
//!   that touch through face neighbours within the node's eight cells.

use std::collections::{BTreeMap, BTreeSet};

use crate::clip::ClippedSolid;
use crate::grid::Grid;

use super::lattice::{component_labels, lattice};
use super::{CellKind, CellMoments, moment_index};

#[path = "fragment_access.rs"]
mod access;
#[path = "fragment_copies.rs"]
mod copies;
#[path = "fragment_exact.rs"]
mod exact;
#[path = "fragment_faces.rs"]
mod faces;
#[path = "fragment_pieces.rs"]
mod pieces;
#[path = "fragment_ray.rs"]
mod ray;
#[path = "fragment_touch.rs"]
mod touch;
#[path = "fragment_wall.rs"]
mod wall;

/// Most dof copies of one grid node; node keys are node · this + copy.
pub const NODE_COPIES: usize = 8;

/// Lattice depth for cells that may split (at least the solid's own):
/// a fragment's moments come from this lattice, so its resolution sets
/// their error.
pub const SPLIT_DEPTH: u8 = 4;

/// A cell split into fragments.
#[derive(Debug, Clone, PartialEq)]
pub struct SplitCell {
    /// Linear cell index.
    pub cell: usize,
    /// Lattice labels (0 outside, k + 1 for fragment k), x fastest.
    pub labels: Vec<u8>,
    /// Material no lattice point sees, fragments after the lattice's.
    pub slivers: Vec<Sliver>,
}

/// Material of a split cell thinner than the lattice spacing, separate
/// from the lattice's fragments (PH-060).
#[derive(Debug, Clone, PartialEq)]
pub struct Sliver {
    /// Indices of its pieces in the clipped solid, ascending.
    pub pieces: Vec<usize>,
    /// Lower corner of its material's bounding box.
    pub lo: [f64; 3],
    /// Upper corner of its material's bounding box.
    pub hi: [f64; 3],
    /// The material box of each of its pieces (the piece's vertices and
    /// those stepped back behind it): their union is its material, so a
    /// strip wrapping round other material does not claim it (PH-062).
    pub boxes: Vec<([f64; 3], [f64; 3])>,
}

/// The fragment label (k + 1 for fragment k, 0 for none) at `p` in a
/// cell with lattice `labels` and `slivers`: a sliver's when `p` lies in
/// one of its material boxes (to 1e-6 h), else the nearest labelled lattice
/// point's. Callers look just inside the material (1e-6 h): a wider
/// margin would hand the resolved material beside a sliver to it.
pub(crate) fn label_at(
    grid: &Grid,
    ijk: [usize; 3],
    depth: u8,
    labels: &[u8],
    slivers: &[Sliver],
    p: [f64; 3],
) -> u8 {
    let lattice_count = labels.iter().copied().max().unwrap_or(0);
    match sliver_at(grid, slivers, p) {
        Some(k) => lattice_count.saturating_add(k as u8 + 1),
        None => nearest_label(grid, ijk, depth, labels, p),
    }
}

/// The sliver whose material boxes hold `p` (to 1e-6 h), if any.
pub(crate) fn sliver_at(grid: &Grid, slivers: &[Sliver], p: [f64; 3]) -> Option<usize> {
    let margin = 1.0e-6 * grid.cell_size();
    slivers.iter().position(|s| s.holds(p, margin))
}

/// Lattice labels and slivers of a labelled cell.
#[derive(Debug, Clone, PartialEq, Default)]
pub(crate) struct CellLabels {
    /// Lattice labels, x fastest.
    pub lattice: Vec<u8>,
    /// Slivers, fragments after the lattice's.
    pub slivers: Vec<Sliver>,
}

impl CellLabels {
    /// Fragments: the lattice's and the slivers.
    pub fn count(&self) -> usize {
        self.lattice.iter().copied().max().unwrap_or(0) as usize + self.slivers.len()
    }
}

/// Why fragments could not be formed.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum FragmentError {
    /// A cell or node holds more fragments than the key encoding allows.
    #[error("cell or node {0} holds more than {NODE_COPIES} separate fragments")]
    TooMany(usize),
}

/// Fragment data of one solid, stored on [`SolidMoments`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Fragments {
    /// Lattice depth used.
    pub depth: u8,
    /// Split cells, by cell index.
    pub split: Vec<SplitCell>,
    /// Fragment of each clipped piece within its cell (parallel to the
    /// pieces).
    pub piece_fragment: Vec<u8>,
    /// Evaluated face pairs of entries (a < b) and whether they touch.
    pub pair_touch: Vec<((usize, usize), bool)>,
    /// Dof copy of each corner node of each entry.
    pub corner_copy: Vec<[u8; 8]>,
    /// Split cells whose fragments keep lattice moments plus a residual
    /// share, because no cell wall's section could be split by fragment
    /// (none on any fixture or gallery part so far).
    pub lattice_moment_cells: usize,
}

fn ipow(x: f64, n: usize) -> f64 {
    let mut r = 1.0;
    for _ in 0..n {
        r *= x;
    }
    r
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Midpoint-rule moments of the lattice points labelled `tag`.
fn component_moments(grid: &Grid, ijk: [usize; 3], depth: u8, labels: &[u8], tag: u8) -> [f64; 27] {
    let h = grid.cell_size();
    let s = h / (1usize << depth) as f64;
    let w = s * s * s;
    let lo = grid.cell_min(ijk);
    let mut m = [0.0; 27];
    for (i, (_, p)) in lattice(grid, ijk, depth).into_iter().enumerate() {
        if labels[i] != tag {
            continue;
        }
        let l = [0, 1, 2].map(|a| 2.0 * (p[a] - lo[a]) / h - 1.0);
        for c in 0..3 {
            for b in 0..3 {
                for a in 0..3 {
                    m[moment_index(a, b, c)] += w * ipow(l[0], a) * ipow(l[1], b) * ipow(l[2], c);
                }
            }
        }
    }
    m
}

/// (centroid, unit normal) of the pieces of `cell`.
fn faces_of(clipped: &ClippedSolid, cell: usize) -> Vec<([f64; 3], [f64; 3])> {
    let pieces = &clipped.pieces;
    let lo = pieces.partition_point(|p| p.cell < cell);
    let hi = pieces.partition_point(|p| p.cell <= cell);
    pieces[lo..hi]
        .iter()
        .map(|p| {
            let v = p.polygon.vertices();
            let n = v.len().max(1) as f64;
            (
                [0, 1, 2].map(|a| v.iter().map(|x| x[a]).sum::<f64>() / n),
                p.normal,
            )
        })
        .collect()
}

/// Whether two pieces face each other across empty space: opposing
/// normals, each in front of the other (a slot's walls). A thin wall's
/// two faces oppose too but point away from each other.
fn facing(a: &[([f64; 3], [f64; 3])], b: &[([f64; 3], [f64; 3])]) -> bool {
    a.iter().any(|(ca, na)| {
        b.iter().any(|(cb, nb)| {
            let d = [cb[0] - ca[0], cb[1] - ca[1], cb[2] - ca[2]];
            dot(*na, *nb) < 0.0 && dot(d, *na) > 0.0 && dot(d, *nb) < 0.0
        })
    })
}

/// The nearest labelled lattice point of a cell to `p`: its label.
pub(crate) fn nearest_label(
    grid: &Grid,
    ijk: [usize; 3],
    depth: u8,
    labels: &[u8],
    p: [f64; 3],
) -> u8 {
    let mut best = (f64::INFINITY, 0u8);
    for (i, (_, x)) in lattice(grid, ijk, depth).into_iter().enumerate() {
        if labels[i] == 0 {
            continue;
        }
        let d = (0..3).map(|a| (x[a] - p[a]) * (x[a] - p[a])).sum::<f64>();
        if d < best.0 {
            best = (d, labels[i]);
        }
    }
    best.1
}

/// Union-find root with path halving.
fn find(p: &mut [usize], mut x: usize) -> usize {
    while p[x] != x {
        p[x] = p[p[x]];
        x = p[x];
    }
    x
}

/// Splits `cells` (one entry per active cell, ascending) into fragment
/// entries and finds the dof copies.
pub(super) fn fragment(
    grid: &Grid,
    clipped: &ClippedSolid,
    cells: Vec<CellMoments>,
    depth: u8,
    exact_moments: bool,
) -> Result<(Vec<CellMoments>, Fragments), FragmentError> {
    let depth = depth.max(SPLIT_DEPTH);
    let n = 1usize << depth;
    let index_of = |cell: usize| cells.binary_search_by_key(&cell, |c| c.cell).ok();
    // Candidates: cut cells with opposing normals inside or across a face.
    let mut candidate: BTreeSet<usize> = BTreeSet::new();
    let neighbours = |cell: usize| -> Vec<usize> {
        let Some(ijk) = grid.cell_ijk(cell) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for axis in 0..3 {
            for step in [-1i64, 1] {
                let v = ijk[axis] as i64 + step;
                if v < 0 {
                    continue;
                }
                let mut o = ijk;
                o[axis] = v as usize;
                if let Some(c) = grid.cell_index(o) {
                    out.push(c);
                }
            }
        }
        out
    };
    let mut internal: BTreeSet<usize> = BTreeSet::new();
    for c in cells.iter().filter(|c| c.kind == CellKind::Cut) {
        let mine = faces_of(clipped, c.cell);
        if facing(&mine, &mine) {
            candidate.insert(c.cell);
            internal.insert(c.cell);
        }
        for d in neighbours(c.cell) {
            if index_of(d).is_some_and(|k| cells[k].kind == CellKind::Cut)
                && facing(&mine, &faces_of(clipped, d))
            {
                candidate.insert(c.cell);
                candidate.insert(d);
            }
        }
    }
    // Lattice labels of cells with opposing pieces inside: only these
    // can hold two fragments.
    let mut labels: BTreeMap<usize, CellLabels> = BTreeMap::new();
    let mut need: BTreeSet<usize> = candidate.clone();
    for &c in &candidate {
        need.extend(neighbours(c).into_iter().filter(|d| index_of(*d).is_some()));
    }
    for &c in &internal {
        let Some(ijk) = grid.cell_ijk(c) else {
            continue;
        };
        let (mut l, _) = component_labels(grid, clipped, ijk, depth);
        pieces::merge_by_surface(grid, clipped, c, ijk, depth, &mut l);
        let slivers = pieces::slivers(grid, clipped, c, ijk, depth, &l);
        let cell = CellLabels {
            lattice: l,
            slivers,
        };
        if cell.count() > NODE_COPIES {
            return Err(FragmentError::TooMany(c));
        }
        labels.insert(c, cell);
    }
    // Entries and fragment moments.
    let face_only = faces::face_only_splits(grid, clipped, &cells, &labels, depth);
    // Their face pairs are evaluated like a candidate's (an unevaluated
    // pair touches, which would weld the entries through a neighbour).
    for &c in face_only.keys() {
        need.insert(c);
        need.extend(neighbours(c).into_iter().filter(|d| index_of(*d).is_some()));
    }
    let mut entries = Vec::with_capacity(cells.len());
    let mut split = Vec::new();
    for c in &cells {
        if let Some(index) = face_only.get(&c.cell) {
            // No volume: the first entry keeps the (round-off) moments.
            let count = index.iter().copied().max().unwrap_or(0) + 1;
            for k in 0..count {
                let m = if k == 0 { c.m } else { [0.0; 27] };
                entries.push(CellMoments {
                    m,
                    fragment: k,
                    ..*c
                });
            }
            split.push(SplitCell {
                cell: c.cell,
                labels: vec![0; n * n * n],
                slivers: Vec::new(),
            });
            continue;
        }
        let count = labels.get(&c.cell).map_or(1, CellLabels::count);
        let (Some(cl), Some(ijk), true) = (labels.get(&c.cell), grid.cell_ijk(c.cell), count > 1)
        else {
            entries.push(CellMoments { fragment: 0, ..*c });
            continue;
        };
        let l = &cl.lattice;
        let mut parts: Vec<[f64; 27]> = (1..=count)
            .map(|t| component_moments(grid, ijk, depth, l, t as u8))
            .collect();
        // The lattice's residual against the exact cell moments, shared
        // by lattice volume, so the fragments sum to the exact moments
        // and each stays close to its own (positive) lattice measure.
        let total: f64 = parts.iter().map(|p| p[0]).sum();
        let mut residual = c.m;
        for part in &parts {
            for (r, x) in residual.iter_mut().zip(part) {
                *r -= x;
            }
        }
        let last = count - 1;
        let mut assigned = [0.0; 27];
        for (k, part) in parts.iter_mut().enumerate() {
            if k == last {
                // The last takes what is left, so the sum is exact.
                for ((x, m), a) in part.iter_mut().zip(&c.m).zip(&assigned) {
                    *x = m - a;
                }
                break;
            }
            let share = if total > 0.0 { part[0] / total } else { 0.0 };
            for ((x, r), a) in part.iter_mut().zip(&residual).zip(assigned.iter_mut()) {
                *x += share * r;
                *a += *x;
            }
        }
        for (k, m) in parts.into_iter().enumerate() {
            entries.push(CellMoments {
                cell: c.cell,
                kind: CellKind::Cut,
                m,
                fragment: k as u8,
            });
        }
        split.push(SplitCell {
            cell: c.cell,
            labels: l.clone(),
            slivers: cl.slivers.clone(),
        });
    }
    let mut out = Fragments {
        depth,
        split,
        piece_fragment: vec![0; clipped.pieces.len()],
        pair_touch: Vec::new(),
        corner_copy: vec![[0; 8]; entries.len()],
        lattice_moment_cells: 0,
    };
    // Pieces of split cells: face-only cells by the material they bound,
    // the others by their surface components (fragment_pieces).
    for (i, piece) in clipped.pieces.iter().enumerate() {
        if let Some(index) = face_only.get(&piece.cell) {
            let first = clipped.pieces.partition_point(|p| p.cell < piece.cell);
            out.piece_fragment[i] = index.get(i - first).copied().unwrap_or(0);
        }
    }
    let face_only_cell = |cell: usize| face_only.contains_key(&cell);
    pieces::assign(
        grid,
        clipped,
        &out.split,
        depth,
        &face_only_cell,
        &mut out.piece_fragment,
    );
    if exact_moments {
        let skip = |cell: usize| face_only.contains_key(&cell);
        exact::apply(grid, clipped, &cells, &mut entries, &mut out, &skip);
    } else {
        out.lattice_moment_cells = out.split.len() - face_only.len();
    }
    // Touching across faces among candidates and their neighbours:
    // material crosses a face exactly where clipped pieces have edges on
    // it, and a triangle crossing the face leaves the same segment in
    // both cells. A face without such edges is all inside or all
    // outside; its centre decides.
    let entry_range = |cell: usize| {
        let lo = entries.partition_point(|e: &CellMoments| e.cell < cell);
        let hi = entries.partition_point(|e: &CellMoments| e.cell <= cell);
        lo..hi
    };
    let split_cell = |cell: usize| entry_range(cell).len() > 1;
    let touch_context = touch::Context {
        grid,
        clipped,
        labels: &labels,
        depth,
        piece_fragment: &out.piece_fragment,
        split: &split_cell,
        face_only: &face_only_cell,
    };
    let mut pairs: BTreeMap<(usize, usize), bool> = BTreeMap::new();
    for &c in &need {
        let Some(ijk) = grid.cell_ijk(c) else {
            continue;
        };
        for axis in 0..3 {
            let mut o = ijk;
            o[axis] += 1;
            let Some(d) = grid.cell_index(o) else {
                continue;
            };
            if !need.contains(&d) || entry_range(d).is_empty() || entry_range(c).is_empty() {
                continue;
            }
            let touching = touch::touching(
                &touch_context,
                touch::Face {
                    lower: (c, ijk),
                    upper: (d, o),
                    axis,
                },
            );
            for a in entry_range(c) {
                for b in entry_range(d) {
                    let t = touching.contains(&(entries[a].fragment, entries[b].fragment));
                    pairs.insert((a.min(b), a.max(b)), t);
                }
            }
        }
    }
    out.pair_touch = pairs.into_iter().collect();
    // Dof copies per node around candidates.
    copies::assign(grid, &need, &entry_range, &mut out)?;
    Ok((entries, out))
}
