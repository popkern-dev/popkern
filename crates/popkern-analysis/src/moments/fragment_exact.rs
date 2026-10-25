//! Exact moments of the fragments of a split cell.
//!
//! Along any axis d the divergence theorem gives, for each fragment,
//!
//! m_e = (h/2) [ Σ_pieces ∫ G(ξ_d) Π ξ_o^(e_o) n_d dA + c_(e_d) · P ],
//!
//! closed on one wall of the cell: the upper wall with
//! G_a(ξ) = (ξ^(a+1) − (−1)^(a+1)) / (a+1), or the lower wall with
//! G_a(ξ) = (ξ^(a+1) − 1) / (a+1); c = (2, 0, 2/3) for either. P is the
//! fragment's section moments on that wall. The cell's exact total
//! moments give the whole wall's section, so when only one fragment has
//! material on the wall it owns P whole, the others none, and every
//! fragment's moments are exact. A fragment has material on a wall when
//! one of its pieces has an edge crossing the wall (not along the
//! wall's boundary) or, for a wall no piece crosses, when its lattice
//! points reach the wall's layer. A wall no fragment reaches (its whole
//! section, from the total, is empty) closes every fragment with P = 0.
//! When every wall is reached by several fragments, one wall's section
//! is split by fragment in the wall's plane (`fragment_wall`).

use crate::clip::ClippedSolid;
use crate::grid::Grid;

use super::super::triangle::monomial_integrals;
use super::super::{CellMoments, consistent, local, moment_index};
use super::wall::{EMPTY, WallInput, split_section};
use super::{Fragments, component_moments};

/// Relative (to h³) disagreement of piece-only moments between axes
/// above which a cell's piece assignment is not trusted.
const PIECE_AGREEMENT: f64 = 1.0e-10;

/// c_a: the wall weight of G_a at the closing wall.
const WALL_WEIGHT: [f64; 3] = [2.0, 0.0, 2.0 / 3.0];

/// Moment index of exponent `a` along `axis` and (b, c) along the next
/// two axes in cyclic order.
fn permuted_index(axis: usize, a: usize, b: usize, c: usize) -> usize {
    let mut e = [0; 3];
    e[axis] = a;
    e[(axis + 1) % 3] = b;
    e[(axis + 2) % 3] = c;
    moment_index(e[0], e[1], e[2])
}

/// Piece integrals ∫ G(ξ_d) Π ξ_o^(e_o) n_d dA per fragment, indexed
/// `[axis][side][fragment]` (side 0 upper, 1 lower).
fn piece_terms(
    grid: &Grid,
    clipped: &ClippedSolid,
    piece_fragment: &[u8],
    cell: usize,
    ijk: [usize; 3],
    count: usize,
) -> Vec<[[[f64; 27]; 2]; 3]> {
    let h = grid.cell_size();
    let lo = grid.cell_min(ijk);
    let pieces = &clipped.pieces;
    let first = pieces.partition_point(|p| p.cell < cell);
    let end = pieces.partition_point(|p| p.cell <= cell);
    let mut out = vec![[[[0.0; 27]; 2]; 3]; count];
    for (i, piece) in pieces[first..end].iter().enumerate() {
        let f = piece_fragment.get(first + i).copied().unwrap_or(0) as usize;
        let Some(q) = out.get_mut(f) else {
            continue;
        };
        let v = piece.polygon.vertices();
        for t in 1..v.len().saturating_sub(1) {
            let tri = [v[0], v[t], v[t + 1]];
            let u = [0, 1, 2].map(|k| tri[1][k] - tri[0][k]);
            let w = [0, 1, 2].map(|k| tri[2][k] - tri[0][k]);
            let x = [
                u[1] * w[2] - u[2] * w[1],
                u[2] * w[0] - u[0] * w[2],
                u[0] * w[1] - u[1] * w[0],
            ];
            let area = 0.5 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
            if area == 0.0 {
                continue;
            }
            let l = tri.map(|p| local(p, lo, h));
            for (axis, qa) in q.iter_mut().enumerate() {
                let nd = piece.normal[axis];
                if nd == 0.0 {
                    continue;
                }
                let perm = l.map(|p| [p[axis], p[(axis + 1) % 3], p[(axis + 2) % 3]]);
                let m = monomial_integrals(perm, area);
                let others = (0..3).flat_map(|c| (0..3).map(move |b| (b, c)));
                for (b, c) in others {
                    for a in 0..3 {
                        let k = (a + 1) as f64;
                        let odd = if a % 2 == 0 { -1.0 } else { 1.0 };
                        let up = (m[a + 1][b][c] - odd * m[0][b][c]) / k;
                        let down = (m[a + 1][b][c] - m[0][b][c]) / k;
                        let idx = permuted_index(axis, a, b, c);
                        qa[0][idx] += nd * up;
                        qa[1][idx] += nd * down;
                    }
                }
            }
        }
    }
    out
}

/// Fragments with material on each wall, indexed `[axis][side]`, as a
/// bit set over fragments.
fn wall_owners(
    grid: &Grid,
    clipped: &ClippedSolid,
    piece_fragment: &[u8],
    s: &SplitInput<'_>,
) -> [[u32; 2]; 3] {
    let (cell, ijk, labels, slivers, depth) = (s.cell, s.ijk, s.labels, s.slivers, s.depth);
    let h = grid.cell_size();
    let tol = 1.0e-9 * h;
    let lo = grid.cell_min(ijk);
    let n = 1usize << depth;
    let pieces = &clipped.pieces;
    let first = pieces.partition_point(|p| p.cell < cell);
    let end = pieces.partition_point(|p| p.cell <= cell);
    let mut owners = [[0u32; 2]; 3];
    for (axis, row) in owners.iter_mut().enumerate() {
        for (side, bits) in row.iter_mut().enumerate() {
            let plane = lo[axis] + if side == 0 { h } else { 0.0 };
            let mut crossed = false;
            for (i, piece) in pieces[first..end].iter().enumerate() {
                let v = piece.polygon.vertices();
                if v.iter().all(|x| (x[axis] - plane).abs() <= tol) {
                    continue;
                }
                for k in 0..v.len() {
                    let (a, b) = (v[k], v[(k + 1) % v.len()]);
                    if (a[axis] - plane).abs() > tol || (b[axis] - plane).abs() > tol {
                        continue;
                    }
                    let along_boundary = (0..3).filter(|&u| u != axis).any(|u| {
                        [lo[u], lo[u] + h]
                            .iter()
                            .any(|&f| (a[u] - f).abs() <= tol && (b[u] - f).abs() <= tol)
                    });
                    if !along_boundary {
                        crossed = true;
                        let f = piece_fragment.get(first + i).copied().unwrap_or(0);
                        *bits |= 1 << f;
                    }
                }
            }
            if crossed {
                continue;
            }
            // No piece crosses the wall: it is all inside or all
            // outside, and the lattice layer next to it says which
            // fragment, if any, reaches it.
            // Slivers have no lattice points: one reaches the wall when
            // one of its material boxes does.
            let lattice_count = labels.iter().copied().max().unwrap_or(0) as usize;
            for (k, sl) in slivers.iter().enumerate() {
                if sl.boxes_on(axis, plane, tol).next().is_some() {
                    *bits |= 1 << (lattice_count + k);
                }
            }
            let layer = if side == 0 { n - 1 } else { 0 };
            for (idx, &l) in labels.iter().enumerate() {
                let at = [idx % n, (idx / n) % n, idx / (n * n)];
                if l > 0 && at[axis] == layer {
                    *bits |= 1 << (l - 1);
                }
            }
        }
    }
    owners
}

/// The inputs of one split cell.
struct SplitInput<'a> {
    /// Linear cell index.
    cell: usize,
    /// Its grid coordinates.
    ijk: [usize; 3],
    /// The cell's exact total moments.
    total: &'a [f64; 27],
    /// Lattice labels of the cell.
    labels: &'a [u8],
    /// Slivers of the cell.
    slivers: &'a [super::Sliver],
    /// Lattice depth.
    depth: u8,
    /// Lattice (midpoint-rule) moments of each fragment.
    lattice: &'a [[f64; 27]],
}

/// Exact moments of every fragment of a split cell, closed on the
/// first wall (x, y, z; upper before lower) that at most one fragment
/// reaches, else on the first whose section splits in its plane;
/// `None` when no wall qualifies or every result leaves some
/// fragment's lattice volume by more than one lattice layer.
fn exact_moments(
    grid: &Grid,
    clipped: &ClippedSolid,
    piece_fragment: &[u8],
    s: &SplitInput<'_>,
) -> Option<Vec<[f64; 27]>> {
    let h = grid.cell_size();
    let count = s.lattice.len();
    let terms = piece_terms(grid, clipped, piece_fragment, s.cell, s.ijk, count);
    let owners = wall_owners(grid, clipped, piece_fragment, s);
    let layer = h / (1usize << s.depth) as f64;
    // Moments with exponent 1 along an axis need no wall term, so two
    // axes give each such moment twice from the pieces alone: they agree
    // only when every piece bounds its own fragment.
    for t in &terms {
        for (e, ex) in (0..27usize).map(|e| (e, [e % 3, (e / 3) % 3, e / 9])) {
            let mut seen: Option<f64> = None;
            for (_, sides) in t.iter().enumerate().filter(|&(d, _)| ex[d] == 1) {
                let v = 0.5 * h * sides[0][e];
                if seen.is_some_and(|w| (w - v).abs() > PIECE_AGREEMENT * h * h * h) {
                    return None;
                }
                seen = Some(v);
            }
        }
    }
    // Closing walls: first those at most one fragment reaches (its
    // owner takes the whole section), then any whose section splits by
    // fragment in the wall's plane (`fragment_wall`).
    for split_in_plane in [false, true] {
        for axis in 0..3 {
            for side in 0..2 {
                // The whole wall's section moments from the exact total.
                let mut wall = [0.0; 9];
                for c in 0..3 {
                    for b in 0..3 {
                        let idx = permuted_index(axis, 0, b, c);
                        let q: f64 = terms.iter().map(|t| t[axis][side][idx]).sum();
                        wall[b + 3 * c] = (2.0 * s.total[idx] / h - q) / WALL_WEIGHT[0];
                    }
                }
                let bits = owners[axis][side];
                let sections: Vec<[f64; 9]> = match (split_in_plane, bits.count_ones()) {
                    (false, 1) => (0..count)
                        .map(|f| if 1 << f == bits { wall } else { [0.0; 9] })
                        .collect(),
                    (false, 0) if wall[0].abs() <= EMPTY * h * h => vec![[0.0; 9]; count],
                    (false, _) => continue,
                    (true, _) => {
                        let input = WallInput {
                            cell: s.cell,
                            ijk: s.ijk,
                            axis,
                            side,
                            total: &wall,
                            labels: s.labels,
                            slivers: s.slivers,
                            depth: s.depth,
                            count,
                        };
                        match split_section(grid, clipped, piece_fragment, &input) {
                            Some(sections) => sections,
                            None => continue,
                        }
                    }
                };
                let parts: Vec<[f64; 27]> = (0..count)
                    .map(|f| {
                        let mut m = [0.0; 27];
                        for c in 0..3 {
                            for b in 0..3 {
                                for (a, w) in WALL_WEIGHT.iter().enumerate() {
                                    let idx = permuted_index(axis, a, b, c);
                                    let p = w * sections[f][b + 3 * c];
                                    m[idx] = 0.5 * h * (terms[f][axis][side][idx] + p);
                                }
                            }
                        }
                        m
                    })
                    .collect();
                let plausible = parts
                    .iter()
                    .zip(s.lattice)
                    .all(|(m, l)| consistent(m, h) && (m[0] - l[0]).abs() <= layer * h * h);
                if plausible {
                    return Some(parts);
                }
            }
        }
    }
    None
}

/// Replaces the fragment moments of every split cell (except `skip`
/// ones, which have exact moments already) by exact ones where a wall
/// qualifies; counts the cells that keep lattice moments.
pub(super) fn apply(
    grid: &Grid,
    clipped: &ClippedSolid,
    cells: &[CellMoments],
    entries: &mut [CellMoments],
    frags: &mut Fragments,
    skip: &dyn Fn(usize) -> bool,
) {
    let depth = frags.depth;
    for sc in frags.split.iter().filter(|sc| !skip(sc.cell)) {
        let Some(ijk) = grid.cell_ijk(sc.cell) else {
            continue;
        };
        let Ok(k) = cells.binary_search_by_key(&sc.cell, |c| c.cell) else {
            continue;
        };
        let lo = entries.partition_point(|e| e.cell < sc.cell);
        let hi = entries.partition_point(|e| e.cell <= sc.cell);
        let lattice: Vec<[f64; 27]> = (1..=hi - lo)
            .map(|t| component_moments(grid, ijk, depth, &sc.labels, t as u8))
            .collect();
        let input = SplitInput {
            cell: sc.cell,
            ijk,
            total: &cells[k].m,
            labels: &sc.labels,
            slivers: &sc.slivers,
            depth,
            lattice: &lattice,
        };
        match exact_moments(grid, clipped, &frags.piece_fragment, &input) {
            Some(parts) => {
                for (e, m) in entries[lo..hi].iter_mut().zip(parts) {
                    e.m = m;
                }
            }
            None => frags.lattice_moment_cells += 1,
        }
    }
}
