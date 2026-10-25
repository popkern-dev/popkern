//! Which fragment of a split cell each clipped piece bounds.
//!
//! A piece's own lattice label (the label just inside its centroid) is
//! wrong where the fragment is thinner than the lattice spacing. The
//! surface is more reliable: two pieces that share an edge inside the
//! cell lie on one sheet of the solid's boundary and bound the same
//! material. So the pieces of a split cell are grouped into edge-
//! connected components, and each component takes the area-weighted
//! majority of its pieces' labels.

use std::collections::BTreeMap;

use crate::clip::ClippedSolid;
use crate::grid::Grid;

use super::{NODE_COPIES, Sliver, SplitCell, find, nearest_label};

/// Quantised end point, for matching shared edges.
type Key = [i64; 3];

fn key(p: [f64; 3], tol: f64) -> Key {
    p.map(|c| (c / tol).round() as i64)
}

/// Union-find parents grouping the pieces of `cell` into edge-connected
/// surface components (edges on a cell wall are shared with the
/// neighbour cell, not inside this one); indices are relative to the
/// cell's first piece.
fn surface_components(
    grid: &Grid,
    clipped: &ClippedSolid,
    cell: usize,
    ijk: [usize; 3],
) -> Vec<usize> {
    let h = grid.cell_size();
    let tol = 1.0e-7 * h;
    let lo = grid.cell_min(ijk);
    let pieces = &clipped.pieces;
    let first = pieces.partition_point(|p| p.cell < cell);
    let end = pieces.partition_point(|p| p.cell <= cell);
    let on_wall = |a: [f64; 3], b: [f64; 3]| {
        (0..3).any(|u| {
            [lo[u], lo[u] + h]
                .iter()
                .any(|&w| (a[u] - w).abs() <= tol && (b[u] - w).abs() <= tol)
        })
    };
    let mut parent: Vec<usize> = (0..end - first).collect();
    let mut edges: BTreeMap<(Key, Key), usize> = BTreeMap::new();
    for (i, piece) in pieces[first..end].iter().enumerate() {
        let v = piece.polygon.vertices();
        for k in 0..v.len() {
            let (a, b) = (v[k], v[(k + 1) % v.len()]);
            if on_wall(a, b) {
                continue;
            }
            let (ka, kb) = (key(a, tol), key(b, tol));
            let e = if ka <= kb { (ka, kb) } else { (kb, ka) };
            match edges.get(&e) {
                Some(&j) => {
                    let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                    if ri != rj {
                        parent[ri.max(rj)] = ri.min(rj);
                    }
                }
                None => {
                    edges.insert(e, i);
                }
            }
        }
    }
    parent
}

/// The lattice label just inside piece `piece`.
fn piece_label(
    grid: &Grid,
    ijk: [usize; 3],
    depth: u8,
    labels: &[u8],
    piece: &crate::clip::FacetPiece,
) -> u8 {
    nearest_label(grid, ijk, depth, labels, probe(grid, depth, piece))
}

/// Merges the lattice components of cell `ijk` that border one surface
/// component: material bounded by one connected sheet of the surface
/// inside the cell is connected there, even where it is thinner than
/// the lattice (a wedge at a tangent join) and the lattice sees it in
/// pieces. Relabels `labels` 1..=count in lattice order; returns count.
pub(super) fn merge_by_surface(
    grid: &Grid,
    clipped: &ClippedSolid,
    cell: usize,
    ijk: [usize; 3],
    depth: u8,
    labels: &mut [u8],
) -> usize {
    let pieces = &clipped.pieces;
    let first = pieces.partition_point(|p| p.cell < cell);
    let mut roots = surface_components(grid, clipped, cell, ijk);
    let mut tag_parent: Vec<usize> = (0..=u8::MAX as usize).collect();
    let mut sheet_tag: BTreeMap<usize, usize> = BTreeMap::new();
    for i in 0..roots.len() {
        let tag = piece_label(grid, ijk, depth, labels, &pieces[first + i]) as usize;
        if tag == 0 {
            continue;
        }
        let root = find(&mut roots, i);
        match sheet_tag.get(&root) {
            Some(&other) => {
                let (a, b) = (find(&mut tag_parent, tag), find(&mut tag_parent, other));
                if a != b {
                    tag_parent[a.max(b)] = a.min(b);
                }
            }
            None => {
                sheet_tag.insert(root, tag);
            }
        }
    }
    let mut renumber = [0u8; u8::MAX as usize + 1];
    let mut count = 0usize;
    for l in labels.iter_mut() {
        if *l == 0 {
            continue;
        }
        let root = find(&mut tag_parent, *l as usize);
        if renumber[root] == 0 {
            count += 1;
            renumber[root] = u8::try_from(count).unwrap_or(u8::MAX);
        }
        *l = renumber[root];
    }
    count
}

/// The probe point just inside piece `piece` (half a lattice spacing
/// along its inward normal from its centroid).
fn probe(grid: &Grid, depth: u8, piece: &crate::clip::FacetPiece) -> [f64; 3] {
    let s = grid.cell_size() / (1usize << depth) as f64;
    let v = piece.polygon.vertices();
    let nv = v.len().max(1) as f64;
    [0, 1, 2].map(|a| v.iter().map(|x| x[a]).sum::<f64>() / nv - 0.5 * s * piece.normal[a])
}

/// A probe just inside a resolved surface finds a labelled lattice
/// point within this many lattice spacings (at most √1.5 ≈ 1.22 for
/// material at least a spacing thick); material across a slot wider
/// than this is farther.
const SUPPORT_SPACINGS: f64 = 1.5;

/// Material of cell `ijk` that no lattice point sees, beside a slot: a
/// surface component most of whose area (PH-060) has its nearest
/// labelled lattice point farther than [`SUPPORT_SPACINGS`] from points
/// just inside its pieces (pieces lying on a cell wall aside) and
/// across empty space (the straight way
/// there leaves the material), so labelling it by that point would weld
/// across the slot. Material connected to its nearest label through
/// the cell (a wedge where a plate crosses a cell corner) keeps it. Pieces whose material lies beyond the cell (a face on a
/// cell wall) do not count. Material thinner than the
/// spacing beside a slot is such a sheet; its probes' nearest labels
/// lie across the slot, so without this it would be welded there. A
/// sheet joined to resolved material (a wedge at a tangent join) has
/// most of its area supported and stays with it. In order of each
/// sheet's first piece.
pub(super) fn slivers(
    grid: &Grid,
    clipped: &ClippedSolid,
    cell: usize,
    ijk: [usize; 3],
    depth: u8,
    labels: &[u8],
) -> Vec<Sliver> {
    let h = grid.cell_size();
    let s = h / (1usize << depth) as f64;
    let lo = grid.cell_min(ijk);
    let pieces = &clipped.pieces;
    let first = pieces.partition_point(|p| p.cell < cell);
    let lattice = super::super::lattice::lattice(grid, ijk, depth);
    let cell_pieces = &pieces[first..pieces.partition_point(|p| p.cell <= cell)];
    let mut roots = surface_components(grid, clipped, cell, ijk);
    let on_wall = |v: &[[f64; 3]]| {
        (0..3).any(|a| {
            [lo[a], lo[a] + h]
                .iter()
                .any(|&w| v.iter().all(|x| (x[a] - w).abs() <= 1.0e-9 * h))
        })
    };
    // Faces of one piece of material thinner than the support reach that
    // face each other across it are one sheet (a fin crossing the cell,
    // its faces meeting only on the cell walls): the way back from one
    // into the material ends on the other.
    let reach = SUPPORT_SPACINGS * s;
    for (i, piece) in cell_pieces.iter().enumerate() {
        let v = piece.polygon.vertices();
        if on_wall(v) {
            continue;
        }
        let nv = v.len().max(1) as f64;
        let c = [0, 1, 2].map(|a| v.iter().map(|x| x[a]).sum::<f64>() / nv);
        let back = piece.normal.map(|x| -x);
        let (_, hit) = super::ray::first_hit(cell_pieces, Some(i), c, back, (1.0e-9 * h, reach));
        if let Some(j) = hit.filter(|&j| !on_wall(cell_pieces[j].polygon.vertices())) {
            let (ri, rj) = (find(&mut roots, i), find(&mut roots, j));
            roots[ri.max(rj)] = ri.min(rj);
        }
    }
    let mut wall_pieces: Vec<usize> = Vec::new();
    // Per sheet: (supported area, unsupported area, pieces), by root.
    let mut sheets: BTreeMap<usize, (f64, f64, Vec<usize>)> = BTreeMap::new();
    for i in 0..roots.len() {
        let piece = &pieces[first + i];
        let v = piece.polygon.vertices();
        // A piece lying on a cell wall has its sheet edges on the wall,
        // where sheets do not join: it does not vote, and joins a sliver
        // it shares an edge with afterwards.
        if on_wall(v) {
            wall_pieces.push(first + i);
            continue;
        }
        let nv = v.len().max(1) as f64;
        let p = [0, 1, 2]
            .map(|a| v.iter().map(|x| x[a]).sum::<f64>() / nv - 1.0e-6 * h * piece.normal[a]);
        if (0..3).any(|a| p[a] < lo[a] || p[a] > lo[a] + h) {
            continue;
        }
        // The nearest labelled lattice point, and whether the straight
        // way to it leaves the material (crosses empty space).
        let (d2, q) = lattice
            .iter()
            .zip(labels)
            .filter(|(_, l)| **l > 0)
            .map(|((_, x), _)| {
                (
                    (0..3).map(|a| (x[a] - p[a]) * (x[a] - p[a])).sum::<f64>(),
                    *x,
                )
            })
            .fold(
                (f64::INFINITY, p),
                |best, x| if x.0 < best.0 { x } else { best },
            );
        let near = d2.sqrt();
        // Exactly, against the cell's surface: a gap thinner than the
        // spacing counts.
        let across = near.is_finite() && near > 0.0 && {
            let dir = [0, 1, 2].map(|a| (q[a] - p[a]) / near);
            let tol = (1.0e-9 * h, near);
            super::ray::first_hit(cell_pieces, Some(i), p, dir, tol).0 < near
        };
        let root = find(&mut roots, i);
        let sheet = sheets.entry(root).or_insert((0.0, 0.0, Vec::new()));
        if !across {
            sheet.0 += piece.area;
        } else {
            sheet.1 += piece.area;
        }
        sheet.2.push(first + i);
    }
    // The box of a piece's material: its vertices and those stepped back
    // behind the surface, up to the support reach, the cell's boundary or
    // the material's own opposite surface.
    let material_box = |k: usize| {
        let (mut blo, mut bhi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        let n = pieces[k].normal;
        for v in pieces[k].polygon.vertices() {
            let exit = (0..3)
                .filter(|&a| n[a] != 0.0)
                .map(|a| {
                    if n[a] < 0.0 {
                        (lo[a] + h - v[a]) / -n[a]
                    } else {
                        (v[a] - lo[a]) / n[a]
                    }
                })
                .fold(reach, f64::min)
                .max(0.0);
            // Not past the material's opposite surface into a second gap.
            let back = n.map(|c| -c);
            let exit =
                super::ray::first_hit(cell_pieces, Some(k - first), *v, back, (1.0e-9 * h, exit)).0;
            let back = [0, 1, 2].map(|a| v[a] - exit * n[a]);
            for p in [*v, back] {
                for a in 0..3 {
                    blo[a] = blo[a].min(p[a]);
                    bhi[a] = bhi[a].max(p[a]);
                }
            }
        }
        (blo, bhi)
    };
    let bounds = |boxes: &[([f64; 3], [f64; 3])]| {
        boxes.iter().fold(
            ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
            |(l, u), (bl, bu)| {
                (
                    std::array::from_fn(|a| l[a].min(bl[a])),
                    std::array::from_fn(|a| u[a].max(bu[a])),
                )
            },
        )
    };
    let out: Vec<Sliver> = sheets
        .into_values()
        // Most of its area unsupported, and at least a lattice face of
        // it: crumbs of tessellation near small features are noise.
        .filter(|(supported, unsupported, _)| unsupported > supported && *unsupported >= s * s)
        .map(|(_, _, members)| {
            let boxes: Vec<_> = members.iter().map(|&k| material_box(k)).collect();
            let (blo, bhi) = bounds(&boxes);
            Sliver {
                pieces: members,
                lo: blo,
                hi: bhi,
                boxes,
            }
        })
        // A sheet with no material behind it inside the cell (its
        // material box has no thickness) is no sliver.
        .filter(|sl: &Sliver| (0..3).all(|a| sl.hi[a] - sl.lo[a] > 1.0e-6 * h))
        .collect();
    // Pieces lying on a cell wall join the sliver they share an edge
    // with, through other wall pieces too: repeated until none joins (a
    // sliver's bottom on a cell wall is tessellated into pieces only
    // some of which touch its side pieces).
    let tol = 1.0e-7 * h;
    let same = |p: [f64; 3], q: [f64; 3]| (0..3).all(|a| (p[a] - q[a]).abs() <= tol);
    let shares_edge = |x: usize, y: usize| {
        let (a, b) = (pieces[x].polygon.vertices(), pieces[y].polygon.vertices());
        (0..a.len()).any(|i| {
            let (p, q) = (a[i], a[(i + 1) % a.len()]);
            (0..b.len()).any(|j| {
                let (r, t) = (b[j], b[(j + 1) % b.len()]);
                (same(p, r) && same(q, t)) || (same(p, t) && same(q, r))
            })
        })
    };
    let mut merged = out;
    let mut pending = wall_pieces;
    loop {
        let before = pending.len();
        pending.retain(|&w| {
            let Some(sl) = merged
                .iter_mut()
                .find(|sl| sl.pieces.iter().any(|&k| shares_edge(k, w)))
            else {
                return true;
            };
            sl.pieces.push(w);
            sl.pieces.sort_unstable();
            sl.boxes.push(material_box(w));
            (sl.lo, sl.hi) = bounds(&sl.boxes);
            false
        });
        if pending.len() == before {
            break;
        }
    }
    // Material whose boxes hold a labelled lattice point is (as far as
    // this test can tell) that resolved material: a thin curved wall
    // whose straight way to its label cuts across the curve
    // (door_hinge's knuckles). Per piece, so a strip wrapping round
    // other material does not hold that material's points (PH-062).
    let mut out: Vec<Sliver> = merged
        .into_iter()
        .filter(|sl: &Sliver| {
            !lattice.iter().zip(labels).any(|((_, x), l)| {
                *l > 0
                    && sl
                        .boxes
                        .iter()
                        .any(|(lo, hi)| (0..3).all(|a| x[a] > lo[a] && x[a] < hi[a]))
            })
        })
        .collect();
    out.sort_by_key(|s| s.pieces.first().copied());
    out
}

/// Sets `piece_fragment` for the pieces of every cell in `split`
/// (except `skip` ones, whose pieces are assigned already); a sliver's
/// pieces take its fragment, after the lattice's.
pub(super) fn assign(
    grid: &Grid,
    clipped: &ClippedSolid,
    split: &[SplitCell],
    depth: u8,
    skip: &dyn Fn(usize) -> bool,
    piece_fragment: &mut [u8],
) {
    let pieces = &clipped.pieces;
    for sc in split.iter().filter(|sc| !skip(sc.cell)) {
        let Some(ijk) = grid.cell_ijk(sc.cell) else {
            continue;
        };
        let first = pieces.partition_point(|p| p.cell < sc.cell);
        let end = pieces.partition_point(|p| p.cell <= sc.cell);
        let n = end - first;
        let lattice_count = sc.labels.iter().copied().max().unwrap_or(0);
        let sliver_of = |k: usize| sc.slivers.iter().position(|s| s.pieces.contains(&k));
        let mut parent = surface_components(grid, clipped, sc.cell, ijk);
        // Each piece's own label, then the component's weighted vote.
        let mut votes: BTreeMap<usize, [f64; NODE_COPIES + 1]> = BTreeMap::new();
        let mut own = vec![0u8; n];
        for (i, piece) in pieces[first..end].iter().enumerate() {
            let tag = piece_label(grid, ijk, depth, &sc.labels, piece);
            own[i] = tag;
            let root = find(&mut parent, i);
            if let Some(slot) = votes
                .entry(root)
                .or_insert([0.0; NODE_COPIES + 1])
                .get_mut(tag as usize)
            {
                *slot += piece.area;
            }
        }
        for (i, own) in own.iter().enumerate() {
            if let Some(k) = sliver_of(first + i) {
                if let Some(slot) = piece_fragment.get_mut(first + i) {
                    *slot = lattice_count.saturating_add(k as u8);
                }
                continue;
            }
            let root = find(&mut parent, i);
            let tally = votes.get(&root).copied().unwrap_or([0.0; NODE_COPIES + 1]);
            // Unlabelled (0) never wins against a real fragment.
            let best = (1..tally.len())
                .filter(|&t| tally[t] > 0.0)
                .max_by(|&a, &b| tally[a].total_cmp(&tally[b]).then(b.cmp(&a)))
                .unwrap_or(*own as usize);
            if let Some(slot) = piece_fragment.get_mut(first + i) {
                *slot = (best as u8).saturating_sub(1);
            }
        }
    }
}
