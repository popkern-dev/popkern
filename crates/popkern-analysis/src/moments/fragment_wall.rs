//! One wall's section of a split cell, split by fragment.
//!
//! When every wall of a split cell is reached by several fragments, the
//! fragments' sections on a wall come from Green's theorem in the wall's
//! plane, the 2-D analogue of the closure in `fragment_exact`. Along
//! in-plane axis u (w the other in-plane axis), closed on one edge line
//! of the wall, each fragment's section moments are
//!
//! P_pq = (h/2) [ Σ_trace ∫ H_p(ξ_u) ξ_w^q ν_u ds + c_p · L_q ],
//!
//! where the trace is the fragment's piece edges on the wall (ν the
//! outward in-plane normal of its section: the piece normal projected
//! into the wall), H_p closes on one edge line and L is the fragment's
//! 1-D moments on the opposite edge line, with c = (2, 0, 2/3) again.
//! The whole wall's section (from the cell's exact total) gives the
//! whole line's L, so when one fragment reaches that line (or none
//! does) every fragment's section is exact. When several do, each
//! fragment's share of the line is its intervals, bounded where its
//! trace meets the line (PH-054), checked against the whole line.

use crate::clip::ClippedSolid;
use crate::grid::Grid;

use super::super::local;

/// c_p: the weight of H_p at the remaining edge line.
const LINE_WEIGHT: [f64; 3] = [2.0, 0.0, 2.0 / 3.0];

/// Three-point Gauss-Legendre on [0, 1] (exact to degree 5; the trace
/// integrands are at most degree 3 + 2).
const GAUSS: [(f64, f64); 3] = [
    (0.5 - 0.5 * 0.774_596_669_241_483_4, 5.0 / 18.0),
    (0.5, 8.0 / 18.0),
    (0.5 + 0.5 * 0.774_596_669_241_483_4, 5.0 / 18.0),
];

/// A section or line measure this small (relative to h² or h) is empty.
pub(super) const EMPTY: f64 = 1.0e-9;

/// One piece edge on the wall.
struct TraceEdge {
    fragment: usize,
    a: [f64; 3],
    b: [f64; 3],
    normal: [f64; 3],
}

/// The wall of a split cell to split, and what is known of it.
pub(super) struct WallInput<'a> {
    /// Linear cell index.
    pub cell: usize,
    /// Its grid coordinates.
    pub ijk: [usize; 3],
    /// Wall normal axis.
    pub axis: usize,
    /// 0 for the upper wall, 1 for the lower.
    pub side: usize,
    /// The whole wall's section moments, index b + 3c with b, c the
    /// exponents along the next two axes in cyclic order.
    pub total: &'a [f64; 9],
    /// Lattice labels of the cell.
    pub labels: &'a [u8],
    /// Slivers of the cell (no lattice points).
    pub slivers: &'a [super::Sliver],
    /// Lattice depth.
    pub depth: u8,
    /// Fragment count.
    pub count: usize,
}

/// An edge line of a wall: along axis `w`, at `position` on axis `u`.
struct Line {
    u: usize,
    w: usize,
    position: f64,
    lo: [f64; 3],
    h: f64,
    tol: f64,
}

/// Each fragment's 1-D moments ∫ ξ_w^q ds on an edge line several
/// fragments reach, from its intervals: a trace edge ending on the line
/// is where its fragment's material starts or ends along +w (ends when
/// the section's outward normal points along +w), and a fragment whose
/// first event ends was inside from the line's start. `None` when the
/// events do not alternate or the intervals miss the whole line's
/// moments `whole`.
fn line_intervals(
    trace: &[TraceEdge],
    line: &Line,
    count: usize,
    whole: &[f64; 3],
) -> Option<Vec<[f64; 3]>> {
    let Line {
        u,
        w,
        position,
        lo,
        h,
        tol,
    } = *line;
    let on = |x: f64| (x - position).abs() <= tol;
    // (ξ_w, fragment, ends here).
    let mut events: Vec<(f64, usize, bool)> = Vec::new();
    for e in trace {
        let p = match (on(e.a[u]), on(e.b[u])) {
            (true, false) => e.a,
            (false, true) => e.b,
            _ => continue,
        };
        let inplane = (e.normal[u] * e.normal[u] + e.normal[w] * e.normal[w]).sqrt();
        let nw = if inplane > 0.0 {
            e.normal[w] / inplane
        } else {
            0.0
        };
        if nw.abs() < 1.0e-9 {
            return None;
        }
        let t = (2.0 * (p[w] - lo[w]) / h - 1.0).clamp(-1.0, 1.0);
        events.push((t, e.fragment, nw > 0.0));
    }
    events.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut out = vec![[0.0; 3]; count];
    for (f, slot) in out.iter_mut().enumerate() {
        let mine: Vec<(f64, bool)> = events
            .iter()
            .filter(|e| e.1 == f)
            .map(|e| (e.0, e.2))
            .collect();
        let Some(&(_, first_ends)) = mine.first() else {
            continue;
        };
        let mut start = first_ends.then_some(-1.0);
        let mut add = |a: f64, b: f64| {
            for (q, m) in slot.iter_mut().enumerate() {
                let k = (q + 1) as f64;
                *m += 0.5 * h * (ipow(b, q + 1) - ipow(a, q + 1)) / k;
            }
        };
        for &(t, ends) in &mine {
            match (start, ends) {
                (Some(a), true) => {
                    add(a, t);
                    start = None;
                }
                (None, false) => start = Some(t),
                _ => return None,
            }
        }
        if let Some(a) = start {
            add(a, 1.0);
        }
    }
    for q in 0..3 {
        let sum: f64 = out.iter().map(|m| m[q]).sum();
        if (sum - whole[q]).abs() > 1.0e-8 * h {
            return None;
        }
    }
    Some(out)
}

/// Index b + 3c of exponents `pu` along `u` and `qw` along the other
/// in-plane axis of the wall normal to `axis`.
fn section_index(axis: usize, u: usize, pu: usize, qw: usize) -> usize {
    if u == (axis + 1) % 3 {
        pu + 3 * qw
    } else {
        qw + 3 * pu
    }
}

fn ipow(x: f64, n: usize) -> f64 {
    let mut r = 1.0;
    for _ in 0..n {
        r *= x;
    }
    r
}

/// Each fragment's section moments on the wall (same indexing as
/// `total`), closed on the first edge line (u = next axis then the one
/// after; upper before lower) that at most one fragment reaches, else
/// on the first whose fragments' intervals account for it; `None` when
/// no edge line qualifies.
pub(super) fn split_section(
    grid: &Grid,
    clipped: &ClippedSolid,
    piece_fragment: &[u8],
    s: &WallInput<'_>,
) -> Option<Vec<[f64; 9]>> {
    let h = grid.cell_size();
    let tol = 1.0e-9 * h;
    let lo = grid.cell_min(s.ijk);
    let (axis, side) = (s.axis, s.side);
    let plane = lo[axis] + if side == 0 { h } else { 0.0 };
    let pieces = &clipped.pieces;
    let first = pieces.partition_point(|p| p.cell < s.cell);
    let end = pieces.partition_point(|p| p.cell <= s.cell);
    let on = |x: f64, at: f64| (x - at).abs() <= tol;
    // The trace: piece edges on the wall, except along its boundary.
    let mut trace: Vec<TraceEdge> = Vec::new();
    for (i, piece) in pieces[first..end].iter().enumerate() {
        let v = piece.polygon.vertices();
        if v.iter().all(|x| on(x[axis], plane)) {
            continue;
        }
        for k in 0..v.len() {
            let (a, b) = (v[k], v[(k + 1) % v.len()]);
            if !on(a[axis], plane) || !on(b[axis], plane) {
                continue;
            }
            let along_boundary = (0..3).filter(|&u| u != axis).any(|u| {
                [lo[u], lo[u] + h]
                    .iter()
                    .any(|&f| on(a[u], f) && on(b[u], f))
            });
            if !along_boundary {
                let fragment = piece_fragment.get(first + i).copied().unwrap_or(0) as usize;
                trace.push(TraceEdge {
                    fragment,
                    a,
                    b,
                    normal: piece.normal,
                });
            }
        }
    }
    let n = 1usize << s.depth;
    let layer_of = |side: usize| if side == 0 { n - 1 } else { 0 };
    // Edge lines at most one fragment reaches first, then lines split
    // into each fragment's intervals.
    for split_line in [false, true] {
        for u in [(axis + 1) % 3, (axis + 2) % 3] {
            let w = 3 - axis - u;
            for line_side in 0..2 {
                let line = lo[u] + if line_side == 0 { h } else { 0.0 };
                // Trace terms per fragment, [p][q].
                let mut terms = vec![[[0.0; 3]; 3]; s.count];
                for e in &trace {
                    let t = terms.get_mut(e.fragment)?;
                    let inplane = (e.normal[u] * e.normal[u] + e.normal[w] * e.normal[w]).sqrt();
                    if inplane == 0.0 {
                        continue;
                    }
                    let nu = e.normal[u] / inplane;
                    let d = [0, 1, 2].map(|k| e.b[k] - e.a[k]);
                    let length = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                    for (g, weight) in GAUSS {
                        let x = local([0, 1, 2].map(|k| e.a[k] + g * d[k]), lo, h);
                        for (p, row) in t.iter_mut().enumerate() {
                            let k = (p + 1) as f64;
                            let closing = if line_side == 0 {
                                // Closed on the lower line: H_p(−1) = 0.
                                if p % 2 == 0 { -1.0 } else { 1.0 }
                            } else {
                                1.0
                            };
                            let hp = (ipow(x[u], p + 1) - closing) / k;
                            for (q, slot) in row.iter_mut().enumerate() {
                                *slot += nu * length * weight * hp * ipow(x[w], q);
                            }
                        }
                    }
                }
                // The whole line's moments from the whole wall's section.
                let mut whole = [0.0; 3];
                for (q, l) in whole.iter_mut().enumerate() {
                    let trace_sum: f64 = terms.iter().map(|t| t[0][q]).sum();
                    let wall = s.total[section_index(axis, u, 0, q)];
                    *l = (2.0 * wall / h - trace_sum) / LINE_WEIGHT[0];
                }
                // Fragments reaching the line: trace edges ending on it, or
                // for a line no edge reaches, the lattice row along it.
                let mut owners = 0u32;
                for e in &trace {
                    if on(e.a[u], line) != on(e.b[u], line) {
                        owners |= 1 << e.fragment;
                    }
                }
                if owners == 0 {
                    // Slivers have no lattice points: one reaches the line
                    // when one of its material boxes touches the wall and
                    // the line.
                    let lattice_count = s.labels.iter().copied().max().unwrap_or(0) as usize;
                    let near = |x: f64, p: f64| (x - p).abs() <= tol;
                    for (k, sl) in s.slivers.iter().enumerate() {
                        if sl
                            .boxes_on(axis, plane, tol)
                            .any(|(lo, hi)| near(lo[u], line) || near(hi[u], line))
                        {
                            owners |= 1 << (lattice_count + k);
                        }
                    }
                    let (la, lu) = (layer_of(side), layer_of(line_side));
                    for (idx, &l) in s.labels.iter().enumerate() {
                        let at = [idx % n, (idx / n) % n, idx / (n * n)];
                        if l > 0 && at[axis] == la && at[u] == lu {
                            owners |= 1 << (l - 1);
                        }
                    }
                }
                let lines: Vec<[f64; 3]> = match (split_line, owners.count_ones()) {
                    (false, 1) => (0..s.count)
                        .map(|f| if 1 << f == owners { whole } else { [0.0; 3] })
                        .collect(),
                    (false, 0) if whole[0].abs() <= EMPTY * h => vec![[0.0; 3]; s.count],
                    (true, c) if c >= 2 => {
                        let at = Line {
                            u,
                            w,
                            position: line,
                            lo,
                            h,
                            tol,
                        };
                        match line_intervals(&trace, &at, s.count, &whole) {
                            Some(lines) => lines,
                            None => continue,
                        }
                    }
                    _ => continue,
                };
                let parts: Vec<[f64; 9]> = (0..s.count)
                    .map(|f| {
                        let mut m = [0.0; 9];
                        for p in 0..3 {
                            for q in 0..3 {
                                let l = LINE_WEIGHT[p] * lines[f][q];
                                m[section_index(axis, u, p, q)] = 0.5 * h * (terms[f][p][q] + l);
                            }
                        }
                        m
                    })
                    .collect();
                if parts.iter().all(|m| m[0] >= -EMPTY * h * h) {
                    return Some(parts);
                }
            }
        }
    }
    None
}
