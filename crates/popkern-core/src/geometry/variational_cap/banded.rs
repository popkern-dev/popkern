//! The stacked least-squares solve `min λ ‖energy rows‖² + ‖A y − b‖²`
//! in the null space's spoke structure (PK-620).
//!
//! Every energy and rim row samples one sector, and a sector's
//! coefficients touch only its interior columns `I_s`, the local modes
//! `M_{s−1}` and `M_s` of the spokes on either side, and the global
//! block `G` (see `nullspace_basis`): the sectors form a ring joined by
//! the spokes' modes. The stacked matrix's QR factorisation runs on
//! that ring:
//!
//! 1. per sector (in parallel), QR on `I_s` over the sector's rows; the
//!    rows left touch `[M_{s−1}, M_s, G]` and are compressed by QR to
//!    at most that many: a segment of one sector;
//! 2. neighbouring segments merge, level by level (each level in
//!    parallel): QR eliminates the modes they share, and the rows left,
//!    on the merged segment's two boundaries and `G`, are compressed;
//!    the last merge closes the ring and eliminates both boundaries;
//! 3. the last rows' QR solves `G`, then back substitution solves each
//!    merge's shared modes (from the last merge back) and each `I_s`.
//!
//! Each step is an orthogonal transformation of a subset of the rows, so
//! the result is the stacked QR's, without forming normal equations.
//! Measured on the octagon h1.2 variable apex's 16-sector degree-9 fit
//! (2754 unknowns, 1792–2688 rim rows): the dense stacked QR took 4.5 s
//! per step; an elimination of the interiors alone, with a dense block
//! over every coupled column (834 wide), 0.3 s.
//!
//! A row the structure does not admit (touching two sectors' interiors,
//! or the modes of spokes that share no sector) makes the solve return
//! `None`, and the caller takes the dense QR.

use faer::Mat;

use super::factor::EnergyFactor;
use super::linalg::{lstsq_with_norms, par_map, qr_eliminate};

/// Relative rank tolerance, as the dense stacked QR's.
const RANK_TOL: f64 = 1e-12;

/// A column's block in the ring.
#[derive(Clone, Copy, PartialEq)]
enum Block {
    Interior(usize),
    /// Spoke p's local modes (spoke p joins sectors p and p + 1).
    Modes(usize),
    Global,
}

/// The ring's column blocks, in the `[x; y; z]` order.
struct Chain {
    n: usize,
    block: Vec<Block>,
    interior: Vec<Vec<usize>>,
    modes: Vec<Vec<usize>>,
    global: Vec<usize>,
}

impl Chain {
    fn new(energy: &EnergyFactor) -> Self {
        let layout = &energy.hs.layout;
        let (n, k) = (layout.interior.len(), energy.k);
        let mut block = vec![Block::Global; 3 * k];
        let mut interior = vec![Vec::new(); n];
        let mut modes = vec![Vec::new(); n];
        let mut global = Vec::new();
        for c in 0..3 {
            for (s, range) in layout.interior.iter().enumerate() {
                for j in range.clone() {
                    block[c * k + j] = Block::Interior(s);
                    interior[s].push(c * k + j);
                }
            }
            for (p, range) in layout.modes.iter().enumerate() {
                for j in range.clone() {
                    block[c * k + j] = Block::Modes(p);
                    modes[p].push(c * k + j);
                }
            }
            global.extend(layout.global.clone().map(|j| c * k + j));
        }
        Chain {
            n,
            block,
            interior,
            modes,
            global,
        }
    }

    /// Sector s's columns: `[I_s | M_{s−1} | M_s | G]`, with the widths
    /// of the four parts.
    fn sector_columns(&self, s: usize) -> (Vec<usize>, [usize; 4]) {
        let previous = &self.modes[(s + self.n - 1) % self.n];
        let parts = [&self.interior[s], previous, &self.modes[s], &self.global];
        let widths = parts.map(|p| p.len());
        (parts.into_iter().flatten().copied().collect(), widths)
    }

    /// The sectors a row touching `col` can belong to, as a bit mask.
    fn admits(&self, col: usize) -> u64 {
        match self.block[col] {
            Block::Interior(s) => 1 << s,
            Block::Modes(p) => (1 << p) | (1 << ((p + 1) % self.n)),
            Block::Global => u64::MAX,
        }
    }
}

/// One sector after step 1: its interior rows (over all its columns and
/// the right-hand side) and its compressed rows over `[M_{s−1}, M_s,
/// G]` and the right-hand side.
struct Eliminated {
    columns: Vec<usize>,
    widths: [usize; 4],
    interior: Mat<f64>,
    rest: Mat<f64>,
}

/// The structured solve; `None` when a row breaks the structure or a
/// factor is numerically rank deficient.
pub(super) fn solve(
    energy: &EnergyFactor,
    lambda: f64,
    rows: &Mat<f64>,
    rhs: &[f64],
) -> Option<[Vec<f64>; 3]> {
    let chain = Chain::new(energy);
    let n = chain.n;
    let k = energy.k;
    if !(2..64).contains(&n) || energy.sectors.len() != n {
        return None;
    }
    let sl = lambda.sqrt();
    // Assign each rim row to a sector, column by column over faer's
    // column-major storage.
    let all = (1u64 << n) - 1;
    let mut mask = vec![all; rows.nrows()];
    for col in 0..3 * k {
        let admits = chain.admits(col);
        if admits == u64::MAX {
            continue;
        }
        for (i, &v) in rows.col_as_slice(col).iter().enumerate() {
            if v != 0.0 {
                mask[i] &= admits;
            }
        }
    }
    let mut groups: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut shared = Vec::new();
    for (i, &m) in mask.iter().enumerate() {
        match m {
            0 => return None,
            m if m == all => shared.push(i),
            m => groups[m.trailing_zeros() as usize].push(i),
        }
    }
    // Every column's norm in the whole stacked system (the rank tests').
    let mut norms2 = vec![0.0; 3 * k];
    for sector in &energy.sectors {
        for (j, &col) in sector.cols.iter().enumerate() {
            norms2[col] += lambda
                * (0..sector.a.nrows())
                    .map(|i| sector.a.read(i, j).powi(2))
                    .sum::<f64>();
        }
    }
    for (col, slot) in norms2.iter_mut().enumerate() {
        *slot += rows.col_as_slice(col).iter().map(|v| v * v).sum::<f64>();
    }
    let norms: Vec<f64> = norms2.iter().map(|v| v.sqrt()).collect();
    // Step 1.
    let eliminate = |s: usize| -> Option<Eliminated> {
        let (columns, widths) = chain.sector_columns(s);
        let mut local = vec![usize::MAX; 3 * k];
        for (t, &col) in columns.iter().enumerate() {
            local[col] = t;
        }
        let sector = &energy.sectors[s];
        let width = columns.len() + 1;
        let total = sector.a.nrows() + groups[s].len();
        let mut m = Mat::<f64>::zeros(total, width);
        // The rim rows first: they carry weights up to 1e6 against the
        // energy's √λ, and Householder QR without pivoting is accurate on
        // such stiff weighted problems only with the heavy rows leading
        // (Powell and Reid). Measured with the energy rows leading: the
        // octagon h0.4 variable apex's corner-cell G1 step landed 8% over
        // the minimum objective, which one refinement step recovered.
        let rim = groups[s].len();
        for (r, &i) in groups[s].iter().enumerate() {
            for (t, &col) in columns.iter().enumerate() {
                m.write(r, t, rows.read(i, col));
            }
            m.write(r, width - 1, rhs[i]);
        }
        for i in 0..sector.a.nrows() {
            for (j, &col) in sector.cols.iter().enumerate() {
                let t = local[col];
                if t == usize::MAX {
                    return None;
                }
                m.write(rim + i, t, sl * sector.a.read(i, j));
            }
            m.write(rim + i, width - 1, -sl * sector.c[i]);
        }
        let ni = widths[0];
        let interior_norms: Vec<f64> = columns[..ni].iter().map(|&c| norms[c]).collect();
        qr_eliminate(&mut m, ni, Some((&interior_norms, RANK_TOL)))?;
        let interior = m.submatrix(0, 0, ni, width).to_owned();
        let rest = compress(m.submatrix(ni, ni, total - ni, width - ni).to_owned())?;
        Some(Eliminated {
            columns,
            widths,
            interior,
            rest,
        })
    };
    let eliminated: Vec<Eliminated> = par_map(n, eliminate).into_iter().collect::<Option<_>>()?;
    // Step 2: merge neighbouring segments, level by level.
    let ng = chain.global.len();
    let mut segments: Vec<Segment> = eliminated
        .iter()
        .map(|e| {
            let [_, wp, wn, _] = e.widths;
            Segment {
                left: e.columns[e.widths[0]..e.widths[0] + wp].to_vec(),
                right: e.columns[e.widths[0] + wp..e.widths[0] + wp + wn].to_vec(),
                rows: e.rest.clone(),
            }
        })
        .collect();
    let mut merges: Vec<Merge> = Vec::new();
    while segments.len() > 1 {
        let pairs = segments.len() / 2;
        let odd = (segments.len() % 2 == 1).then(|| segments.pop()).flatten();
        let merged = par_map(pairs, |i| {
            merge(
                &segments[2 * i],
                &segments[2 * i + 1],
                &chain.global,
                &norms,
            )
        });
        segments = Vec::with_capacity(pairs + 1);
        for result in merged {
            let (segment, record) = result?;
            segments.push(segment);
            merges.push(record);
        }
        segments.extend(odd);
    }
    let top = segments.pop()?;
    if !(top.left.is_empty() && top.right.is_empty()) {
        return None;
    }
    // Step 3: G, then back substitution.
    let carry = &top.rows;
    let total = carry.nrows() + shared.len();
    let mut a = Mat::<f64>::zeros(total, ng);
    let mut b = Mat::<f64>::zeros(total, 1);
    for i in 0..carry.nrows() {
        for g in 0..ng {
            a.write(i, g, carry.read(i, g));
        }
        b.write(i, 0, carry.read(i, ng));
    }
    for (r, &i) in shared.iter().enumerate() {
        for (g, &col) in chain.global.iter().enumerate() {
            a.write(carry.nrows() + r, g, rows.read(i, col));
        }
        b.write(carry.nrows() + r, 0, rhs[i]);
    }
    let global_norms: Vec<f64> = chain.global.iter().map(|&c| norms[c]).collect();
    let yg = lstsq_with_norms(a, b, RANK_TOL, &global_norms)?;
    let mut y = vec![0.0; 3 * k];
    for (g, &col) in chain.global.iter().enumerate() {
        y[col] = yg.read(g, 0);
    }
    // Each shared boundary, from the last merge back: a merge's known
    // columns are its segment's boundaries, solved by a later merge.
    for record in merges.iter().rev() {
        let solved = back_substitute(&record.rows, record.shared.len(), &record.known, &y)?;
        for (t, &col) in record.shared.iter().enumerate() {
            y[col] = solved[t];
        }
    }
    // Each interior.
    let interiors = par_map(n, |s| {
        let e = &eliminated[s];
        back_substitute(&e.interior, e.widths[0], &e.columns[e.widths[0]..], &y)
    });
    for (s, solved) in interiors.into_iter().enumerate() {
        let solved = solved?;
        for (t, &col) in chain.interior[s].iter().enumerate() {
            y[col] = solved[t];
        }
    }
    if y.iter().any(|v| !v.is_finite()) {
        return None;
    }
    Some(std::array::from_fn(|c| y[c * k..(c + 1) * k].to_vec()))
}

/// Consecutive sectors' rows after their interiors and inner spokes'
/// modes are eliminated: over `[left | right | G | rhs]`, `left` the
/// modes of the spoke before the first sector, `right` of the spoke
/// after the last (both empty once the ring is closed).
struct Segment {
    left: Vec<usize>,
    right: Vec<usize>,
    rows: Mat<f64>,
}

/// One merge's eliminated rows: `shared` (the modes between the two
/// segments) over `[shared | known | rhs]`.
struct Merge {
    shared: Vec<usize>,
    known: Vec<usize>,
    rows: Mat<f64>,
}

/// Merge two neighbouring segments: QR eliminates the modes between
/// them (`a.right`, which is `b.left`, and `a.left` too when the merge
/// closes the ring), and the rows left are compressed onto the merged
/// segment's boundaries and `G`.
fn merge(a: &Segment, b: &Segment, global: &[usize], norms: &[f64]) -> Option<(Segment, Merge)> {
    if a.right != b.left {
        return None;
    }
    let closes = !a.left.is_empty() && a.left == b.right;
    let shared: Vec<usize> = if closes {
        a.right.iter().chain(&a.left).copied().collect()
    } else {
        a.right.clone()
    };
    let (left, right) = if closes {
        (Vec::new(), Vec::new())
    } else {
        (a.left.clone(), b.right.clone())
    };
    let ng = global.len();
    let (ws, wl, wr) = (shared.len(), left.len(), right.len());
    let width = ws + wl + wr + ng + 1;
    // Where each part of a segment's rows lands: its left boundary, its
    // right boundary, then G and the right-hand side.
    let (a_left, a_right) = if closes { (a.right.len(), 0) } else { (ws, 0) };
    let (b_left, b_right) = if closes {
        (0, a.right.len())
    } else {
        (0, ws + wl)
    };
    let place =
        |m: &mut Mat<f64>, at: usize, seg: &Segment, (to_left, to_right): (usize, usize)| {
            let (sl, sr) = (seg.left.len(), seg.right.len());
            for i in 0..seg.rows.nrows() {
                for t in 0..sl {
                    m.write(at + i, to_left + t, seg.rows.read(i, t));
                }
                for t in 0..sr {
                    m.write(at + i, to_right + t, seg.rows.read(i, sl + t));
                }
                for g in 0..=ng {
                    m.write(at + i, ws + wl + wr + g, seg.rows.read(i, sl + sr + g));
                }
            }
        };
    let total = a.rows.nrows() + b.rows.nrows();
    let mut m = Mat::<f64>::zeros(total, width);
    place(&mut m, 0, a, (a_left, a_right));
    place(&mut m, a.rows.nrows(), b, (b_left, b_right));
    let shared_norms: Vec<f64> = shared.iter().map(|&c| norms[c]).collect();
    qr_eliminate(&mut m, ws, Some((&shared_norms, RANK_TOL)))?;
    let known: Vec<usize> = left.iter().chain(&right).chain(global).copied().collect();
    let record = Merge {
        shared,
        known,
        rows: m.submatrix(0, 0, ws, width).to_owned(),
    };
    let rows = compress(m.submatrix(ws, ws, total - ws, width - ws).to_owned())?;
    Some((Segment { left, right, rows }, record))
}

/// Solve the upper-triangular leading `w` columns of `rows` (whose last
/// column is the right-hand side) with the columns after them, named by
/// `known`, already in `y`.
fn back_substitute(rows: &Mat<f64>, w: usize, known: &[usize], y: &[f64]) -> Option<Vec<f64>> {
    let rhs_col = rows.ncols() - 1;
    if rows.nrows() < w || w + known.len() != rhs_col {
        return None;
    }
    let mut x = Mat::from_fn(w, 1, |i, _| {
        rows.read(i, rhs_col)
            - known
                .iter()
                .enumerate()
                .map(|(t, &col)| rows.read(i, w + t) * y[col])
                .sum::<f64>()
    });
    let r = popkern_math::linalg::DenseMat::from_fn(w, w, |i, j| rows.read(i, j));
    let mut xd = crate::geometry::dense::dense(&x);
    popkern_math::linalg::solve_upper_in_place(&r, &mut xd);
    x = crate::geometry::dense::mat(&xd);
    Some((0..w).map(|i| x.read(i, 0)).collect())
}

/// The rows of `m` (its last column the right-hand side) compressed by
/// QR to at most its column count less one: the rows dropped hold only
/// residual.
fn compress(mut m: Mat<f64>) -> Option<Mat<f64>> {
    let left = m.ncols() - 1;
    if m.nrows() <= left {
        return Some(m);
    }
    qr_eliminate(&mut m, left, None)?;
    Some(m.submatrix(0, 0, left, m.ncols()).to_owned())
}
