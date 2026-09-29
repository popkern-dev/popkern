//! The KP20 bi-9 cap's free coefficients and the null space of the
//! homogeneous relation rows over them (PK-600): internal C² across the
//! lattice seams, sector-seam C⁰, and eq (8)/(9) collocated at bi-9.
//! Rim bands and o-blocks are held. The rows are homogeneous: a start
//! that satisfies every relation keeps satisfying them under any `D`
//! in their null space.

// Coefficient nets are addressed by (i, j) index pairs throughout.
#![allow(clippy::needless_range_loop)]

use super::kp20_bi9_seams::{Side, b9, samples, slot};

const P1: usize = 0;
const P2: usize = 1;
const P3: usize = 2;
const P4: usize = 3;

/// Internal seams of a sector: (minus piece, minus side, plus piece, plus side).
pub(super) const INTERNAL: [(usize, Side, usize, Side); 4] = [
    (P1, Side::I9, P2, Side::I0),
    (P4, Side::I9, P3, Side::I0),
    (P1, Side::J9, P4, Side::J0),
    (P2, Side::J9, P3, Side::J0),
];

/// Whether coefficient `(i, j)` of piece `k` is held: a rim band
/// (placement-owned) or the o-block (G2 at the extraordinary point).
pub(super) fn held(k: usize, i: usize, j: usize) -> bool {
    match k {
        P1 => i <= 2 && j <= 2,
        P2 => i >= 7,
        P3 => i >= 7 || j >= 7,
        P4 => j >= 7,
        _ => false,
    }
}

/// Free coefficients and the homogeneous relation rows over them.
pub(super) struct RelationSpace {
    /// Unknown index per `[sector][piece][i][j]` (`usize::MAX` when held).
    pub index: Vec<[[[usize; 10]; 10]; 4]>,
    /// Number of free coefficients.
    pub count: usize,
    /// Sparse relation rows `(unknown, weight)`, each scaled to unit norm.
    pub rows: Vec<Vec<(usize, f64)>>,
}

/// Build the free-coefficient index and the relation rows for an
/// `n`-sided cap with eq (10) coefficients `b` (near o, near m).
pub(super) fn relation_space(n: usize, b: [[f64; 3]; 2]) -> Option<RelationSpace> {
    // Unknown index per free coefficient.
    let mut index = vec![[[[usize::MAX; 10]; 10]; 4]; n];
    let mut count = 0usize;
    for s in 0..n {
        for k in 0..4 {
            for i in 0..10 {
                for j in 0..10 {
                    if !held(k, i, j) {
                        index[s][k][i][j] = count;
                        count += 1;
                    }
                }
            }
        }
    }
    // Homogeneous relation rows over the unknowns.
    let mut rows: Vec<Vec<(usize, f64)>> = Vec::new();
    let term = |row: &mut Vec<(usize, f64)>, s: usize, k: usize, (i, j): (usize, usize), w: f64| {
        let x = index[s][k][i][j];
        if x != usize::MAX && w != 0.0 {
            row.push((x, w));
        }
    };
    for s in 0..n {
        for &(m, sm, p, sp) in &INTERNAL {
            for along in 0..10 {
                for q in 0..3 {
                    let mut row = Vec::new();
                    // plus side difference of order q minus minus side's.
                    let coeffs: &[f64] = match q {
                        0 => &[1.0],
                        1 => &[-1.0, 1.0],
                        _ => &[1.0, -2.0, 1.0],
                    };
                    for (d, &c) in coeffs.iter().enumerate() {
                        term(&mut row, s, p, slot(sp, d, along), c);
                    }
                    // minus side: M0 - M1 (q=1), M0 - 2M1 + M2 (q=2): depth d weight
                    let minus: &[f64] = match q {
                        0 => &[1.0],
                        1 => &[1.0, -1.0],
                        _ => &[1.0, -2.0, 1.0],
                    };
                    for (d, &c) in minus.iter().enumerate() {
                        term(&mut row, s, m, slot(sm, d, along), -c);
                    }
                    if !row.is_empty() {
                        rows.push(row);
                    }
                }
            }
        }
        let t = (s + 1) % n;
        for (pair, (a, c)) in [(P1, P1), (P2, P4)].into_iter().enumerate() {
            let bb = b[pair];
            for along in 0..10 {
                let mut row = Vec::new();
                term(&mut row, s, a, slot(Side::J0, 0, along), 1.0);
                term(&mut row, t, c, slot(Side::I0, 0, along), -1.0);
                if !row.is_empty() {
                    rows.push(row);
                }
            }
            for u in samples(16) {
                let bu =
                    bb[0] * (1.0 - u) * (1.0 - u) + 2.0 * bb[1] * u * (1.0 - u) + bb[2] * u * u;
                let bp = 2.0 * ((bb[1] - bb[0]) * (1.0 - u) + (bb[2] - bb[1]) * u);
                let (mut e8, mut e9) = (Vec::new(), Vec::new());
                for along in 0..10 {
                    let [w, w1, w2] = b9(along, u);
                    let pr = |d: usize| slot(Side::J0, d, along);
                    let cr = |d: usize| slot(Side::I0, d, along);
                    // eq (8): 9 Fc + 9 Fp − b V'
                    term(&mut e8, t, c, cr(1), 9.0 * w);
                    term(&mut e8, t, c, cr(0), -9.0 * w);
                    term(&mut e8, s, a, pr(1), 9.0 * w);
                    term(&mut e8, s, a, pr(0), -9.0 * w - bu * w1);
                    // eq (9): 72 Sc − 72 Sp + 18 b Fp' − b² V'' − b b' V'
                    term(&mut e9, t, c, cr(2), 72.0 * w);
                    term(&mut e9, t, c, cr(1), -144.0 * w);
                    term(&mut e9, t, c, cr(0), 72.0 * w);
                    term(&mut e9, s, a, pr(2), -72.0 * w);
                    term(&mut e9, s, a, pr(1), 144.0 * w + 18.0 * bu * w1);
                    term(
                        &mut e9,
                        s,
                        a,
                        pr(0),
                        -72.0 * w - 18.0 * bu * w1 - bu * bu * w2 - bu * bp * w1,
                    );
                }
                rows.push(e8);
                rows.push(e9);
            }
        }
    }
    // Merge duplicate terms and scale every row to unit norm.
    let rows: Vec<Vec<(usize, f64)>> = rows
        .into_iter()
        .filter_map(|row| {
            let mut merged: Vec<(usize, f64)> = Vec::with_capacity(row.len());
            let mut sorted = row;
            sorted.sort_by_key(|&(x, _)| x);
            for (x, w) in sorted {
                match merged.last_mut() {
                    Some(last) if last.0 == x => last.1 += w,
                    _ => merged.push((x, w)),
                }
            }
            let norm = merged.iter().map(|(_, w)| w * w).sum::<f64>().sqrt();
            (norm > 0.0).then(|| merged.into_iter().map(|(x, w)| (x, w / norm)).collect())
        })
        .collect();
    Some(RelationSpace { index, count, rows })
}
