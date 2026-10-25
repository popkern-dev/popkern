use super::csr::{Block, BlockCsr, BlockTriplet, CsrError};
use super::spd_test_helpers::{Lcg, bits_hash, frobenius};

fn random_block(rng: &mut Lcg) -> Block {
    [0, 1, 2].map(|_| [0, 1, 2].map(|_| rng.next_signed()))
}

fn build(rows: usize, cols: usize, triplets: &[BlockTriplet]) -> BlockCsr {
    BlockCsr::from_triplets(rows, cols, triplets).unwrap_or_else(|e| panic!("{e}"))
}

fn apply(m: &BlockCsr, x: &[f64]) -> Vec<f64> {
    let mut y = vec![f64::NAN; 3 * m.block_rows()];
    m.mul_vec(x, &mut y).unwrap_or_else(|e| panic!("{e}"));
    y
}

/// The dense scalar matrix of the triplets, duplicates added in input
/// order.
fn dense(rows: usize, cols: usize, triplets: &[BlockTriplet]) -> Vec<f64> {
    let width = 3 * cols;
    let mut d = vec![0.0; 3 * rows * width];
    for &(r, c, ref block) in triplets {
        for a in 0..3 {
            for b in 0..3 {
                d[(3 * r + a) * width + 3 * c + b] += block[a][b];
            }
        }
    }
    d
}

#[test]
fn the_product_matches_a_dense_reference() {
    let mut rng = Lcg(31);
    for case in 0..20u64 {
        let rows = 1 + (case as usize * 7) % 13;
        let cols = 1 + (case as usize * 5) % 11;
        let count = (case as usize * 3) % 40;
        let mut triplets: Vec<BlockTriplet> = (0..count)
            .map(|_| {
                let r = ((rng.next_signed() + 1.0) * 0.5 * rows as f64) as usize % rows;
                let c = ((rng.next_signed() + 1.0) * 0.5 * cols as f64) as usize % cols;
                (r, c, random_block(&mut rng))
            })
            .collect();
        // a guaranteed duplicate and, with rows > 1, a guaranteed empty
        // last row (no triplet uses it after the filter)
        if let Some(&(r, c, _)) = triplets.first() {
            triplets.push((r, c, random_block(&mut rng)));
        }
        if rows > 1 {
            triplets.retain(|t| t.0 != rows - 1);
        }
        let m = build(rows, cols, &triplets);
        let x: Vec<f64> = (0..3 * cols).map(|_| rng.next_signed()).collect();
        let y = apply(&m, &x);
        let d = dense(rows, cols, &triplets);
        let reference: Vec<f64> = (0..3 * rows)
            .map(|i| (0..3 * cols).map(|j| d[i * 3 * cols + j] * x[j]).sum())
            .collect();
        let diff: Vec<f64> = y.iter().zip(&reference).map(|(p, q)| p - q).collect();
        let scale = frobenius(&reference).max(f64::MIN_POSITIVE);
        assert!(
            frobenius(&diff) <= 1e-15 * scale,
            "case {case}: {} relative",
            frobenius(&diff) / scale
        );
        if rows > 1 {
            assert_eq!(m.row(rows - 1).count(), 0);
            assert_eq!(&y[3 * (rows - 1)..], &[0.0; 3]);
        }
    }
}

#[test]
fn duplicates_sum_in_input_order() {
    // (1e16 + 1) - 1e16 = 0 in f64, but (1e16 - 1e16) + 1 = 1: the sum
    // depends on order, so this pins input order.
    let block = |v: f64| -> Block { [[v, 0.0, 0.0], [0.0; 3], [0.0; 3]] };
    let other_row: BlockTriplet = (1, 0, block(5.0));
    let forward = [
        (0, 0, block(1e16)),
        other_row,
        (0, 0, block(1.0)),
        (0, 0, block(-1e16)),
    ];
    let m = build(2, 1, &forward);
    assert_eq!(m.stored_blocks(), 2);
    let (_, sum) = m.row(0).next().unwrap_or_else(|| panic!("row 0 empty"));
    assert_eq!(sum[0][0], 0.0);
    let reordered = [
        (0, 0, block(1e16)),
        (0, 0, block(-1e16)),
        (0, 0, block(1.0)),
    ];
    let m = build(1, 1, &reordered);
    let (_, sum) = m.row(0).next().unwrap_or_else(|| panic!("row 0 empty"));
    assert_eq!(sum[0][0], 1.0);
}

#[test]
fn rows_are_sorted_by_column() {
    let b = [[1.0; 3]; 3];
    let m = build(1, 4, &[(0, 3, b), (0, 0, b), (0, 2, b)]);
    let cols: Vec<usize> = m.row(0).map(|(c, _)| c).collect();
    assert_eq!(cols, [0, 2, 3]);
    assert_eq!(m.row(9).count(), 0);
    assert_eq!((m.block_rows(), m.block_cols()), (1, 4));
}

#[test]
fn malformed_input_is_typed() {
    let b = [[1.0; 3]; 3];
    assert_eq!(
        BlockCsr::from_triplets(2, 2, &[(0, 0, b), (2, 1, b)]),
        Err(CsrError::OutOfRange {
            index: 1,
            row: 2,
            col: 1,
            rows: 2,
            cols: 2
        })
    );
    let mut bad = b;
    bad[1][2] = f64::INFINITY;
    assert_eq!(
        BlockCsr::from_triplets(2, 2, &[(1, 1, bad)]),
        Err(CsrError::NonFinite(0))
    );
    let m = build(2, 3, &[(0, 0, b)]);
    let mut y = vec![0.0; 6];
    assert_eq!(
        m.mul_vec(&[0.0; 6], &mut y),
        Err(CsrError::VectorLength {
            len: 6,
            expected: 9
        })
    );
    assert_eq!(
        m.mul_vec(&[0.0; 9], &mut [0.0; 3]),
        Err(CsrError::VectorLength {
            len: 3,
            expected: 6
        })
    );
}

const PIN_ROWS: usize = 1000;

/// A banded block matrix like a coarse 1-D stencil: block row r couples
/// to columns r−3..=r+3.
fn pin_problem() -> (BlockCsr, Vec<f64>) {
    let mut rng = Lcg(0x5eed);
    let mut triplets = Vec::new();
    for r in 0..PIN_ROWS {
        for c in r.saturating_sub(3)..(r + 4).min(PIN_ROWS) {
            triplets.push((r, c, random_block(&mut rng)));
        }
    }
    let x = (0..3 * PIN_ROWS).map(|_| rng.next_signed()).collect();
    (build(PIN_ROWS, PIN_ROWS, &triplets), x)
}

/// The same product accumulated in reverse column order.
fn reversed_product(m: &BlockCsr, x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; 3 * m.block_rows()];
    for r in 0..m.block_rows() {
        let row: Vec<(usize, &Block)> = m.row(r).collect();
        for &(c, block) in row.iter().rev() {
            for a in 0..3 {
                for b in 0..3 {
                    y[3 * r + a] += block[a][b] * x[3 * c + b];
                }
            }
        }
    }
    y
}

/// FNV-1a of the product's bits on the pin problem. A change means the
/// SpMV's arithmetic order changed.
const PINNED_SPMV_HASH: u64 = 6885687272492572238;

#[test]
fn the_spmv_result_bits_are_pinned() {
    let (m, x) = pin_problem();
    let y = apply(&m, &x);
    assert_eq!(bits_hash(&apply(&m, &x)), bits_hash(&y));
    assert_eq!(bits_hash(&y), PINNED_SPMV_HASH);
}

#[test]
fn the_pin_sees_a_change_of_summation_order() {
    let (m, x) = pin_problem();
    let forward = apply(&m, &x);
    let reversed = reversed_product(&m, &x);
    // Same value to round-off, different bits: the pin is sensitive.
    let diff: Vec<f64> = forward.iter().zip(&reversed).map(|(p, q)| p - q).collect();
    assert!(frobenius(&diff) <= 1e-14 * frobenius(&forward));
    assert_ne!(bits_hash(&forward), bits_hash(&reversed));
}

#[test]
fn block_pattern_assembly_matches_triplets_bit_for_bit() {
    // Three overlapping "cells" on 12 nodes, 3 dofs each, with values
    // that do not sum exactly; the pattern path must sum each entry in
    // the same (cell) order as the triplet sort.
    let elements: [[usize; 8]; 3] = [
        [0, 1, 2, 3, 4, 5, 6, 7],
        [2, 3, 8, 9, 6, 7, 10, 11],
        [4, 5, 6, 7, 8, 9, 10, 11],
    ];
    let value = |e: usize, r: usize, c: usize| 0.1 * (e + 1) as f64 + 1.0 / (1 + r + 3 * c) as f64;
    let mut triplets = Vec::new();
    use super::sparse::Csr;
    let mut a = Csr::block_pattern(12, 3, &elements).unwrap_or_else(|e| panic!("{e}"));
    for (e, nodes) in elements.iter().enumerate() {
        for r in 0..24 {
            for c in 0..24 {
                let (gr, gc) = (3 * nodes[r / 3] + r % 3, 3 * nodes[c / 3] + c % 3);
                triplets.push((gr, gc, value(e, r, c)));
                a.add_at(gr, gc, value(e, r, c))
                    .unwrap_or_else(|e| panic!("{e}"));
            }
        }
    }
    let b = Csr::from_triplets(36, 36, &triplets).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(a, b);
    // An entry outside the pattern is refused.
    assert!(a.add_at(0, 3 * 11, 1.0).is_err());
}

#[test]
fn block_diagonal_matches_offset_triplets() {
    use super::sparse::Csr;
    let a = Csr::from_triplets(2, 2, &[(0, 0, 1.5), (0, 1, 0.0), (1, 1, -2.0)])
        .unwrap_or_else(|e| panic!("{e}"));
    let b = Csr::from_triplets(3, 3, &[(0, 2, 4.0), (1, 0, 0.25), (2, 1, 7.0), (2, 2, 1.0)])
        .unwrap_or_else(|e| panic!("{e}"));
    let mut triplets = Vec::new();
    for (m, off) in [(&a, 0), (&b, 2)] {
        for r in 0..m.rows() {
            triplets.extend(m.row(r).map(|(c, v)| (r + off, c + off, v)));
        }
    }
    let want = Csr::from_triplets(5, 5, &triplets).unwrap_or_else(|e| panic!("{e}"));
    let got = Csr::block_diagonal(vec![a.clone(), b]).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(got, want);
    assert_eq!(Csr::block_diagonal(vec![a.clone()]), Ok(a));
    // A non-finite entry is refused at its global position.
    let big = Csr::from_triplets(1, 1, &[(0, 0, f64::MAX)]).unwrap_or_else(|e| panic!("{e}"));
    let two = Csr::from_triplets(1, 1, &[(0, 0, 2.0)]).unwrap_or_else(|e| panic!("{e}"));
    let bad = big.mul(&two).unwrap_or_else(|e| panic!("{e}"));
    let refused = Csr::block_diagonal(vec![want, bad]);
    assert!(matches!(
        refused,
        Err(super::sparse::SparseError::NonFinite { row: 5, col: 5 })
    ));
}
