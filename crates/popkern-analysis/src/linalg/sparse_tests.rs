use super::sparse::{Csr, PARALLEL_MIN_ENTRIES};
use super::spd_test_helpers::Lcg;

/// A random matrix with `rows` rows and a varying number of entries per
/// row, so blocks of rows carry uneven work.
fn random(rows: usize, cols: usize, rng: &mut Lcg) -> Csr {
    let mut triplets = Vec::new();
    for r in 0..rows {
        for k in 0..(r % 13) {
            let c = (r * 31 + k * 977) % cols;
            triplets.push((r, c, rng.next_signed()));
        }
    }
    Csr::from_triplets(rows, cols, &triplets).unwrap_or_else(|e| panic!("{e}"))
}

/// The threaded product equals the sequential one bit for bit, at
/// every worker count, including a last block shorter than the rest
/// and more workers than rows.
#[test]
fn threaded_products_match_the_sequential_one_bit_for_bit() {
    let mut rng = Lcg(7);
    for (rows, cols) in [(1001, 997), (5, 5)] {
        let a = random(rows, cols, &mut rng);
        let x: Vec<f64> = (0..cols).map(|_| rng.next_signed()).collect();
        let mut sequential = vec![0.0; rows];
        a.mul_vec_on(&x, &mut sequential, 1)
            .unwrap_or_else(|e| panic!("{e}"));
        for workers in [2, 3, 7, 8, 64] {
            let mut y = vec![f64::NAN; rows];
            a.mul_vec_on(&x, &mut y, workers)
                .unwrap_or_else(|e| panic!("{e}"));
            let same = y
                .iter()
                .zip(&sequential)
                .all(|(p, q)| p.to_bits() == q.to_bits());
            assert!(same, "{rows} rows on {workers} workers");
        }
    }
}

/// A product large enough to use threads gives the sequential bits,
/// and a shape mismatch is still refused.
#[test]
fn a_large_product_is_threaded_and_identical() {
    let mut rng = Lcg(11);
    let rows = PARALLEL_MIN_ENTRIES / 5;
    let a = random(rows, rows, &mut rng);
    let x: Vec<f64> = (0..rows).map(|_| rng.next_signed()).collect();
    let (mut y, mut z) = (vec![0.0; rows], vec![0.0; rows]);
    a.mul_vec(&x, &mut y).unwrap_or_else(|e| panic!("{e}"));
    a.mul_vec_on(&x, &mut z, 1)
        .unwrap_or_else(|e| panic!("{e}"));
    assert!(y.iter().zip(&z).all(|(p, q)| p.to_bits() == q.to_bits()));
    assert!(a.mul_vec(&x[1..], &mut y).is_err());
}

/// A 3×3 block added at once equals its nine entries added one by one,
/// and a block outside the stored pattern is refused.
#[test]
fn a_block_add_equals_its_entries_added_singly() {
    let mut rng = Lcg(3);
    let elements = [[0, 1, 2, 3, 4, 5, 6, 7], [4, 5, 6, 7, 8, 9, 10, 11]];
    let pattern = Csr::block_pattern(12, 3, &elements).unwrap_or_else(|e| panic!("{e}"));
    let (mut by_block, mut by_entry) = (pattern.clone(), pattern.clone());
    for element in &elements {
        for &a in element {
            for &b in element {
                let block: [[f64; 3]; 3] =
                    std::array::from_fn(|_| std::array::from_fn(|_| rng.next_signed()));
                by_block
                    .add_block3(3 * a, 3 * b, &block)
                    .unwrap_or_else(|e| panic!("{e}"));
                for (i, row) in block.iter().enumerate() {
                    for (j, v) in row.iter().enumerate() {
                        by_entry
                            .add_at(3 * a + i, 3 * b + j, *v)
                            .unwrap_or_else(|e| panic!("{e}"));
                    }
                }
            }
        }
    }
    assert_eq!(by_block, by_entry);
    // Nodes 0 and 8 share no element.
    assert!(by_block.add_block3(0, 24, &[[1.0; 3]; 3]).is_err());
    assert!(by_block.add_block3(0, 0, &[[f64::NAN; 3]; 3]).is_err());
}

/// The threaded sparse product equals the sequential one bit for bit,
/// at every worker count, and the automatic split is the same matrix.
#[test]
fn threaded_sparse_products_match_the_sequential_one_bit_for_bit() {
    let mut rng = Lcg(5);
    let a = random(1001, 997, &mut rng);
    let b = random(997, 1003, &mut rng);
    let sequential = a.mul_on(&b, 1).unwrap_or_else(|e| panic!("{e}"));
    for workers in [2, 3, 7, 8, 2000] {
        let threaded = a.mul_on(&b, workers).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(threaded, sequential, "{workers} workers");
    }
    assert_eq!(a.mul(&b), Ok(sequential));
    assert!(a.mul(&a).is_err());
}

/// A rescaled copy has the same graph; one with an entry zeroed or a
/// different shape does not.
#[test]
fn the_graph_counts_zeros() {
    let mut rng = Lcg(9);
    let a = random(40, 40, &mut rng);
    let mut scaled = a.clone();
    for v in scaled.values_mut() {
        *v *= -3.0;
    }
    assert!(a.same_graph(&scaled));
    let mut zeroed = a.clone();
    zeroed.values_mut()[7] = 0.0;
    assert!(!a.same_graph(&zeroed));
    assert!(!a.same_graph(&random(41, 40, &mut rng)));
    let r = a.position(5, a.row(5).nth(1).map_or(0, |e| e.0));
    assert!(r.is_some());
}
