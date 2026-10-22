use super::*;

fn sample(rows: usize, cols: usize, seed: f64) -> DenseMat {
    use crate::PortableMath;
    DenseMat::from_fn(rows, cols, |i, j| {
        ((i * 7 + j * 3) as f64 * seed).psin() + if i == j { 2.0 } else { 0.0 }
    })
}

fn max_diff(a: &DenseMat, b: &DenseMat) -> f64 {
    a.as_slice()
        .iter()
        .zip(b.as_slice())
        .fold(0.0, |m, (x, y)| m.max((x - y).abs()))
}

#[test]
fn matmul_matches_the_definition_and_accumulates() {
    let (a, b) = (sample(4, 3, 0.3), sample(3, 5, 0.7));
    let mut c = DenseMat::from_fn(4, 5, |i, j| (i + j) as f64);
    let before = c.clone();
    matmul(&mut c, Some(2.0), &a, &b, -1.0);
    for i in 0..4 {
        for j in 0..5 {
            let s: f64 = (0..3).map(|k| a.get(i, k) * b.get(k, j)).sum();
            assert!((c.get(i, j) - (2.0 * before.get(i, j) - s)).abs() < 1e-14);
        }
    }
    assert_eq!(matmul_tn(&a, &a), {
        let mut d = DenseMat::zeros(3, 3);
        matmul(&mut d, None, &a.transpose(), &a, 1.0);
        d
    });
}

#[test]
fn syrk_upper_is_the_upper_triangle_of_the_product() {
    let q = sample(4, 6, 0.53);
    let mut full = DenseMat::zeros(4, 4);
    matmul(&mut full, None, &q, &q.transpose(), 1.0);
    let upper = syrk_upper(&q);
    for i in 0..4 {
        for j in 0..4 {
            let want = if i <= j { full.get(i, j) } else { 0.0 };
            assert!((upper.get(i, j) - want).abs() < 1e-14);
        }
    }
}

#[test]
fn lu_cholesky_and_qr_solve_what_they_factor() {
    let a = sample(6, 6, 0.41);
    let x = DenseMat::from_fn(6, 2, |i, j| (i as f64) - 2.0 * j as f64);
    let mut b = DenseMat::zeros(6, 2);
    matmul(&mut b, None, &a, &x, 1.0);
    let lu = lu_solve(a.clone(), b.clone()).expect("lu");
    assert!(max_diff(&lu, &x) < 1e-12);

    let spd = matmul_tn(&a, &a);
    let l = cholesky(&spd).expect("spd");
    let mut y = DenseMat::zeros(6, 2);
    matmul(&mut y, None, &spd, &x, 1.0);
    solve_lower_in_place(&l, &mut y);
    solve_lower_transposed_in_place(&l, &mut y);
    assert!(max_diff(&y, &x) < 1e-10);

    let tall = sample(9, 4, 0.23);
    let xt = DenseMat::from_fn(4, 1, |i, _| 1.0 + i as f64);
    let mut bt = DenseMat::zeros(9, 1);
    matmul(&mut bt, None, &tall, &xt, 1.0);
    let qr = HouseholderQr::factor_left(tall, 4);
    assert!(max_diff(&qr.solve(bt), &xt) < 1e-12);
}

#[test]
fn column_pivoted_qr_solves_least_squares_in_the_original_order() {
    let a = sample(8, 4, 0.61);
    let x = DenseMat::from_fn(4, 1, |i, _| 0.5 - i as f64);
    let mut b = DenseMat::zeros(8, 1);
    matmul(&mut b, None, &a, &x, 1.0);
    let solved = ColPivQr::new(a).solve_lstsq(b);
    assert!(max_diff(&solved, &x) < 1e-12);
}

#[test]
fn singular_and_indefinite_inputs_decline() {
    let singular = DenseMat::from_fn(3, 3, |i, _| i as f64);
    assert!(lu_solve(singular, DenseMat::zeros(3, 1)).is_none());
    let indefinite = DenseMat::from_fn(2, 2, |i, j| if i == j { -1.0 } else { 0.0 });
    assert!(cholesky(&indefinite).is_none());
}

#[test]
fn svd_reproduces_a_rank_deficient_matrix() {
    let left = sample(7, 3, 0.37);
    let right = sample(3, 5, 0.71);
    let mut a = DenseMat::zeros(7, 5);
    matmul(&mut a, None, &left, &right, 1.0);
    let d = jacobi_svd(&a).expect("svd");
    let worst = (0..7)
        .flat_map(|i| (0..5).map(move |j| (i, j)))
        .map(|(i, j)| {
            let v: f64 = (0..5)
                .map(|t| d.u.get(i, t) * d.sigma[t] * d.v.get(j, t))
                .sum();
            (v - a.get(i, j)).abs()
        })
        .fold(0.0, f64::max);
    assert!(worst < 1e-12, "{worst:e}");
    assert_eq!(d.sigma.iter().filter(|&&s| s > 1e-10).count(), 3);
}

/// Bits pinned: the kernels are part of the kernel's numeric contract
/// (one order of operations everywhere); a change here moves every
/// corner cap and must be deliberate.
#[test]
fn bits_are_pinned() {
    let a = sample(5, 5, 0.29);
    let b = DenseMat::from_fn(5, 1, |i, _| 1.0 / (1.0 + i as f64));
    let x = lu_solve(a, b).expect("lu");
    let bits: Vec<u64> = x.as_slice().iter().map(|v| v.to_bits()).collect();
    assert_eq!(
        bits,
        [
            0x3fca_918c_39a0_8adf,
            0x3fcb_2cb9_51be_3ae5,
            0x3fdf_8568_e005_985c,
            0xbfc0_25fb_8158_d1e7,
            0xbf6a_4409_9869_3bb8,
        ]
    );
}

/// A quasi-definite `[[H, Aᵀ], [A, −εI]]` with sparse blocks, as lower
/// triangle triplets in columns.
fn quasi_definite(nh: usize, na: usize, eps: f64) -> (usize, Vec<usize>, Vec<usize>, Vec<f64>) {
    use crate::PortableMath;
    let n = nh + na;
    let mut cols: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for (j, column) in cols.iter_mut().enumerate().take(nh) {
        column.push((j, 4.0 + (j as f64 * 0.3).psin()));
        if j + 1 < nh {
            column.push((j + 1, -1.0));
        }
    }
    for r in 0..na {
        for j in [r % nh, (3 * r + 1) % nh] {
            let w = 0.5 + 0.1 * (r as f64);
            if !cols[j].iter().any(|&(row, _)| row == nh + r) {
                cols[j].push((nh + r, w));
            }
        }
        cols[nh + r].push((nh + r, -eps));
    }
    let mut col_ptr = vec![0usize];
    let (mut rows, mut vals) = (Vec::new(), Vec::new());
    for column in &cols {
        for &(r, v) in column {
            rows.push(r);
            vals.push(v);
        }
        col_ptr.push(rows.len());
    }
    (n, col_ptr, rows, vals)
}

#[test]
fn sparse_ldlt_solves_a_quasi_definite_system_under_any_ordering() {
    let (n, col_ptr, rows, vals) = quasi_definite(9, 4, 1e-3);
    let dense = DenseMat::from_fn(n, n, |i, j| {
        let (r, c) = (i.max(j), i.min(j));
        (col_ptr[c]..col_ptr[c + 1])
            .find(|&p| rows[p] == r)
            .map_or(0.0, |p| vals[p])
    });
    let b: Vec<f64> = (0..n).map(|i| 1.0 + i as f64 * 0.25).collect();
    let want = lu_solve(dense, DenseMat::from_fn(n, 1, |i, _| b[i])).expect("lu");
    for perm in [
        (0..n).collect::<Vec<_>>(),
        (0..n).rev().collect(),
        (0..n).map(|k| (k * 5) % n).collect(),
    ] {
        let sym = LdltSymbolic::new(n, &col_ptr, &rows, perm).expect("symbolic");
        let factor = sym.factor(&vals).expect("factor");
        let mut x = b.clone();
        sym.solve(&factor, &mut x);
        for (i, xi) in x.iter().enumerate() {
            assert!(
                (xi - want.get(i, 0)).abs() < 1e-10,
                "x[{i}] {} vs {}",
                xi,
                want.get(i, 0)
            );
        }
    }
}

/// PK-600: dense blocks coupled by relation rows (the KP20 curvature
/// system's shape) factor as wide supernodal panels and solve like LU
/// under every ordering, including ones that interleave the blocks.
#[test]
fn sparse_ldlt_blocks_factor_as_wide_panels() {
    use crate::PortableMath;
    let (blocks, width, relations) = (3usize, 24usize, 5usize);
    let nh = blocks * width;
    let n = nh + relations;
    let mut cols: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for blk in 0..blocks {
        for j in 0..width {
            let c = blk * width + j;
            for i in j..width {
                let r = blk * width + i;
                let v = if i == j {
                    width as f64 + 1.0
                } else {
                    ((r * 3 + c) as f64 * 0.37).psin() * 0.5
                };
                cols[c].push((r, v));
            }
        }
    }
    // Each relation ties the last unknown of one block to the first of
    // the next: the seams.
    for r in 0..relations {
        let (a, b) = ((r % blocks) * width + width - 1, ((r + 1) % blocks) * width);
        let (lo, hi) = (a.min(b), a.max(b));
        cols[lo].push((nh + r, 1.0));
        cols[hi].push((nh + r, -1.0 - 0.1 * r as f64));
        // −G well away from zero: a no-pivot order that eliminates a
        // relation before its block grows elements by 1/G.
        cols[nh + r].push((nh + r, -1e-3));
    }
    let mut col_ptr = vec![0usize];
    let (mut rows, mut vals) = (Vec::new(), Vec::new());
    for column in &mut cols {
        column.sort_by_key(|&(r, _)| r);
        for &(r, v) in column.iter() {
            rows.push(r);
            vals.push(v);
        }
        col_ptr.push(rows.len());
    }
    let dense = DenseMat::from_fn(n, n, |i, j| {
        let (r, c) = (i.max(j), i.min(j));
        (col_ptr[c]..col_ptr[c + 1])
            .find(|&p| rows[p] == r)
            .map_or(0.0, |p| vals[p])
    });
    let b: Vec<f64> = (0..n).map(|i| (i as f64 * 0.21).pcos()).collect();
    let want = lu_solve(dense, DenseMat::from_fn(n, 1, |i, _| b[i])).expect("lu");
    let natural: Vec<usize> = (0..n).collect();
    let interleaved: Vec<usize> = (0..n).map(|k| (k * 5) % n).collect();
    for (name, perm) in [("natural", natural), ("interleaved", interleaved)] {
        let sym = LdltSymbolic::new(n, &col_ptr, &rows, perm).expect("symbolic");
        if name == "natural" {
            // Each dense block is one panel (plus the relation tail).
            assert!(
                sym.supernodes() <= blocks + relations,
                "{}",
                sym.supernodes()
            );
        }
        let factor = sym.factor(&vals).expect("factor");
        let mut x = b.clone();
        sym.solve(&factor, &mut x);
        for (i, xi) in x.iter().enumerate() {
            assert!(
                (xi - want.get(i, 0)).abs() < 1e-9 * (1.0 + want.get(i, 0).abs()),
                "{name} x[{i}] {} vs {}",
                xi,
                want.get(i, 0)
            );
        }
    }
    // A zero pivot declines.
    let sym = LdltSymbolic::new(1, &[0, 1], &[0], vec![0]).expect("symbolic");
    assert!(sym.factor(&[0.0]).is_none());
}

#[test]
fn sparse_ldlt_refuses_a_malformed_pattern() {
    // An upper-triangle entry where the lower triangle is expected.
    assert!(LdltSymbolic::new(2, &[0, 1, 2], &[1, 0], vec![0, 1]).is_none());
    // A permutation that repeats a row.
    assert!(LdltSymbolic::new(2, &[0, 1, 2], &[0, 1], vec![0, 0]).is_none());
}

#[test]
fn symmetric_eigenvalues_are_the_singular_values_up_to_sign() {
    for (n, seed) in [(1, 0.3), (2, 0.4), (7, 0.45), (30, 0.17)] {
        let b = sample(n, n, seed);
        // A symmetric indefinite matrix: B + Bᵀ minus a shift.
        let a = DenseMat::from_fn(n, n, |i, j| {
            b.get(i, j) + b.get(j, i) - if i == j { 3.0 } else { 0.0 }
        });
        let mut eig: Vec<f64> = symmetric_eigenvalues(&a)
            .expect("eig")
            .iter()
            .map(|v| v.abs())
            .collect();
        eig.sort_by(|x, y| y.total_cmp(x));
        let sv = jacobi_svd(&a).expect("svd").sigma;
        for (e, s) in eig.iter().zip(&sv) {
            assert!((e - s).abs() <= 1e-12 * (1.0 + s), "n={n}: {e} vs {s}");
        }
    }
}

#[test]
fn symmetric_eigenvalues_resolve_a_badly_conditioned_kkt_matrix() {
    // [[H, Cᵀ], [C, −1e-12 I]] where constraint row 2 repeats row 0:
    // only the −εI block resolves the redundant direction.
    let (nh, nc) = (6, 3);
    let n = nh + nc;
    let c = |r: usize, k: usize| -> f64 {
        let row = if r == 2 { 0 } else { r };
        (row + k).is_multiple_of(3) as u8 as f64
    };
    let a = DenseMat::from_fn(n, n, |i, j| match (i < nh, j < nh) {
        (true, true) if i == j => 2.0,
        (true, true) if i.abs_diff(j) == 1 => -1.0,
        (true, true) => 0.0,
        (false, true) => c(i - nh, j),
        (true, false) => c(j - nh, i),
        (false, false) if i == j => -1e-12,
        (false, false) => 0.0,
    });
    let eig = symmetric_eigenvalues(&a).expect("converges");
    assert_eq!(eig.len(), n);
    let smallest = eig.iter().map(|v| v.abs()).fold(f64::INFINITY, f64::min);
    assert!(smallest > 0.0 && smallest < 1e-9, "{smallest:e}");
}
