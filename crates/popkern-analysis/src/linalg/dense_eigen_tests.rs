use super::dense_eigen::{max_generalized_eigenvalue_dense, symmetric_eigen};

#[test]
fn the_largest_generalized_eigenvalue() {
    let flat = |m: &[[f64; 3]; 3]| m.iter().flatten().copied().collect::<Vec<f64>>();
    let b = [[2.0, 0.5, 0.0], [0.5, 1.0, 0.2], [0.0, 0.2, 3.0]];
    let a = [[1.0, 0.3, -0.2], [0.3, 4.0, 0.1], [-0.2, 0.1, 0.5]];
    let got = max_generalized_eigenvalue_dense(&flat(&a), &flat(&b), 3, 300).unwrap_or(f64::NAN);
    // Against the eigenvalues of b⁻¹ᐟ² a b⁻¹ᐟ² by Jacobi.
    let (bv, bq) = symmetric_eigen(b);
    let half: [[f64; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| (0..3).map(|k| bq[i][k] * bq[j][k] / bv[k].sqrt()).sum())
    });
    let m: [[f64; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            (0..3)
                .map(|p| {
                    (0..3)
                        .map(|q| half[i][p] * a[p][q] * half[q][j])
                        .sum::<f64>()
                })
                .sum()
        })
    });
    let (values, _) = symmetric_eigen(m);
    assert!(
        (got - values[2]).abs() <= 1e-10 * values[2],
        "{got} vs {}",
        values[2]
    );
    let scaled: Vec<f64> = flat(&b).iter().map(|x| 3.0 * x).collect();
    let three = max_generalized_eigenvalue_dense(&scaled, &flat(&b), 3, 300).unwrap_or(f64::NAN);
    assert!((three - 3.0).abs() <= 1e-12, "{three}");
    // Not positive definite: no answer.
    let indefinite = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]];
    assert!(max_generalized_eigenvalue_dense(&flat(&a), &flat(&indefinite), 3, 300).is_none());
}

#[test]
fn the_dense_jacobi_matches_the_fixed_size_one() {
    use super::dense_eigen::symmetric_eigen_dense;
    let a = [[4.0, 1.0, 0.5], [1.0, 3.0, -0.2], [0.5, -0.2, 1.0]];
    let (fixed, _) = symmetric_eigen(a);
    let flat: Vec<f64> = a.iter().flatten().copied().collect();
    let (values, vectors) = symmetric_eigen_dense(flat.clone(), 3);
    for (x, y) in values.iter().zip(fixed) {
        assert!((x - y).abs() <= 1e-12, "{values:?} vs {fixed:?}");
    }
    // A v = λ v for each column.
    for k in 0..3 {
        for r in 0..3 {
            let av: f64 = (0..3).map(|c| flat[r * 3 + c] * vectors[c * 3 + k]).sum();
            assert!((av - values[k] * vectors[r * 3 + k]).abs() <= 1e-12);
        }
    }
}
