//! The blocked kernels against the scalar loops they replaced (PK-1223):
//! bit for bit, not to a tolerance, since a ranking or certified bound
//! downstream may hinge on the last ulp.

use super::*;
use crate::PortableMath;

/// An awkward matrix: mixed signs and magnitudes, scattered exact zeros.
fn awkward(rows: usize, cols: usize, seed: f64) -> DenseMat {
    DenseMat::from_fn(rows, cols, |i, j| {
        if (i * 5 + j * 3).is_multiple_of(7) {
            return 0.0;
        }
        let x = ((i * 13 + j * 7) as f64 * seed).psin();
        x * 10f64.powi(((i + 2 * j) % 5) as i32 - 2)
    })
}

fn bits(m: &DenseMat) -> Vec<u64> {
    m.as_slice().iter().map(|x| x.to_bits()).collect()
}

/// The scalar reflector loop `reflect_column` ran before the blocking.
fn reflect_column_scalar(m: &mut DenseMat, k: usize) -> f64 {
    let rows = m.nrows();
    let alpha = m.get(k, k);
    let mut tail = 0.0;
    for i in k + 1..rows {
        let v = m.get(i, k);
        tail += v * v;
    }
    if tail == 0.0 {
        return 0.0;
    }
    let norm = alpha.phypot(tail.sqrt());
    let beta = if alpha >= 0.0 { -norm } else { norm };
    let t = (beta - alpha) / beta;
    let scale = 1.0 / (alpha - beta);
    for i in k + 1..rows {
        m.set(i, k, m.get(i, k) * scale);
    }
    m.set(k, k, beta);
    for j in k + 1..m.ncols() {
        let mut w = m.get(k, j);
        for i in k + 1..rows {
            w += m.get(i, k) * m.get(i, j);
        }
        w *= t;
        m.set(k, j, m.get(k, j) - w);
        for i in k + 1..rows {
            m.set(i, j, m.get(i, j) - w * m.get(i, k));
        }
    }
    t
}

#[test]
fn matmul_tn_is_the_scalar_dot_product_bit_for_bit() {
    for (rows, left, right) in [(13, 7, 9), (40, 17, 5), (3, 1, 1), (64, 33, 31)] {
        let a = awkward(rows, left, 0.37);
        let b = awkward(rows, right, 0.91);
        let scalar = DenseMat::from_fn(left, right, |i, j| dot(a.col(i), b.col(j)));
        assert_eq!(
            bits(&matmul_tn(&a, &b)),
            bits(&scalar),
            "{rows}x{left} · {rows}x{right}"
        );
    }
}

#[test]
fn householder_reflections_are_the_scalar_loop_bit_for_bit() {
    for (rows, cols, left) in [
        (17, 11, 11),
        (30, 9, 9),
        (12, 12, 7),
        (25, 14, 3),
        (6, 5, 5),
    ] {
        let input = awkward(rows, cols, 0.53);
        let mut scalar = input.clone();
        let tau: Vec<f64> = (0..left.min(rows))
            .map(|k| reflect_column_scalar(&mut scalar, k))
            .collect();
        let blocked = HouseholderQr::factor_left(input, left);
        assert_eq!(
            bits(&blocked.qr),
            bits(&scalar),
            "{rows}x{cols} left {left}"
        );
        let t = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(t(&blocked.tau), t(&tau));
    }
}
