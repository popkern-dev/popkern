//! Deterministic test matrices with a known spectrum.

use popkern_math::PortableMath;

/// A 64-bit linear congruential generator (Knuth's MMIX constants):
/// the same sequence on every target.
pub(crate) struct Lcg(pub(crate) u64);

impl Lcg {
    /// The next value, uniform in [−1, 1).
    pub(crate) fn next_signed(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 52) as f64) - 1.0
    }
}

/// Row-major n×n product a·b.
pub(crate) fn matmul(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..n {
            let aik = a[i * n + k];
            for j in 0..n {
                c[i * n + j] += aik * b[k * n + j];
            }
        }
    }
    c
}

/// An orthogonal matrix: a product of n Householder reflectors with
/// generated vectors.
pub(crate) fn orthogonal(n: usize, rng: &mut Lcg) -> Vec<f64> {
    let mut q = vec![0.0; n * n];
    for i in 0..n {
        q[i * n + i] = 1.0;
    }
    for _ in 0..n {
        let v: Vec<f64> = (0..n).map(|_| rng.next_signed()).collect();
        let vv: f64 = v.iter().map(|x| x * x).sum();
        if vv == 0.0 {
            continue;
        }
        let mut h = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                let identity = if i == j { 1.0 } else { 0.0 };
                h[i * n + j] = identity - 2.0 * v[i] * v[j] / vv;
            }
        }
        q = matmul(&q, &h, n);
    }
    q
}

/// Eigenvalues spaced geometrically from 1 down to 1/κ.
pub(crate) fn spectrum(n: usize, condition: f64) -> Vec<f64> {
    if n == 1 {
        return vec![1.0];
    }
    (0..n)
        .map(|k| condition.ppowf(-(k as f64) / ((n - 1) as f64)))
        .collect()
}

/// Q·diag(λ)·Qᵀ, symmetrised exactly.
pub(crate) fn with_spectrum(n: usize, lambda: &[f64], seed: u64) -> Vec<f64> {
    let mut rng = Lcg(seed);
    let q = orthogonal(n, &mut rng);
    let mut ql = q.clone();
    for i in 0..n {
        for j in 0..n {
            ql[i * n + j] *= lambda[j];
        }
    }
    let mut qt = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            qt[j * n + i] = q[i * n + j];
        }
    }
    let mut a = matmul(&ql, &qt, n);
    for i in 0..n {
        for j in 0..i {
            let m = 0.5 * (a[i * n + j] + a[j * n + i]);
            a[i * n + j] = m;
            a[j * n + i] = m;
        }
    }
    a
}

/// Frobenius norm.
pub(crate) fn frobenius(a: &[f64]) -> f64 {
    a.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// FNV-1a over the bit patterns of `values`.
pub(crate) fn bits_hash(values: &[f64]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for v in values {
        for byte in v.to_bits().to_le_bytes() {
            h ^= u64::from(byte);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    h
}
