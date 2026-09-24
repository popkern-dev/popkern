//! Dense symmetric positive-definite linear solve (Cholesky) with
//! three right-hand sides, shared by the setback rim reconciles.
//! Extracted from `setback_rim_curvature` (structural line budget).

use popkern_core::geometry::Vec3;

/// Dense symmetric positive-definite solve (Cholesky, in place),
/// three right-hand sides; `None` on a non-SPD system.
pub(super) fn solve_dense_spd(ata: &mut [Vec<f64>], atb: &[[f64; 3]]) -> Option<Vec<Vec3>> {
    let n = ata.len();
    for (i, row) in ata.iter_mut().enumerate() {
        row[i] += 1e-12 * (1.0 + row[i].abs());
    }
    for j in 0..n {
        let head: f64 = ata[j].iter().take(j).map(|v| v * v).sum();
        let d = ata[j][j] - head;
        if d <= 0.0 || !d.is_finite() {
            return None;
        }
        let d = d.sqrt();
        ata[j][j] = d;
        for i in (j + 1)..n {
            let cross: f64 = ata[i]
                .iter()
                .zip(ata[j].iter())
                .take(j)
                .map(|(a, b)| a * b)
                .sum();
            ata[i][j] = (ata[i][j] - cross) / d;
        }
    }
    // All three axes at once with Vec3 arithmetic.
    let b: Vec<Vec3> = atb.iter().map(|r| Vec3::new(r[0], r[1], r[2])).collect();
    let mut y = vec![Vec3::zero(); n];
    for i in 0..n {
        let mut s = b[i];
        for (k, yk) in y.iter().enumerate().take(i) {
            s = s - *yk * ata[i][k];
        }
        y[i] = s * (1.0 / ata[i][i]);
    }
    let mut x = vec![Vec3::zero(); n];
    for i in (0..n).rev() {
        let mut s = y[i];
        for (k, xk) in x.iter().enumerate().skip(i + 1) {
            s = s - *xk * ata[k][i];
        }
        x[i] = s * (1.0 / ata[i][i]);
    }
    Some(x)
}
