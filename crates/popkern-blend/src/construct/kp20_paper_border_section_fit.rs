//! The exact strip-section rim rows for the setback border (PK-398):
//! one two-piece degree-6 spline per arc, C⁰/C¹/C² at the spoke point
//! `m` by construction, with the arc's corner points interpolated
//! exactly (adjacent arcs share corners bit-for-bit). Fitting the
//! halves independently was the measured author of the earlier
//! prototype's skew-corner regression (octagon8_mixed spoke transverse
//! G2 2.6e-1): independent fits agree at `m` in value only.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{bernstein6, eval_sextic_row};

/// Number of section samples per arc (uniform in the span parameter).
const SAMPLES: usize = 65;

/// One arc's fitted section spline.
#[derive(Debug, Clone, Copy)]
pub(super) struct SectionSplineFit {
    /// `[corner_prev → m]`, `[m → corner_next]`.
    pub halves: [[Point3; 7]; 2],
    /// Worst sample distance between the spline and the section.
    pub residual: f64,
}

/// Fit the section `sample(t)`, `t ∈ [0, 1]`, as a C² two-piece
/// degree-6 spline with `sample(0) = c_start` and `sample(1) = c_end`
/// pinned. Returns `None` when the KKT system is singular.
pub(super) fn fit_section_spline(
    sample: &dyn Fn(f64) -> Point3,
    c_start: Point3,
    c_end: Point3,
) -> Option<SectionSplineFit> {
    // Unknowns per axis: A1..A6 (slots 0..6) and B0..B5 (slots 6..12).
    const VARS: usize = 12;
    const CONS: usize = 3;
    let mut gram = [[0.0_f64; VARS]; VARS];
    let mut moment = [[0.0_f64; 3]; VARS];
    let samples: Vec<(f64, Point3)> = (0..SAMPLES)
        .map(|k| {
            let t = k as f64 / (SAMPLES - 1) as f64;
            (t, sample(t))
        })
        .collect();
    for &(t, target) in &samples {
        let (basis, offset, pinned, pinned_weight) = if t <= 0.5 {
            let b = bernstein6(2.0 * t);
            (b, 0_usize, c_start, b[0])
        } else {
            let b = bernstein6(2.0 * t - 1.0);
            (b, 6_usize, c_end, b[6])
        };
        let r = [
            target.x - pinned_weight * pinned.x,
            target.y - pinned_weight * pinned.y,
            target.z - pinned_weight * pinned.z,
        ];
        // Half A uses basis indices 1..=6 → slots 0..6; half B uses
        // basis indices 0..=5 → slots 6..12.
        let coeff = |slot: usize| -> f64 {
            if offset == 0 {
                basis[slot + 1]
            } else {
                basis[slot - 6]
            }
        };
        let range = if offset == 0 { 0..6 } else { 6..12 };
        for i in range.clone() {
            for j in range.clone() {
                gram[i][j] += coeff(i) * coeff(j);
            }
            for c in 0..3 {
                moment[i][c] += coeff(i) * r[c];
            }
        }
    }
    // C⁰: A6 − B0 = 0; C¹: A6 − A5 − B1 + B0 = 0;
    // C²: A4 − 2A5 + A6 − B0 + 2B1 − B2 = 0 (equal half lengths).
    let constraints: [[f64; VARS]; CONS] = {
        let mut rows = [[0.0_f64; VARS]; CONS];
        rows[0][5] = 1.0;
        rows[0][6] = -1.0;
        rows[1][5] = 1.0;
        rows[1][4] = -1.0;
        rows[1][7] = -1.0;
        rows[1][6] = 1.0;
        rows[2][3] = 1.0;
        rows[2][4] = -2.0;
        rows[2][5] = 1.0;
        rows[2][6] = -1.0;
        rows[2][7] = 2.0;
        rows[2][8] = -1.0;
        rows
    };
    const SIZE: usize = VARS + CONS;
    let mut matrix = [[0.0_f64; SIZE]; SIZE];
    for i in 0..VARS {
        matrix[i][..VARS].copy_from_slice(&gram[i]);
        for (c, row) in constraints.iter().enumerate() {
            matrix[i][VARS + c] = row[i];
            matrix[VARS + c][i] = row[i];
        }
    }
    let mut rhs = [[0.0_f64; 3]; SIZE];
    rhs[..VARS].copy_from_slice(&moment[..VARS]);
    let solution = solve_dense_15(&matrix, &rhs)?;
    let at = |slot: usize| Point3::new(solution[slot][0], solution[slot][1], solution[slot][2]);
    let first: [Point3; 7] = [c_start, at(0), at(1), at(2), at(3), at(4), at(5)];
    let second: [Point3; 7] = [at(6), at(7), at(8), at(9), at(10), at(11), c_end];
    let mut residual = 0.0_f64;
    for &(t, target) in &samples {
        let p = if t <= 0.5 {
            eval_sextic_row(&first, 2.0 * t)
        } else {
            eval_sextic_row(&second, 2.0 * t - 1.0)
        };
        residual = residual.max(p.distance_to(&target));
    }
    if !residual.is_finite() {
        return None;
    }
    Some(SectionSplineFit {
        halves: [first, second],
        residual,
    })
}

/// Worst sample distance between a cubic row and the section (the
/// interpolant's sagitta), for the fit acceptance rule.
pub(super) fn cubic_section_residual(cubic: &[Point3; 4], sample: &dyn Fn(f64) -> Point3) -> f64 {
    let mut worst = 0.0_f64;
    for k in 0..SAMPLES {
        let t = k as f64 / (SAMPLES - 1) as f64;
        let u = 1.0 - t;
        let b = [u * u * u, 3.0 * u * u * t, 3.0 * u * t * t, t * t * t];
        let p = Point3::new(
            b[0] * cubic[0].x + b[1] * cubic[1].x + b[2] * cubic[2].x + b[3] * cubic[3].x,
            b[0] * cubic[0].y + b[1] * cubic[1].y + b[2] * cubic[2].y + b[3] * cubic[3].y,
            b[0] * cubic[0].z + b[1] * cubic[1].z + b[2] * cubic[2].z + b[3] * cubic[3].z,
        );
        worst = worst.max(p.distance_to(&sample(t)));
    }
    worst
}

/// Gaussian elimination with partial pivoting on the 15×15 KKT
/// system (three right-hand sides). The KKT matrix is symmetric
/// indefinite; pivoting keeps the elimination stable on the
/// well-conditioned Bernstein Gram block.
fn solve_dense_15(mat: &[[f64; 15]; 15], rhs: &[[f64; 3]; 15]) -> Option<[[f64; 3]; 15]> {
    let mut a = *mat;
    let mut b = *rhs;
    for col in 0..15 {
        let mut piv = col;
        for r in (col + 1)..15 {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-14 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let pivot = a[col];
        let pivot_b = b[col];
        for r in (col + 1)..15 {
            let f = a[r][col] / pivot[col];
            if f == 0.0 {
                continue;
            }
            for (c, pv) in pivot.iter().enumerate().skip(col) {
                a[r][c] -= f * pv;
            }
            for (c, pv) in pivot_b.iter().enumerate() {
                b[r][c] -= f * pv;
            }
        }
    }
    let mut x = [[0.0_f64; 3]; 15];
    for row in (0..15).rev() {
        for c in 0..3 {
            let mut acc = b[row][c];
            for col in (row + 1)..15 {
                acc -= a[row][col] * x[col][c];
            }
            x[row][c] = acc / a[row][row];
        }
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn circle(t: f64) -> Point3 {
        let a = core::f64::consts::FRAC_PI_2 * t;
        Point3::new(0.1 * a.cos(), 0.1 * a.sin(), 0.3 * t)
    }

    #[test]
    fn quarter_helix_section_fits_below_the_cubic_sagitta_with_c2_junction() {
        let c0 = circle(0.0);
        let c1 = circle(1.0);
        let fit = fit_section_spline(&circle, c0, c1).expect("fit");
        let cubic = {
            let q1 = circle(1.0 / 3.0);
            let q2 = circle(2.0 / 3.0);
            let inner = |a: f64, b: f64, c: f64, d: f64| {
                Point3::new(
                    (a * c0.x + b * q1.x + c * q2.x + d * c1.x) / 6.0,
                    (a * c0.y + b * q1.y + c * q2.y + d * c1.y) / 6.0,
                    (a * c0.z + b * q1.z + c * q2.z + d * c1.z) / 6.0,
                )
            };
            [
                c0,
                inner(-5.0, 18.0, -9.0, 2.0),
                inner(2.0, -9.0, 18.0, -5.0),
                c1,
            ]
        };
        let cubic_residual = cubic_section_residual(&cubic, &circle);
        assert!(
            fit.residual < 1e-3 * cubic_residual,
            "fit {:e} vs cubic {:e}",
            fit.residual,
            cubic_residual
        );
        assert!(fit.halves[0][0].distance_to(&c0) == 0.0);
        assert!(fit.halves[1][6].distance_to(&c1) == 0.0);
        let [a, b] = fit.halves;
        assert!(a[6].distance_to(&b[0]) < 1e-12);
        assert!(((a[6] - a[5]) - (b[1] - b[0])).length() < 1e-12);
        let sa = (a[4] - a[5]) - (a[5] - a[6]);
        let sb = (b[0] - b[1]) - (b[1] - b[2]);
        assert!((sa - sb).length() < 1e-11);
    }
}
