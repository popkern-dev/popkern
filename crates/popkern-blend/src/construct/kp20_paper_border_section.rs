//! Exact strip-section rim rows for the setback border (PK-398).
//!
//! The setback fan's rim row used to be the 4-point cubic interpolant
//! of the strip's corner-to-corner section; its sagitta (7.775e-4 at
//! r = 0.1 against the degree-5 G2 section) was the measured source
//! of the served cap↔strip tangent jump (6.85e-2, PK-392). This
//! module fits the section as ONE two-piece degree-6 spline per arc:
//! both halves share the spoke point `m` with C⁰/C¹/C² agreement by
//! construction (the spoke strip rows consume exactly those three end
//! coefficients on both sides), and the arc's corner points are
//! interpolated exactly so adjacent arcs share corners bit-for-bit.
//! Fitting the halves independently was the measured author of the
//! prototype's skew-corner regression (octagon8_mixed spoke transverse
//! G2 2.6e-1): independent fits agree at `m` in value only.
//!
//! The section spline fit itself lives in `kp20_paper_border_section_fit`;
//! the ring-wide cross-field fit that carries the junction data two
//! pieces must share (spoke value/step/second difference, the corner
//! block) lives in `kp20_paper_border_ring_fit` (PK-400). This module
//! keeps the per-piece pinned fit (the ring fit's fallback) and the
//! corner twist step used by the n-gon path.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::bernstein6;

/// Fit a vector field sampled along `t ∈ [0, 1]` as a degree-6 BB row
/// with the end coefficients pinned to the samples at `t = 0` and
/// `t = 1` (which must be present) and the five interior coefficients
/// in least squares.
#[cfg(test)]
pub(super) fn fit_sextic_field(samples: &[(f64, Vec3)]) -> [Vec3; 7] {
    fit_sextic_field_pinned(samples, [None; 7])
}

/// [`fit_sextic_field`] with explicit coefficient pins: every `Some`
/// coefficient is held, the others are fitted in least squares. The
/// end coefficients default to the samples at `t = 0` / `t = 1`.
/// Junction data shared between two pieces (a corner twist step, the
/// spoke value/slope/second difference) enters here so each field is
/// one globally smooth polynomial that carries the shared data by
/// construction — a post-hoc nudge of a single coefficient was
/// measured to leak its wiggle through the de Casteljau split into
/// the neighbouring quarter (rim-jet tilt 5e-10 → 1.4e-3).
pub(super) fn fit_sextic_field_pinned(
    samples: &[(f64, Vec3)],
    pins: [Option<Vec3>; 7],
) -> [Vec3; 7] {
    let at = |target: f64| -> Vec3 {
        samples
            .iter()
            .min_by(|a, b| {
                (a.0 - target)
                    .abs()
                    .partial_cmp(&(b.0 - target).abs())
                    .unwrap_or(core::cmp::Ordering::Equal)
            })
            .map(|s| s.1)
            .unwrap_or(Vec3::zero())
    };
    let mut fixed: [Option<Vec3>; 7] = pins;
    if fixed[0].is_none() {
        fixed[0] = Some(at(0.0));
    }
    if fixed[6].is_none() {
        fixed[6] = Some(at(1.0));
    }
    let free: Vec<usize> = (0..7).filter(|i| fixed[*i].is_none()).collect();
    let mut out = [Vec3::zero(); 7];
    for (i, slot) in out.iter_mut().enumerate() {
        if let Some(v) = fixed[i] {
            *slot = v;
        }
    }
    if free.is_empty() {
        return out;
    }
    let m = free.len();
    let mut gram = vec![vec![0.0_f64; m]; m];
    let mut moment = vec![Vec3::zero(); m];
    for &(t, value) in samples {
        let b = bernstein6(t);
        let mut r = value;
        for (i, v) in fixed.iter().enumerate() {
            if let Some(v) = v {
                r = r - *v * b[i];
            }
        }
        for (fi, &i) in free.iter().enumerate() {
            for (fj, &j) in free.iter().enumerate() {
                gram[fi][fj] += b[i] * b[j];
            }
            moment[fi] = moment[fi] + r * b[i];
        }
    }
    match solve_small(&gram, &moment) {
        Some(solution) => {
            for (fi, &i) in free.iter().enumerate() {
                out[i] = solution[fi];
            }
        }
        None => {
            // Degenerate sampling: linear interpolation of the pins.
            let (v0, v6) = (out[0], out[6]);
            for &i in &free {
                out[i] = v0 + (v6 - v0) * (i as f64 / 6.0);
            }
        }
    }
    out
}

/// Shared cross-jet twist step at a sector corner `c` between piece
/// `a` (ending at `c`) and piece `b` (starting at `c`).
///
/// The KP20 corner block's `(5, 5)` slot is claimed by both bands as
/// `c + L_a + L_b + Δ/4` — `L` the c-quarter rim legs from `c`, `Δ`
/// the step of each piece's cross-jet field next to `c` (its mixed
/// derivative). The two claims agree only for one shared `Δ`, and the
/// slot's legs must stay in both strips' tangent planes at the
/// adjacent stations (`n_a`, `n_b`): `n_a·(L_b + Δ/4) = 0`,
/// `n_b·(L_a + Δ/4) = 0`. `Δ` is the minimum-norm particular solution
/// of those two rows plus the component of the pieces' natural mean
/// step along the planes' intersection line (the homogeneous
/// freedom). Perpendicular flat corners (the cube) give a purely
/// in-plane twist; skew corners get the compensating normal part.
pub(super) fn corner_twist_step(
    natural_a: Vec3,
    natural_b: Vec3,
    leg_a: Vec3,
    leg_b: Vec3,
    n_a: Vec3,
    n_b: Vec3,
) -> Vec3 {
    let mean = (natural_a + natural_b) * 0.5;
    let rhs = [-4.0 * n_a.dot(&leg_b), -4.0 * n_b.dot(&leg_a)];
    let line = n_a.cross(&n_b);
    if line.length() > 1.0e-6 {
        let g = [
            [n_a.dot(&n_a), n_a.dot(&n_b)],
            [n_b.dot(&n_a), n_b.dot(&n_b)],
        ];
        let det = g[0][0] * g[1][1] - g[0][1] * g[1][0];
        let alpha = (rhs[0] * g[1][1] - rhs[1] * g[0][1]) / det;
        let beta = (g[0][0] * rhs[1] - g[1][0] * rhs[0]) / det;
        let particular = n_a * alpha + n_b * beta;
        let line = line * (1.0 / line.length());
        particular + line * mean.dot(&line)
    } else {
        // Coincident planes: one normal row, the rest in-plane.
        match n_a.normalized() {
            Some(n) => {
                let normal_part = 0.5 * (rhs[0] + rhs[1]) / n_a.length().max(1e-300);
                mean - n * mean.dot(&n) + n * normal_part
            }
            None => mean,
        }
    }
}

fn solve_small(mat: &[Vec<f64>], rhs: &[Vec3]) -> Option<Vec<Vec3>> {
    let n = rhs.len();
    let mut a: Vec<Vec<f64>> = mat.to_vec();
    let mut b: Vec<Vec3> = rhs.to_vec();
    for col in 0..n {
        let mut piv = col;
        for r in (col + 1)..n {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-14 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let pivot = a[col].clone();
        let pivot_b = b[col];
        for r in (col + 1)..n {
            let f = a[r][col] / pivot[col];
            for (c, pv) in pivot.iter().enumerate().skip(col) {
                a[r][c] -= f * pv;
            }
            b[r] = b[r] - pivot_b * f;
        }
    }
    let mut x = vec![Vec3::zero(); n];
    for row in (0..n).rev() {
        let mut acc = b[row];
        for col in (row + 1)..n {
            acc = acc - x[col] * a[row][col];
        }
        x[row] = acc * (1.0 / a[row][row]);
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pinned_fit_holds_every_pin_and_reproduces_polynomial_content() {
        let truth: [Vec3; 7] = core::array::from_fn(|i| {
            let x = i as f64;
            Vec3::new(0.1 * x, (x * x) * 0.02 - 0.3, 0.05 * (x - 3.0).powi(3))
        });
        let samples: Vec<(f64, Vec3)> = (0..25)
            .map(|k| {
                let t = k as f64 / 24.0;
                let b = bernstein6(t);
                let mut acc = Vec3::zero();
                for (w, c) in b.iter().zip(&truth) {
                    acc = acc + *c * *w;
                }
                (t, acc)
            })
            .collect();
        let free = fit_sextic_field(&samples);
        for (a, b) in free.iter().zip(&truth) {
            assert!((*a - *b).length() < 1e-10);
        }
        let mut pins = [None; 7];
        pins[4] = Some(truth[4]);
        pins[5] = Some(truth[5]);
        pins[6] = Some(truth[6]);
        pins[1] = Some(truth[1]);
        let pinned = fit_sextic_field_pinned(&samples, pins);
        for (a, b) in pinned.iter().zip(&truth) {
            assert!((*a - *b).length() < 1e-10);
        }
        assert_eq!(pinned[6], truth[6]);
        let moved = Vec3::new(0.0, 0.0, 0.3);
        let mut pins = [None; 7];
        pins[0] = Some(truth[0] + moved);
        pins[1] = Some(truth[1] + moved);
        let pinned = fit_sextic_field_pinned(&samples, pins);
        assert_eq!(pinned[0], truth[0] + moved);
        assert_eq!(pinned[1], truth[1] + moved);
        assert!((pinned[6] - truth[6]).length() < 1e-12);
    }

    #[test]
    fn corner_twist_step_keeps_both_corner_legs_in_plane() {
        let tilt = (7.5_f64).to_radians();
        let n_a = Vec3::new(0.0, tilt.cos(), tilt.sin());
        let n_b = Vec3::new(tilt.sin(), tilt.cos(), 0.0);
        let leg_a = Vec3::new(0.0, 0.0, 0.01);
        let leg_b = Vec3::new(0.01, 0.0, 0.0);
        let step = corner_twist_step(
            Vec3::new(0.002, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.002),
            leg_a,
            leg_b,
            n_a,
            n_b,
        );
        assert!(n_a.dot(&(leg_b + step * 0.25)).abs() < 1e-12);
        assert!(n_b.dot(&(leg_a + step * 0.25)).abs() < 1e-12);
        let skew_leg_b = Vec3::new(0.01, 0.0, 0.004);
        let step = corner_twist_step(Vec3::zero(), Vec3::zero(), leg_a, skew_leg_b, n_a, n_b);
        assert!(n_a.dot(&(skew_leg_b + step * 0.25)).abs() < 1e-12);
        assert!(n_b.dot(&(leg_a + step * 0.25)).abs() < 1e-12);
    }
}
