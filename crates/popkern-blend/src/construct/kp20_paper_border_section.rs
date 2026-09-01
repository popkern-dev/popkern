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
//! this module owns the cross-field fits and the junction data the two
//! pieces of a sector must share (spoke value/step/second difference,
//! the corner block).

use popkern_core::geometry::gregory::bernstein6;
use popkern_core::geometry::{Point3, Vec3};

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

/// Junction data of a field row at one end, oriented ALONG the arc:
/// `(value, step, second)` with `step = b₁ − b₀` and
/// `second = b₂ − 2b₁ + b₀` at the start, mirrored at the end.
pub(super) fn end_junction(row: &[Vec3; 7], at_end: bool) -> (Vec3, Vec3, Vec3) {
    if at_end {
        (
            row[6],
            row[6] - row[5],
            (row[6] - row[5]) - (row[5] - row[4]),
        )
    } else {
        (
            row[0],
            row[1] - row[0],
            (row[2] - row[1]) - (row[1] - row[0]),
        )
    }
}

/// Mean of two ends' junction data (a C² junction shared by both).
pub(super) fn shared_junction(a: (Vec3, Vec3, Vec3), b: (Vec3, Vec3, Vec3)) -> (Vec3, Vec3, Vec3) {
    ((a.0 + b.0) * 0.5, (a.1 + b.1) * 0.5, (a.2 + b.2) * 0.5)
}

/// Coefficient pins that make a row END with the shared junction data.
pub(super) fn pins_ending_with(j: (Vec3, Vec3, Vec3), pins: &mut [Option<Vec3>; 7]) {
    let (value, step, second) = j;
    pins[6] = Some(value);
    pins[5] = Some(value - step);
    pins[4] = Some(value - step * 2.0 + second);
}

/// Coefficient pins that make a row START with the shared junction data.
pub(super) fn pins_starting_with(j: (Vec3, Vec3, Vec3), pins: &mut [Option<Vec3>; 7]) {
    let (value, step, second) = j;
    pins[0] = Some(value);
    pins[1] = Some(value + step);
    pins[2] = Some(value + step * 2.0 + second);
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

/// De Casteljau halves of a degree-6 vector row at `t = ½`.
pub(super) fn split_sextic_vec_row(row: &[Vec3; 7]) -> ([Vec3; 7], [Vec3; 7]) {
    let mut levels = [[Vec3::zero(); 7]; 7];
    levels[0] = *row;
    for level in 1..7 {
        for i in 0..(7 - level) {
            levels[level][i] = (levels[level - 1][i] + levels[level - 1][i + 1]) * 0.5;
        }
    }
    (
        core::array::from_fn(|i| levels[i][0]),
        core::array::from_fn(|i| levels[6 - i][i]),
    )
}

/// The KP20 corner block shared by the two pieces meeting at `c`
/// (PK-398). The band placements claim the corner piece's 3×3 block at
/// `c` from BOTH pieces — A's rim row with its depth rows in the
/// B-direction, B's rim row with its depth rows in the A-direction —
/// and average the four interior slots. Independently constructed
/// pieces disagree there by tangential terms living in slightly
/// different tangent planes, which the average turns into normal
/// cross-curvature two stations from `c` (measured 1.4e-1 on the cube
/// rim+verticals class). This solver takes the MEAN block and returns
/// the c-quarter jet/curvature coefficients (indices 4..6 for A, 0..2
/// for B) each piece must carry to reproduce it, so the assembly's
/// average is a no-op by construction. Inputs are the c-QUARTER rows
/// (A: rim `a0`, jets `a1`, curvatures `a2`; B likewise): the band's
/// depth-1 row is `r0 + ½d1`, its depth-2 row `r0 + d1 + 0.2·d2`.
/// Returns `(a1[4..=6], a2[4..=6], b1[0..=2], b2[0..=2])` (a
/// [`CornerBlockPins`]).
///
/// `n_a` / `n_b` are the strips' unit normals at the pieces' quarter
/// stations `[4, 5]` (A) and `[1, 2]` (B). The jet pins are projected
/// into those tangent planes and the curvature pins re-solved around
/// them: on G2 input the mean block is in-plane already (no-op), but
/// on G1 rolling-ball input the two strips disagree in curvature at
/// `c` (each section carries `1/r` where the other is straight) and
/// the mean would push half of that disagreement into the JETS —
/// measured as rim tangency 1.96e-1 on the blunt pentagon G1 control.
/// The curvature rows are the honest home of a G2 incompatibility.
/// Per-piece c-quarter pins reproducing a shared corner block:
/// `(a1[4..=6], a2[4..=6], b1[0..=2], b2[0..=2])`.
pub(super) type CornerBlockPins = ([Vec3; 3], [Vec3; 3], [Vec3; 3], [Vec3; 3]);

#[allow(clippy::too_many_arguments)]
pub(super) fn shared_corner_block(
    a0: &[Point3; 7],
    a1: &[Vec3; 7],
    a2: &[Vec3; 7],
    b0: &[Point3; 7],
    b1: &[Vec3; 7],
    b2: &[Vec3; 7],
    n_a: [Vec3; 2],
    n_b: [Vec3; 2],
) -> CornerBlockPins {
    let in_plane = |v: Vec3, n: Vec3| match n.normalized() {
        Some(n) => v - n * v.dot(&n),
        None => v,
    };
    let depth1 = |r0: Point3, d1: Vec3| r0 + d1 * 0.5;
    let depth2 = |r0: Point3, d1: Vec3, d2: Vec3| r0 + d1 + d2 * 0.2;
    let mean =
        |p: Point3, q: Point3| Point3::new(0.5 * (p.x + q.x), 0.5 * (p.y + q.y), 0.5 * (p.z + q.z));
    // Slot (i, j): i along A (6 at c), j along B (6 at c).
    let m55 = mean(depth1(a0[5], a1[5]), depth1(b0[1], b1[1]));
    let m45 = mean(depth1(a0[4], a1[4]), depth2(b0[1], b1[1], b2[1]));
    let m54 = mean(depth2(a0[5], a1[5], a2[5]), depth1(b0[2], b1[2]));
    let m44 = mean(depth2(a0[4], a1[4], a2[4]), depth2(b0[2], b1[2], b2[2]));
    let a1_5 = in_plane((m55 - a0[5]) * 2.0, n_a[1]);
    let a1_4 = in_plane((m45 - a0[4]) * 2.0, n_a[0]);
    let a2_5 = (m54 - a0[5] - a1_5) * 5.0;
    let a2_4 = (m44 - a0[4] - a1_4) * 5.0;
    let b1_1 = in_plane((m55 - b0[1]) * 2.0, n_b[0]);
    let b1_2 = in_plane((m54 - b0[2]) * 2.0, n_b[1]);
    let b2_1 = (m45 - b0[1] - b1_1) * 5.0;
    let b2_2 = (m44 - b0[2] - b1_2) * 5.0;
    (
        [a1_4, a1_5, a1[6]],
        [a2_4, a2_5, a2[6]],
        [b1[0], b1_1, b1_2],
        [b2[0], b2_1, b2_2],
    )
}

/// Full-row pins reproducing given c-quarter end coefficients: for a
/// row ENDING at `c` (`q[0..3]` = quarter indices 4, 5, 6) the last
/// three coefficients; for a row STARTING at `c` (`q` = quarter
/// indices 0, 1, 2) the first three. Inverse of the de Casteljau
/// half: `q6 = b6`, `q5 = (b5 + b6)/2`, `q4 = (b4 + 2b5 + b6)/4`.
pub(super) fn pins_from_quarter(q: [Vec3; 3], at_end: bool, pins: &mut [Option<Vec3>; 7]) {
    if at_end {
        let b6 = q[2];
        let b5 = q[1] * 2.0 - b6;
        let b4 = q[0] * 4.0 - b5 * 2.0 - b6;
        pins[4] = Some(b4);
        pins[5] = Some(b5);
        pins[6] = Some(b6);
    } else {
        let b0 = q[0];
        let b1 = q[1] * 2.0 - b0;
        let b2 = q[2] * 4.0 - b1 * 2.0 - b0;
        pins[0] = Some(b0);
        pins[1] = Some(b1);
        pins[2] = Some(b2);
    }
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
        let junction = end_junction(&truth, true);
        pins_ending_with(junction, &mut pins);
        pins[1] = Some(truth[1]);
        let pinned = fit_sextic_field_pinned(&samples, pins);
        for (a, b) in pinned.iter().zip(&truth) {
            assert!((*a - *b).length() < 1e-10);
        }
        assert_eq!(end_junction(&pinned, true), junction);
        let moved = Vec3::new(0.0, 0.0, 0.3);
        let mut pins = [None; 7];
        pins_starting_with(
            (truth[0] + moved, truth[1] - truth[0], Vec3::zero()),
            &mut pins,
        );
        let pinned = fit_sextic_field_pinned(&samples, pins);
        assert_eq!(pinned[0], truth[0] + moved);
        assert_eq!(pinned[1], truth[0] + moved + (truth[1] - truth[0]));
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
