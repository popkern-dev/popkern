//! Stack-linear folded-G2 Jacobian block for the seam polish
//! (PK-282).
//!
//! The mixed Jacobian's folded-G2 block was per-column forward FD:
//! rebuild the ring and re-evaluate every seam sample's surface
//! derivatives for each of the D parameter columns — measured as the
//! binding term of the extreme-ladder repair once the Gregory SI
//! gate was cut (D = 504, 75 samples, ~10 bidegree-9 evaluations
//! per sample per column per iteration).
//!
//! The bi-nonic patch is a plain tensor-product Bézier, so every
//! derivative the G2 formula reads is LINEAR in the control points:
//! perturbing one CP coordinate by ε shifts the sample's derivative
//! stack by exactly ε·(basis product). This block therefore gathers
//! the stacks ONCE, and per column applies the closed-form stack
//! delta and re-runs only the pure formula (branch gates included) —
//! the same forward-difference values as the rebuild-based FD up to
//! floating-point associativity, at zero surface evaluations per
//! column. Columns whose basis sensitivity is identically zero for a
//! sample keep the exact 0.0 the full FD would have produced.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::{
    RectangularGregoryBiNonic, bernstein_binonic_all, bernstein_binonic_deriv_all,
    bernstein_binonic_second_deriv_all,
};

use super::three_quad_chart_orientation_metrics::{
    G2SampleStack, g2_residual_from_stack, g2_sample_stack,
};
use super::three_quad_chart_orientation_polish::seam_row_cells;

/// One control-point write a packed parameter column performs.
struct CpWrite {
    patch: usize,
    i: usize,
    j: usize,
    axis: usize,
}

/// Per-column CP writes, mirroring `pack_seam_rows`' order exactly:
/// per-patch seam cells first, then the shared per-spoke blocks —
/// each spoke column writes `grid[9][j]` of patch `k` AND
/// `grid[9−j][0]` of patch `k+1`.
fn column_cp_writes(n: usize) -> Vec<Vec<CpWrite>> {
    let cells = seam_row_cells();
    let mut cols = Vec::with_capacity(n * cells.len() * 3 + n * 8 * 3);
    for p in 0..n {
        for &(i, j) in &cells {
            for axis in 0..3 {
                cols.push(vec![CpWrite {
                    patch: p,
                    i,
                    j,
                    axis,
                }]);
            }
        }
    }
    for p in 0..n {
        for j in 1..=8 {
            for axis in 0..3 {
                cols.push(vec![
                    CpWrite {
                        patch: p,
                        i: 9,
                        j,
                        axis,
                    },
                    CpWrite {
                        patch: (p + 1) % n,
                        i: 9 - j,
                        j: 0,
                        axis,
                    },
                ]);
            }
        }
    }
    cols
}

fn add_axis(v: &mut Vec3, axis: usize, delta: f64) {
    match axis {
        0 => v.x += delta,
        1 => v.y += delta,
        _ => v.z += delta,
    }
}

/// `[su, sv, suu, suv, svv]` basis sensitivities of one CP `(i, j)`
/// at a sample point, from value/derivative/second-derivative basis
/// rows in each direction.
fn stack_coefficients(
    (bu0, bu1, bu2): (&[f64; 10], &[f64; 10], &[f64; 10]),
    (bv0, bv1, bv2): (&[f64; 10], &[f64; 10], &[f64; 10]),
    i: usize,
    j: usize,
) -> [f64; 5] {
    [
        bu1[i] * bv0[j],
        bu0[i] * bv1[j],
        bu2[i] * bv0[j],
        bu1[i] * bv1[j],
        bu0[i] * bv2[j],
    ]
}

/// The folded-G2 Jacobian block: `n·25` rows by `d` columns of
/// forward differences at step `eps`, weighted by `g2_weight` —
/// value-compatible with FD over the polish's `g2_block`.
pub(super) fn g2_jacobian_stack_linear(
    ring: &[RectangularGregoryBiNonic],
    d: usize,
    eps: f64,
    g2_weight: f64,
) -> Vec<Vec<f64>> {
    let n = ring.len();
    let samples: Vec<f64> = (0..25).map(|s| 0.05 + 0.90 * (s as f64) / 24.0).collect();

    // Base stacks and (weighted, zero-collapsed) residuals.
    let mut stacks: Vec<G2SampleStack> = Vec::with_capacity(n * 25);
    let mut base = Vec::with_capacity(n * 25);
    for k in 0..n {
        let a = &ring[k];
        let b = &ring[(k + 1) % n];
        for &t in &samples {
            let stack = g2_sample_stack(a, b, t);
            base.push(g2_weight * g2_residual_from_stack(&stack).unwrap_or(0.0));
            stacks.push(stack);
        }
    }

    // Basis rows. Patch `a` evaluates at (1, t): the u-direction
    // rows are constants at 1, the v-direction rows vary per
    // sample. Patch `b` evaluates at (1−t, 0): u varies, v is
    // constant at 0.
    let bu_at_1 = (
        bernstein_binonic_all(1.0),
        bernstein_binonic_deriv_all(1.0),
        bernstein_binonic_second_deriv_all(1.0),
    );
    let bv_at_0 = (
        bernstein_binonic_all(0.0),
        bernstein_binonic_deriv_all(0.0),
        bernstein_binonic_second_deriv_all(0.0),
    );
    let bv_at_t: Vec<([f64; 10], [f64; 10], [f64; 10])> = samples
        .iter()
        .map(|&t| {
            (
                bernstein_binonic_all(t),
                bernstein_binonic_deriv_all(t),
                bernstein_binonic_second_deriv_all(t),
            )
        })
        .collect();
    let bu_at_tb: Vec<([f64; 10], [f64; 10], [f64; 10])> = samples
        .iter()
        .map(|&t| {
            let tb = 1.0 - t;
            (
                bernstein_binonic_all(tb),
                bernstein_binonic_deriv_all(tb),
                bernstein_binonic_second_deriv_all(tb),
            )
        })
        .collect();

    let columns = column_cp_writes(n);
    debug_assert_eq!(columns.len(), d, "column table must mirror pack_seam_rows");

    let mut rows = vec![vec![0.0; d]; n * 25];
    for (col, writes) in columns.iter().enumerate() {
        // A patch participates in pair `p` as side `a` and in pair
        // `p − 1` as side `b`; a column's writes may land on both
        // sides of one pair (spoke blocks), so deltas accumulate
        // into one stack before the formula runs — matching the
        // simultaneous perturbation the rebuild-based FD applies.
        for pair in 0..n {
            let a_patch = pair;
            let b_patch = (pair + 1) % n;
            for (s, _) in samples.iter().enumerate() {
                let row = pair * 25 + s;
                let mut stack = stacks[row];
                let mut touched = false;
                for w in writes {
                    if w.patch == a_patch {
                        let c = stack_coefficients(
                            (&bu_at_1.0, &bu_at_1.1, &bu_at_1.2),
                            (&bv_at_t[s].0, &bv_at_t[s].1, &bv_at_t[s].2),
                            w.i,
                            w.j,
                        );
                        for (entry, coef) in stack.a.iter_mut().zip(c) {
                            if coef != 0.0 {
                                add_axis(entry, w.axis, eps * coef);
                                touched = true;
                            }
                        }
                    }
                    if w.patch == b_patch {
                        let c = stack_coefficients(
                            (&bu_at_tb[s].0, &bu_at_tb[s].1, &bu_at_tb[s].2),
                            (&bv_at_0.0, &bv_at_0.1, &bv_at_0.2),
                            w.i,
                            w.j,
                        );
                        for (entry, coef) in stack.b.iter_mut().zip(c) {
                            if coef != 0.0 {
                                add_axis(entry, w.axis, eps * coef);
                                touched = true;
                            }
                        }
                    }
                }
                if !touched {
                    continue;
                }
                let pert = g2_weight * g2_residual_from_stack(&stack).unwrap_or(0.0);
                rows[row][col] = (pert - base[row]) / eps;
            }
        }
    }
    rows
}
