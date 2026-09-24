//! Fairing rows for the KP18/KP20 paper guide, and the null-space
//! selector that resolves the guide where the paper construction
//! leaves it underdetermined (PK-518).
//!
//! The guide fit (border corner jets of `g ∘ χ̃`) plus the KP18 join
//! constraints and the centre pin determine every guide coefficient
//! at n >= 5. At n = 3 and at the regular valence n = 4 they do not:
//! the KKT system is singular (condition `inf`; ~1e6..1e7 at every
//! n >= 5), and the solve's absolute 1e-12 ridge then picked the
//! null-space component by distance to the WORLD ORIGIN — a guide
//! whose interior depends on where the model sits. At n = 4 that was
//! wild even on the characteristic border (the paper's own validation
//! universe, where the guide must be the identity): the strip
//! projection had to move the sampled cap by 2.7.
//!
//! The selector is lexicographic in spirit: fit first, then — inside
//! the null space — the fairest guide, measured by second differences
//! along the three lattice directions of each sector's total-degree-6
//! net. Those rows vanish EXACTLY on affine maps (including at the
//! net's boundary, unlike the one-sided Laplacian of the recorded
//! shape-conditioning fairing), so an exactly representable fit —
//! the identity on characteristic data — is reproduced exactly, and
//! only the undetermined directions are steered. Well-conditioned
//! systems (n >= 5) never see the selector and keep their
//! paper-faithful solution bit for bit.

use super::kp20_bi6_guide_lsq::{LinearRow, var_index};

/// Condition number above which the paper guide's KKT system is
/// numerically singular (measured: n = 3, 4 `inf`; n >= 5 at most
/// ~1e7).
const SINGULAR_GUIDE_KKT_CONDITION: f64 = 1.0e12;

/// Weight of the null-space selector rows against the unit fit rows.
/// Exact on affine fits, so its only range-space effect is on fits
/// the guide cannot represent exactly (order `weight²`); large enough
/// that `weight²` out-votes the 1e-12 KKT ridge in the null space.
pub(super) const NULL_SPACE_SELECTOR_WEIGHT: f64 = 1.0e-4;

/// Whether an unregularized guide solve's KKT system was numerically
/// singular.
pub(super) fn guide_kkt_is_singular(kkt_condition_number: f64) -> bool {
    !kkt_condition_number.is_finite() || kkt_condition_number > SINGULAR_GUIDE_KKT_CONDITION
}

/// Append the null-space selector rows (target zero): for every
/// sector, the second difference of the total-degree-6 control net
/// along each lattice direction `(1, 0)`, `(0, 1)` and `(-1, 1)`
/// wherever it fits, weighted by `weight`; no-op at zero weight.
pub(super) fn push_null_space_selector_rows(
    rows: &mut Vec<LinearRow>,
    targets: &mut [Vec<f64>; 3],
    n: usize,
    weight: f64,
) {
    if weight <= 0.0 {
        return;
    }
    let at = |sector: usize, i: usize, j: usize| var_index(sector, super::tri6_flat_index(6, i, j));
    for sector in 0..n {
        for i in 0..=6_usize {
            for j in 0..=(6 - i) {
                let mut triples: Vec<[(usize, usize); 3]> = Vec::with_capacity(3);
                if i + 2 + j <= 6 {
                    triples.push([(i, j), (i + 1, j), (i + 2, j)]);
                    triples.push([(i, j), (i, j + 1), (i, j + 2)]);
                    triples.push([(i + 2, j), (i + 1, j + 1), (i, j + 2)]);
                }
                for [a, b, c] in triples {
                    rows.push(LinearRow {
                        terms: vec![
                            (at(sector, a.0, a.1), weight),
                            (at(sector, b.0, b.1), -2.0 * weight),
                            (at(sector, c.0, c.1), weight),
                        ],
                    });
                    for axis_targets in targets.iter_mut() {
                        axis_targets.push(0.0);
                    }
                }
            }
        }
    }
}

/// Append the discrete-Laplacian fairing rows (target zero) for every
/// sector's total-degree-6 control net, weighted by `weight`; no-op
/// at zero weight.
pub(super) fn push_guide_fairing_rows(
    rows: &mut Vec<LinearRow>,
    targets: &mut [Vec<f64>; 3],
    n: usize,
    weight: f64,
) {
    if weight <= 0.0 {
        return;
    }
    for sector in 0..n {
        for i in 0..=6_usize {
            for j in 0..=(6 - i) {
                let neighbors: [(i64, i64); 6] =
                    [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)];
                let mut terms = vec![(var_index(sector, super::tri6_flat_index(6, i, j)), 0.0_f64)];
                let mut count = 0.0;
                for (di, dj) in neighbors {
                    let ni = i as i64 + di;
                    let nj = j as i64 + dj;
                    if ni < 0 || nj < 0 || ni + nj > 6 {
                        continue;
                    }
                    terms.push((
                        var_index(sector, super::tri6_flat_index(6, ni as usize, nj as usize)),
                        -weight,
                    ));
                    count += 1.0;
                }
                if count < 2.0 {
                    continue;
                }
                terms[0].1 = weight * count;
                rows.push(LinearRow { terms });
                for axis_targets in targets.iter_mut() {
                    axis_targets.push(0.0);
                }
            }
        }
    }
}
