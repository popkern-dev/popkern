//! De Casteljau subdivision helpers for KP20 bi-6 patches.

#[path = "kp20_polar_support_cell_boundary_projection.rs"]
pub(super) mod support_cell_boundary_projection;

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::Kp20Bi6Patch;

pub(super) fn split_bi6_patch_u(patch: &Kp20Bi6Patch, t: f64) -> (Kp20Bi6Patch, Kp20Bi6Patch) {
    let mut left = Kp20Bi6Patch::zero();
    let mut right = Kp20Bi6Patch::zero();
    for column in 0..7 {
        let mut curve = [Point3::new(0.0, 0.0, 0.0); 7];
        for (row, point) in curve.iter_mut().enumerate() {
            *point = patch.cps[row][column];
        }
        let (left_curve, right_curve) = split_bi6_curve(curve, t);
        for row in 0..7 {
            left.cps[row][column] = left_curve[row];
            right.cps[row][column] = right_curve[row];
        }
    }
    (left, right)
}

pub(super) fn split_bi6_patch_v(patch: &Kp20Bi6Patch, t: f64) -> (Kp20Bi6Patch, Kp20Bi6Patch) {
    let mut lower = Kp20Bi6Patch::zero();
    let mut upper = Kp20Bi6Patch::zero();
    for row in 0..7 {
        let (lower_curve, upper_curve) = split_bi6_curve(patch.cps[row], t);
        lower.cps[row] = lower_curve;
        upper.cps[row] = upper_curve;
    }
    (lower, upper)
}

fn split_bi6_curve(mut levels: [Point3; 7], t: f64) -> ([Point3; 7], [Point3; 7]) {
    let mut left = [Point3::new(0.0, 0.0, 0.0); 7];
    let mut right = [Point3::new(0.0, 0.0, 0.0); 7];
    left[0] = levels[0];
    right[6] = levels[6];
    for degree in 1..=6 {
        for idx in 0..=(6 - degree) {
            levels[idx] = lerp(levels[idx], levels[idx + 1], t);
        }
        left[degree] = levels[0];
        right[6 - degree] = levels[6 - degree];
    }
    (left, right)
}

fn lerp(a: Point3, b: Point3, t: f64) -> Point3 {
    a + (b - a) * t
}

/// Splits a bi-6 patch into `count` sub-patches uniformly along the U axis.
///
/// Body retained verbatim from the pre-refactor implementation so the
/// existing uniform-path kernel output stays bit-identical (a chained
/// FP-division reformulation drifts by 1 ULP at `count = 3`). Stage A3
/// of the cap-vertex-aware split design routes this through
/// [`split_by_params`] only when the resulting sequence agrees with
/// `i / count` to within `UNIFORM_EQUIVALENCE_TOL` — at which point the
/// drift question is moot because the legacy formula is used directly.
pub(super) fn split_uniform_u(base: &Kp20Bi6Patch, count: usize) -> Vec<Kp20Bi6Patch> {
    split_uniform_legacy(base, count, split_bi6_patch_u)
}

/// Splits a bi-6 patch into `count` sub-patches uniformly along the V axis.
/// See [`split_uniform_u`] for the bit-identity rationale.
pub(super) fn split_uniform_v(base: &Kp20Bi6Patch, count: usize) -> Vec<Kp20Bi6Patch> {
    split_uniform_legacy(base, count, split_bi6_patch_v)
}

fn split_uniform_legacy(
    base: &Kp20Bi6Patch,
    count: usize,
    split: fn(&Kp20Bi6Patch, f64) -> (Kp20Bi6Patch, Kp20Bi6Patch),
) -> Vec<Kp20Bi6Patch> {
    let mut patches = Vec::with_capacity(count);
    let mut remainder = *base;
    for index in 0..count {
        if index + 1 == count {
            patches.push(remainder);
        } else {
            let local_t = 1.0 / ((count - index) as f64);
            let (head, tail) = split(&remainder, local_t);
            patches.push(head);
            remainder = tail;
        }
    }
    patches
}

/// Splits `base` at the explicit global parameter sequence `params`
/// (each in (0, 1), strictly ascending). Returns `params.len() + 1`
/// sub-patches; the i-th spans the global parameter range
/// `[params[i-1], params[i]]` (with `params[-1] = 0` and `params[len] = 1`).
///
/// Used by the cap-vertex-aware support-cell split path
/// (Stage A3 of `docs/decisions/2026-06-21-kp20-cap-vertex-aware-
/// support-cell-split.md`). Not used by `split_uniform_{u,v}` — those
/// keep the legacy formula to avoid the 1-ULP drift of chained
/// FP-division reformulations at `count = 3` and similar.
#[allow(dead_code)] // first caller lands in Stage A3
pub(super) fn split_by_params(
    base: &Kp20Bi6Patch,
    params: &[f64],
    split: fn(&Kp20Bi6Patch, f64) -> (Kp20Bi6Patch, Kp20Bi6Patch),
) -> Vec<Kp20Bi6Patch> {
    let mut patches = Vec::with_capacity(params.len() + 1);
    let mut remainder = *base;
    let mut last_global = 0.0;
    for &global_t in params {
        let local_t = (global_t - last_global) / (1.0 - last_global);
        let (head, tail) = split(&remainder, local_t);
        patches.push(head);
        remainder = tail;
        last_global = global_t;
    }
    patches.push(remainder);
    patches
}

#[cfg(test)]
mod tests {
    use super::*;

    fn synthetic_bi6_patch() -> Kp20Bi6Patch {
        let mut patch = Kp20Bi6Patch::zero();
        for row in 0..7 {
            for column in 0..7 {
                patch.cps[row][column] = Point3::new(
                    row as f64 / 6.0,
                    column as f64 / 6.0,
                    (row * column) as f64 / 36.0,
                );
            }
        }
        patch
    }

    fn max_cp_distance(a: &Kp20Bi6Patch, b: &Kp20Bi6Patch) -> f64 {
        let mut worst = 0.0_f64;
        for row in 0..7 {
            for column in 0..7 {
                worst = worst.max(a.cps[row][column].distance_to(&b.cps[row][column]));
            }
        }
        worst
    }

    /// `split_by_params` with the uniform global sequence
    /// `[1/N, 2/N, ..., (N-1)/N]` partitions the parameter space the
    /// same way as `split_uniform_*`, within floating-point precision.
    /// 1e-12 is two orders above the worst observed drift (1 ULP ≈
    /// 1.1e-16 at unit-extent CPs); confirms structural equivalence
    /// without locking the test to a particular FP-rounding regime.
    #[test]
    fn split_by_params_uniform_matches_split_uniform_within_fp_precision() {
        let base = synthetic_bi6_patch();
        for &count in &[1, 2, 3, 4] {
            let params: Vec<f64> = (1..count).map(|i| i as f64 / count as f64).collect();

            let by_params_u = split_by_params(&base, &params, split_bi6_patch_u);
            let uniform_u = split_uniform_u(&base, count);
            assert_eq!(by_params_u.len(), uniform_u.len(), "u count={count}");
            for (sub, (a, b)) in by_params_u.iter().zip(uniform_u.iter()).enumerate() {
                let drift = max_cp_distance(a, b);
                assert!(
                    drift < 1.0e-12,
                    "u count={count} sub={sub} drift={drift:.3e}"
                );
            }

            let by_params_v = split_by_params(&base, &params, split_bi6_patch_v);
            let uniform_v = split_uniform_v(&base, count);
            assert_eq!(by_params_v.len(), uniform_v.len(), "v count={count}");
            for (sub, (a, b)) in by_params_v.iter().zip(uniform_v.iter()).enumerate() {
                let drift = max_cp_distance(a, b);
                assert!(
                    drift < 1.0e-12,
                    "v count={count} sub={sub} drift={drift:.3e}"
                );
            }
        }
    }
}
