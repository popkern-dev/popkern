use super::*;

#[test]
fn adaptive_sign_clear() {
    assert_eq!(adaptive_sign(1.0, 1e-10), Some(std::cmp::Ordering::Greater));
    assert_eq!(adaptive_sign(-1.0, 1e-10), Some(std::cmp::Ordering::Less));
}

#[test]
fn adaptive_sign_ambiguous() {
    assert_eq!(adaptive_sign(1e-15, 1e-10), None);
}

#[test]
fn double_precision_add() {
    let a = DoublePrecision::from_f64(1.0);
    let b = DoublePrecision::from_f64(1e-16);
    let sum = a.dd_add(b);
    // Should preserve the small value that f64 would lose
    assert!((sum.to_f64() - 1.0).abs() < 1e-14);
}

#[test]
fn double_precision_mul() {
    let a = DoublePrecision::from_f64(3.0);
    let b = DoublePrecision::from_f64(7.0);
    let product = a.dd_mul(b);
    assert!((product.to_f64() - 21.0).abs() < 1e-14);
}

#[test]
fn double_precision_sign() {
    let val = DoublePrecision::from_f64(0.5);
    assert_eq!(val.sign(1e-10), Some(std::cmp::Ordering::Greater));

    let neg = DoublePrecision::from_f64(-0.5);
    assert_eq!(neg.sign(1e-10), Some(std::cmp::Ordering::Less));

    let tiny = DoublePrecision::from_f64(1e-15);
    assert_eq!(tiny.sign(1e-10), None);
}

#[test]
fn double_precision_from_sub_matches_naive_for_well_conditioned() {
    // For well-conditioned subtraction (no catastrophic cancellation
    // relative to the input magnitudes), from_sub should agree with
    // naive a - b to full f64 precision.
    for (a, b) in [(5.0, 3.0), (1e10, 1e5), (-1.0, 2.0), (42.5, 42.5)] {
        let dp = DoublePrecision::from_sub(a, b);
        assert_eq!(
            dp.to_f64(),
            a - b,
            "from_sub({a}, {b}).to_f64() must equal a-b"
        );
    }
}

#[test]
fn double_precision_from_sub_captures_error_term() {
    // When a and b land on specific bit patterns the naive
    // subtraction `a - b` can have a small f64-representable error.
    // Construct such a case: `a = 1.0, b = 1.0 - 1e-16` — the f64
    // value of `1.0 - 1e-16` is NOT exactly representable, so
    // `a - (computed b)` introduces a tiny error that `from_sub`'s
    // `lo` term should make explicit. We verify the error term is
    // finite and non-negative in magnitude.
    let a = 1.0f64;
    let b: f64 = 1.0 - 1e-16;
    let dp = DoublePrecision::from_sub(a, b);
    assert!(dp.hi.is_finite());
    assert!(dp.lo.is_finite());
    // The reconstructed value should agree with the naive
    // subtraction (both are f64-level accurate for this case).
    let naive = a - b;
    assert!(
        (dp.to_f64() - naive).abs() <= f64::EPSILON * naive.abs().max(1.0),
        "from_sub and naive diverge: dp={}, naive={}",
        dp.to_f64(),
        naive
    );
}

#[test]
fn double_precision_from_mul_exact() {
    let dp = DoublePrecision::from_mul(3.0, 7.0);
    assert_eq!(dp.hi, 21.0);
    assert_eq!(dp.lo, 0.0);
}

#[test]
fn double_precision_dd_sub() {
    let a = DoublePrecision::from_f64(5.0);
    let b = DoublePrecision::from_f64(3.0);
    let result = a.dd_sub(b);
    assert_eq!(result.to_f64(), 2.0);
}

#[test]
fn double_precision_dd_neg() {
    let a = DoublePrecision::from_f64(2.5);
    let neg = a.dd_neg();
    assert_eq!(neg.to_f64(), -2.5);
}

#[test]
fn adaptive_dist_sq_lt_clear_inside() {
    // Two points very close together — clearly inside the 1.0
    // threshold squared = 1.0.
    let p1 = (0.0, 0.0, 0.0);
    let p2 = (0.1, 0.1, 0.1);
    assert!(adaptive_dist_sq_lt(p1, p2, 1.0));
}

#[test]
fn adaptive_dist_sq_lt_clear_outside() {
    let p1 = (0.0, 0.0, 0.0);
    let p2 = (2.0, 0.0, 0.0);
    // Distance² = 4.0, threshold² = 1.0 → false.
    assert!(!adaptive_dist_sq_lt(p1, p2, 1.0));
}

#[test]
fn adaptive_dist_sq_lt_exact_threshold() {
    // Exactly at the threshold: |p1 - p2|² = 1.0 vs threshold_sq = 1.0.
    // The f64 margin puts us in the ambiguous zone; the DP recomputation
    // should produce a deterministic answer (NOT less-than since exact
    // equality is on the boundary).
    let p1 = (0.0, 0.0, 0.0);
    let p2 = (1.0, 0.0, 0.0);
    assert!(!adaptive_dist_sq_lt(p1, p2, 1.0));
}

#[test]
fn adaptive_dist_sq_lt_near_threshold_from_below() {
    // Distance² slightly less than threshold² — the adaptive
    // predicate must return true even when f64 is on the edge.
    let p1 = (0.0, 0.0, 0.0);
    let p2 = (0.999_999_999, 0.0, 0.0);
    assert!(adaptive_dist_sq_lt(p1, p2, 1.0));
}

#[test]
fn adaptive_dist_sq_lt_degenerate_identical_points() {
    // Identical points — distance 0, clearly below any positive
    // threshold.
    let p = (42.0, -17.0, 3.5);
    assert!(adaptive_dist_sq_lt(p, p, 1e-12));
}

#[test]
fn adaptive_cross_ratio_sq_lt_parallel_vectors_are_tangent() {
    // Exactly parallel vectors: cross product is zero, sine is
    // zero, clearly below any positive threshold.
    let n1 = (1.0, 0.0, 0.0);
    let n2 = (2.0, 0.0, 0.0);
    assert!(adaptive_cross_ratio_sq_lt(n1, n2, 1e-4));
}

#[test]
fn adaptive_cross_ratio_sq_lt_perpendicular_vectors_are_not_tangent() {
    // 90° apart: sin = 1, clearly above any tiny threshold.
    let n1 = (1.0, 0.0, 0.0);
    let n2 = (0.0, 1.0, 0.0);
    assert!(!adaptive_cross_ratio_sq_lt(n1, n2, 0.1));
}

#[test]
fn adaptive_cross_ratio_sq_lt_slightly_offset_is_tangent_at_loose_threshold() {
    // Two vectors 1e-6 radians apart: sin ≈ 1e-6. A 1e-4
    // threshold accepts them as tangent; a 1e-9 threshold rejects.
    let n1 = (1.0, 0.0, 0.0);
    let n2 = (1.0, 1e-6, 0.0);
    assert!(adaptive_cross_ratio_sq_lt(n1, n2, 1e-4));
    assert!(!adaptive_cross_ratio_sq_lt(n1, n2, 1e-9));
}

#[test]
fn adaptive_cross_ratio_sq_lt_zero_vector_is_not_tangent() {
    // Degenerate: one of the normals is zero. The predicate must
    // NOT report tangent (no defined angle).
    let n1 = (1.0, 0.0, 0.0);
    let n_zero = (0.0, 0.0, 0.0);
    assert!(!adaptive_cross_ratio_sq_lt(n1, n_zero, 1.0));
    assert!(!adaptive_cross_ratio_sq_lt(n_zero, n1, 1.0));
}

#[test]
fn adaptive_cross_ratio_sq_lt_antiparallel_vectors_are_tangent() {
    // Antiparallel: cross product is zero, sine is zero — the
    // predicate should report tangent (this test matches the
    // existing intersect.rs semantics, which doesn't distinguish
    // parallel from antiparallel).
    let n1 = (1.0, 0.0, 0.0);
    let n2 = (-1.0, 0.0, 0.0);
    assert!(adaptive_cross_ratio_sq_lt(n1, n2, 1e-4));
}

#[test]
fn adaptive_abs_plane_distance_lt_on_plane_is_below_threshold() {
    // Point on the plane: signed distance is zero, clearly below
    // any positive threshold.
    let point = (0.5, 0.5, 0.0);
    let origin = (0.0, 0.0, 0.0);
    let normal = (0.0, 0.0, 1.0); // XY plane
    assert!(adaptive_abs_plane_distance_lt(point, origin, normal, 1e-6));
}

#[test]
fn adaptive_abs_plane_distance_lt_far_from_plane() {
    let point = (0.0, 0.0, 5.0);
    let origin = (0.0, 0.0, 0.0);
    let normal = (0.0, 0.0, 1.0);
    // Signed distance = 5 > 0.1 threshold.
    assert!(!adaptive_abs_plane_distance_lt(point, origin, normal, 0.1));
}

#[test]
fn adaptive_abs_plane_distance_lt_both_signs_respected() {
    let origin = (0.0, 0.0, 0.0);
    let normal = (0.0, 0.0, 1.0);
    // Point below the plane: signed distance is negative, magnitude is 2.
    assert!(!adaptive_abs_plane_distance_lt(
        (0.0, 0.0, -2.0),
        origin,
        normal,
        1.0
    ));
}

#[test]
fn adaptive_abs_plane_distance_lt_at_threshold_boundary() {
    // Point exactly at threshold distance — the strict < makes
    // this false, and the adaptive recomputation must agree.
    let point = (0.0, 0.0, 1.0);
    let origin = (0.0, 0.0, 0.0);
    let normal = (0.0, 0.0, 1.0);
    assert!(!adaptive_abs_plane_distance_lt(point, origin, normal, 1.0));
}

// ── DP slow-path coverage ──────────────────────────────────
//
// The next group of tests crafts inputs that fall INSIDE the f64
// ambiguity margin so the double-precision recomputation branch
// is forced to fire. Without these, the slow path is dead code
// from the test suite's perspective and a regression that breaks
// DP escalation could go unnoticed.

#[test]
fn adaptive_dist_sq_lt_slow_path_at_exact_equality_returns_false() {
    // dist_sq exactly equals threshold_sq → fast path's two
    // branches both fail (`+ margin < threshold` is false because
    // the margin is positive; `> threshold + margin` is also
    // false). Falls through to the DP recomputation, which
    // computes 0.0 - 0.0 = 0.0 and returns `0.0 < 0.0 = false`.
    let p1 = (0.0, 0.0, 0.0);
    let p2 = (1.0, 0.0, 0.0);
    let threshold_sq = 1.0;
    assert!(!adaptive_dist_sq_lt(p1, p2, threshold_sq));
}

#[test]
fn adaptive_dist_sq_lt_slow_path_just_below_returns_true() {
    // Construct a case where dist_sq is just below threshold_sq
    // by an amount inside the f64 margin. The fast path can't
    // decide; the DP slow path computes the difference exactly.
    //
    // dist_sq = 1.0, threshold_sq = 1.0 + ulp.
    let p1 = (0.0, 0.0, 0.0);
    let p2 = (1.0, 0.0, 0.0);
    // f64::EPSILON is the spacing at 1.0; threshold_sq = 1 + ε
    // is one ulp above 1.0.
    let threshold_sq = 1.0 + f64::EPSILON;
    assert!(
        adaptive_dist_sq_lt(p1, p2, threshold_sq),
        "dist² = 1.0 vs threshold² = 1.0 + ε should report below"
    );
}

#[test]
fn adaptive_dist_sq_lt_slow_path_just_above_returns_false() {
    // Mirror case: threshold_sq = 1.0 - ulp, dist_sq = 1.0.
    let p1 = (0.0, 0.0, 0.0);
    let p2 = (1.0, 0.0, 0.0);
    let threshold_sq = 1.0 - f64::EPSILON;
    assert!(
        !adaptive_dist_sq_lt(p1, p2, threshold_sq),
        "dist² = 1.0 vs threshold² = 1.0 - ε should report not-below"
    );
}

#[test]
fn adaptive_cross_ratio_sq_lt_slow_path_exact_perpendicular() {
    // Vectors at exactly 90° produce cross_len_sq = 1 and
    // norms_sq = 1, with threshold = 1 → rhs = 1. The fast path
    // can't separate them; the DP path computes 0.0 - 0.0 = 0.0
    // and the strict `< 0.0` returns false.
    let n1 = (1.0, 0.0, 0.0);
    let n2 = (0.0, 1.0, 0.0);
    assert!(!adaptive_cross_ratio_sq_lt(n1, n2, 1.0));
}

#[test]
fn adaptive_cross_ratio_sq_lt_slow_path_just_above_threshold() {
    // sin(89.99°) ≈ 0.99998; with threshold = 0.99 we expect
    // sin² > threshold² → cross_len_sq > rhs → return false.
    // The angle is far enough from the threshold that DP
    // recomputation just confirms the f64 answer, but the test
    // pins the behaviour.
    let n1 = (1.0, 0.0, 0.0);
    let n2 = (
        (89.99_f64).to_radians().cos(),
        (89.99_f64).to_radians().sin(),
        0.0,
    );
    assert!(!adaptive_cross_ratio_sq_lt(n1, n2, 0.99));
}

#[test]
fn adaptive_abs_plane_distance_lt_slow_path_exact_threshold() {
    // Point exactly at distance 1.0 from a plane through the
    // origin with normal +Z, threshold = 1.0. Falls into the
    // ambiguity zone; DP computes |1.0| < 1.0 = false.
    let point = (0.0, 0.0, 1.0);
    let origin = (0.0, 0.0, 0.0);
    let normal = (0.0, 0.0, 1.0);
    assert!(!adaptive_abs_plane_distance_lt(point, origin, normal, 1.0));
}

#[test]
fn adaptive_abs_plane_distance_lt_slow_path_just_below() {
    // Threshold one ulp above the point's distance — should
    // return true. The DP recomputation is required because the
    // f64 fast path's margin can't distinguish the two values.
    let point = (0.0, 0.0, 1.0);
    let origin = (0.0, 0.0, 0.0);
    let normal = (0.0, 0.0, 1.0);
    let threshold = 1.0 + f64::EPSILON;
    assert!(adaptive_abs_plane_distance_lt(
        point, origin, normal, threshold
    ));
}

// ── Microbenchmarks (run with `--ignored --nocapture`) ────────────

#[test]
fn adaptive_dist_sq_lt_at_extreme_scale_collapses_to_origin() {
    // Documents the fundamental limit: at coordinate scale 1e10
    // the f64 ulp is ~1e-6, so any tolerance below that resolves
    // points that *appear* close to be exactly identical. The
    // adaptive predicate respects this — DP cannot recover
    // precision that f64 lost at INPUT time.
    let p1 = (1e10, 0.0, 0.0);
    let p2 = (1e10 + 1e-10, 0.0, 0.0); // collapses to p1 in f64
    let threshold_sq = 1e-20;
    // The points are literally the same in f64 → distance is 0
    // → 0 < threshold² → returns true. This is the CORRECT
    // answer for inputs that f64 cannot distinguish.
    assert!(adaptive_dist_sq_lt(p1, p2, threshold_sq));
}
