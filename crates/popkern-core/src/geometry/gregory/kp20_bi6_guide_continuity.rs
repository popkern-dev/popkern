//! KP18 total-degree-6 guide sector continuity formulas.
//!
//! KP20 reuses KP18's degree-6 guide. These helpers encode KP18
//! equations (1)-(5) for adjacent triangular guide sectors in the
//! same side-by-side KP20 namespace as the cap assembler.
//!
//! # Convention mismatch with `solve_kp20_guide_lsq`
//!
//! The helpers below address BB CPs at `flat_index(6, 6 − i, j)`
//! (i.e., the `(i_bary, j_bary)` pair with `k_bary` implicit).
//! `solve_kp20_guide_lsq::constraint_rows` encodes KP18 (1)–(5)
//! using `flat_index(6, j_bary, k_bary)` — a DIFFERENT address.
//! The two encodings constrain different CPs.
//!
//! Audit (2026-06-23): no production caller invokes these
//! helpers. The LSQ encoding is geometrically correct (it
//! identifies CPs along the shared spoke between adjacent
//! sectors); these helpers do not. They remain in the codebase
//! only because their unit tests pin the formulas they happen
//! to encode. Treat them as INFORMATIONAL: do not call them
//! from production paths. See
//! `docs/decisions/2026-06-23-kp20-bi6-paper-grade-guide-from-kp18.md`
//! Step 4 of the investigation plan.

use super::bernstein_tri6::flat_index;
use super::kp20_bi6_guide_sample::Kp20Tri6GuideSector;
use crate::geometry::Point3;

/// KP18 equation (1) join weights for adjacent sector wedges.
pub fn kp18_guide_join_weights(n: usize) -> Result<(f64, f64, f64), String> {
    if n < 3 {
        return Err(format!(
            "KP18/KP20 guide join weights require n >= 3; got {n}"
        ));
    }
    let c = (core::f64::consts::TAU / n as f64).cos();
    Ok((-1.0, 2.0 * c, 2.0 * (1.0 - c)))
}

/// Apply KP18 equation (2), C0 edge identification.
///
/// **Convention warning**: this helper addresses CPs differently
/// than `solve_kp20_guide_lsq::constraint_rows`. The LSQ
/// constraint (which is the geometrically-correct one) writes
/// `r.cps[flat_index(6, i, 0)] = s.cps[flat_index(6, 0, i)]` —
/// this helper writes a DIFFERENT pair. Do not use it in
/// production; the LSQ path is authoritative.
#[deprecated(
    note = "indexing-convention mismatch with solve_kp20_guide_lsq; see module-level note"
)]
pub fn kp18_apply_c0_join_degree6(left: &Kp20Tri6GuideSector, right: &mut Kp20Tri6GuideSector) {
    for i in 0..=6 {
        set(right, 6 - i, i, 0, get(left, 6 - i, 0, i));
    }
}

/// Apply KP18 equation (3), C1 cross-edge row propagation.
///
/// See module-level note: the LSQ encoding is authoritative.
#[deprecated(
    note = "indexing-convention mismatch with solve_kp20_guide_lsq; see module-level note"
)]
pub fn kp18_apply_c1_join_degree6(
    left: &Kp20Tri6GuideSector,
    right: &mut Kp20Tri6GuideSector,
    n: usize,
) -> Result<(), String> {
    let (w0, w1, w2) = kp18_guide_join_weights(n)?;
    for i in 1..=6 {
        let value = point_sum(&[
            (w0, get(left, 6 - i, 1, i - 1)),
            (w1, get(left, 6 - i, 0, i)),
            (w2, get(left, 6 - i + 1, 0, i - 1)),
        ]);
        set(right, 6 - i, i - 1, 1, value);
    }
    Ok(())
}

/// Apply KP18 equation (4), C2 cross-edge row propagation.
///
/// See module-level note: the LSQ encoding is authoritative.
#[deprecated(
    note = "indexing-convention mismatch with solve_kp20_guide_lsq; see module-level note"
)]
pub fn kp18_apply_c2_join_degree6(
    left: &Kp20Tri6GuideSector,
    right: &mut Kp20Tri6GuideSector,
    n: usize,
) -> Result<(), String> {
    let (w0, w1, w2) = kp18_guide_join_weights(n)?;
    for i in 2..=6 {
        let value = point_sum(&[
            (w0 * w0, get(left, 6 - i, 2, i - 2)),
            (2.0 * w0 * w1, get(left, 6 - i, 1, i - 1)),
            (w1 * w1, get(left, 6 - i, 0, i)),
            (2.0 * w0 * w2, get(left, 6 - i + 1, 1, i - 2)),
            (2.0 * w1 * w2, get(left, 6 - i + 1, 0, i - 1)),
            (w2 * w2, get(left, 6 - i + 2, 0, i - 2)),
        ]);
        set(right, 6 - i, i - 2, 2, value);
    }
    Ok(())
}

/// Apply KP18 equation (5), C3 cross-edge row propagation.
///
/// See module-level note: the LSQ encoding is authoritative.
#[deprecated(
    note = "indexing-convention mismatch with solve_kp20_guide_lsq; see module-level note"
)]
pub fn kp18_apply_c3_join_degree6(
    left: &Kp20Tri6GuideSector,
    right: &mut Kp20Tri6GuideSector,
    n: usize,
) -> Result<(), String> {
    let (w0, w1, w2) = kp18_guide_join_weights(n)?;
    for i in 3..=6 {
        let value = point_sum(&[
            (w0.powi(3), get(left, 6 - i, 3, i - 3)),
            (3.0 * w0 * w0 * w1, get(left, 6 - i, 2, i - 2)),
            (3.0 * w0 * w1 * w1, get(left, 6 - i, 1, i - 1)),
            (w1.powi(3), get(left, 6 - i, 0, i)),
            (3.0 * w0 * w0 * w2, get(left, 6 - i + 1, 2, i - 3)),
            (6.0 * w0 * w1 * w2, get(left, 6 - i + 1, 1, i - 2)),
            (3.0 * w1 * w1 * w2, get(left, 6 - i + 1, 0, i - 1)),
            (3.0 * w0 * w2 * w2, get(left, 6 - i + 2, 1, i - 3)),
            (3.0 * w1 * w2 * w2, get(left, 6 - i + 2, 0, i - 2)),
            (w2.powi(3), get(left, 6 - i + 3, 0, i - 3)),
        ]);
        set(right, 6 - i, i - 3, 3, value);
    }
    Ok(())
}

/// Apply KP18 equations (2)-(5) through third order.
///
/// See module-level note: the LSQ encoding is authoritative.
#[deprecated(
    note = "indexing-convention mismatch with solve_kp20_guide_lsq; see module-level note"
)]
#[allow(deprecated)]
pub fn kp18_apply_cubic_join_degree6(
    left: &Kp20Tri6GuideSector,
    right: &mut Kp20Tri6GuideSector,
    n: usize,
) -> Result<(), String> {
    kp18_apply_c0_join_degree6(left, right);
    kp18_apply_c1_join_degree6(left, right, n)?;
    kp18_apply_c2_join_degree6(left, right, n)?;
    kp18_apply_c3_join_degree6(left, right, n)
}

fn get(sector: &Kp20Tri6GuideSector, i: usize, j: usize, k: usize) -> Point3 {
    debug_assert_eq!(i + j + k, 6);
    sector.cps[flat_index(6, i, j)]
}

fn set(sector: &mut Kp20Tri6GuideSector, i: usize, j: usize, k: usize, value: Point3) {
    debug_assert_eq!(i + j + k, 6);
    sector.cps[flat_index(6, i, j)] = value;
}

fn point_sum(terms: &[(f64, Point3)]) -> Point3 {
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for (weight, point) in terms {
        out.x += weight * point.x;
        out.y += weight * point.y;
        out.z += weight * point.z;
    }
    out
}

#[cfg(test)]
#[path = "kp20_bi6_guide_continuity_tests.rs"]
mod tests;
