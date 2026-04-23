//! Phase IV PnS — Central-coefficient formulas for the n-gon cap.
//!
//! Two closed-form rules sit at the centre of every sector of a
//! PnS n-gon cap:
//!
//! 1. **Eq (4) — `p_{22}` symmetric average.** Pins down the
//!    bi-4 cap's centre coefficient from its own cross-pattern
//!    neighbourhood. Based on the degree-reduction identity
//!    `(2/3)(q_1 + q_3) − (1/6)(q_0 + q_4)` that projects a
//!    degree-4 curve's middle CP to the corresponding degree-3
//!    curve's middle — applied symmetrically in u and v and
//!    averaged. Paper: `14bi234.pdf` eq (4).
//!
//! 2. **Eq (5) — Central point `eop`.** The shared `p^k_{00}`
//!    across all N sectors is a two-weight convex combination of
//!    the DS-net inputs:
//!    ```text
//!      eop = e · avg_k(d^k_1) + (1 − e) · avg_k((d^k_2 + d^k_4) / 2)
//!    ```
//!    The weight `e = e^def` is valence-dependent, set by
//!    experiment against the paper's obstacle course. Variants:
//!    - **Bi-4** (§4): `e^def = (1.5 + ε)/(n + 1.6)` with
//!      `ε = 0` for `n ∈ {5, 6, 7}`, `ε = n/100` for
//!      `n ∈ {8, 9, 10}`.
//!    - **Bi-3 n ∈ {3, 5}** (§4.2): `e^def := 0.22`.
//!    - **2×2 bi-3 macro-patch** (Appendix C):
//!      `e^def := (0.02n + 1.5)/(n + 1.75)`.

#![allow(dead_code)]

use crate::geometry::Point3;

/// Compute `p_{22}` of a bi-4 sector by eq (4)'s symmetric
/// degree-reduction average.
///
/// Inputs (all from the sector's own 5×5 BB-grid):
/// - `p02, p12, p32, p42`: column 2, rows 0, 1, 3, 4.
/// - `p20, p21, p23, p24`: row 2, columns 0, 1, 3, 4.
///
/// ```text
///   p_{22} := (1/2)·((2/3)(p_{12} + p_{32}) − (1/6)(p_{02} + p_{42}))
///           + (1/2)·((2/3)(p_{21} + p_{23}) − (1/6)(p_{20} + p_{24}))
/// ```
#[allow(clippy::too_many_arguments)]
pub fn p22_symmetric_average(
    p02: Point3,
    p12: Point3,
    p32: Point3,
    p42: Point3,
    p20: Point3,
    p21: Point3,
    p23: Point3,
    p24: Point3,
) -> Point3 {
    // Horizontal reduction of column 2:
    //   h = (2/3)(p12 + p32) − (1/6)(p02 + p42)
    let h = Point3::new(
        (2.0 / 3.0) * (p12.x + p32.x) - (1.0 / 6.0) * (p02.x + p42.x),
        (2.0 / 3.0) * (p12.y + p32.y) - (1.0 / 6.0) * (p02.y + p42.y),
        (2.0 / 3.0) * (p12.z + p32.z) - (1.0 / 6.0) * (p02.z + p42.z),
    );
    // Vertical reduction of row 2:
    //   v = (2/3)(p21 + p23) − (1/6)(p20 + p24)
    let v = Point3::new(
        (2.0 / 3.0) * (p21.x + p23.x) - (1.0 / 6.0) * (p20.x + p24.x),
        (2.0 / 3.0) * (p21.y + p23.y) - (1.0 / 6.0) * (p20.y + p24.y),
        (2.0 / 3.0) * (p21.z + p23.z) - (1.0 / 6.0) * (p20.z + p24.z),
    );
    // Average the two reductions.
    Point3::new(0.5 * (h.x + v.x), 0.5 * (h.y + v.y), 0.5 * (h.z + v.z))
}

/// Compute the central point `eop = p^k_{00}` via eq (5) as a
/// weighted average of the DS-net inputs.
///
/// `d1`, `d2_plus_d4_half` are the per-sector inputs:
/// - `d1[k] = d^k_1`: layer-1 DS-net point for sector k.
/// - `d2_plus_d4_half[k] = (d^k_2 + d^k_4) / 2`: averaged pair for
///   sector k.
///
/// `e` is the weight. Pass one of the paper defaults from
/// [`e_def_bi4`], [`E_DEF_BI3_N3_N5`], or [`e_def_bi3_2x2`].
///
/// ```text
///   eop = e · avg_k(d^k_1) + (1 − e) · avg_k((d^k_2 + d^k_4) / 2)
/// ```
pub fn central_point(d1: &[Point3], d2_plus_d4_half: &[Point3], e: f64) -> Point3 {
    assert_eq!(
        d1.len(),
        d2_plus_d4_half.len(),
        "central_point: d1 and d2_plus_d4_half length mismatch"
    );
    assert!(!d1.is_empty(), "central_point: need N ≥ 1 sectors");
    let n = d1.len() as f64;
    let mut avg1 = Point3::new(0.0, 0.0, 0.0);
    let mut avg24 = Point3::new(0.0, 0.0, 0.0);
    for i in 0..d1.len() {
        avg1.x += d1[i].x;
        avg1.y += d1[i].y;
        avg1.z += d1[i].z;
        avg24.x += d2_plus_d4_half[i].x;
        avg24.y += d2_plus_d4_half[i].y;
        avg24.z += d2_plus_d4_half[i].z;
    }
    avg1 = Point3::new(avg1.x / n, avg1.y / n, avg1.z / n);
    avg24 = Point3::new(avg24.x / n, avg24.y / n, avg24.z / n);
    Point3::new(
        e * avg1.x + (1.0 - e) * avg24.x,
        e * avg1.y + (1.0 - e) * avg24.y,
        e * avg1.z + (1.0 - e) * avg24.z,
    )
}

/// Paper-default weight `e^def` for the **bi-4 cap** construction
/// (§4, eq 5 context):
/// ```text
///   e^def := (1.5 + ε) / (n + 1.6)
///   ε := 0          for n ∈ {5, 6, 7}
///   ε := n / 100    for n ∈ {8, 9, 10}
/// ```
/// Outside `n ∈ [5, 10]`, we use `ε = n/100` — the authors' tested
/// range is 5-10 but the formula extends monotonically.
pub fn e_def_bi4(n: usize) -> f64 {
    let epsilon = if (5..=7).contains(&n) {
        0.0
    } else {
        (n as f64) / 100.0
    };
    (1.5 + epsilon) / (n as f64 + 1.6)
}

/// Paper default `e^def := 0.22` for the bi-3 n ∈ {3, 5}
/// construction (§4.2).
pub const E_DEF_BI3_N3_N5: f64 = 0.22;

/// Paper default for the 2×2 bi-3 macro-patch cap (Appendix C):
/// ```text
///   e^def := (0.02n + 1.5) / (n + 1.75)
/// ```
pub fn e_def_bi3_2x2(n: usize) -> f64 {
    (0.02 * (n as f64) + 1.5) / (n as f64 + 1.75)
}

#[cfg(test)]
#[path = "pns_center_tests.rs"]
mod tests;
