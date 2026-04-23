//! Phase IV PnS — Bi-3 twist solver (eq 7) for `n ∈ {3, 5}`.
//!
//! In the PnS n ∈ {3, 5} bi-3 cap construction (§4.2 of
//! `docs/papers/14bi234.pdf`), the G¹ row-0 is only partially pinned
//! by the local eq (2³) — that gives `p_{20}` alone. The remaining
//! G¹ rows need the pair of twist CPs `(p_{11}^k, p_{11}^{k+1})` for
//! each shared boundary to satisfy the **vertex-enclosure** identity
//! at the central point. That identity is eq (7):
//!
//! ```text
//!   (3 − c)·[3c·(p_{11}^k + p_{11}^{k+1})
//!            − (3 − 2c)·((p_{01}^k + p_{01}^{k+1}) − 2·(1 − c)·p_{00})]
//!     = −2c²·[c·(p_{31}^k + p_{31}^{k+1}) − 3·(p_{21}^k + p_{21}^{k+1})]
//! ```
//!
//! where `c := cos(2π/n)` and `p_{00}` is the shared central point.
//! This gives `n` equations (one per shared boundary) in `n`
//! unknowns (one `p_{11}^k` per sector), with the structure
//! `x_k + x_{k+1} = rhs_k` — a **circulant pair-sum system**.
//!
//! For odd `n` (which is all that the PnS bi-3 path covers,
//! n ∈ {3, 5}), this system is of full rank and has the
//! closed-form solution
//! ```text
//!   x_k = (1/2) · Σ_{j=0}^{n-1} (−1)^j · rhs_{(k+j) mod n}
//! ```
//! derived by alternating-sum elimination on the pair-sum rows.
//! No iterative solve needed.
//!
//! The paper notes "When n is odd, the system is of full rank" —
//! this module only supports `n ∈ {3, 5}` accordingly.

#![allow(dead_code)]

use crate::geometry::Point3;

/// Solve eq (7) of `14bi234.pdf` for the N twist CPs
/// `p_{11}^k, k = 0..n` of a bi-3 n ∈ {3, 5} cap.
///
/// # Arguments
///
/// - `p00`: shared central point (`p_{00}` across all sectors).
/// - `p01`: `p_{01}^k` for `k ∈ 0..n` — first-interior CP on each
///   sector's lower edge.
/// - `p21`: `p_{21}^k` for `k ∈ 0..n`.
/// - `p31`: `p_{31}^k` for `k ∈ 0..n`.
/// - `n`: valence (must be 3 or 5).
///
/// # Returns
///
/// `Vec<Point3>` of length `n`, where entry `k` is the solved
/// `p_{11}^k`. The indexing convention is that boundary `k` is
/// shared between sectors `k` and `(k+1) mod n`.
///
/// # Panics
///
/// - `n ∉ {3, 5}` (the paper gives eq 7 only for odd n ∈ {3, 5};
///   even n is rank-deficient and not covered).
/// - `p01.len() != n`, `p21.len() != n`, or `p31.len() != n`.
pub fn solve_bi3_twists(
    p00: Point3,
    p01: &[Point3],
    p21: &[Point3],
    p31: &[Point3],
    n: usize,
) -> Vec<Point3> {
    assert!(
        n == 3 || n == 5,
        "PnS bi-3 twist solver: only n ∈ {{3, 5}} supported, got {n}"
    );
    assert_eq!(p01.len(), n, "p01 length mismatch");
    assert_eq!(p21.len(), n, "p21 length mismatch");
    assert_eq!(p31.len(), n, "p31 length mismatch");

    let c = (2.0 * std::f64::consts::PI / (n as f64)).cos();

    // Per-boundary RHS: eq (7) rearranged as
    //   x_k + x_{k+1} = rhs_k
    // where
    //   rhs_k = [(3 − 2c)·(p_{01}^k + p_{01}^{k+1} − 2(1-c)·p_{00})
    //           − (2c² / (3 − c))·(c·(p_{31}^k + p_{31}^{k+1})
    //                               − 3·(p_{21}^k + p_{21}^{k+1}))] / (3c)
    //
    // Derivation: eq (7) is
    //   (3-c)·[3c·(x_k+x_{k+1}) − (3-2c)·A_k] = −2c²·B_k
    // where A_k = (p_{01}^k + p_{01}^{k+1}) − 2(1-c)·p_{00}
    //       B_k = c·(p_{31}^k + p_{31}^{k+1}) − 3·(p_{21}^k + p_{21}^{k+1})
    // Solving for x_k + x_{k+1}:
    //   3c·(x_k+x_{k+1}) = (3-2c)·A_k − 2c²·B_k/(3-c)
    //   x_k+x_{k+1} = [(3-2c)·A_k − 2c²·B_k/(3-c)] / (3c)
    let three_minus_2c = 3.0 - 2.0 * c;
    let two_c2_over_3mc = 2.0 * c * c / (3.0 - c);
    let inv_3c = 1.0 / (3.0 * c);
    let two_one_minus_c_p00 = scale_pt(2.0 * (1.0 - c), p00);

    let mut rhs: Vec<Point3> = Vec::with_capacity(n);
    for k in 0..n {
        let kp1 = (k + 1) % n;
        // A_k = (p01[k] + p01[kp1]) − 2(1-c)·p00
        let a_k = sub_pt(add_pt(p01[k], p01[kp1]), two_one_minus_c_p00);
        // B_k = c·(p31[k] + p31[kp1]) − 3·(p21[k] + p21[kp1])
        let b_k = sub_pt(
            scale_pt(c, add_pt(p31[k], p31[kp1])),
            scale_pt(3.0, add_pt(p21[k], p21[kp1])),
        );
        // Combined RHS:
        let r = scale_pt(
            inv_3c,
            sub_pt(
                scale_pt(three_minus_2c, a_k),
                scale_pt(two_c2_over_3mc, b_k),
            ),
        );
        rhs.push(r);
    }

    // Solve the pair-sum circulant system `x_k + x_{k+1} = rhs_k`
    // for odd n via the closed-form alternating sum:
    //   x_k = (1/2) · Σ_{j=0}^{n-1} (−1)^j · rhs_{(k+j) mod n}
    let mut out: Vec<Point3> = Vec::with_capacity(n);
    for k in 0..n {
        let mut acc = Point3::new(0.0, 0.0, 0.0);
        let mut sign = 1.0;
        for j in 0..n {
            let idx = (k + j) % n;
            let r = rhs[idx];
            acc.x += sign * r.x;
            acc.y += sign * r.y;
            acc.z += sign * r.z;
            sign = -sign;
        }
        out.push(Point3::new(0.5 * acc.x, 0.5 * acc.y, 0.5 * acc.z));
    }
    out
}

#[inline]
fn add_pt(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}

#[inline]
fn sub_pt(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

#[inline]
fn scale_pt(s: f64, p: Point3) -> Point3 {
    Point3::new(s * p.x, s * p.y, s * p.z)
}

#[cfg(test)]
#[path = "pns_bi3_twists_tests.rs"]
mod tests;
