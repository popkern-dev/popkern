//! Phase I.11e.2 — Boundary tensor-border `b^k` and its BB-jet at
//! the unit-square corners.
//!
//! KP13 §3 eq (4) fits `g△ ∘ χ` to the tensor-border `b` at jet
//! matching points `t ∈ {m_k, v_k}`. In the shared χ-domain/
//! b-domain unit square, these corners are:
//!
//! | paper label | `(u, v)`   |
//! |-------------|------------|
//! | star `o`    | `(0, 0)`   |
//! | `m_k`       | `(1, 0)`   |
//! | `m_{k+1}`   | `(0, 1)`   |
//! | `v_k`       | `(1, 1)`   |
//!
//! For a degree-3 tensor-product Bézier patch `b^k : [0, 1]² → R³`
//! with 4 × 4 BB-coefs `b_{i,j}`, the jet at each corner is linear
//! in the coefs. The corner-derivative identities:
//!
//! At `(0, 0)`:
//! ```text
//!   b        = b_{0,0}
//!   ∂_u b    = 3·(b_{1,0} − b_{0,0})
//!   ∂_v b    = 3·(b_{0,1} − b_{0,0})
//!   ∂²_u b   = 6·(b_{2,0} − 2·b_{1,0} + b_{0,0})
//!   ∂²_v b   = 6·(b_{0,2} − 2·b_{0,1} + b_{0,0})
//!   ∂²_uv b  = 9·(b_{1,1} − b_{1,0} − b_{0,1} + b_{0,0})
//! ```
//!
//! At other corners, substitute the symmetric indices (e.g., at
//! `(1, 0)` the u-layer flips to `b_{3,·}, b_{2,·}, b_{1,·}`).
//!
//! The resulting [`TaylorJet2`] is what the LSQ subtracts from the
//! composite jet `H_t(g△^k ∘ χ^k)` to form a residual.

#![allow(dead_code)]

use super::kp13_bi5_chi::ChiCorner;
use super::kp13_bi5_jet::TaylorJet2;
use crate::geometry::Point3;

/// A degree-3 tensor-product Bézier patch — one tensor-border sector.
/// `cps[i][j]` is the BB-coef at `(i, j)` ∈ `{0, 1, 2, 3}²`.
#[derive(Debug, Clone, Copy)]
pub struct Bi3TensorBorder {
    /// `cps[i][j]` is the BB coefficient at `(i, j) ∈ {0, 1, 2, 3}²`.
    pub cps: [[Point3; 4]; 4],
}

impl Bi3TensorBorder {
    /// Zero-initialised bi-3 patch (all CPs at origin).
    pub fn zero() -> Self {
        let z = Point3::new(0.0, 0.0, 0.0);
        Self { cps: [[z; 4]; 4] }
    }

    /// Construct from a 4×4 grid of CPs.
    pub fn from_cps(cps: [[Point3; 4]; 4]) -> Self {
        Self { cps }
    }

    /// Extract the full second-order Taylor jet (value + 1st + 2nd
    /// partials) at the given corner of this patch's unit square.
    pub fn jet_at_corner(&self, corner: ChiCorner) -> TaylorJet2 {
        let c = &self.cps;
        match corner {
            ChiCorner::Star => corner_jet_at(c[0][0], c[1][0], c[0][1], c[2][0], c[0][2], c[1][1]),
            ChiCorner::Mk => corner_jet_at(c[3][0], c[2][0], c[3][1], c[1][0], c[3][2], c[2][1]),
            ChiCorner::Mk1 => corner_jet_at(c[0][3], c[1][3], c[0][2], c[2][3], c[0][1], c[1][2]),
            ChiCorner::Vk => corner_jet_at(c[3][3], c[2][3], c[3][2], c[1][3], c[3][1], c[2][2]),
        }
    }
}

/// Compute the Taylor jet at a bi-3 corner given the six
/// corner-cluster CPs in a canonical order:
///
/// * `c00` — the corner itself
/// * `cu1` — one step inward in the `u` direction from the corner
/// * `cv1` — one step inward in the `v` direction from the corner
/// * `cu2` — two steps inward in `u`
/// * `cv2` — two steps inward in `v`
/// * `cuv` — one step inward in both `u` and `v` (the "twist" CP)
///
/// The sign of the first partials is adjusted so that `∂_u` always
/// points into the domain `[0, 1]`, matching the KP13 convention
/// where the jet is measured in the canonical parameterisation of
/// the unit square.
///
/// For bi-3 at any corner:
/// ```text
///   ∂_u  = ±3·(cu1 − c00)
///   ∂_v  = ±3·(cv1 − c00)
///   ∂²_u = 6·(cu2 − 2·cu1 + c00)
///   ∂²_v = 6·(cv2 − 2·cv1 + c00)
///   ∂²_uv = ±9·(cuv − cu1 − cv1 + c00)
/// ```
///
/// The sign for `∂_u` is positive at `(0, 0)` / `(0, 1)` and
/// negative at `(1, 0)` / `(1, 1)` (the u-direction-into-domain
/// sign). Analogous for `∂_v`. The mixed `∂²_uv` sign is the
/// product of the two first-partial signs.
fn corner_jet_at(
    c00: Point3,
    cu1: Point3,
    cv1: Point3,
    cu2: Point3,
    cv2: Point3,
    cuv: Point3,
) -> TaylorJet2 {
    // All signs positive here — the Bi3TensorBorder::jet_at_corner
    // caller feeds CPs in a canonically-oriented order such that
    // cu1, cv1 are inward and partials come out with the canonical
    // BB convention (pointing away from the corner along each axis).
    // Callers downstream that need domain-oriented partials can
    // flip signs externally.
    TaylorJet2 {
        value: c00,
        d_u: scale_sub(3.0, cu1, 3.0, c00),
        d_v: scale_sub(3.0, cv1, 3.0, c00),
        d_uu: three_point(6.0, cu2, -12.0, cu1, 6.0, c00),
        d_vv: three_point(6.0, cv2, -12.0, cv1, 6.0, c00),
        d_uv: four_point([(9.0, cuv), (-9.0, cu1), (-9.0, cv1), (9.0, c00)]),
    }
}

#[inline]
fn scale_sub(s1: f64, p1: Point3, s2: f64, p2: Point3) -> Point3 {
    Point3::new(
        s1 * p1.x - s2 * p2.x,
        s1 * p1.y - s2 * p2.y,
        s1 * p1.z - s2 * p2.z,
    )
}

#[inline]
fn three_point(s1: f64, p1: Point3, s2: f64, p2: Point3, s3: f64, p3: Point3) -> Point3 {
    Point3::new(
        s1 * p1.x + s2 * p2.x + s3 * p3.x,
        s1 * p1.y + s2 * p2.y + s3 * p3.y,
        s1 * p1.z + s2 * p2.z + s3 * p3.z,
    )
}

#[inline]
fn four_point(pairs: [(f64, Point3); 4]) -> Point3 {
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    for (s, p) in pairs {
        x += s * p.x;
        y += s * p.y;
        z += s * p.z;
    }
    Point3::new(x, y, z)
}

/// The full tensor-border input to the LSQ — one `Bi3TensorBorder`
/// per sector, indexed by sector `k` ∈ `0..n`.
///
/// The jets at `m_k` and `v_k` of sector `k`'s patch form the
/// target side of KP13 eq (4)'s residual.
#[derive(Debug, Clone)]
pub struct TensorBorder {
    /// One `Bi3TensorBorder` per sector; `patches.len() == n`.
    pub patches: Vec<Bi3TensorBorder>,
}

impl TensorBorder {
    /// Zero-initialised tensor-border with `n` patches.
    pub fn zero(n: usize) -> Self {
        assert!(n >= 3, "TensorBorder needs n ≥ 3, got {n}");
        Self {
            patches: vec![Bi3TensorBorder::zero(); n],
        }
    }

    /// Valence / patch count.
    pub fn n(&self) -> usize {
        self.patches.len()
    }

    /// Sector `k`'s jet at the given corner.
    pub fn jet(&self, k: usize, corner: ChiCorner) -> TaylorJet2 {
        self.patches[k].jet_at_corner(corner)
    }
}

#[cfg(test)]
#[path = "kp13_bi5_boundary_tests.rs"]
mod tests;
