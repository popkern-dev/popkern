//! Karčiauskas-Peters 2009 C¹/C² inter-sector
//! connectivity for `g△`.
//!
//! Paper: [KP09a] K. Karčiauskas & J. Peters (2009), "Assembling
//! curvature continuous surfaces from triangular patches",
//! *Computers and Graphics* 33(3):204-210. `docs/papers/09tri8c2_smi.pdf`.
//! Math summary: `docs/journal/polar-kp-math.md` §3.2a.
//!
//! # Two-index paper labelling (KP09a §2.1)
//!
//! KP09a uses 2-index labels `(x, y)` with the implicit third
//! index `d − x − y` (here `d = 5` for tri-5). The star `o` is at
//! `(0, 0)`, `v_i` at `(d, 0)`, `v_{i+1}` at `(0, d)`.
//!
//! This maps **directly** onto our `flat_index(5, i, j)` indexing
//! (paper `x` is our `i`, paper `y` is our `j`). No permutation
//! needed — so these helpers work naturally in our basis.
//!
//! # Shared edge
//!
//! The boundary between sector `i` and sector `i+1` is the ray
//! `o → v_{i+1}` (in sector i's coords), equivalently `o → v_i`
//! (in sector i+1's coords). The paper identity across this edge
//! [KP09a §2.1] is
//!
//! ```text
//!   p^{i+1}_{j, 0} = p^i_{0, j}     for j = 0..d.
//! ```
//!
//! Both sides run from the star (`j = 0`) to the shared vertex
//! (`j = d = 5`).
//!
//! # C¹ and C² join equations
//!
//! [KP09a eq (3), C¹]:
//!
//! ```text
//!   p^{i+1}_{j, 1} = ν₀·p^i_{0, j} + ν₁·p^i_{0, j+1} + ν₂·p^i_{1, j}
//! ```
//!
//! [KP09a eq (4), C²]:
//!
//! ```text
//!   p^{i+1}_{j, 2} = ν₀²·p^i_{0,j}     + 2ν₀ν₁·p^i_{0,j+1} + ν₁²·p^i_{0,j+2}
//!                  + 2ν₀ν₂·p^i_{1,j}   + 2ν₁ν₂·p^i_{1,j+1} + ν₂²·p^i_{2,j}
//! ```
//!
//! [KP09a eq (5), weights]: with `v_{i+1} = ν₀·o + ν₁·v_i + ν₂·v_{i-1}`,
//!
//! ```text
//!   ν₀ := 2c̄     where c̄ := 1 − c
//!   ν₁ := 2c     where c := cos(2π/n)
//!   ν₂ := −1
//! ```
//!
//! The three weights sum to `1`, so the equations preserve affine
//! combinations — in particular a constant cap stays constant
//! under either join.

#![allow(dead_code)]

use super::bernstein_tri5::flat_index;
use super::kp13_bi5_gtri5::Gtri5Sector;
use crate::geometry::Point3;
use std::f64::consts::PI;

/// KP09a eq (5) weights `(ν₀, ν₁, ν₂)` for valence `n`.
/// Panics if `n < 3`.
#[derive(Debug, Clone, Copy)]
pub struct ConnectWeights {
    pub nu0: f64,
    pub nu1: f64,
    pub nu2: f64,
}

impl ConnectWeights {
    /// Default weights from [KP09a eq (5)].
    pub fn for_valence(n: usize) -> Self {
        assert!(n >= 3, "KP09a connectivity needs n ≥ 3; got {n}");
        let c = (2.0 * PI / n as f64).cos();
        let c_bar = 1.0 - c;
        Self {
            nu0: 2.0 * c_bar,
            nu1: 2.0 * c,
            nu2: -1.0,
        }
    }

    /// Sum of the three weights — should always be `1` (affine
    /// combination). Useful as a sanity check in property tests.
    #[inline]
    pub fn sum(&self) -> f64 {
        self.nu0 + self.nu1 + self.nu2
    }
}

/// Paper-index → flat-index for tri-5.
/// Paper `(x, y)` with implicit w-exp = `5 − x − y` maps directly
/// onto our `flat_index(5, x, y)`.
#[inline]
pub fn paper_flat(x: usize, y: usize) -> usize {
    flat_index(5, x, y)
}

/// Compute `p^{i+1}_{j, 1}` via KP09a eq (3).
///
/// `sector_i` provides the three CPs `p^i_{0, j}`, `p^i_{0, j+1}`,
/// `p^i_{1, j}`; valid for `j ∈ 0..=4` (so that `j + 1 ≤ 5`).
///
/// Panics if `j > 4`.
pub fn c1_join(weights: &ConnectWeights, sector_i: &Gtri5Sector, j: usize) -> Point3 {
    assert!(j <= 4, "c1_join requires j ≤ 4; got j = {j}");
    let p0j = sector_i.cps[paper_flat(0, j)];
    let p0jp1 = sector_i.cps[paper_flat(0, j + 1)];
    let p1j = sector_i.cps[paper_flat(1, j)];
    linear3(weights.nu0, p0j, weights.nu1, p0jp1, weights.nu2, p1j)
}

/// Compute `p^{i+1}_{j, 2}` via KP09a eq (4).
///
/// `sector_i` provides the six CPs at paper positions `(0, j)`,
/// `(0, j+1)`, `(0, j+2)`, `(1, j)`, `(1, j+1)`, `(2, j)`. Valid
/// for `j ∈ 0..=3` (so that `j + 2 ≤ 5`).
///
/// Panics if `j > 3`.
pub fn c2_join(weights: &ConnectWeights, sector_i: &Gtri5Sector, j: usize) -> Point3 {
    assert!(j <= 3, "c2_join requires j ≤ 3; got j = {j}");
    let nu0 = weights.nu0;
    let nu1 = weights.nu1;
    let nu2 = weights.nu2;

    let p0j = sector_i.cps[paper_flat(0, j)];
    let p0jp1 = sector_i.cps[paper_flat(0, j + 1)];
    let p0jp2 = sector_i.cps[paper_flat(0, j + 2)];
    let p1j = sector_i.cps[paper_flat(1, j)];
    let p1jp1 = sector_i.cps[paper_flat(1, j + 1)];
    let p2j = sector_i.cps[paper_flat(2, j)];

    let pairs: [(f64, Point3); 6] = [
        (nu0 * nu0, p0j),
        (2.0 * nu0 * nu1, p0jp1),
        (nu1 * nu1, p0jp2),
        (2.0 * nu0 * nu2, p1j),
        (2.0 * nu1 * nu2, p1jp1),
        (nu2 * nu2, p2j),
    ];
    combine6(pairs)
}

#[inline]
fn linear3(a: f64, pa: Point3, b: f64, pb: Point3, c: f64, pc: Point3) -> Point3 {
    Point3::new(
        a * pa.x + b * pb.x + c * pc.x,
        a * pa.y + b * pb.y + c * pc.y,
        a * pa.z + b * pb.z + c * pc.z,
    )
}

#[inline]
fn combine6(pairs: [(f64, Point3); 6]) -> Point3 {
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for (w, p) in pairs {
        out.x += w * p.x;
        out.y += w * p.y;
        out.z += w * p.z;
    }
    out
}

#[cfg(test)]
#[path = "kp13_bi5_connect_tests.rs"]
mod tests;
