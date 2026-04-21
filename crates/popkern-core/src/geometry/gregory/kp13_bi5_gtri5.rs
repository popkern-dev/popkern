//! Phase I.11 — Karčiauskas-Peters 2013 auxiliary map `g△`.
//!
//! `g△` is an `n`-piece C² piecewise polynomial of **total degree
//! 5** on an n-gon domain Ω. Each sector `g△^k` is a tri-5 patch
//! on the wedge `Δ_k = {(0,0), (c_k, s_k), (c_{k+1}, s_{k+1})}`
//! where `c_k := cos(k·2π/n)`, `s_k := sin(k·2π/n)`.
//!
//! Paper: [KP13] §3, Fig. 4 - 7.
//!
//! # Barycentric labelling
//!
//! Per-sector BB-coefficients use barycentric indices `(i, j, k)`
//! with `i + j + k = 5`. In our convention the three corners of
//! `Δ_k` map to:
//!
//! ```text
//!   (5, 0, 0) := star point `o = (0, 0)`
//!   (0, 5, 0) := `m_k    = (c_k,     s_k)`        [edge mid-point]
//!   (0, 0, 5) := `m_{k+1} = (c_{k+1}, s_{k+1})`   [shared with sector k+1]
//! ```
//!
//! So `g△^k_{050}` is the mid-edge vertex that the paper labels
//! `m_k`, and `g△^k_{005}` is the mid-edge vertex of the *next*
//! sector boundary (shared with `g△^{k+1}_{050}`).
//!
//! Free BB-coefficients per-sector from [KP13 Fig. 4] are marked
//! with disks / circles — in this module, see [`free_bb_indices`].
//!
//! # Star point and global DoF count
//!
//! The star-point coefficient `g△^k_{500}` is shared across every
//! sector (all sectors meet at the origin) and is pinned to the
//! Catmull-Clark limit of the cap-net (math summary §1.5). This
//! leaves the remaining free BB-coefs determined by the LSQ match
//! against tensor-border jets at `m_k` and `v_k` (eq 4).
//!
//! Paper reports the resulting linear system as `(6n + 5) × (6n + 5)`
//! [KP13 §5.1].
//!
//! # This file
//!
//! I.11a — data structures, star-point pinning, per-sector tri-5
//! evaluation. C² inter-sector fill-in (I.11b) and LSQ assembly
//! (I.11c-d) follow in later commits, which is why several
//! public items below are not yet consumed.

#![allow(dead_code)]

use crate::geometry::Point3;
use std::f64::consts::PI;

use super::bernstein_tri5::{
    TRI5_CP_COUNT, bernstein_tri5_all, bernstein_tri5_partial_u, bernstein_tri5_partial_v,
    bernstein_tri5_second_partial_uu, bernstein_tri5_second_partial_uv,
    bernstein_tri5_second_partial_vv, flat_index,
};

/// Flat index of the star-point BB-coef `(5, 0, 0)` ↔ `b(500)`.
pub const IDX_500: usize = 20;
/// Flat index of `m_k` BB-coef `(0, 5, 0)`.
pub const IDX_050: usize = 5;
/// Flat index of `m_{k+1}` BB-coef `(0, 0, 5)`.
pub const IDX_005: usize = 0;
/// Flat index of `(2, 3, 0)` — interior free coef per Fig. 4.
pub const IDX_230: usize = 14;
/// Flat index of `(1, 2, 2)` — interior free coef per Fig. 4.
pub const IDX_122: usize = 8;
/// Flat index of `(0, 3, 2)` — interior free coef per Fig. 4.
pub const IDX_032: usize = 3;
/// Flat index of `(0, 2, 3)` — interior free coef per Fig. 4.
pub const IDX_023: usize = 2;

/// Free BB-coef flat-indices per sector per [KP13 Fig. 4].
/// Seven positions: `{005, 023, 032, 050, 122, 230, 500}`.
///
/// Note `500` is shared across every sector (star point) and
/// `050` / `005` form pairs between adjacent sectors — these
/// sharings reduce the global free-parameter count from `7n`
/// to `5n + 1` (before the star-point pinning).
pub const fn free_bb_indices() -> [usize; 7] {
    [
        IDX_005, IDX_023, IDX_032, IDX_050, IDX_122, IDX_230, IDX_500,
    ]
}

/// Sector-corner barycentric coordinates in the ambient `(x, y)`
/// plane — `(c_k, s_k)` in paper notation.
#[inline]
pub fn sector_corner(n: usize, k: usize) -> (f64, f64) {
    let t = (k as f64) * 2.0 * PI / (n as f64);
    (t.cos(), t.sin())
}

/// `n`-piece auxiliary tri-5 map `g△ : Ω → R³`. One sector per
/// wedge, 21 BB-coefs per sector.
///
/// Invariants (maintained by constructors / setters, **not**
/// re-validated by the LSQ solver until it lands):
///
/// * `sectors.len() == n`.
/// * All sectors' `cps[IDX_500]` are equal — they share the
///   star-point.
/// * For each sector pair `(k, k+1 mod n)`: sector `k`'s
///   `cps[IDX_005]` equals sector `k+1`'s `cps[IDX_050]` (shared
///   mid-edge vertex).
///
/// The second invariant (mid-edge sharing) is **not** enforced at
/// construction — individual sector mutation can break it. Call
/// [`Gtri5::resync_mid_edges`] after free-form setters to restore.
#[derive(Debug, Clone)]
pub struct Gtri5 {
    /// Valence `n` — number of sectors.
    pub n: usize,
    /// Per-sector BB-coefs, `n` sectors × 21 coefs each, `Point3`-
    /// valued.
    pub sectors: Vec<Gtri5Sector>,
}

/// A single tri-5 sector of `g△`. `cps[flat_index(5, i, j)]` is
/// the BB-coefficient at barycentric position `(i, j, 5-i-j)`.
#[derive(Debug, Clone, Copy)]
pub struct Gtri5Sector {
    /// 21 BB-coefficients in the row-major `(i, j)` layout
    /// documented in [`super::bernstein_tri5`].
    pub cps: [Point3; TRI5_CP_COUNT],
}

impl Gtri5Sector {
    /// All-zero sector (placeholder — not a valid guide).
    pub fn zero() -> Self {
        Self {
            cps: [Point3::new(0.0, 0.0, 0.0); TRI5_CP_COUNT],
        }
    }

    /// Get the BB-coef at barycentric `(i, j, 5-i-j)`.
    /// Panics if `i + j > 5`.
    #[inline]
    pub fn cp(&self, i: usize, j: usize) -> Point3 {
        self.cps[flat_index(5, i, j)]
    }

    /// Set the BB-coef at barycentric `(i, j, 5-i-j)`.
    /// Panics if `i + j > 5`.
    #[inline]
    pub fn set_cp(&mut self, i: usize, j: usize, p: Point3) {
        self.cps[flat_index(5, i, j)] = p;
    }

    /// Evaluate this sector's tri-5 polynomial at barycentric
    /// `(u, v, w = 1 - u - v)`.
    pub fn evaluate(&self, u: f64, v: f64) -> Point3 {
        let basis = bernstein_tri5_all(u, v);
        combine(&basis, &self.cps)
    }

    /// First partial ∂/∂u at `(u, v)`.
    pub fn partial_u(&self, u: f64, v: f64) -> Point3 {
        let b = bernstein_tri5_partial_u(u, v);
        combine(&b, &self.cps)
    }

    /// First partial ∂/∂v at `(u, v)`.
    pub fn partial_v(&self, u: f64, v: f64) -> Point3 {
        let b = bernstein_tri5_partial_v(u, v);
        combine(&b, &self.cps)
    }

    /// Second partial ∂²/∂u² at `(u, v)`.
    pub fn second_partial_uu(&self, u: f64, v: f64) -> Point3 {
        let b = bernstein_tri5_second_partial_uu(u, v);
        combine(&b, &self.cps)
    }

    /// Second partial ∂²/∂v² at `(u, v)`.
    pub fn second_partial_vv(&self, u: f64, v: f64) -> Point3 {
        let b = bernstein_tri5_second_partial_vv(u, v);
        combine(&b, &self.cps)
    }

    /// Second mixed partial ∂²/∂u∂v at `(u, v)`.
    pub fn second_partial_uv(&self, u: f64, v: f64) -> Point3 {
        let b = bernstein_tri5_second_partial_uv(u, v);
        combine(&b, &self.cps)
    }
}

impl Gtri5 {
    /// New `g△` with `n` sectors, all BB-coefs at origin.
    /// Panics if `n < 3` (bi-5 construction needs valence ≥ 3).
    pub fn new(n: usize) -> Self {
        assert!(n >= 3, "KP13 bi-5 g△ needs valence ≥ 3; got {n}");
        let sectors = vec![Gtri5Sector::zero(); n];
        Self { n, sectors }
    }

    /// New `g△` with `n` sectors and the star-point pinned to
    /// `star` (propagated to every sector's `IDX_500` coef).
    pub fn new_with_star_point(n: usize, star: Point3) -> Self {
        let mut g = Self::new(n);
        g.set_star_point(star);
        g
    }

    /// Set the shared star-point coefficient across every sector.
    /// Paper: star-point = Catmull-Clark limit of the cap-net
    /// (math summary §1.5, [KP13 §3]).
    pub fn set_star_point(&mut self, star: Point3) {
        for s in &mut self.sectors {
            s.cps[IDX_500] = star;
        }
    }

    /// Read the current star-point coefficient. Panics if sectors
    /// disagree (invariant violation).
    pub fn star_point(&self) -> Point3 {
        let expected = self.sectors[0].cps[IDX_500];
        for (k, s) in self.sectors.iter().enumerate().skip(1) {
            let got = s.cps[IDX_500];
            assert!(
                points_equal(got, expected),
                "star-point mismatch at sector {k}"
            );
        }
        expected
    }

    /// Mutable access to sector `k`. Bypasses invariant
    /// maintenance — call [`Self::resync_mid_edges`] and
    /// [`Self::set_star_point`] afterwards if shared coefs were
    /// touched.
    #[inline]
    pub fn sector_mut(&mut self, k: usize) -> &mut Gtri5Sector {
        &mut self.sectors[k]
    }

    /// Immutable access to sector `k`.
    #[inline]
    pub fn sector(&self, k: usize) -> &Gtri5Sector {
        &self.sectors[k]
    }

    /// Restore the shared mid-edge invariant: for each sector
    /// pair `(k, k+1 mod n)`, set sector `k`'s `IDX_005` equal
    /// to sector `k+1`'s `IDX_050`.
    ///
    /// Uses sector `k+1`'s `IDX_050` as authoritative (the
    /// LSQ solver owns m_k's position on each sector edge).
    pub fn resync_mid_edges(&mut self) {
        for k in 0..self.n {
            let next = (k + 1) % self.n;
            let authoritative = self.sectors[next].cps[IDX_050];
            self.sectors[k].cps[IDX_005] = authoritative;
        }
    }

    /// Set the mid-edge vertex `m_k` (the barycentric `(0, 5, 0)`
    /// of sector `k`, equivalently `(0, 0, 5)` of sector `k-1`).
    /// Propagates to both sectors simultaneously, preserving the
    /// shared-mid-edge invariant.
    pub fn set_mid_edge(&mut self, k: usize, p: Point3) {
        let kk = k % self.n;
        let prev = (kk + self.n - 1) % self.n;
        self.sectors[kk].cps[IDX_050] = p;
        self.sectors[prev].cps[IDX_005] = p;
    }

    /// Read the mid-edge vertex `m_k`. Panics if sectors disagree
    /// (invariant violation — caller forgot [`Self::resync_mid_edges`]).
    pub fn mid_edge(&self, k: usize) -> Point3 {
        let kk = k % self.n;
        let prev = (kk + self.n - 1) % self.n;
        let left = self.sectors[kk].cps[IDX_050];
        let right = self.sectors[prev].cps[IDX_005];
        assert!(
            points_equal(left, right),
            "mid-edge m_k mismatch at k = {kk}: lhs = {left:?}, rhs = {right:?}"
        );
        left
    }

    /// Evaluate `g△` on sector `k` at barycentric `(u, v)`.
    #[inline]
    pub fn evaluate_sector(&self, k: usize, u: f64, v: f64) -> Point3 {
        self.sectors[k].evaluate(u, v)
    }

    /// Validate structural invariants. Returns `Err` with the
    /// first violation found.
    pub fn validate(&self) -> Result<(), String> {
        if self.n < 3 {
            return Err(format!("g△ valence n = {} < 3", self.n));
        }
        if self.sectors.len() != self.n {
            return Err(format!(
                "g△ sector count {} != valence n = {}",
                self.sectors.len(),
                self.n
            ));
        }
        let star = self.sectors[0].cps[IDX_500];
        for (k, s) in self.sectors.iter().enumerate().skip(1) {
            if !points_equal(s.cps[IDX_500], star) {
                return Err(format!("g△ star-point (IDX_500) mismatch at sector {k}"));
            }
        }
        for k in 0..self.n {
            let next = (k + 1) % self.n;
            let left = self.sectors[k].cps[IDX_005];
            let right = self.sectors[next].cps[IDX_050];
            if !points_equal(left, right) {
                return Err(format!(
                    "g△ mid-edge share mismatch: sector {k} IDX_005 vs sector {next} IDX_050"
                ));
            }
        }
        Ok(())
    }
}

#[inline]
fn combine(basis: &[f64; TRI5_CP_COUNT], cps: &[Point3; TRI5_CP_COUNT]) -> Point3 {
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..TRI5_CP_COUNT {
        let b = basis[i];
        out.x += b * cps[i].x;
        out.y += b * cps[i].y;
        out.z += b * cps[i].z;
    }
    out
}

#[inline]
fn points_equal(a: Point3, b: Point3) -> bool {
    const EQ_EPS: f64 = 1e-12;
    (a.x - b.x).abs() < EQ_EPS && (a.y - b.y).abs() < EQ_EPS && (a.z - b.z).abs() < EQ_EPS
}

#[cfg(test)]
#[path = "kp13_bi5_gtri5_tests.rs"]
mod tests;
