//! Phase I.14 — KP13 § 5 bi-5 cap assembly `g → a → p`.
//!
//! This module holds the data structures and invariants of the
//! final bi-5 cap surface `p` and its auxiliary predecessor `a`.
//! The cap covers the extraordinary vertex with `n` sectors, each
//! subdivided into a `2 × 2` layout of bi-5 sub-patches (see KP13
//! Fig. 10):
//!
//! ```text
//!   ┌───────┬───────┐
//!   │   m   │   r   │   ← "far" sub-patches (at v_k corner)
//!   ├───────┼───────┤
//!   │   ℓ   │   o   │   ← "near-star" sub-patches (at o corner)
//!   └───────┴───────┘
//!        |       |
//!     m_{k+1}   m_k (inter-sector boundary, shared with sector k+1)
//! ```
//!
//! Each sub-patch lives on its own `[0, 1]²` parameter domain
//! (piece-local coords). The four sub-patches share edges
//! internally (`C²`-connected within a sector, paper § 4.2),
//! and adjacent sectors share their `r`↔`ℓ` edges (which become
//! the `G²` strips in paper § 5).
//!
//! # Phase layout for the cap assembly
//!
//! KP13 § 5 builds the cap in two stages:
//!
//! 1. **Auxiliary cap `a`** (not yet `G²` across sector
//!    boundaries):
//!    - Outer rings of `a^{s,k}` for `s ∈ {ℓ, m, r}` match the
//!      tensor-border `b ∘ β` (paper § 4.2).
//!    - Corner jets of `a^{o,k}` from eq (17):
//!      `H_t a^{o,k} := H_t (C(g)^k · (1/λ) · L⁻¹ ∘ σ^o)`.
//!    - Internal `C²` connections between `a^{o,k}` and the
//!      three other pieces.
//!
//! 2. **Final cap `p`** with `G²` continuity across sector
//!    boundaries (paper § 5 eqs (18)–(21)):
//!    - Transfer the four `r^{k-1}` ↔ `ℓ^k` strip CPs from `a`.
//!    - Solve a global `F_1` minimisation (eq 21) over the free
//!      DoFs `G_3 ∪ {p^{o,k}_{i,2}}` (`5n + 6` per axis, or
//!      `5n + 7` at `n = 6`) against the `a`-target.
//!
//! # Scope of this module
//!
//! Scaffolding only — data types, zero-initialisation, piece-
//! wise `eval`, basic structural validation. Per-phase wiring
//! (eq (17) jet matching from `Bi5Guide`, β-reparameterised
//! tensor-border transfer, `C²` internal connections, eq (21)
//! `F_1` minimisation) lives in sibling modules.

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use super::bernstein_bi5::bernstein_bi5_all;
use crate::geometry::Point3;

/// Number of BB CPs per sub-patch (bi-5 tensor-product).
pub const BI5_PIECE_CP_COUNT: usize = 36;

/// Number of sub-patches per sector (`2 × 2` layout).
pub const PIECES_PER_SECTOR: usize = 4;

/// Identifier for one of the four sub-patches within a sector.
///
/// Layout reference — with `(u, v) ∈ [0, 1]²` the sector's own
/// coords, `(0, 0)` = star, `(1, 0)` = `m_k`, `(0, 1)` = `m_{k+1}`,
/// `(1, 1)` = `v_k`:
///
/// * [`CapPieceId::O`] — origin / star quadrant, sector-local
///   `(u, v) ∈ [0, 1/2]²`.
/// * [`CapPieceId::L`] — left (toward `m_{k+1}`) quadrant,
///   `(u, v) ∈ [0, 1/2] × [1/2, 1]`.
/// * [`CapPieceId::R`] — right (toward `m_k`) quadrant,
///   `(u, v) ∈ [1/2, 1] × [0, 1/2]`.
/// * [`CapPieceId::M`] — middle (toward `v_k`) quadrant,
///   `(u, v) ∈ [1/2, 1]²`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CapPieceId {
    /// Origin / star sub-patch.
    O,
    /// Left sub-patch (toward `m_{k+1}` along the `v`-axis).
    L,
    /// Right sub-patch (toward `m_k` along the `u`-axis).
    R,
    /// Middle sub-patch (toward `v_k`).
    M,
}

impl CapPieceId {
    /// Return all four piece identifiers in canonical order.
    pub const fn all() -> [CapPieceId; PIECES_PER_SECTOR] {
        [CapPieceId::O, CapPieceId::L, CapPieceId::R, CapPieceId::M]
    }

    /// Array index into a sector's `pieces` array.
    #[inline]
    pub const fn index(self) -> usize {
        match self {
            CapPieceId::O => 0,
            CapPieceId::L => 1,
            CapPieceId::R => 2,
            CapPieceId::M => 3,
        }
    }
}

/// One bi-5 sub-patch of the cap — 6 × 6 BB-grid on its own
/// `[0, 1]²` piece-local parameter domain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bi5CapPiece {
    pub cps: [[Point3; 6]; 6],
}

impl Bi5CapPiece {
    /// Zero piece (all CPs at origin).
    pub fn zero() -> Self {
        Self {
            cps: [[Point3::new(0.0, 0.0, 0.0); 6]; 6],
        }
    }

    /// Evaluate the sub-patch at piece-local `(u, v)` via tensor-
    /// product Bernstein basis.
    pub fn eval(&self, u: f64, v: f64) -> Point3 {
        let bu = bernstein_bi5_all(u);
        let bv = bernstein_bi5_all(v);
        let mut out = Point3::new(0.0, 0.0, 0.0);
        for p in 0..6 {
            for q in 0..6 {
                let w = bu[p] * bv[q];
                let cp = self.cps[p][q];
                out.x += w * cp.x;
                out.y += w * cp.y;
                out.z += w * cp.z;
            }
        }
        out
    }
}

/// One sector's four sub-patches (`O`, `L`, `R`, `M`) in canonical
/// order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bi5CapSector {
    pub pieces: [Bi5CapPiece; PIECES_PER_SECTOR],
}

impl Bi5CapSector {
    /// All pieces zero.
    pub fn zero() -> Self {
        Self {
            pieces: [Bi5CapPiece::zero(); PIECES_PER_SECTOR],
        }
    }

    /// Access a piece by [`CapPieceId`].
    #[inline]
    pub fn piece(&self, id: CapPieceId) -> &Bi5CapPiece {
        &self.pieces[id.index()]
    }

    /// Mutable access to a piece.
    #[inline]
    pub fn piece_mut(&mut self, id: CapPieceId) -> &mut Bi5CapPiece {
        &mut self.pieces[id.index()]
    }

    /// Evaluate the sector at sector-local `(u, v) ∈ [0, 1]²`
    /// by dispatching to the correct sub-patch and mapping
    /// into piece-local coords.
    ///
    /// The sector-local parameter domain is the unit square with
    /// `(0, 0) = o`, `(1, 0) = m_k`, `(0, 1) = m_{k+1}`,
    /// `(1, 1) = v_k`. Each sub-patch spans a `1/2 × 1/2` block
    /// of this square.
    pub fn eval(&self, u: f64, v: f64) -> Point3 {
        let (id, u_local, v_local) = dispatch_piece(u, v);
        self.piece(id).eval(u_local, v_local)
    }
}

/// The bi-5 cap: `n` sectors of 4 sub-patches each, with
/// `C⁰`/`C²`/`G²` invariants enforced by the assembly pipeline
/// in sibling modules (I.14.b–I.14.e).
#[derive(Debug, Clone, PartialEq)]
pub struct Bi5Cap {
    pub n: usize,
    pub sectors: Vec<Bi5CapSector>,
}

impl Bi5Cap {
    /// New zero cap with `n ≥ 3` sectors.
    pub fn zero(n: usize) -> Self {
        assert!(n >= 3, "KP13 bi-5 cap needs n ≥ 3, got {n}");
        Self {
            n,
            sectors: vec![Bi5CapSector::zero(); n],
        }
    }

    /// Access sector `k` (mod `n`).
    #[inline]
    pub fn sector(&self, k: usize) -> &Bi5CapSector {
        &self.sectors[k % self.n]
    }

    /// Mutable sector access (indexed mod `n`).
    #[inline]
    pub fn sector_mut(&mut self, k: usize) -> &mut Bi5CapSector {
        let n = self.n;
        &mut self.sectors[k % n]
    }

    /// Structural sanity: `n ≥ 3`, sector count matches `n`, all
    /// CP components are finite.
    pub fn validate(&self) -> Result<(), String> {
        if self.n < 3 {
            return Err(format!("Bi5Cap::validate: n = {} < 3", self.n));
        }
        if self.sectors.len() != self.n {
            return Err(format!(
                "Bi5Cap::validate: sectors.len() = {} ≠ n = {}",
                self.sectors.len(),
                self.n
            ));
        }
        for (k, sector) in self.sectors.iter().enumerate() {
            for id in CapPieceId::all() {
                let piece = sector.piece(id);
                for p in 0..6 {
                    for q in 0..6 {
                        let cp = piece.cps[p][q];
                        if !(cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite()) {
                            return Err(format!(
                                "Bi5Cap::validate: non-finite CP at sector {k}, piece {id:?}, ({p}, {q}): {cp:?}"
                            ));
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

/// Given sector-local `(u, v) ∈ [0, 1]²`, identify which sub-patch
/// the query falls into and map to piece-local `(u', v') ∈ [0, 1]²`.
///
/// Boundary convention: queries on the internal `u = 1/2` or
/// `v = 1/2` lines resolve to the **lower-index** piece (e.g.,
/// `u = 1/2` → O or L, not R or M). Adjacent sub-patches are
/// `C²`-connected after assembly so either choice agrees on the
/// shared boundary up to floating-point precision.
pub fn dispatch_piece(u: f64, v: f64) -> (CapPieceId, f64, f64) {
    let right = u > 0.5;
    let upper = v > 0.5;
    let id = match (right, upper) {
        (false, false) => CapPieceId::O,
        (false, true) => CapPieceId::L,
        (true, false) => CapPieceId::R,
        (true, true) => CapPieceId::M,
    };
    let u_local = if right { (u - 0.5) * 2.0 } else { u * 2.0 };
    let v_local = if upper { (v - 0.5) * 2.0 } else { v * 2.0 };
    (id, u_local, v_local)
}

#[cfg(test)]
#[path = "kp13_bi5_cap_tests.rs"]
mod tests;
