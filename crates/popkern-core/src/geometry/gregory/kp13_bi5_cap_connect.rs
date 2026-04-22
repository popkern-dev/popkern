//! KP13 § 5 I.14.d — internal `C²` propagation between pieces
//! within each sector.
//!
//! Paper § 5: "We complete `a` by connecting `a^{o,k}` `C²` to
//! `a^{s,k}_{ij}`, `s ∈ {ℓ, m, r}`." This module implements that
//! `C²` join across each of the four internal shared boundaries
//! per sector:
//!
//! ```text
//!   ┌─── M ───┬─── R ───┐
//!   │    │    │         │
//!   ├────┼────o────┼────┤     ← sector-local (u, v) = (0.5, 0.5)
//!   │         │    │    │        is the 4-way crossing.
//!   └─── L ───┴─── O ───┘
//! ```
//!
//! Shared edges (in sector-local coords):
//! - `O ↔ R`: `u = 0.5` line, `v ∈ [0, 0.5]`.
//! - `O ↔ L`: `v = 0.5` line, `u ∈ [0, 0.5]`.
//! - `R ↔ M`: `v = 0.5` line, `u ∈ [0.5, 1]`.
//! - `L ↔ M`: `u = 0.5` line, `v ∈ [0.5, 1]`.
//!
//! # Bi-5 `C²` boundary formulas
//!
//! For two bi-5 pieces `A` (at piece-local `u ∈ [0, 1]`
//! corresponding to sector `u ∈ [u*, u* + 1/2]`) and `B` (at
//! piece-local `u ∈ [0, 1]` corresponding to sector
//! `u ∈ [u* − 1/2, u*]`), sharing the boundary at A's `u = 0` edge
//! = B's `u = 1` edge:
//!
//! Because both pieces scale piece-local parameter to sector
//! parameter by the same factor (`2`), the scalar `C²` boundary
//! constraints in sector coords become the following constraints
//! on the bi-5 CPs in piece-local coords:
//!
//! - `C⁰`:  `A[0][j] = B[5][j]`
//! - `C¹`:  `A[1][j] = 2 B[5][j] − B[4][j]`
//! - `C²`:  `A[2][j] = 4 B[5][j] − 4 B[4][j] + B[3][j]`
//!
//! (for every column index `j`).
//!
//! Analogous formulas for the `v`-direction shared boundary:
//! `A[i][0] = B[i][5]`, `A[i][1] = 2 B[i][5] − B[i][4]`,
//! `A[i][2] = 4 B[i][5] − 4 B[i][4] + B[i][3]`.
//!
//! # 4-way crossing at sector center
//!
//! The sector-local point `(0.5, 0.5)` is a 4-way meeting of `O`,
//! `R`, `L`, `M`. Propagating from `O` to `M` via `R` and via `L`
//! produces different values for `M`'s CPs near the crossing —
//! the propagation paths disagree at overlap because the chain of
//! `C²` constraints overdetermines the `M` interior. This module
//! applies the propagation in a deterministic order:
//!
//! 1. `O → R`, `O → L` (fills R's first 3 columns / L's first 3 rows).
//! 2. `R → M` (fills M's first 3 rows).
//! 3. `L → M` (fills M's first 3 columns; **overwrites R's
//!    contribution at the 3×3 overlap**).
//!
//! The overwrite choice is arbitrary — paper does not prescribe a
//! specific priority. A future refinement could average the two
//! paths or use F_1 minimisation (anticipating I.14.e's eq (21)).
//!
//! # Scope limitation: `O` interior
//!
//! `O`'s 12 non-corner-block CPs are currently zero (I.14.c's MVP
//! only sets corner-jet CPs). The propagation formulas use `O`'s
//! columns/rows at indices 3, 4, 5 — several of which *are*
//! still zero (`O[4][2], O[4][3], O[3][1..=4], ...`). Those
//! propagate as zero into the neighbors' derived CPs, leaving
//! holes in `R`, `L`, `M` interiors. Filling `O`'s interior via
//! de Casteljau subdivision of the guide at `t = 1/(2λ)` is a
//! paper-faithful follow-up that plugs these holes directly.

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use super::kp13_bi5_cap::{Bi5Cap, Bi5CapPiece, Bi5CapSector, CapPieceId};
use crate::geometry::Point3;

/// Apply the full I.14.d `C²` internal connection to every sector
/// of `cap`. Equivalent to calling [`connect_sector`] on each
/// sector in turn.
pub fn connect_all_sectors(cap: &mut Bi5Cap) {
    for k in 0..cap.n {
        connect_sector(cap.sector_mut(k));
    }
}

/// Apply the `C²` internal connections within one sector.
///
/// Propagation order:
/// 1. `O → R` (O's `u = 1` edge → R's `u = 0, 1, 2` columns).
/// 2. `O → L` (O's `v = 1` edge → L's `v = 0, 1, 2` rows).
/// 3. `R → M` (R's `v = 1` edge → M's `v = 0, 1, 2` rows).
/// 4. `L → M` (L's `u = 1` edge → M's `u = 0, 1, 2` columns;
///    overwrites R's contribution at the 3 × 3 overlap).
pub fn connect_sector(sector: &mut Bi5CapSector) {
    // Snapshot O before any propagation — other pieces read from
    // O, and the propagation formulas need all three of O's
    // `u = 5, 4, 3` columns (and analogous rows).
    let o = *sector.piece(CapPieceId::O);

    // Step 1: O → R. R's piece-local u-axis matches sector-local
    // u; the shared boundary is O's u = 1 (piece-local `i = 5`)
    // ↔ R's u = 0 (piece-local `i = 0`).
    {
        let r = sector.piece_mut(CapPieceId::R);
        for j in 0..6 {
            r.cps[0][j] = o.cps[5][j];
            r.cps[1][j] = sub(scale(o.cps[5][j], 2.0), o.cps[4][j]);
            r.cps[2][j] = add3(
                scale(o.cps[5][j], 4.0),
                scale(o.cps[4][j], -4.0),
                o.cps[3][j],
            );
        }
    }

    // Step 2: O → L. L's piece-local v-axis matches sector-local
    // v; the shared boundary is O's v = 1 ↔ L's v = 0.
    {
        let l = sector.piece_mut(CapPieceId::L);
        for i in 0..6 {
            l.cps[i][0] = o.cps[i][5];
            l.cps[i][1] = sub(scale(o.cps[i][5], 2.0), o.cps[i][4]);
            l.cps[i][2] = add3(
                scale(o.cps[i][5], 4.0),
                scale(o.cps[i][4], -4.0),
                o.cps[i][3],
            );
        }
    }

    // Step 3: R → M. R's v = 1 ↔ M's v = 0.
    let r_snapshot = *sector.piece(CapPieceId::R);
    {
        let m = sector.piece_mut(CapPieceId::M);
        for i in 0..6 {
            m.cps[i][0] = r_snapshot.cps[i][5];
            m.cps[i][1] = sub(scale(r_snapshot.cps[i][5], 2.0), r_snapshot.cps[i][4]);
            m.cps[i][2] = add3(
                scale(r_snapshot.cps[i][5], 4.0),
                scale(r_snapshot.cps[i][4], -4.0),
                r_snapshot.cps[i][3],
            );
        }
    }

    // Step 4: L → M. L's u = 1 ↔ M's u = 0. Overwrites R's
    // contribution to M at the 3 × 3 overlap `(0..=2, 0..=2)`.
    let l_snapshot = *sector.piece(CapPieceId::L);
    {
        let m = sector.piece_mut(CapPieceId::M);
        for j in 0..6 {
            m.cps[0][j] = l_snapshot.cps[5][j];
            m.cps[1][j] = sub(scale(l_snapshot.cps[5][j], 2.0), l_snapshot.cps[4][j]);
            m.cps[2][j] = add3(
                scale(l_snapshot.cps[5][j], 4.0),
                scale(l_snapshot.cps[4][j], -4.0),
                l_snapshot.cps[3][j],
            );
        }
    }
}

// ---- small vector helpers ----

#[inline]
fn add(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}

#[inline]
fn sub(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

#[inline]
fn scale(p: Point3, s: f64) -> Point3 {
    Point3::new(s * p.x, s * p.y, s * p.z)
}

#[inline]
fn add3(a: Point3, b: Point3, c: Point3) -> Point3 {
    Point3::new(a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z)
}

// Public for tests
#[doc(hidden)]
pub fn __test_reach_piece_cps(p: &Bi5CapPiece) -> &[[Point3; 6]; 6] {
    &p.cps
}

#[cfg(test)]
#[path = "kp13_bi5_cap_connect_tests.rs"]
mod tests;
