//! Dual-jet corner cells of the fit's sectors (PK-595).
//!
//! The setback corner theorem (`docs/architecture/variational-corner-fit.md`)
//! rules out a cap that is C³ at a setback corner and G2 to both strips
//! there. KP20 met the same obstruction with dual cross-jet rational
//! corner cells (PK-519); the fit adopts them at each sector's rim corner
//! `S_s(1, 1)`.
//!
//! The four cells both rims' 2-jets read at the corner, control indices
//! `(m−2, m−2)`, `(m−2, m−3)`, `(m−3, m−2)`, `(m−3, m−3)`, carry a
//! u-claim `P^u` (the rim `u = 1`) and a v-claim `P^v` (the rim `v = 1`).
//! The cell's coefficient is `a·P^u + b·P^v` with `w_u = 1 − u`,
//! `w_v = 1 − v`, `D = w_u² + w_v²`, `a = w_v² / D`, `b = 1 − a`.
//!
//! On `u = 1`, `a = 1` and `∂a/∂u = ∂a/∂v = ∂²a/∂u∂v = 0`; `∂²a/∂u²`
//! is unbounded but multiplies the cell's basis, which vanishes on the
//! boundary. So the rim `u = 1` reads `P^u` alone through second order,
//! and `v = 1` reads `P^v` alone: each rim owns its corner jet. At the
//! corner itself the blend is 0/0, and a rim sample takes its own rim's
//! claim ([`CornerOwner`]).
//!
//! A corner-cell sector is degree 9, one span: it is then exactly a
//! [`RectangularGregoryBiNonicDualJet`] with its `(9, 9)` corner
//! activated (the same blend, the same cell positions and the same claim
//! ownership), so it serves as that patch ([`dual_jet_patch`]) and
//! inherits its exact NURBS export and the certifier's polynomial
//! boundary surrogate (PK-541): the rim seams certify against the
//! sector's plain u-claim or v-claim net, which carries the 2-jet on
//! that rim exactly.

use crate::geometry::Point3;
use crate::geometry::gregory::{DualJetCorner, RectangularGregoryBiNonicDualJet};

/// Which rim a sample at the sector corner belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum CornerOwner {
    /// An interior sample (the corner itself is never one).
    Interior,
    /// A sample of the rim `u = 1`.
    URim,
    /// A sample of the rim `v = 1`.
    VRim,
}

/// Dual cells per sector corner.
pub(super) const CELLS: usize = 4;

/// The degree of a corner-cell sector: the dual-jet patch's.
pub(super) const CORNER_CELL_DEGREE: usize = 9;

/// The dual-jet patch's activated corner, `(9, 9)` in its
/// `DUAL_JET_CORNER_POSITIONS` order.
const ACTIVE_CORNER: usize = 2;

/// The dual-cell slot of control `(i, j)` in an `m × m` sector.
pub(super) fn cell_slot(m: usize, i: usize, j: usize) -> Option<usize> {
    if m < 4 {
        return None;
    }
    let (near, far) = (m - 2, m - 3);
    match (i, j) {
        _ if i == near && j == near => Some(0),
        _ if i == near && j == far => Some(1),
        _ if i == far && j == near => Some(2),
        _ if i == far && j == far => Some(3),
        _ => None,
    }
}

/// The control `(i, j)` of dual-cell `slot` in an `m × m` sector.
#[cfg(test)]
pub(super) fn cell_control(m: usize, slot: usize) -> (usize, usize) {
    let (near, far) = (m - 2, m - 3);
    match slot {
        0 => (near, near),
        1 => (near, far),
        2 => (far, near),
        _ => (far, far),
    }
}

/// The u-claim weight `a` and its derivatives `[a, a_u, a_v, a_uu,
/// a_uv, a_vv]` at `(u, v)`.
pub(super) fn claim_weight(u: f64, v: f64, owner: CornerOwner) -> [f64; 6] {
    let (wu, wv) = (1.0 - u, 1.0 - v);
    let d = wu * wu + wv * wv;
    if d <= 0.0 {
        // The corner: the owning rim's claim, constant along it.
        return match owner {
            CornerOwner::VRim => [0.0; 6],
            _ => [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        };
    }
    let (d2, d3) = (d * d, d * d * d);
    [
        wv * wv / d,
        2.0 * wu * wv * wv / d2,
        -2.0 * wv * wu * wu / d2,
        2.0 * wv * wv * (3.0 * wu * wu - wv * wv) / d3,
        4.0 * wu * wv * (wv * wv - wu * wu) / d3,
        -2.0 * wu * wu * (3.0 * wv * wv - wu * wu) / d3,
    ]
}

/// The derivatives of `basis × weight`, both given as `[f, f_u, f_v,
/// f_uu, f_uv, f_vv]`.
pub(super) fn product(b: [f64; 6], w: [f64; 6]) -> [f64; 6] {
    [
        b[0] * w[0],
        b[1] * w[0] + b[0] * w[1],
        b[2] * w[0] + b[0] * w[2],
        b[3] * w[0] + 2.0 * b[1] * w[1] + b[0] * w[3],
        b[4] * w[0] + b[1] * w[2] + b[2] * w[1] + b[0] * w[4],
        b[5] * w[0] + 2.0 * b[2] * w[2] + b[0] * w[5],
    ]
}

/// `1 − a` from `a`'s derivatives.
pub(super) fn complement(w: [f64; 6]) -> [f64; 6] {
    [1.0 - w[0], -w[1], -w[2], -w[3], -w[4], -w[5]]
}

/// The degree-9 corner-cell sector as a dual-jet patch. `plain(i, j)`
/// is control `(i, j)` (the u-claim at a dual cell) and `claim_v(slot)`
/// the dual cell's v-claim.
pub(super) fn dual_jet_patch(
    plain: &dyn Fn(usize, usize) -> Point3,
    claim_v: &dyn Fn(usize) -> Point3,
) -> RectangularGregoryBiNonicDualJet {
    let grid: [[Point3; 10]; 10] = std::array::from_fn(|i| std::array::from_fn(|j| plain(i, j)));
    let pair = |slot: usize| {
        let (i, j) = slot_control(slot);
        (plain(i, j), claim_v(slot))
    };
    let mut corners = [None; 4];
    corners[ACTIVE_CORNER] = Some(DualJetCorner {
        cell_11: pair(0),
        cell_12: pair(1),
        cell_21: pair(2),
        cell_22: pair(3),
    });
    RectangularGregoryBiNonicDualJet { grid, corners }
}

/// The control of dual-cell `slot` in a degree-9 sector.
fn slot_control(slot: usize) -> (usize, usize) {
    let m = CORNER_CELL_DEGREE + 1;
    let (near, far) = (m - 2, m - 3);
    match slot {
        0 => (near, near),
        1 => (near, far),
        2 => (far, near),
        _ => (far, far),
    }
}

#[cfg(test)]
#[path = "corner_cells_tests.rs"]
mod tests;
