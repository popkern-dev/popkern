//! Hermite-form [`KpGuidePatch`] constructor.
//!
//! Turns per-corner Hermite data (position + u/v tangent + mixed
//! twist) into a bicubic Bezier guide patch. This is the **low-
//! level primitive** the K-P 2004 guided-surface construction
//! builds on; the higher-level orchestration (N-way Shirman-
//! Séquin enclosure at the central vertex, K-P 2004 §4.2 outer-
//! corner formula) computes the Hermite data and then delegates
//! here.
//!
//! # Hermite → CPs formulas
//!
//! Given corner positions `P_{00}`, `P_{30}`, `P_{33}`, `P_{03}`
//! and per-corner tangents `∂G/∂u` (`Pu`), `∂G/∂v` (`Pv`), and
//! mixed twist `∂²G/∂u∂v` (`Puv`), the 16 CPs of a bicubic Bezier
//! patch follow from the classical boundary-derivative identities
//! `B'(0) = n·(P_1 − P_0)` and `∂²B/∂u∂v` at each corner. Result:
//!
//! ```text
//!   Corners:      grid[0][0], grid[3][0], grid[3][3], grid[0][3]
//!                 = corners[0..=3].pos (CCW: 0 = u=0 v=0, 1 = u=1 v=0,
//!                                       2 = u=1 v=1, 3 = u=0 v=1)
//!
//!   Edge CPs:     grid[1][0] = c0.pos + (1/3)·c0.Pu
//!                 grid[2][0] = c1.pos − (1/3)·c1.Pu
//!                 grid[2][3] = c2.pos − (1/3)·c2.Pu
//!                 grid[1][3] = c3.pos + (1/3)·c3.Pu
//!                 grid[0][1] = c0.pos + (1/3)·c0.Pv
//!                 grid[3][1] = c1.pos + (1/3)·c1.Pv
//!                 grid[3][2] = c2.pos − (1/3)·c2.Pv
//!                 grid[0][2] = c3.pos − (1/3)·c3.Pv
//!
//!   Inner CPs:    grid[1][1] = c0.pos + (1/3)·(Pu+Pv) + (1/9)·Puv
//!                 grid[2][1] = c1.pos + (1/3)·(−Pu+Pv) − (1/9)·Puv
//!                 grid[2][2] = c2.pos + (1/3)·(−Pu−Pv) + (1/9)·Puv
//!                 grid[1][2] = c3.pos + (1/3)·(Pu−Pv) − (1/9)·Puv
//! ```
//!
//! The sign pattern is `u_sign · v_sign` on the twist term, and
//! individual `u_sign` / `v_sign` on the tangent terms, where
//! `u_sign = +1` at u-corners with u = 0 and `−1` with u = 1
//! (likewise `v_sign`). Symbolically:
//!
//! ```text
//!   inner(c) = c.pos + u_sign(c)·(1/3)·c.Pu
//!                    + v_sign(c)·(1/3)·c.Pv
//!                    + u_sign(c)·v_sign(c)·(1/9)·c.Puv
//! ```
//!
//! # Reconstruction invariant
//!
//! By construction, the returned [`KpGuidePatch`] satisfies
//!
//! ```text
//!   G.point_at(u_c, v_c)                = c.pos
//!   G.partial_u(u_c, v_c)               = c.Pu
//!   G.partial_v(u_c, v_c)               = c.Pv
//!   G.second_partial_uv(u_c, v_c)       = c.Puv
//! ```
//!
//! at each of the 4 corners `(u_c, v_c) ∈ {0, 1}²`. This is the
//! load-bearing contract tested in
//! [`super::kp_guide_build_tests`].

use super::kp_guide_patch::KpGuidePatch;
use crate::geometry::{Point3, Vec3};

/// Per-corner Hermite data for [`build_guide_from_hermite`].
///
/// Position + two tangent directions (u, v) + one twist (mixed
/// partial ∂²G/∂u∂v). The twist can be zero for a "no-twist"
/// guide; this corresponds to the parallelogram rule for inner
/// CPs. Non-zero twist is required for the K-P 2004 construction
/// at the central vertex, where the Shirman-Séquin enclosure
/// dictates a specific mixed-partial value.
#[derive(Debug, Clone, Copy)]
pub struct GuideCornerData {
    /// Corner position.
    pub position: Point3,
    /// Tangent vector `∂G/∂u` at the corner.
    pub tangent_u: Vec3,
    /// Tangent vector `∂G/∂v` at the corner.
    pub tangent_v: Vec3,
    /// Mixed partial `∂²G/∂u∂v` at the corner.
    pub twist: Vec3,
}

impl GuideCornerData {
    /// Construct a corner with zero twist. Use for "no-twist"
    /// guides where the mixed partial is not constrained by the
    /// caller.
    pub fn no_twist(position: Point3, tangent_u: Vec3, tangent_v: Vec3) -> Self {
        Self {
            position,
            tangent_u,
            tangent_v,
            twist: Vec3::zero(),
        }
    }
}

/// Construct a [`KpGuidePatch`] from 4 corners of Hermite data.
///
/// Corner order is counter-clockwise starting at `(u, v) = (0, 0)`:
///
/// ```text
///   corners[0] ↔ (u=0, v=0)     corners[3] ↔ (u=0, v=1)
///   corners[1] ↔ (u=1, v=0)     corners[2] ↔ (u=1, v=1)
/// ```
///
/// The returned guide satisfies `point_at(u_c, v_c) = c.position`,
/// `partial_u(u_c, v_c) = c.tangent_u`,
/// `partial_v(u_c, v_c) = c.tangent_v`, and
/// `second_partial_uv(u_c, v_c) = c.twist` at each of the four
/// corners by construction (not by solve).
pub fn build_guide_from_hermite(corners: [GuideCornerData; 4]) -> KpGuidePatch {
    let c0 = corners[0];
    let c1 = corners[1];
    let c2 = corners[2];
    let c3 = corners[3];

    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];

    // Corners.
    grid[0][0] = c0.position;
    grid[3][0] = c1.position;
    grid[3][3] = c2.position;
    grid[0][3] = c3.position;

    // Edge CPs — u direction at the v = 0 and v = 1 boundaries.
    grid[1][0] = c0.position + c0.tangent_u * (1.0 / 3.0);
    grid[2][0] = c1.position + c1.tangent_u * (-1.0 / 3.0);
    grid[2][3] = c2.position + c2.tangent_u * (-1.0 / 3.0);
    grid[1][3] = c3.position + c3.tangent_u * (1.0 / 3.0);

    // Edge CPs — v direction at the u = 0 and u = 1 boundaries.
    grid[0][1] = c0.position + c0.tangent_v * (1.0 / 3.0);
    grid[3][1] = c1.position + c1.tangent_v * (1.0 / 3.0);
    grid[3][2] = c2.position + c2.tangent_v * (-1.0 / 3.0);
    grid[0][2] = c3.position + c3.tangent_v * (-1.0 / 3.0);

    // Inner CPs — Hermite reconstruction with per-corner sign
    // pattern. `u_sign`, `v_sign` track which side of the unit
    // square the corner is on.
    grid[1][1] = inner_cp(&c0, 1.0, 1.0);
    grid[2][1] = inner_cp(&c1, -1.0, 1.0);
    grid[2][2] = inner_cp(&c2, -1.0, -1.0);
    grid[1][2] = inner_cp(&c3, 1.0, -1.0);

    KpGuidePatch { grid }
}

/// Hermite inner-CP formula. `u_sign` and `v_sign` are `±1`
/// depending on which corner of the unit square the data is
/// attached to (see [`build_guide_from_hermite`] for the
/// convention).
#[inline]
fn inner_cp(c: &GuideCornerData, u_sign: f64, v_sign: f64) -> Point3 {
    c.position
        + c.tangent_u * (u_sign / 3.0)
        + c.tangent_v * (v_sign / 3.0)
        + c.twist * (u_sign * v_sign / 9.0)
}

#[cfg(test)]
#[path = "kp_guide_build_tests.rs"]
mod tests;
