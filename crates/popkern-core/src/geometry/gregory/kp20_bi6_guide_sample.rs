//! KP20 guide-Hermite sampling into bi-6 patch grids.
//!
//! This module turns a total-degree-6 guide sector into bi-6
//! tensor-product corner blocks by sampling the guide's cubic Hermite
//! jet at the four corners of a patch-local domain. The maps here are
//! affine in the guide domain; the later KP20 τ wiring can supply the
//! same corner-jet interface after applying the full cubic chain rule.

#![allow(clippy::needless_range_loop)]

use super::bernstein_tri6::{
    TRI6_CP_COUNT, bernstein_tri6_all, bernstein_tri6_partial_u, bernstein_tri6_partial_v,
    bernstein_tri6_second_partial_uu, bernstein_tri6_second_partial_uv,
    bernstein_tri6_second_partial_vv, bernstein_tri6_third_partial_uuu,
    bernstein_tri6_third_partial_uuv, bernstein_tri6_third_partial_uvv,
    bernstein_tri6_third_partial_vvv,
};
use super::kp_polar_shared::weighted_sum;
use super::kp20_bi6_cap::{Kp20Bi6Patch, Kp20MacroPatch, Kp20PieceId};
use super::kp20_bi6_math::HERMITE_TO_BI6_W;
use crate::geometry::Point3;

/// Total-degree-6 KP20 guide sector in triangular BB form.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20Tri6GuideSector {
    /// Triangular BB coefficients in [`super::tri6_flat_index`] order.
    pub cps: [Point3; TRI6_CP_COUNT],
}

impl Kp20Tri6GuideSector {
    /// Construct a guide sector with all coefficients at the origin.
    pub fn zero() -> Self {
        Self {
            cps: [Point3::new(0.0, 0.0, 0.0); TRI6_CP_COUNT],
        }
    }

    /// Evaluate the guide value and derivatives through order three.
    pub fn jet3(&self, u: f64, v: f64) -> Kp20GuideJet3 {
        Kp20GuideJet3 {
            h: [
                [
                    accumulate(&self.cps, &bernstein_tri6_all(u, v)),
                    accumulate(&self.cps, &bernstein_tri6_partial_v(u, v)),
                    accumulate(&self.cps, &bernstein_tri6_second_partial_vv(u, v)),
                    accumulate(&self.cps, &bernstein_tri6_third_partial_vvv(u, v)),
                ],
                [
                    accumulate(&self.cps, &bernstein_tri6_partial_u(u, v)),
                    accumulate(&self.cps, &bernstein_tri6_second_partial_uv(u, v)),
                    accumulate(&self.cps, &bernstein_tri6_third_partial_uvv(u, v)),
                    Point3::new(0.0, 0.0, 0.0),
                ],
                [
                    accumulate(&self.cps, &bernstein_tri6_second_partial_uu(u, v)),
                    accumulate(&self.cps, &bernstein_tri6_third_partial_uuv(u, v)),
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(0.0, 0.0, 0.0),
                ],
                [
                    accumulate(&self.cps, &bernstein_tri6_third_partial_uuu(u, v)),
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(0.0, 0.0, 0.0),
                ],
            ],
        }
    }
}

/// Cubic Hermite jet `h[p][q] = ∂u^p ∂v^q f` for `p + q <= 3`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20GuideJet3 {
    /// Coordinate triples for the value and derivatives.
    pub h: [[Point3; 4]; 4],
}

/// Affine map from patch-local inward corner coordinates to guide
/// coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20AffineGuideMap {
    /// Guide-domain coordinate at the patch corner.
    pub origin: [f64; 2],
    /// Guide-domain delta for increasing the local first coordinate.
    pub du: [f64; 2],
    /// Guide-domain delta for increasing the local second coordinate.
    pub dv: [f64; 2],
}

impl Kp20AffineGuideMap {
    /// Construct a map `origin + s * du + t * dv`.
    pub const fn new(origin: [f64; 2], du: [f64; 2], dv: [f64; 2]) -> Self {
        Self { origin, du, dv }
    }
}

/// Corner of one bi-6 tensor-product patch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kp20PatchCorner {
    /// Lower-left corner `(0, 0)`.
    U0V0,
    /// Lower-right corner `(1, 0)`.
    U1V0,
    /// Upper-left corner `(0, 1)`.
    U0V1,
    /// Upper-right corner `(1, 1)`.
    U1V1,
}

impl Kp20PatchCorner {
    /// All corners in stable sampling order.
    pub const fn all() -> [Self; 4] {
        [Self::U0V0, Self::U1V0, Self::U0V1, Self::U1V1]
    }
}

/// Sample a 4x4 bi-6 corner block from the guide.
pub fn kp20_sample_corner_block(
    guide: &Kp20Tri6GuideSector,
    map: Kp20AffineGuideMap,
) -> [[Point3; 4]; 4] {
    let jet = compose_affine_jet(guide.jet3(map.origin[0], map.origin[1]), map.du, map.dv);
    kp20_corner_block_from_jet(jet)
}

/// Convert an already-composed cubic Hermite jet into a bi-6 corner block.
pub fn kp20_corner_block_from_jet(jet: Kp20GuideJet3) -> [[Point3; 4]; 4] {
    hermite_to_bi6_point_block(jet.h)
}

/// Populate a full bi-6 patch by sampling all four corner blocks and
/// averaging BB coefficients reached from multiple corners.
pub fn kp20_sample_patch_from_guide(
    guide: &Kp20Tri6GuideSector,
    corner_maps: [(Kp20PatchCorner, Kp20AffineGuideMap); 4],
) -> Kp20Bi6Patch {
    let zero = Point3::new(0.0, 0.0, 0.0);
    let mut sums = [[zero; 7]; 7];
    let mut counts = [[0_u8; 7]; 7];

    for (corner, map) in corner_maps {
        let block = kp20_sample_corner_block(guide, map);
        for i in 0..4 {
            for j in 0..4 {
                let (ii, jj) = oriented_index(corner, i, j);
                sums[ii][jj].x += block[i][j].x;
                sums[ii][jj].y += block[i][j].y;
                sums[ii][jj].z += block[i][j].z;
                counts[ii][jj] += 1;
            }
        }
    }

    let mut out = Kp20Bi6Patch::zero();
    for i in 0..7 {
        for j in 0..7 {
            debug_assert!(
                counts[i][j] > 0,
                "corner blocks must cover the full bi-6 grid"
            );
            let inv = 1.0 / counts[i][j] as f64;
            out.cps[i][j] = Point3::new(sums[i][j].x * inv, sums[i][j].y * inv, sums[i][j].z * inv);
        }
    }
    out
}

/// Populate all four quadrants of one KP20 macro-patch from guide
/// corner maps.
pub fn kp20_sample_macro_patch_from_guide(
    guide: &Kp20Tri6GuideSector,
    piece_maps: [(Kp20PieceId, [(Kp20PatchCorner, Kp20AffineGuideMap); 4]); 4],
) -> Kp20MacroPatch {
    let mut out = Kp20MacroPatch::zero();
    for (piece_id, maps) in piece_maps {
        *out.piece_mut(piece_id) = kp20_sample_patch_from_guide(guide, maps);
    }
    out
}

fn compose_affine_jet(jet: Kp20GuideJet3, du: [f64; 2], dv: [f64; 2]) -> Kp20GuideJet3 {
    let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for a in 0..4 {
        for b in 0..(4 - a) {
            let mut terms = Vec::new();
            for i in 0..=a {
                for j in 0..=b {
                    let x_order = i + j;
                    let y_order = a + b - x_order;
                    let weight = binom(a, i) as f64
                        * binom(b, j) as f64
                        * du[0].powi(i as i32)
                        * du[1].powi((a - i) as i32)
                        * dv[0].powi(j as i32)
                        * dv[1].powi((b - j) as i32);
                    terms.push((weight, jet.h[x_order][y_order]));
                }
            }
            out[a][b] = weighted_sum(&terms);
        }
    }
    Kp20GuideJet3 { h: out }
}

fn hermite_to_bi6_point_block(h: [[Point3; 4]; 4]) -> [[Point3; 4]; 4] {
    let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            let mut terms = Vec::new();
            for p in 0..=i {
                for q in 0..=j {
                    terms.push((HERMITE_TO_BI6_W[i][p] * HERMITE_TO_BI6_W[j][q], h[p][q]));
                }
            }
            out[i][j] = weighted_sum(&terms);
        }
    }
    out
}

fn oriented_index(corner: Kp20PatchCorner, i: usize, j: usize) -> (usize, usize) {
    match corner {
        Kp20PatchCorner::U0V0 => (i, j),
        Kp20PatchCorner::U1V0 => (6 - i, j),
        Kp20PatchCorner::U0V1 => (i, 6 - j),
        Kp20PatchCorner::U1V1 => (6 - i, 6 - j),
    }
}

fn accumulate(cps: &[Point3; TRI6_CP_COUNT], weights: &[f64; TRI6_CP_COUNT]) -> Point3 {
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..TRI6_CP_COUNT {
        out.x += weights[i] * cps[i].x;
        out.y += weights[i] * cps[i].y;
        out.z += weights[i] * cps[i].z;
    }
    out
}

const fn binom(n: usize, k: usize) -> usize {
    match (n, k) {
        (_, 0) => 1,
        (1, 1) => 1,
        (2, 1) => 2,
        (2, 2) => 1,
        (3, 1) | (3, 2) => 3,
        (3, 3) => 1,
        _ => 0,
    }
}

#[cfg(test)]
#[path = "kp20_bi6_guide_sample_tests.rs"]
mod tests;
