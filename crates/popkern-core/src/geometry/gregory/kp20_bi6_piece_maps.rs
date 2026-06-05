//! KP20 sector piece-to-guide map wiring.
//!
//! The KP20 macro-patch has four bi-6 pieces per sector. These helpers
//! provide a stable default subdivision of the triangular guide domain
//! into the four sector-local pieces used by the side-by-side assembler.

use super::kp20_bi6_cap::{Kp20MacroPatch, Kp20PieceId};
use super::kp20_bi6_guide_sample::{
    Kp20AffineGuideMap, Kp20PatchCorner, Kp20Tri6GuideSector, kp20_sample_macro_patch_from_guide,
};
use super::kp20_bi6_tau::{Kp20TauMap, kp20_sample_corner_block_with_tau_affine};
use crate::geometry::Point3;

/// Corner maps for one KP20 bi-6 piece.
pub type Kp20PieceCornerMaps = [(Kp20PatchCorner, Kp20AffineGuideMap); 4];

/// All four piece maps for one sector.
pub type Kp20SectorPieceMaps = [(Kp20PieceId, Kp20PieceCornerMaps); 4];

/// Default piece-to-guide maps for the current KP20 macro-patch layout.
pub fn kp20_default_sector_piece_maps() -> Kp20SectorPieceMaps {
    [
        (
            Kp20PieceId::NearOriginLeft,
            square_maps([0.0, 0.0], [0.5, 0.0], [0.0, 0.5]),
        ),
        (
            Kp20PieceId::NearOriginRight,
            square_maps([0.0, 0.0], [0.0, 0.5], [0.5, 0.0]),
        ),
        (
            Kp20PieceId::NearMidpointRight,
            square_maps([0.5, 0.5], [0.5, -0.5], [-0.5, 0.0]),
        ),
        (
            Kp20PieceId::NearMidpointLeft,
            square_maps([0.0, 0.5], [0.5, 0.0], [0.0, -0.5]),
        ),
    ]
}

/// Orientation-valid piece-to-guide maps for comparison against the
/// current bounded-affine scaffold.
pub fn kp20_orientation_valid_sector_piece_maps() -> Kp20SectorPieceMaps {
    [
        (
            Kp20PieceId::NearOriginLeft,
            square_maps([0.0, 0.0], [0.25, 0.0], [0.0, 0.25]),
        ),
        (
            Kp20PieceId::NearOriginRight,
            square_maps([0.25, 0.0], [0.25, 0.0], [0.0, 0.25]),
        ),
        (
            Kp20PieceId::NearMidpointRight,
            square_maps([0.5, 0.0], [0.25, 0.0], [0.0, 0.25]),
        ),
        (
            Kp20PieceId::NearMidpointLeft,
            square_maps([0.0, 0.25], [0.25, 0.0], [0.0, 0.25]),
        ),
    ]
}

/// Default-footprint piece maps with every base affine map oriented
/// positively in the sector-local guide domain.
pub fn kp20_paper_oriented_sector_piece_maps() -> Kp20SectorPieceMaps {
    let mut maps = kp20_default_sector_piece_maps();
    for (_, corner_maps) in &mut maps {
        let base = corner_maps[0].1;
        let det = base.du[0] * base.dv[1] - base.du[1] * base.dv[0];
        if det < 0.0 {
            *corner_maps = square_maps(base.origin, base.dv, base.du);
        }
    }
    maps
}

/// Positive-orientation maps whose support-to-origin guide-image edges
/// match the C2 row ownership used by the KP20 assembler.
pub fn kp20_c2_stitched_sector_piece_maps() -> Kp20SectorPieceMaps {
    [
        (
            Kp20PieceId::NearOriginLeft,
            square_maps([0.25, 0.25], [0.25, 0.0], [0.0, 0.25]),
        ),
        (
            Kp20PieceId::NearOriginRight,
            square_maps([0.5, 0.0], [0.25, 0.0], [0.0, 0.25]),
        ),
        (
            Kp20PieceId::NearMidpointRight,
            square_maps([0.25, 0.0], [0.25, 0.0], [0.0, 0.25]),
        ),
        (
            Kp20PieceId::NearMidpointLeft,
            square_maps([0.0, 0.25], [0.25, 0.0], [0.0, 0.25]),
        ),
    ]
}

/// Domain-coverage report for the current affine piece-map scaffold.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20PieceMapDomainReport {
    /// Number of piece corners whose guide-domain coordinate lies outside
    /// the local triangular sector domain `u >= 0`, `v >= 0`, `u + v <= 1`.
    pub outside_triangle_corners: usize,
    /// Largest positive violation of `u + v <= 1`.
    pub max_simplex_excess: f64,
    /// Stable labels for offending piece corners.
    pub offending_corners: Vec<String>,
    /// Number of base piece maps whose affine orientation is not positive.
    pub non_positive_orientation_maps: usize,
    /// Stable labels for non-positive base piece maps.
    pub orientation_offenders: Vec<String>,
}

/// Guide-domain compatibility report for map edges consumed by assembly.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20PieceMapStitchReport {
    /// Largest mismatch between lower `u=1` support rows and upper `u=0`
    /// origin rows used by KP20 equation (12).
    pub max_internal_c2_edge_mismatch: f64,
    /// Largest mismatch between origin-piece `v=0` rows used by KP20
    /// equation (13)'s origin strip.
    pub max_origin_strip_edge_mismatch: f64,
}

impl Kp20PieceMapDomainReport {
    /// Whether every affine mapped corner stays inside the sector triangle.
    pub fn stays_inside_triangle(&self) -> bool {
        self.outside_triangle_corners == 0
    }
}

/// Measure whether the current affine piece maps extrapolate outside the
/// sector-local triangular guide domain.
pub fn kp20_default_piece_map_domain_report() -> Kp20PieceMapDomainReport {
    kp20_piece_map_domain_report(kp20_default_sector_piece_maps())
}

/// Measure whether caller-supplied affine piece maps stay inside the
/// sector-local triangular guide domain and preserve positive orientation.
pub fn kp20_piece_map_domain_report(maps: Kp20SectorPieceMaps) -> Kp20PieceMapDomainReport {
    let mut outside_triangle_corners = 0;
    let mut max_simplex_excess = 0.0_f64;
    let mut offending_corners = Vec::new();
    let mut non_positive_orientation_maps = 0;
    let mut orientation_offenders = Vec::new();
    for (piece, maps) in maps {
        let base_map = maps[0].1;
        let det = base_map.du[0] * base_map.dv[1] - base_map.du[1] * base_map.dv[0];
        if det <= 1.0e-12 {
            non_positive_orientation_maps += 1;
            orientation_offenders.push(format!("{piece:?}:det={det:.3e}"));
        }
        for (corner, map) in maps {
            let u = map.origin[0];
            let v = map.origin[1];
            let simplex_excess = (u + v - 1.0).max(0.0);
            let outside = u < -1e-12 || v < -1e-12 || simplex_excess > 1e-12;
            max_simplex_excess = max_simplex_excess.max(simplex_excess);
            if outside {
                outside_triangle_corners += 1;
                offending_corners.push(format!("{piece:?}.{corner:?}@({u:.3},{v:.3})"));
            }
        }
    }
    Kp20PieceMapDomainReport {
        outside_triangle_corners,
        max_simplex_excess,
        offending_corners,
        non_positive_orientation_maps,
        orientation_offenders,
    }
}

/// Measure guide-domain edge compatibility for caller-supplied piece maps.
pub fn kp20_piece_map_stitch_report(maps: Kp20SectorPieceMaps) -> Kp20PieceMapStitchReport {
    let Some(near_origin_left) = base_map(maps, Kp20PieceId::NearOriginLeft) else {
        return invalid_stitch_report();
    };
    let Some(near_origin_right) = base_map(maps, Kp20PieceId::NearOriginRight) else {
        return invalid_stitch_report();
    };
    let Some(near_midpoint_right) = base_map(maps, Kp20PieceId::NearMidpointRight) else {
        return invalid_stitch_report();
    };
    let Some(near_midpoint_left) = base_map(maps, Kp20PieceId::NearMidpointLeft) else {
        return invalid_stitch_report();
    };
    let max_internal_c2_edge_mismatch = max_edge_mismatch(
        near_midpoint_left,
        EdgeCoord::U1,
        near_origin_left,
        EdgeCoord::U0,
    )
    .max(max_edge_mismatch(
        near_midpoint_right,
        EdgeCoord::U1,
        near_origin_right,
        EdgeCoord::U0,
    ));
    let max_origin_strip_edge_mismatch = max_edge_mismatch(
        near_origin_left,
        EdgeCoord::V0,
        near_origin_right,
        EdgeCoord::V0,
    );
    Kp20PieceMapStitchReport {
        max_internal_c2_edge_mismatch,
        max_origin_strip_edge_mismatch,
    }
}

/// Sample one guide sector into one KP20 macro-patch using the
/// default piece maps.
pub fn sample_kp20_macro_patch_from_default_maps(guide: &Kp20Tri6GuideSector) -> Kp20MacroPatch {
    kp20_sample_macro_patch_from_guide(guide, kp20_default_sector_piece_maps())
}

/// Sample one guide sector into one KP20 macro-patch using caller-supplied
/// sector-local piece maps.
pub fn sample_kp20_macro_patch_from_piece_maps(
    guide: &Kp20Tri6GuideSector,
    maps: Kp20SectorPieceMaps,
) -> Kp20MacroPatch {
    kp20_sample_macro_patch_from_guide(guide, maps)
}

/// Sample one guide sector into one KP20 macro-patch using the
/// orientation-valid comparison maps.
pub fn sample_kp20_macro_patch_from_orientation_valid_maps(
    guide: &Kp20Tri6GuideSector,
) -> Kp20MacroPatch {
    kp20_sample_macro_patch_from_guide(guide, kp20_orientation_valid_sector_piece_maps())
}

/// Sample one guide sector into one KP20 macro-patch using the
/// paper-oriented default-footprint comparison maps.
pub fn sample_kp20_macro_patch_from_paper_oriented_maps(
    guide: &Kp20Tri6GuideSector,
) -> Kp20MacroPatch {
    kp20_sample_macro_patch_from_guide(guide, kp20_paper_oriented_sector_piece_maps())
}

/// Sample one guide sector into one KP20 macro-patch using the
/// C2-stitched comparison maps.
pub fn sample_kp20_macro_patch_from_c2_stitched_maps(
    guide: &Kp20Tri6GuideSector,
) -> Kp20MacroPatch {
    kp20_sample_macro_patch_from_guide(guide, kp20_c2_stitched_sector_piece_maps())
}

/// Sample one guide sector into one KP20 macro-patch using the
/// sector-local KP20 τ map before each default piece map.
pub fn sample_kp20_macro_patch_from_tau_maps(
    guide: &Kp20Tri6GuideSector,
    valence: usize,
    sector: usize,
) -> Result<Kp20MacroPatch, String> {
    let tau = Kp20TauMap::for_sector_local(valence, sector)?;
    let mut out = Kp20MacroPatch::zero();
    for (piece_id, maps) in kp20_default_sector_piece_maps() {
        *out.piece_mut(piece_id) = kp20_sample_patch_from_tau_maps(guide, &tau, maps);
    }
    Ok(out)
}

fn kp20_sample_patch_from_tau_maps(
    guide: &Kp20Tri6GuideSector,
    tau: &Kp20TauMap,
    corner_maps: Kp20PieceCornerMaps,
) -> super::kp20_bi6_cap::Kp20Bi6Patch {
    let zero = Point3::new(0.0, 0.0, 0.0);
    let mut sums = [[zero; 7]; 7];
    let mut counts = [[0_u8; 7]; 7];

    for (corner, map) in corner_maps {
        let block = kp20_sample_corner_block_with_tau_affine(guide, tau, map);
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

    let mut out = super::kp20_bi6_cap::Kp20Bi6Patch::zero();
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

fn square_maps(origin: [f64; 2], du: [f64; 2], dv: [f64; 2]) -> Kp20PieceCornerMaps {
    [
        (
            Kp20PatchCorner::U0V0,
            Kp20AffineGuideMap::new(origin, du, dv),
        ),
        (
            Kp20PatchCorner::U1V0,
            Kp20AffineGuideMap::new(add2(origin, du), scale2(du, -1.0), dv),
        ),
        (
            Kp20PatchCorner::U0V1,
            Kp20AffineGuideMap::new(add2(origin, dv), du, scale2(dv, -1.0)),
        ),
        (
            Kp20PatchCorner::U1V1,
            Kp20AffineGuideMap::new(
                add2(add2(origin, du), dv),
                scale2(du, -1.0),
                scale2(dv, -1.0),
            ),
        ),
    ]
}

const fn add2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] + b[0], a[1] + b[1]]
}

const fn scale2(a: [f64; 2], s: f64) -> [f64; 2] {
    [a[0] * s, a[1] * s]
}

fn invalid_stitch_report() -> Kp20PieceMapStitchReport {
    Kp20PieceMapStitchReport {
        max_internal_c2_edge_mismatch: f64::INFINITY,
        max_origin_strip_edge_mismatch: f64::INFINITY,
    }
}

fn base_map(maps: Kp20SectorPieceMaps, piece: Kp20PieceId) -> Option<Kp20AffineGuideMap> {
    for (id, corner_maps) in maps {
        if id == piece {
            return Some(corner_maps[0].1);
        }
    }
    None
}

fn max_edge_mismatch(
    first: Kp20AffineGuideMap,
    first_edge: EdgeCoord,
    second: Kp20AffineGuideMap,
    second_edge: EdgeCoord,
) -> f64 {
    [0.0, 0.5, 1.0]
        .into_iter()
        .map(|t| {
            distance2(
                edge_point(first, first_edge, t),
                edge_point(second, second_edge, t),
            )
        })
        .fold(0.0, f64::max)
}

#[derive(Debug, Clone, Copy)]
enum EdgeCoord {
    U0,
    U1,
    V0,
}

fn edge_point(map: Kp20AffineGuideMap, edge: EdgeCoord, t: f64) -> [f64; 2] {
    match edge {
        EdgeCoord::U0 => add2(map.origin, scale2(map.dv, t)),
        EdgeCoord::U1 => add2(add2(map.origin, map.du), scale2(map.dv, t)),
        EdgeCoord::V0 => add2(map.origin, scale2(map.du, t)),
    }
}

fn distance2(a: [f64; 2], b: [f64; 2]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    (dx * dx + dy * dy).sqrt()
}

const fn oriented_index(corner: Kp20PatchCorner, i: usize, j: usize) -> (usize, usize) {
    match corner {
        Kp20PatchCorner::U0V0 => (i, j),
        Kp20PatchCorner::U1V0 => (6 - i, j),
        Kp20PatchCorner::U0V1 => (i, 6 - j),
        Kp20PatchCorner::U1V1 => (6 - i, 6 - j),
    }
}

#[cfg(test)]
#[path = "kp20_bi6_piece_maps_tests.rs"]
mod tests;
