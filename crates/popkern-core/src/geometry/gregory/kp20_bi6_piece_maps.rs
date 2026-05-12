//! KP20 sector piece-to-guide map wiring.
//!
//! The KP20 macro-patch has four bi-6 pieces per sector. These helpers
//! provide a stable default subdivision of the triangular guide domain
//! into the four sector-local pieces used by the side-by-side assembler.

use super::kp20_bi6_cap::{Kp20MacroPatch, Kp20PieceId};
use super::kp20_bi6_guide_sample::{
    Kp20AffineGuideMap, Kp20PatchCorner, Kp20Tri6GuideSector, kp20_sample_macro_patch_from_guide,
};

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
            square_maps([0.0, 0.0], [0.5, 0.0], [0.5, 0.5]),
        ),
        (
            Kp20PieceId::NearMidpointRight,
            square_maps([0.5, 0.5], [0.5, 0.0], [0.0, -0.5]),
        ),
        (
            Kp20PieceId::NearMidpointLeft,
            square_maps([0.0, 0.5], [0.5, 0.0], [0.0, -0.5]),
        ),
    ]
}

/// Sample one guide sector into one KP20 macro-patch using the
/// default piece maps.
pub fn sample_kp20_macro_patch_from_default_maps(guide: &Kp20Tri6GuideSector) -> Kp20MacroPatch {
    kp20_sample_macro_patch_from_guide(guide, kp20_default_sector_piece_maps())
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

#[cfg(test)]
#[path = "kp20_bi6_piece_maps_tests.rs"]
mod tests;
