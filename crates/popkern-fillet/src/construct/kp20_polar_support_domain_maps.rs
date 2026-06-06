//! Support-domain map parameters for KP20 support-owned fold repair.

use popkern_core::geometry::gregory::{
    Kp20AffineGuideMap, Kp20PatchCorner, Kp20PieceCornerMaps, Kp20PieceId, Kp20SectorPieceMaps,
    kp20_positive_nonoverlap_c2_sector_piece_maps,
};

use super::TARGET_PIECE;

#[derive(Debug, Clone, Copy)]
pub(crate) struct Kp20SupportDomainParams {
    pub(super) origin: [f64; 2],
    pub(super) span_u: f64,
    pub(super) span_v: f64,
    pub(super) skew: f64,
}

impl Kp20SupportDomainParams {
    pub(super) const fn family() -> [Self; 8] {
        [
            Self::new([0.00, 0.20], 0.30, 0.30, 0.00),
            Self::new([0.00, 0.22], 0.30, 0.28, 0.02),
            Self::new([0.02, 0.20], 0.28, 0.30, -0.02),
            Self::new([0.00, 0.18], 0.32, 0.30, 0.00),
            Self::new([0.00, 0.24], 0.28, 0.28, 0.04),
            Self::new([0.04, 0.20], 0.26, 0.30, -0.04),
            Self::new([0.00, 0.16], 0.34, 0.30, 0.02),
            Self::new([0.02, 0.22], 0.30, 0.26, -0.02),
        ]
    }

    const fn new(origin: [f64; 2], span_u: f64, span_v: f64, skew: f64) -> Self {
        Self {
            origin,
            span_u,
            span_v,
            skew,
        }
    }

    pub(super) fn refinement_family(self) -> [Self; 10] {
        [
            self.shifted(0.00, -0.02, 0.00, 0.02, 0.00),
            self.shifted(0.00, 0.02, 0.00, -0.02, 0.00),
            self.shifted(0.02, 0.00, -0.02, 0.00, 0.00),
            self.shifted(-0.02, 0.00, 0.02, 0.00, 0.00),
            self.shifted(0.00, 0.00, 0.02, -0.02, 0.02),
            self.shifted(0.00, 0.00, -0.02, 0.02, -0.02),
            self.shifted(0.01, -0.01, 0.01, 0.00, 0.02),
            self.shifted(-0.01, 0.01, 0.00, 0.01, -0.02),
            self.shifted(0.00, -0.01, 0.02, 0.02, 0.00),
            self.shifted(0.02, 0.00, 0.00, 0.02, 0.02),
        ]
    }

    pub(super) fn shifted(
        self,
        origin_u: f64,
        origin_v: f64,
        span_u: f64,
        span_v: f64,
        skew: f64,
    ) -> Self {
        Self::new(
            [
                (self.origin[0] + origin_u).clamp(0.0, 0.10),
                (self.origin[1] + origin_v).clamp(0.14, 0.30),
            ],
            (self.span_u + span_u).clamp(0.24, 0.36),
            (self.span_v + span_v).clamp(0.22, 0.34),
            (self.skew + skew).clamp(-0.06, 0.06),
        )
    }
}

pub(super) fn support_domain_piece_maps(params: Kp20SupportDomainParams) -> Kp20SectorPieceMaps {
    support_domain_piece_maps_for_target(TARGET_PIECE, params)
}

pub(super) fn support_domain_piece_maps_for_target(
    target_piece: Kp20PieceId,
    params: Kp20SupportDomainParams,
) -> Kp20SectorPieceMaps {
    support_domain_piece_maps_for_targets(&[target_piece], params)
}

pub(super) fn support_domain_piece_maps_for_targets(
    target_pieces: &[Kp20PieceId],
    params: Kp20SupportDomainParams,
) -> Kp20SectorPieceMaps {
    let mut maps = kp20_positive_nonoverlap_c2_sector_piece_maps();
    for (piece, corner_maps) in &mut maps {
        if target_pieces.contains(piece) {
            *corner_maps = support_domain_corner_maps(params);
        }
    }
    maps
}

pub(super) fn support_domain_corner_maps(params: Kp20SupportDomainParams) -> Kp20PieceCornerMaps {
    square_maps(
        params.origin,
        [params.span_u, params.skew],
        [params.skew, params.span_v],
    )
}

pub(super) fn changed_target_pieces(n: usize) -> Vec<[bool; 4]> {
    changed_pieces_for_target(n, TARGET_PIECE)
}

pub(super) fn changed_pieces_for_target(n: usize, target_piece: Kp20PieceId) -> Vec<[bool; 4]> {
    changed_pieces_for_targets(n, &[target_piece])
}

pub(super) fn changed_pieces_for_targets(
    n: usize,
    target_pieces: &[Kp20PieceId],
) -> Vec<[bool; 4]> {
    let mut changed = vec![[false; 4]; n];
    for pieces in &mut changed {
        for target_piece in target_pieces {
            pieces[target_piece.index()] = true;
        }
    }
    changed
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
