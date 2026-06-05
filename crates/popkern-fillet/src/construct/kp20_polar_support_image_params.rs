//! Fold-site-aware support-image affine map parameters.

use popkern_core::geometry::gregory::Kp20PieceId;
use popkern_core::geometry::gregory::{Kp20AffineGuideMap, Kp20PatchCorner, Kp20PieceCornerMaps};

use crate::construct::kp20_polar_fold_summary::Kp20FoldSite;

#[derive(Debug, Clone, Copy)]
pub(in crate::construct) struct Kp20SupportImageParams {
    uniform: f64,
    du_scale: f64,
    dv_scale: f64,
    shear: f64,
    fold_center_shift: f64,
    frame_kind: Kp20SupportImageFrameKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kp20SupportImageFrameKind {
    DefaultAffine,
    SectorLocalExactBoundary,
}

impl Kp20SupportImageParams {
    pub(super) fn identity() -> Self {
        Self::default_affine(1.0, 1.0, 1.0, 0.0, 0.0)
    }

    pub(super) fn coordinate_search_family() -> Vec<Self> {
        vec![
            Self::default_affine(0.92, 1.0, 1.0, 0.04, 0.0),
            Self::default_affine(0.89, 1.0, 1.0, 0.04, 0.0),
            Self::default_affine(0.95, 1.0, 1.0, 0.04, 0.0),
            Self::default_affine(0.92, 0.96, 1.0, 0.04, 0.0),
            Self::default_affine(0.92, 1.04, 1.0, 0.04, 0.0),
            Self::default_affine(0.92, 1.0, 0.96, 0.04, 0.0),
            Self::default_affine(0.92, 1.0, 1.04, 0.04, 0.0),
            Self::default_affine(0.92, 1.0, 1.0, 0.02, 0.0),
            Self::default_affine(0.92, 1.0, 1.0, 0.06, 0.0),
            Self::default_affine(0.92, 1.0, 1.0, 0.04, 0.08),
            Self::default_affine(0.92, 1.0, 1.0, 0.04, -0.08),
            Self::sector_local_exact_boundary(0.28, 0.0),
            Self::sector_local_exact_boundary(0.24, 0.08),
            Self::sector_local_exact_boundary(0.24, -0.08),
        ]
    }

    pub(super) fn lift_informed_family(anchor: Self) -> Vec<Self> {
        vec![
            Self::with_kind(
                (anchor.uniform - 0.04).max(0.84),
                anchor.du_scale,
                anchor.dv_scale,
                anchor.shear,
                anchor.fold_center_shift,
                anchor.frame_kind,
            ),
            Self::with_kind(
                (anchor.uniform + 0.04).min(1.02),
                anchor.du_scale,
                anchor.dv_scale,
                anchor.shear,
                anchor.fold_center_shift,
                anchor.frame_kind,
            ),
            Self::with_kind(
                anchor.uniform,
                (anchor.du_scale + 0.04).min(1.08),
                (anchor.dv_scale - 0.04).max(0.92),
                anchor.shear,
                anchor.fold_center_shift,
                anchor.frame_kind,
            ),
            Self::with_kind(
                anchor.uniform,
                (anchor.du_scale - 0.04).max(0.92),
                (anchor.dv_scale + 0.04).min(1.08),
                anchor.shear,
                anchor.fold_center_shift,
                anchor.frame_kind,
            ),
            Self::with_kind(
                anchor.uniform,
                anchor.du_scale,
                anchor.dv_scale,
                anchor.shear,
                (anchor.fold_center_shift + 0.06).clamp(-0.14, 0.14),
                anchor.frame_kind,
            ),
            Self::with_kind(
                anchor.uniform,
                anchor.du_scale,
                anchor.dv_scale,
                anchor.shear,
                (anchor.fold_center_shift - 0.06).clamp(-0.14, 0.14),
                anchor.frame_kind,
            ),
        ]
    }

    pub(super) fn support_side_reparameterization_family() -> Vec<Self> {
        vec![
            Self::sector_local_exact_boundary(0.30, 0.00),
            Self::sector_local_exact_boundary(0.28, 0.08),
            Self::sector_local_exact_boundary(0.28, -0.08),
            Self::sector_local_exact_boundary(0.24, 0.12),
            Self::sector_local_exact_boundary(0.24, -0.12),
            Self::default_affine(0.88, 0.96, 1.04, 0.02, 0.10),
            Self::default_affine(0.88, 1.04, 0.96, 0.02, -0.10),
            Self::default_affine(0.84, 1.0, 1.0, 0.04, 0.0),
        ]
    }

    const fn default_affine(
        uniform: f64,
        du_scale: f64,
        dv_scale: f64,
        shear: f64,
        fold_center_shift: f64,
    ) -> Self {
        Self::with_kind(
            uniform,
            du_scale,
            dv_scale,
            shear,
            fold_center_shift,
            Kp20SupportImageFrameKind::DefaultAffine,
        )
    }

    const fn sector_local_exact_boundary(uniform: f64, fold_center_shift: f64) -> Self {
        Self::with_kind(
            uniform,
            1.0,
            1.0,
            0.0,
            fold_center_shift,
            Kp20SupportImageFrameKind::SectorLocalExactBoundary,
        )
    }

    const fn with_kind(
        uniform: f64,
        du_scale: f64,
        dv_scale: f64,
        shear: f64,
        fold_center_shift: f64,
        frame_kind: Kp20SupportImageFrameKind,
    ) -> Self {
        Self {
            uniform,
            du_scale,
            dv_scale,
            shear,
            fold_center_shift,
            frame_kind,
        }
    }

    pub(super) fn reparameterized_maps(
        self,
        corner_maps: Kp20PieceCornerMaps,
        fold_site: Kp20FoldSite,
    ) -> Kp20PieceCornerMaps {
        if self.frame_kind == Kp20SupportImageFrameKind::SectorLocalExactBoundary {
            return self.sector_local_exact_boundary_maps(fold_site);
        }
        let base = corner_maps[0].1;
        let mut center = [
            base.origin[0] + 0.5 * base.du[0] + 0.5 * base.dv[0],
            base.origin[1] + 0.5 * base.du[1] + 0.5 * base.dv[1],
        ];
        let site_offset = [fold_site.u - 0.5, fold_site.v - 0.5];
        center[0] +=
            self.fold_center_shift * (site_offset[0] * base.du[0] + site_offset[1] * base.dv[0]);
        center[1] +=
            self.fold_center_shift * (site_offset[0] * base.du[1] + site_offset[1] * base.dv[1]);
        let du = [
            self.uniform * self.du_scale * base.du[0] + self.shear * base.dv[0],
            self.uniform * self.du_scale * base.du[1] + self.shear * base.dv[1],
        ];
        let dv = [
            self.shear * base.du[0] + self.uniform * self.dv_scale * base.dv[0],
            self.shear * base.du[1] + self.uniform * self.dv_scale * base.dv[1],
        ];
        let origin = [
            center[0] - 0.5 * du[0] - 0.5 * dv[0],
            center[1] - 0.5 * du[1] - 0.5 * dv[1],
        ];
        square_maps(origin, du, dv)
    }

    fn sector_local_exact_boundary_maps(self, fold_site: Kp20FoldSite) -> Kp20PieceCornerMaps {
        let center = sector_piece_center(fold_site.piece);
        let site_offset = [fold_site.u - 0.5, fold_site.v - 0.5];
        let center = [
            center[0] + self.fold_center_shift * site_offset[0] * self.uniform,
            center[1] + self.fold_center_shift * site_offset[1] * self.uniform,
        ];
        let half = 0.5 * self.uniform;
        let center = clamp_simplex_square_center(center, half);
        square_maps(
            [center[0] - half, center[1] - half],
            [self.uniform, 0.0],
            [0.0, self.uniform],
        )
    }

    pub(super) fn label(self) -> String {
        if self.frame_kind == Kp20SupportImageFrameKind::SectorLocalExactBoundary {
            return format!(
                "frame:sector_local_exact_boundary,uniform:{:.3},fold_center_shift:{:.3}",
                self.uniform, self.fold_center_shift
            );
        }
        if self.fold_center_shift.abs() <= 1.0e-12 {
            format!(
                "uniform:{:.3},du:{:.3},dv:{:.3},shear:{:.3}",
                self.uniform, self.du_scale, self.dv_scale, self.shear
            )
        } else {
            format!(
                "uniform:{:.3},du:{:.3},dv:{:.3},shear:{:.3},fold_center_shift:{:.3}",
                self.uniform, self.du_scale, self.dv_scale, self.shear, self.fold_center_shift
            )
        }
    }
}

fn sector_piece_center(piece: Kp20PieceId) -> [f64; 2] {
    match piece {
        Kp20PieceId::NearOriginLeft => [0.18, 0.18],
        Kp20PieceId::NearOriginRight => [0.46, 0.18],
        Kp20PieceId::NearMidpointRight => [0.56, 0.18],
        Kp20PieceId::NearMidpointLeft => [0.18, 0.46],
    }
}

fn clamp_simplex_square_center(mut center: [f64; 2], half: f64) -> [f64; 2] {
    center[0] = center[0].clamp(half, 1.0 - half);
    center[1] = center[1].clamp(half, 1.0 - half);
    let excess = center[0] + center[1] + 2.0 * half - 1.0;
    if excess > 0.0 {
        center[0] -= 0.5 * excess;
        center[1] -= 0.5 * excess;
    }
    center[0] = center[0].max(half);
    center[1] = center[1].max(half);
    center
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
mod tests {
    use popkern_core::geometry::gregory::Kp20PieceId;

    use super::*;

    #[test]
    fn unshifted_label_matches_existing_diagnostic_shape() {
        let params = Kp20SupportImageParams::coordinate_search_family()[0];

        assert_eq!(
            params.label(),
            "uniform:0.920,du:1.000,dv:1.000,shear:0.040"
        );
    }

    #[test]
    fn fold_center_shift_moves_map_center_toward_measured_site() {
        let base = square_maps([0.0, 0.0], [1.0, 0.0], [0.0, 1.0]);
        let params = Kp20SupportImageParams::default_affine(1.0, 1.0, 1.0, 0.0, 0.1);
        let shifted = params.reparameterized_maps(base, fold_site(0.8, 0.2));
        let shifted_center = map_center(shifted);

        assert!((shifted_center[0] - 0.53).abs() <= 1.0e-12);
        assert!((shifted_center[1] - 0.47).abs() <= 1.0e-12);
        assert_eq!(shifted[0].1.du, [1.0, 0.0]);
        assert_eq!(shifted[0].1.dv, [0.0, 1.0]);
    }

    #[test]
    fn sector_local_exact_boundary_frame_stays_inside_simplex() {
        let base = square_maps([0.0, 0.0], [0.5, 0.0], [0.0, 0.5]);
        let params = Kp20SupportImageParams::sector_local_exact_boundary(0.24, 0.08);
        let maps = params.reparameterized_maps(base, fold_site(0.9, 0.1));

        assert_eq!(
            params.label(),
            "frame:sector_local_exact_boundary,uniform:0.240,fold_center_shift:0.080"
        );
        assert!(base_map_determinant(maps) > 0.0);
        for (_, map) in maps {
            let u = map.origin[0];
            let v = map.origin[1];
            assert!(u >= -1.0e-12);
            assert!(v >= -1.0e-12);
            assert!(u + v <= 1.0 + 1.0e-12);
        }
    }

    fn fold_site(u: f64, v: f64) -> Kp20FoldSite {
        Kp20FoldSite {
            sector: 0,
            piece: Kp20PieceId::NearOriginLeft,
            u,
            v,
            min_jacobian_ratio: 0.0,
            min_jacobian_normal_cos: -1.0,
        }
    }

    fn map_center(maps: Kp20PieceCornerMaps) -> [f64; 2] {
        let map = maps[0].1;
        [
            map.origin[0] + 0.5 * map.du[0] + 0.5 * map.dv[0],
            map.origin[1] + 0.5 * map.du[1] + 0.5 * map.dv[1],
        ]
    }

    fn base_map_determinant(maps: Kp20PieceCornerMaps) -> f64 {
        let map = maps[0].1;
        map.du[0] * map.dv[1] - map.du[1] * map.dv[0]
    }
}
