//! KP20 sector piece-to-guide map wiring.
//!
//! The KP20 macro-patch has four bi-6 pieces per sector. These helpers
//! provide a stable default subdivision of the triangular guide domain
//! into the four sector-local pieces used by the side-by-side assembler.

use super::bernstein_tri6::flat_index as tri6_flat_index;
use super::kp20_bi6_cap::{Kp20MacroPatch, Kp20PieceId};
use super::kp20_bi6_guide_sample::{
    Kp20AffineGuideMap, Kp20PatchCorner, Kp20Tri6GuideSector, kp20_corner_block_from_jet,
    kp20_sample_macro_patch_from_guide, kp20_sample_patch_from_guide,
};
use super::kp20_bi6_tau::{
    Kp20TauJet3, Kp20TauMap, kp20_compose_source_jet_with_tau_jet,
    kp20_sample_corner_block_with_tau_affine,
};
use super::{bernstein_bi6_all, bernstein_bi6_deriv_all};
use crate::geometry::Point3;

/// Corner maps for one KP20 bi-6 piece.
pub type Kp20PieceCornerMaps = [(Kp20PatchCorner, Kp20AffineGuideMap); 4];

/// All four piece maps for one sector.
pub type Kp20SectorPieceMaps = [(Kp20PieceId, Kp20PieceCornerMaps); 4];

/// Which KP20 τ construction is composed with the guide sampler.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kp20TauFamily {
    /// Existing kernel τ used by the first tau-based diagnostic sampler.
    Kernel,
    /// KP20 Appendix A / equation (22) τ in sector-local guide coordinates.
    Paper,
    /// KP20 Appendix A refined central-expansion map `τref = 2τ(u/2,v/2)`.
    PaperRefinedApex,
}

/// Sampler family for one macro-patch piece.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kp20PieceSamplerFamily {
    /// Sample the guide directly through the piece's affine map.
    AffineGuide,
    /// Compose the guide with a named τ family before Hermite sampling.
    Tau(Kp20TauFamily),
}

/// Per-piece sampler assignment.
pub type Kp20PieceSamplerFamilies = [(Kp20PieceId, Kp20PieceSamplerFamily); 4];

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

/// Positive-orientation maps whose C2 and origin-strip guide-image edges
/// both match the row ownership consumed by the KP20 assembler.
pub fn kp20_coupled_stitched_sector_piece_maps() -> Kp20SectorPieceMaps {
    [
        (
            Kp20PieceId::NearOriginLeft,
            square_maps([0.25, 0.25], [0.25, 0.0], [0.0, 0.25]),
        ),
        (
            Kp20PieceId::NearOriginRight,
            square_maps([0.25, 0.25], [0.25, 0.0], [0.0, 0.25]),
        ),
        (
            Kp20PieceId::NearMidpointRight,
            square_maps([0.0, 0.25], [0.25, 0.0], [0.0, 0.25]),
        ),
        (
            Kp20PieceId::NearMidpointLeft,
            square_maps([0.0, 0.25], [0.25, 0.0], [0.0, 0.25]),
        ),
    ]
}

/// Positive-orientation adjacent maps with C2-compatible support-origin
/// edges and a reversed-parameter origin-strip seam.
pub fn kp20_nonoverlap_coupled_sector_piece_maps() -> Kp20SectorPieceMaps {
    [
        (
            Kp20PieceId::NearOriginLeft,
            square_maps([0.25, 0.25], [0.25, 0.0], [0.0, 0.25]),
        ),
        (
            Kp20PieceId::NearOriginRight,
            square_maps([0.5, 0.25], [-0.25, 0.0], [0.0, -0.25]),
        ),
        (
            Kp20PieceId::NearMidpointRight,
            square_maps([0.75, 0.25], [-0.25, 0.0], [0.0, -0.25]),
        ),
        (
            Kp20PieceId::NearMidpointLeft,
            square_maps([0.0, 0.25], [0.25, 0.0], [0.0, 0.25]),
        ),
    ]
}

/// Positive-orientation nonoverlapping maps with C2-compatible support-origin
/// edges and reversed origin-strip compatibility.
pub fn kp20_positive_nonoverlap_c2_sector_piece_maps() -> Kp20SectorPieceMaps {
    kp20_nonoverlap_coupled_sector_piece_maps()
}

/// Paper-domain quadrant layout inferred from KP20 Fig. 8 and the predecessor
/// KP bi-5 quadrant convention: quadrant 1 is the lower-left origin quadrant,
/// then 2, 3, 4 run counter-clockwise around the 2x2 sector parameterization.
///
/// This is a domain-layout candidate for the full four-piece τ sampler. It is
/// intentionally separate from the C2 strip-ownership map families.
pub fn kp20_paper_quadrant_sector_piece_maps() -> Kp20SectorPieceMaps {
    [
        (
            Kp20PieceId::NearOriginLeft,
            square_maps([0.0, 0.0], [0.5, 0.0], [0.0, 0.5]),
        ),
        (
            Kp20PieceId::NearOriginRight,
            square_maps([0.5, 0.0], [0.5, 0.0], [0.0, 0.5]),
        ),
        (
            Kp20PieceId::NearMidpointRight,
            square_maps([0.5, 0.5], [0.5, 0.0], [0.0, 0.5]),
        ),
        (
            Kp20PieceId::NearMidpointLeft,
            square_maps([0.0, 0.5], [0.5, 0.0], [0.0, 0.5]),
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

/// Domain report for the composed map `τ ∘ q`, where `q` is a sector-local
/// affine piece map and τ is one of the KP20 cubic reparameterisations.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20TauPieceMapDomainReport {
    /// Tau construction used for the report.
    pub tau_family: Kp20TauFamily,
    /// Valence of the extraordinary point.
    pub valence: usize,
    /// Sector index.
    pub sector: usize,
    /// Number of composed domain samples.
    pub sample_count: usize,
    /// Number of samples outside `u >= 0`, `v >= 0`, `u + v <= 1`.
    pub outside_triangle_samples: usize,
    /// Minimum sampled guide-local `u`.
    pub min_u: f64,
    /// Minimum sampled guide-local `v`.
    pub min_v: f64,
    /// Maximum sampled guide-local `u`.
    pub max_u: f64,
    /// Maximum sampled guide-local `v`.
    pub max_v: f64,
    /// Maximum sampled `u + v`.
    pub max_u_plus_v: f64,
    /// Largest positive violation of the local guide triangle.
    pub max_simplex_excess: f64,
    /// Stable labels for offending samples.
    pub offending_samples: Vec<String>,
}

/// Combined sampled validity report for `τ ∘ q`, where `q` is a
/// sector-local affine piece map.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20TauPieceMapValidityReport {
    /// Sector-domain boundedness of `τ ∘ q`.
    pub domain: Kp20TauPieceMapDomainReport,
    /// Orientation and nonsingularity of `τ ∘ q`.
    pub orientation: Kp20PlanarTauOrientationReport,
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
    /// Whether the origin-strip edge match uses reversed parameter order.
    pub origin_strip_reversed: bool,
    /// Smallest center-to-center distance between semantic piece maps.
    pub min_piece_center_separation: f64,
}

/// Sampled orientation report for a planar bi-6 τ macro-patch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20PlanarTauOrientationReport {
    /// Number of sampled Jacobian determinants.
    pub sample_count: usize,
    /// Number of samples with negative determinant below tolerance.
    pub negative_determinants: usize,
    /// Number of samples whose determinant magnitude is at or below tolerance.
    pub near_zero_determinants: usize,
    /// Smallest sampled determinant.
    pub min_determinant: f64,
    /// Largest sampled determinant.
    pub max_determinant: f64,
    /// Piece containing the smallest sampled determinant.
    pub worst_piece: Kp20PieceId,
    /// Piece-local `u` coordinate of the smallest sampled determinant.
    pub worst_u: f64,
    /// Piece-local `v` coordinate of the smallest sampled determinant.
    pub worst_v: f64,
}

/// Sampled sector-domain report for a planar bi-6 τ macro-patch.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20PlanarTauDomainReport {
    /// Number of sampled image points.
    pub sample_count: usize,
    /// Number of samples outside `u >= 0`, `v >= 0`, `u + v <= 1`.
    pub outside_triangle_samples: usize,
    /// Minimum sampled guide-local `u`.
    pub min_u: f64,
    /// Minimum sampled guide-local `v`.
    pub min_v: f64,
    /// Maximum sampled guide-local `u`.
    pub max_u: f64,
    /// Maximum sampled guide-local `v`.
    pub max_v: f64,
    /// Maximum sampled `u + v`.
    pub max_u_plus_v: f64,
    /// Largest positive violation of the local guide triangle.
    pub max_simplex_excess: f64,
    /// Stable label for the sample with the largest triangle violation.
    pub worst_sample: String,
}

/// Combined sampled validity gate for a planar bi-6 τ reparameterization.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20PlanarTauValidityReport {
    /// Sector-domain boundedness of τ.
    pub domain: Kp20PlanarTauDomainReport,
    /// Orientation and nonsingularity of τ.
    pub orientation: Kp20PlanarTauOrientationReport,
}

/// Rejection returned when a planar τ does not satisfy the sampled paper gate.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20PlanarTauSamplingError {
    /// Full sampled validity report for diagnostics and candidate ranking.
    pub validity: Box<Kp20PlanarTauValidityReport>,
}

impl Kp20PieceMapDomainReport {
    /// Whether every affine mapped corner stays inside the sector triangle.
    pub fn stays_inside_triangle(&self) -> bool {
        self.outside_triangle_corners == 0
    }
}

impl Kp20PlanarTauOrientationReport {
    /// Whether every sampled determinant is strictly positive at tolerance.
    pub fn is_orientation_preserving(&self) -> bool {
        self.negative_determinants == 0 && self.near_zero_determinants == 0
    }
}

impl Kp20TauPieceMapDomainReport {
    /// Whether all sampled composed points stay inside the sector triangle.
    pub fn stays_inside_triangle(&self) -> bool {
        self.outside_triangle_samples == 0
    }
}

impl Kp20PlanarTauDomainReport {
    /// Whether all sampled τ values stay inside the sector triangle.
    pub fn stays_inside_triangle(&self) -> bool {
        self.outside_triangle_samples == 0
    }
}

impl Kp20PlanarTauValidityReport {
    /// Whether the sampled τ is orientation preserving and nonsingular.
    ///
    /// KP20 Appendix A's explicit τ is a reparameterization in this sense; it
    /// is not generally bounded by the unit sector triangle.
    pub fn is_paper_valid_reparameterization(&self) -> bool {
        self.orientation.is_orientation_preserving()
    }

    /// Whether the sampled τ also stays inside the sector triangle.
    ///
    /// This is a strict sampler-containment gate, not the paper's definition
    /// of a valid reparameterization.
    pub fn is_strict_sector_bounded_reparameterization(&self) -> bool {
        self.domain.stays_inside_triangle() && self.orientation.is_orientation_preserving()
    }
}

impl Kp20TauPieceMapValidityReport {
    /// Whether all sampled composed maps are orientation preserving and
    /// nonsingular.
    pub fn is_paper_valid_reparameterization(&self) -> bool {
        self.orientation.is_orientation_preserving()
    }

    /// Whether all sampled composed maps also stay inside the sector triangle.
    pub fn is_strict_sector_bounded_reparameterization(&self) -> bool {
        self.domain.stays_inside_triangle() && self.orientation.is_orientation_preserving()
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

/// Measure whether `τ ∘ q` stays inside the sector-local guide triangle for
/// caller-supplied piece maps.
pub fn kp20_tau_piece_map_domain_report(
    maps: Kp20SectorPieceMaps,
    valence: usize,
    sector: usize,
    tau_family: Kp20TauFamily,
) -> Result<Kp20TauPieceMapDomainReport, String> {
    kp20_tau_piece_map_domain_report_for_pieces(maps, valence, sector, tau_family, |_| true)
}

/// Measure whether `τ ∘ q` stays inside the sector-local guide triangle for a
/// selected subset of semantic pieces.
pub fn kp20_tau_piece_map_domain_report_for_pieces(
    maps: Kp20SectorPieceMaps,
    valence: usize,
    sector: usize,
    tau_family: Kp20TauFamily,
    include_piece: impl Fn(Kp20PieceId) -> bool,
) -> Result<Kp20TauPieceMapDomainReport, String> {
    let tau = tau_for_family(tau_family, valence, sector)?;
    let mut sample_count = 0_usize;
    let mut outside_triangle_samples = 0_usize;
    let mut min_u = f64::INFINITY;
    let mut min_v = f64::INFINITY;
    let mut max_u = f64::NEG_INFINITY;
    let mut max_v = f64::NEG_INFINITY;
    let mut max_u_plus_v = f64::NEG_INFINITY;
    let mut max_simplex_excess = 0.0_f64;
    let mut offending_samples = Vec::new();

    for (piece, corner_maps) in maps {
        if !include_piece(piece) {
            continue;
        }
        for (corner, map) in corner_maps {
            for (su, sv) in COMPOSED_DOMAIN_SAMPLES {
                let q = add2(add2(map.origin, scale2(map.du, su)), scale2(map.dv, sv));
                let p = tau.jet3(q[0], q[1]).d[0][0];
                let u = p[0];
                let v = p[1];
                let u_plus_v = u + v;
                min_u = min_u.min(u);
                min_v = min_v.min(v);
                max_u = max_u.max(u);
                max_v = max_v.max(v);
                max_u_plus_v = max_u_plus_v.max(u_plus_v);
                let simplex_excess = (u_plus_v - 1.0).max(0.0);
                max_simplex_excess = max_simplex_excess.max(simplex_excess);
                let outside = u < -1.0e-12 || v < -1.0e-12 || simplex_excess > 1.0e-12;
                if outside {
                    outside_triangle_samples += 1;
                    if offending_samples.len() < 16 {
                        offending_samples.push(format!(
                            "{piece:?}.{corner:?}[{su:.2},{sv:.2}]->({u:.6},{v:.6})"
                        ));
                    }
                }
                sample_count += 1;
            }
        }
    }

    Ok(Kp20TauPieceMapDomainReport {
        tau_family,
        valence,
        sector,
        sample_count,
        outside_triangle_samples,
        min_u,
        min_v,
        max_u,
        max_v,
        max_u_plus_v,
        max_simplex_excess,
        offending_samples,
    })
}

/// Measure whether `τ ∘ q` is sampled as an oriented, nonsingular map for
/// caller-supplied sector-local affine piece maps.
pub fn kp20_tau_piece_map_orientation_report(
    maps: Kp20SectorPieceMaps,
    valence: usize,
    sector: usize,
    tau_family: Kp20TauFamily,
    determinant_tolerance: f64,
) -> Result<Kp20PlanarTauOrientationReport, String> {
    kp20_tau_piece_map_orientation_report_for_pieces(
        maps,
        valence,
        sector,
        tau_family,
        determinant_tolerance,
        |_| true,
    )
}

/// Measure whether `τ ∘ q` is sampled as an oriented, nonsingular map for a
/// selected subset of semantic pieces.
pub fn kp20_tau_piece_map_orientation_report_for_pieces(
    maps: Kp20SectorPieceMaps,
    valence: usize,
    sector: usize,
    tau_family: Kp20TauFamily,
    determinant_tolerance: f64,
    include_piece: impl Fn(Kp20PieceId) -> bool,
) -> Result<Kp20PlanarTauOrientationReport, String> {
    let tau = tau_for_family(tau_family, valence, sector)?;
    let mut report = Kp20PlanarTauOrientationReport {
        sample_count: 0,
        negative_determinants: 0,
        near_zero_determinants: 0,
        min_determinant: f64::INFINITY,
        max_determinant: f64::NEG_INFINITY,
        worst_piece: Kp20PieceId::NearOriginLeft,
        worst_u: 0.0,
        worst_v: 0.0,
    };

    for (piece, corner_maps) in maps {
        if !include_piece(piece) {
            continue;
        }
        for (corner, map) in corner_maps {
            let (sign_u, sign_v) = corner_derivative_signs(corner);
            for (su, sv) in COMPOSED_DOMAIN_SAMPLES {
                let q = add2(add2(map.origin, scale2(map.du, su)), scale2(map.dv, sv));
                let jet = tau.jet3(q[0], q[1]);
                let dx_dsu = sign_u * (jet.d[1][0][0] * map.du[0] + jet.d[0][1][0] * map.du[1]);
                let dy_dsu = sign_u * (jet.d[1][0][1] * map.du[0] + jet.d[0][1][1] * map.du[1]);
                let dx_dsv = sign_v * (jet.d[1][0][0] * map.dv[0] + jet.d[0][1][0] * map.dv[1]);
                let dy_dsv = sign_v * (jet.d[1][0][1] * map.dv[0] + jet.d[0][1][1] * map.dv[1]);
                let determinant = dx_dsu * dy_dsv - dx_dsv * dy_dsu;

                report.sample_count += 1;
                report.negative_determinants += usize::from(determinant < -determinant_tolerance);
                report.near_zero_determinants +=
                    usize::from(determinant.abs() <= determinant_tolerance);
                if determinant < report.min_determinant {
                    report.min_determinant = determinant;
                    report.worst_piece = piece;
                    report.worst_u = su;
                    report.worst_v = sv;
                }
                report.max_determinant = report.max_determinant.max(determinant);
            }
        }
    }

    Ok(report)
}

fn corner_derivative_signs(corner: Kp20PatchCorner) -> (f64, f64) {
    let sign_u = match corner {
        Kp20PatchCorner::U0V0 | Kp20PatchCorner::U0V1 => 1.0,
        Kp20PatchCorner::U1V0 | Kp20PatchCorner::U1V1 => -1.0,
    };
    let sign_v = match corner {
        Kp20PatchCorner::U0V0 | Kp20PatchCorner::U1V0 => 1.0,
        Kp20PatchCorner::U0V1 | Kp20PatchCorner::U1V1 => -1.0,
    };
    (sign_u, sign_v)
}

/// Measure the sampled paper-validity gate for explicit τ-family piece maps.
pub fn kp20_tau_piece_map_validity_report(
    maps: Kp20SectorPieceMaps,
    valence: usize,
    sector: usize,
    tau_family: Kp20TauFamily,
    determinant_tolerance: f64,
) -> Result<Kp20TauPieceMapValidityReport, String> {
    Ok(Kp20TauPieceMapValidityReport {
        domain: kp20_tau_piece_map_domain_report(maps, valence, sector, tau_family)?,
        orientation: kp20_tau_piece_map_orientation_report(
            maps,
            valence,
            sector,
            tau_family,
            determinant_tolerance,
        )?,
    })
}

/// Measure whether a planar bi-6 τ macro-patch is sampled as an oriented,
/// nonsingular map on each of its four square pieces.
pub fn kp20_planar_bi6_tau_orientation_report(
    tau_macro_patch: &Kp20MacroPatch,
    samples: &[f64],
    determinant_tolerance: f64,
) -> Kp20PlanarTauOrientationReport {
    let mut report = Kp20PlanarTauOrientationReport {
        sample_count: 0,
        negative_determinants: 0,
        near_zero_determinants: 0,
        min_determinant: f64::INFINITY,
        max_determinant: f64::NEG_INFINITY,
        worst_piece: Kp20PieceId::NearOriginLeft,
        worst_u: 0.0,
        worst_v: 0.0,
    };

    for piece_id in Kp20PieceId::all() {
        let patch = tau_macro_patch.piece(piece_id);
        for &u in samples {
            let bu = bernstein_bi6_all(u);
            let du = bernstein_bi6_deriv_all(u);
            for &v in samples {
                let bv = bernstein_bi6_all(v);
                let dv = bernstein_bi6_deriv_all(v);
                let mut dx_du = 0.0_f64;
                let mut dy_du = 0.0_f64;
                let mut dx_dv = 0.0_f64;
                let mut dy_dv = 0.0_f64;
                for i in 0..=6 {
                    for j in 0..=6 {
                        let cp = patch.cps[i][j];
                        let weight_du = du[i] * bv[j];
                        let weight_dv = bu[i] * dv[j];
                        dx_du += weight_du * cp.x;
                        dy_du += weight_du * cp.y;
                        dx_dv += weight_dv * cp.x;
                        dy_dv += weight_dv * cp.y;
                    }
                }
                let determinant = dx_du * dy_dv - dx_dv * dy_du;
                report.sample_count += 1;
                report.negative_determinants += usize::from(determinant < -determinant_tolerance);
                report.near_zero_determinants +=
                    usize::from(determinant.abs() <= determinant_tolerance);
                if determinant < report.min_determinant {
                    report.min_determinant = determinant;
                    report.worst_piece = piece_id;
                    report.worst_u = u;
                    report.worst_v = v;
                }
                report.max_determinant = report.max_determinant.max(determinant);
            }
        }
    }

    report
}

/// Measure whether a planar bi-6 τ macro-patch maps the sampled square pieces
/// into the sector-local guide triangle.
pub fn kp20_planar_bi6_tau_domain_report(
    tau_macro_patch: &Kp20MacroPatch,
    samples: &[f64],
    domain_tolerance: f64,
) -> Kp20PlanarTauDomainReport {
    let mut report = Kp20PlanarTauDomainReport {
        sample_count: 0,
        outside_triangle_samples: 0,
        min_u: f64::INFINITY,
        min_v: f64::INFINITY,
        max_u: f64::NEG_INFINITY,
        max_v: f64::NEG_INFINITY,
        max_u_plus_v: f64::NEG_INFINITY,
        max_simplex_excess: 0.0,
        worst_sample: String::new(),
    };

    for piece_id in Kp20PieceId::all() {
        let patch = tau_macro_patch.piece(piece_id);
        for &u in samples {
            let bu = bernstein_bi6_all(u);
            for &v in samples {
                let bv = bernstein_bi6_all(v);
                let mut x = 0.0_f64;
                let mut y = 0.0_f64;
                for (i, &bu_i) in bu.iter().enumerate() {
                    for (j, &bv_j) in bv.iter().enumerate() {
                        let cp = patch.cps[i][j];
                        let weight = bu_i * bv_j;
                        x += weight * cp.x;
                        y += weight * cp.y;
                    }
                }
                let u_plus_v = x + y;
                let simplex_excess = (u_plus_v - 1.0).max(0.0);
                let lower_excess = (-x).max(0.0).max((-y).max(0.0));
                let total_excess = simplex_excess.max(lower_excess);
                let outside = x < -domain_tolerance
                    || y < -domain_tolerance
                    || simplex_excess > domain_tolerance;

                report.sample_count += 1;
                report.outside_triangle_samples += usize::from(outside);
                report.min_u = report.min_u.min(x);
                report.min_v = report.min_v.min(y);
                report.max_u = report.max_u.max(x);
                report.max_v = report.max_v.max(y);
                report.max_u_plus_v = report.max_u_plus_v.max(u_plus_v);
                if total_excess > report.max_simplex_excess {
                    report.max_simplex_excess = total_excess;
                    report.worst_sample = format!("{piece_id:?}@({u:.2},{v:.2})=({x:.6},{y:.6})");
                }
            }
        }
    }

    report
}

/// Measure the sampled paper-validity gate for a planar bi-6 τ macro-patch.
pub fn kp20_planar_bi6_tau_validity_report(
    tau_macro_patch: &Kp20MacroPatch,
    samples: &[f64],
    domain_tolerance: f64,
    determinant_tolerance: f64,
) -> Kp20PlanarTauValidityReport {
    Kp20PlanarTauValidityReport {
        domain: kp20_planar_bi6_tau_domain_report(tau_macro_patch, samples, domain_tolerance),
        orientation: kp20_planar_bi6_tau_orientation_report(
            tau_macro_patch,
            samples,
            determinant_tolerance,
        ),
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
        false,
    )
    .max(max_edge_mismatch(
        near_midpoint_right,
        EdgeCoord::U1,
        near_origin_right,
        EdgeCoord::U0,
        false,
    ));
    let origin_same = max_edge_mismatch(
        near_origin_left,
        EdgeCoord::V0,
        near_origin_right,
        EdgeCoord::V0,
        false,
    );
    let origin_reversed = max_edge_mismatch(
        near_origin_left,
        EdgeCoord::V0,
        near_origin_right,
        EdgeCoord::V0,
        true,
    );
    let (max_origin_strip_edge_mismatch, origin_strip_reversed) = if origin_reversed < origin_same {
        (origin_reversed, true)
    } else {
        (origin_same, false)
    };
    let min_piece_center_separation = min_center_separation([
        near_origin_left,
        near_origin_right,
        near_midpoint_right,
        near_midpoint_left,
    ]);
    Kp20PieceMapStitchReport {
        max_internal_c2_edge_mismatch,
        max_origin_strip_edge_mismatch,
        origin_strip_reversed,
        min_piece_center_separation,
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
/// sector-local KP20 τ map before each default piece map.
pub fn sample_kp20_macro_patch_from_tau_maps(
    guide: &Kp20Tri6GuideSector,
    valence: usize,
    sector: usize,
) -> Result<Kp20MacroPatch, String> {
    sample_kp20_macro_patch_from_tau_family_maps(
        guide,
        valence,
        sector,
        Kp20TauFamily::Kernel,
        kp20_default_sector_piece_maps(),
    )
}

/// Sample one guide sector using the paper Appendix A τ and the current
/// positive nonoverlap map family. This is an explicit diagnostic/acceptance
/// entry point; production assembly keeps its existing sampler until the
/// paper-sampler campaign satisfies the cap-level gates.
pub fn sample_kp20_macro_patch_from_paper_tau_maps(
    guide: &Kp20Tri6GuideSector,
    valence: usize,
    sector: usize,
) -> Result<Kp20MacroPatch, String> {
    sample_kp20_macro_patch_from_tau_family_maps(
        guide,
        valence,
        sector,
        Kp20TauFamily::Paper,
        kp20_positive_nonoverlap_c2_sector_piece_maps(),
    )
}

/// Sample one guide sector with an explicit τ family and piece-map family.
pub fn sample_kp20_macro_patch_from_tau_family_maps(
    guide: &Kp20Tri6GuideSector,
    valence: usize,
    sector: usize,
    tau_family: Kp20TauFamily,
    maps: Kp20SectorPieceMaps,
) -> Result<Kp20MacroPatch, String> {
    let tau = tau_for_family(tau_family, valence, sector)?;
    let mut out = Kp20MacroPatch::zero();
    for (piece_id, maps) in maps {
        *out.piece_mut(piece_id) = kp20_sample_patch_from_tau_maps(guide, &tau, maps);
    }
    Ok(out)
}

/// Sample one guide sector with explicit per-piece sampler choices.
pub fn sample_kp20_macro_patch_from_piece_sampler_maps(
    guide: &Kp20Tri6GuideSector,
    valence: usize,
    sector: usize,
    sampler_families: Kp20PieceSamplerFamilies,
    maps: Kp20SectorPieceMaps,
) -> Result<Kp20MacroPatch, String> {
    let mut out = Kp20MacroPatch::zero();
    for (piece_id, corner_maps) in maps {
        let sampler = sampler_family_for_piece(sampler_families, piece_id)
            .ok_or_else(|| format!("missing sampler family for {piece_id:?}"))?;
        *out.piece_mut(piece_id) = match sampler {
            Kp20PieceSamplerFamily::AffineGuide => kp20_sample_patch_from_guide(guide, corner_maps),
            Kp20PieceSamplerFamily::Tau(tau_family) => {
                let tau = tau_for_family(tau_family, valence, sector)?;
                kp20_sample_patch_from_tau_maps(guide, &tau, corner_maps)
            }
        };
    }
    Ok(out)
}

/// Sample one guide sector through a full piecewise bi-6 planar τ macro-patch.
///
/// The supplied macro-patch is interpreted as a planar map into the guide
/// sector: each CP's `(x, y)` coordinates are guide-local coordinates sampled
/// by `g o τ`; `z` is ignored.
pub fn sample_kp20_macro_patch_from_planar_bi6_tau_macro_patch(
    guide: &Kp20Tri6GuideSector,
    tau_macro_patch: &Kp20MacroPatch,
) -> Kp20MacroPatch {
    let mut out = Kp20MacroPatch::zero();
    for piece_id in Kp20PieceId::all() {
        *out.piece_mut(piece_id) =
            kp20_sample_patch_from_planar_bi6_tau(guide, tau_macro_patch.piece(piece_id));
    }
    out
}

/// Sample a square-domain bi-6 source surface through a piecewise planar bi-6
/// τ macro-patch.
///
/// This is the KP20 Fig. 8/§4.4-shaped counterpart to the triangular guide
/// sampler above: τ is interpreted as a reparameterization of the tensor source
/// square, so the only source-domain containment gate is `0 <= u,v <= 1`.
/// Like the guide sampler, this constructs the output from total third-order
/// corner jets; it is exact for affine source/τ cases and is the local Hermite
/// macro-patch approximation for higher-degree source/τ composition.
pub fn sample_kp20_macro_patch_from_square_bi6_source_tau_macro_patch(
    source: &super::kp20_bi6_cap::Kp20Bi6Patch,
    tau_macro_patch: &Kp20MacroPatch,
) -> Kp20MacroPatch {
    let mut out = Kp20MacroPatch::zero();
    for piece_id in Kp20PieceId::all() {
        *out.piece_mut(piece_id) =
            kp20_sample_patch_from_square_bi6_source_tau(source, tau_macro_patch.piece(piece_id));
    }
    out
}

/// Sample one guide sector through a planar bi-6 τ only after τ satisfies the
/// strict sampled-source containment gate.
pub fn sample_kp20_macro_patch_from_valid_planar_bi6_tau_macro_patch(
    guide: &Kp20Tri6GuideSector,
    tau_macro_patch: &Kp20MacroPatch,
    validity_samples: &[f64],
    domain_tolerance: f64,
    determinant_tolerance: f64,
) -> Result<Kp20MacroPatch, Kp20PlanarTauSamplingError> {
    let validity = kp20_planar_bi6_tau_validity_report(
        tau_macro_patch,
        validity_samples,
        domain_tolerance,
        determinant_tolerance,
    );
    if !validity.is_strict_sector_bounded_reparameterization() {
        return Err(Kp20PlanarTauSamplingError {
            validity: Box::new(validity),
        });
    }
    Ok(sample_kp20_macro_patch_from_planar_bi6_tau_macro_patch(
        guide,
        tau_macro_patch,
    ))
}

/// Sample one guide sector through a planar bi-6 τ after checking the paper
/// reparameterization gate: orientation preserving and nonsingular.
///
/// Unlike [`sample_kp20_macro_patch_from_valid_planar_bi6_tau_macro_patch`],
/// this permits polynomial guide evaluation outside the unit sector triangle.
/// KP20 Appendix A's printed τ has this behavior, so callers that intentionally
/// follow the paper path should use this entry point and inspect the returned
/// validity report when strict containment is also desired.
pub fn sample_kp20_macro_patch_from_paper_valid_planar_bi6_tau_macro_patch(
    guide: &Kp20Tri6GuideSector,
    tau_macro_patch: &Kp20MacroPatch,
    validity_samples: &[f64],
    domain_tolerance: f64,
    determinant_tolerance: f64,
) -> Result<Kp20MacroPatch, Kp20PlanarTauSamplingError> {
    let validity = kp20_planar_bi6_tau_validity_report(
        tau_macro_patch,
        validity_samples,
        domain_tolerance,
        determinant_tolerance,
    );
    if !validity.is_paper_valid_reparameterization() {
        return Err(Kp20PlanarTauSamplingError {
            validity: Box::new(validity),
        });
    }
    Ok(sample_kp20_macro_patch_from_planar_bi6_tau_macro_patch(
        guide,
        tau_macro_patch,
    ))
}

/// KP20 §4.3 `F3` energy of a planar bi-6 τ macro-patch.
///
/// The functional is evaluated over all four pieces as
/// `∫∫ Σ_{a+b=3} 3!/(a!b!) ||∂u^a ∂v^b τ||² du dv`. Each CP's
/// `(x, y)` coordinates are used; `z` is ignored. The fixed
/// 7-point Gauss-Legendre tensor quadrature is exact for the squared
/// third-derivative polynomials of a bi-6 patch.
pub fn kp20_planar_bi6_tau_f3_energy(tau_macro_patch: &Kp20MacroPatch) -> f64 {
    Kp20PieceId::all()
        .into_iter()
        .map(|piece_id| kp20_planar_bi6_tau_patch_f3_energy(tau_macro_patch.piece(piece_id)))
        .sum()
}

fn sampler_family_for_piece(
    sampler_families: Kp20PieceSamplerFamilies,
    piece_id: Kp20PieceId,
) -> Option<Kp20PieceSamplerFamily> {
    sampler_families
        .into_iter()
        .find_map(|(id, family)| (id == piece_id).then_some(family))
}

fn tau_for_family(
    tau_family: Kp20TauFamily,
    valence: usize,
    sector: usize,
) -> Result<Kp20TauMap, String> {
    match tau_family {
        Kp20TauFamily::Kernel => Kp20TauMap::for_sector_local(valence, sector),
        Kp20TauFamily::Paper => Kp20TauMap::for_sector_local_paper(valence, sector),
        Kp20TauFamily::PaperRefinedApex => {
            Kp20TauMap::for_sector_local_paper_refined_apex(valence, sector)
        }
    }
}

const COMPOSED_DOMAIN_SAMPLES: [(f64, f64); 6] = [
    (0.0, 0.0),
    (1.0, 0.0),
    (0.0, 1.0),
    (1.0, 1.0),
    (0.5, 0.0),
    (0.0, 0.5),
];

const GAUSS7_01: &[(f64, f64); 7] = &[
    (0.025446043828620737, 0.06474248308443485),
    (0.12923440720030277, 0.13985269574463832),
    (0.2970774243113014, 0.19091502525255946),
    (0.5, 0.2089795918367347),
    (0.7029225756886985, 0.19091502525255946),
    (0.8707655927996972, 0.13985269574463832),
    (0.9745539561713793, 0.06474248308443485),
];

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
        for (i, row) in block.iter().enumerate() {
            for (j, point) in row.iter().enumerate() {
                let (ii, jj) = oriented_index(corner, i, j);
                sums[ii][jj].x += point.x;
                sums[ii][jj].y += point.y;
                sums[ii][jj].z += point.z;
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

fn kp20_sample_patch_from_planar_bi6_tau(
    guide: &Kp20Tri6GuideSector,
    tau: &super::kp20_bi6_cap::Kp20Bi6Patch,
) -> super::kp20_bi6_cap::Kp20Bi6Patch {
    let zero = Point3::new(0.0, 0.0, 0.0);
    let mut sums = [[zero; 7]; 7];
    let mut counts = [[0_u8; 7]; 7];

    for corner in [
        Kp20PatchCorner::U0V0,
        Kp20PatchCorner::U1V0,
        Kp20PatchCorner::U0V1,
        Kp20PatchCorner::U1V1,
    ] {
        let (u, v) = match corner {
            Kp20PatchCorner::U0V0 => (0.0, 0.0),
            Kp20PatchCorner::U1V0 => (1.0, 0.0),
            Kp20PatchCorner::U0V1 => (0.0, 1.0),
            Kp20PatchCorner::U1V1 => (1.0, 1.0),
        };
        let sign_u = match corner {
            Kp20PatchCorner::U0V0 | Kp20PatchCorner::U0V1 => 1.0,
            Kp20PatchCorner::U1V0 | Kp20PatchCorner::U1V1 => -1.0,
        };
        let sign_v = match corner {
            Kp20PatchCorner::U0V0 | Kp20PatchCorner::U1V0 => 1.0,
            Kp20PatchCorner::U0V1 | Kp20PatchCorner::U1V1 => -1.0,
        };
        let block = kp20_corner_block_from_jet(compose_tri6_guide_planar_bi6_tau_rect_jet(
            guide, tau, u, v, sign_u, sign_v,
        ));
        for (i, row) in block.iter().enumerate() {
            for (j, point) in row.iter().enumerate() {
                let (ii, jj) = oriented_index(corner, i, j);
                sums[ii][jj].x += point.x;
                sums[ii][jj].y += point.y;
                sums[ii][jj].z += point.z;
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

fn kp20_sample_patch_from_square_bi6_source_tau(
    source: &super::kp20_bi6_cap::Kp20Bi6Patch,
    tau: &super::kp20_bi6_cap::Kp20Bi6Patch,
) -> super::kp20_bi6_cap::Kp20Bi6Patch {
    let zero = Point3::new(0.0, 0.0, 0.0);
    let mut sums = [[zero; 7]; 7];
    let mut counts = [[0_u8; 7]; 7];

    for corner in Kp20PatchCorner::all() {
        let (u, v) = patch_corner_uv(corner);
        let (sign_u, sign_v) = patch_corner_derivative_signs(corner);
        let tau_jet = orient_tau_jet(planar_bi6_tau_jet3(tau, u, v), sign_u, sign_v);
        let source_xy = tau_jet.d[0][0];
        let source_jet = square_bi6_source_jet3(source, source_xy[0], source_xy[1]);
        let block =
            kp20_corner_block_from_jet(kp20_compose_source_jet_with_tau_jet(source_jet, tau_jet));
        for (i, row) in block.iter().enumerate() {
            for (j, point) in row.iter().enumerate() {
                let (ii, jj) = oriented_index(corner, i, j);
                sums[ii][jj].x += point.x;
                sums[ii][jj].y += point.y;
                sums[ii][jj].z += point.z;
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

fn patch_corner_uv(corner: Kp20PatchCorner) -> (f64, f64) {
    match corner {
        Kp20PatchCorner::U0V0 => (0.0, 0.0),
        Kp20PatchCorner::U1V0 => (1.0, 0.0),
        Kp20PatchCorner::U0V1 => (0.0, 1.0),
        Kp20PatchCorner::U1V1 => (1.0, 1.0),
    }
}

fn patch_corner_derivative_signs(corner: Kp20PatchCorner) -> (f64, f64) {
    let sign_u = match corner {
        Kp20PatchCorner::U0V0 | Kp20PatchCorner::U0V1 => 1.0,
        Kp20PatchCorner::U1V0 | Kp20PatchCorner::U1V1 => -1.0,
    };
    let sign_v = match corner {
        Kp20PatchCorner::U0V0 | Kp20PatchCorner::U1V0 => 1.0,
        Kp20PatchCorner::U0V1 | Kp20PatchCorner::U1V1 => -1.0,
    };
    (sign_u, sign_v)
}

fn square_bi6_source_jet3(
    source: &super::kp20_bi6_cap::Kp20Bi6Patch,
    u: f64,
    v: f64,
) -> super::kp20_bi6_guide_sample::Kp20GuideJet3 {
    let bu = [
        bernstein1d_derivative_basis(6, 0, u),
        bernstein1d_derivative_basis(6, 1, u),
        bernstein1d_derivative_basis(6, 2, u),
        bernstein1d_derivative_basis(6, 3, u),
    ];
    let bv = [
        bernstein1d_derivative_basis(6, 0, v),
        bernstein1d_derivative_basis(6, 1, v),
        bernstein1d_derivative_basis(6, 2, v),
        bernstein1d_derivative_basis(6, 3, v),
    ];
    let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for a in 0..4 {
        for b in 0..(4 - a) {
            for (i, &bu_i) in bu[a].iter().enumerate() {
                for (j, &bv_j) in bv[b].iter().enumerate() {
                    let cp = source.cps[i][j];
                    let w = bu_i * bv_j;
                    out[a][b].x += w * cp.x;
                    out[a][b].y += w * cp.y;
                    out[a][b].z += w * cp.z;
                }
            }
        }
    }
    super::kp20_bi6_guide_sample::Kp20GuideJet3 { h: out }
}

#[derive(Debug, Clone, Copy)]
struct RectTaylorJet {
    c: [[f64; 4]; 4],
}

impl RectTaylorJet {
    fn zero() -> Self {
        Self { c: [[0.0; 4]; 4] }
    }

    fn constant(value: f64) -> Self {
        let mut out = Self::zero();
        out.c[0][0] = value;
        out
    }

    fn sub(self, rhs: Self) -> Self {
        let mut out = Self::zero();
        for a in 0..4 {
            for b in 0..4 {
                out.c[a][b] = self.c[a][b] - rhs.c[a][b];
            }
        }
        out
    }

    fn scale(self, scale: f64) -> Self {
        let mut out = Self::zero();
        for a in 0..4 {
            for b in 0..4 {
                out.c[a][b] = scale * self.c[a][b];
            }
        }
        out
    }

    fn mul(self, rhs: Self) -> Self {
        let mut out = Self::zero();
        for a in 0..4 {
            for b in 0..4 {
                let mut value = 0.0;
                for i in 0..=a {
                    for j in 0..=b {
                        value += self.c[i][j] * rhs.c[a - i][b - j];
                    }
                }
                out.c[a][b] = value;
            }
        }
        out
    }

    fn powi(self, exponent: usize) -> Self {
        let mut out = Self::constant(1.0);
        for _ in 0..exponent {
            out = out.mul(self);
        }
        out
    }
}

fn compose_tri6_guide_planar_bi6_tau_rect_jet(
    guide: &Kp20Tri6GuideSector,
    tau: &super::kp20_bi6_cap::Kp20Bi6Patch,
    u: f64,
    v: f64,
    sign_u: f64,
    sign_v: f64,
) -> super::kp20_bi6_guide_sample::Kp20GuideJet3 {
    let [x, y] = planar_bi6_tau_rect_taylor(tau, u, v, sign_u, sign_v);
    let one_minus_x_minus_y = RectTaylorJet::constant(1.0).sub(x).sub(y);
    let mut coeffs = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];

    for i in 0..=6 {
        for j in 0..=(6 - i) {
            let k = 6 - i - j;
            let basis = x
                .powi(i)
                .mul(y.powi(j))
                .mul(one_minus_x_minus_y.powi(k))
                .scale(multinomial3(6, i, j, k));
            let cp = guide.cps[tri6_flat_index(6, i, j)];
            for a in 0..4 {
                for b in 0..4 {
                    let weight = basis.c[a][b];
                    coeffs[a][b].x += weight * cp.x;
                    coeffs[a][b].y += weight * cp.y;
                    coeffs[a][b].z += weight * cp.z;
                }
            }
        }
    }

    let mut h = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for a in 0..4 {
        for b in 0..4 {
            let scale = small_factorial(a) * small_factorial(b);
            h[a][b] = Point3::new(
                scale * coeffs[a][b].x,
                scale * coeffs[a][b].y,
                scale * coeffs[a][b].z,
            );
        }
    }
    super::kp20_bi6_guide_sample::Kp20GuideJet3 { h }
}

fn planar_bi6_tau_rect_taylor(
    tau: &super::kp20_bi6_cap::Kp20Bi6Patch,
    u: f64,
    v: f64,
    sign_u: f64,
    sign_v: f64,
) -> [RectTaylorJet; 2] {
    let bu = [
        bernstein1d_derivative_basis(6, 0, u),
        bernstein1d_derivative_basis(6, 1, u),
        bernstein1d_derivative_basis(6, 2, u),
        bernstein1d_derivative_basis(6, 3, u),
    ];
    let bv = [
        bernstein1d_derivative_basis(6, 0, v),
        bernstein1d_derivative_basis(6, 1, v),
        bernstein1d_derivative_basis(6, 2, v),
        bernstein1d_derivative_basis(6, 3, v),
    ];
    let mut out = [RectTaylorJet::zero(), RectTaylorJet::zero()];
    for a in 0..4 {
        for b in 0..4 {
            let derivative_sign = sign_u.powi(a as i32) * sign_v.powi(b as i32);
            let taylor_scale = derivative_sign / (small_factorial(a) * small_factorial(b));
            for (i, &bu_i) in bu[a].iter().enumerate() {
                for (j, &bv_j) in bv[b].iter().enumerate() {
                    let cp = tau.cps[i][j];
                    let weight = taylor_scale * bu_i * bv_j;
                    out[0].c[a][b] += weight * cp.x;
                    out[1].c[a][b] += weight * cp.y;
                }
            }
        }
    }
    out
}

fn kp20_planar_bi6_tau_patch_f3_energy(tau: &super::kp20_bi6_cap::Kp20Bi6Patch) -> f64 {
    let mut energy = 0.0_f64;
    for &(u, wu) in GAUSS7_01 {
        for &(v, wv) in GAUSS7_01 {
            let jet = planar_bi6_tau_jet3(tau, u, v);
            for a in 0..=3 {
                let b = 3 - a;
                let d = jet.d[a][b];
                energy += wu * wv * binom(3, a) * (d[0] * d[0] + d[1] * d[1]);
            }
        }
    }
    energy
}

fn orient_tau_jet(tau: Kp20TauJet3, sign_u: f64, sign_v: f64) -> Kp20TauJet3 {
    let mut out = tau;
    for a in 0..4 {
        for b in 0..(4 - a) {
            let sign = sign_u.powi(a as i32) * sign_v.powi(b as i32);
            out.d[a][b][0] *= sign;
            out.d[a][b][1] *= sign;
        }
    }
    out
}

fn planar_bi6_tau_jet3(tau: &super::kp20_bi6_cap::Kp20Bi6Patch, u: f64, v: f64) -> Kp20TauJet3 {
    let bu = [
        bernstein1d_derivative_basis(6, 0, u),
        bernstein1d_derivative_basis(6, 1, u),
        bernstein1d_derivative_basis(6, 2, u),
        bernstein1d_derivative_basis(6, 3, u),
    ];
    let bv = [
        bernstein1d_derivative_basis(6, 0, v),
        bernstein1d_derivative_basis(6, 1, v),
        bernstein1d_derivative_basis(6, 2, v),
        bernstein1d_derivative_basis(6, 3, v),
    ];
    let mut out = [[[0.0_f64; 2]; 4]; 4];
    for a in 0..4 {
        for b in 0..(4 - a) {
            for (i, &bu_i) in bu[a].iter().enumerate() {
                for (j, &bv_j) in bv[b].iter().enumerate() {
                    let w = bu_i * bv_j;
                    out[a][b][0] += w * tau.cps[i][j].x;
                    out[a][b][1] += w * tau.cps[i][j].y;
                }
            }
        }
    }
    Kp20TauJet3 { d: out }
}

fn bernstein1d_derivative_basis(degree: usize, order: usize, t: f64) -> [f64; 7] {
    let mut out = [0.0_f64; 7];
    if order > degree {
        return out;
    }
    let scale = falling_factorial(degree, order);
    let lower_degree = degree - order;
    for (i, out_i) in out.iter_mut().enumerate().take(degree + 1) {
        let start = i.saturating_sub(order);
        let end = i.min(lower_degree);
        let mut value = 0.0_f64;
        for r in start..=end {
            let sign = if ((order - (i - r)) & 1) == 0 {
                1.0
            } else {
                -1.0
            };
            value += sign * binom(order, i - r) * bernstein1d(lower_degree, r, t);
        }
        *out_i = scale * value;
    }
    out
}

fn bernstein1d(degree: usize, index: usize, t: f64) -> f64 {
    binom(degree, index) * t.powi(index as i32) * (1.0 - t).powi((degree - index) as i32)
}

fn binom(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    (0..k).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64)
}

fn multinomial3(n: usize, i: usize, j: usize, k: usize) -> f64 {
    small_factorial(n) / (small_factorial(i) * small_factorial(j) * small_factorial(k))
}

fn small_factorial(n: usize) -> f64 {
    (1..=n).fold(1.0, |acc, value| acc * value as f64)
}

fn falling_factorial(n: usize, order: usize) -> f64 {
    (0..order).fold(1.0, |acc, k| acc * (n - k) as f64)
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
        origin_strip_reversed: false,
        min_piece_center_separation: 0.0,
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
    reverse_second: bool,
) -> f64 {
    [0.0, 0.5, 1.0]
        .into_iter()
        .map(|t| {
            distance2(
                edge_point(first, first_edge, t),
                edge_point(
                    second,
                    second_edge,
                    if reverse_second { 1.0 - t } else { t },
                ),
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

fn min_center_separation(maps: [Kp20AffineGuideMap; 4]) -> f64 {
    let mut min_separation = f64::INFINITY;
    for a in 0..maps.len() {
        for b in (a + 1)..maps.len() {
            min_separation =
                min_separation.min(distance2(map_center(maps[a]), map_center(maps[b])));
        }
    }
    min_separation
}

fn map_center(map: Kp20AffineGuideMap) -> [f64; 2] {
    add2(add2(map.origin, scale2(map.du, 0.5)), scale2(map.dv, 0.5))
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
