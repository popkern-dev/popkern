//! Primary KP20 polar macro candidates.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20AffineGuideMap, Kp20Bi6CapInput, Kp20Bi6MacroCapInput,
    Kp20Bi6ValidationReport, Kp20GuideBuildInput, Kp20PatchCorner, Kp20PieceCornerMaps,
    Kp20PieceId, Kp20PieceMapDomainReport, Kp20PieceMapStitchReport, Kp20SectorPieceMaps,
    RectangularGregoryBiNonic, TensorBorderBi3, assemble_kp20_bi6_cap, assemble_kp20_bi6_macro_cap,
    build_kp20_guide_ring, elevate_bi6_grid_to_bi9, kp20_piece_map_domain_report,
    kp20_piece_map_stitch_report, kp20_positive_nonoverlap_c2_sector_piece_maps,
    measure_kp20_bi6_macro_cap, sample_kp20_macro_patch_from_piece_maps,
};

use super::super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, candidate_summary_from_cap,
    fast_candidate_summary_from_changed_patches, gregory_patches_from_macro_patches,
    is_better_summary, validation_is_clean,
};
use super::labels::{fold_ownership_label, piece_counts_label, validation_label};
use crate::types::ApexConstructionDiagnostic;

const TARGET_PIECE: Kp20PieceId = Kp20PieceId::NearMidpointLeft;

pub(super) fn evaluate_positive_primary_candidates(
    n: usize,
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> (Vec<ApexConstructionDiagnostic>, Option<Vec<GregoryPatch>>) {
    let mut diagnostics = Vec::new();
    let mut baseline = None;
    match try_positive_nonoverlap_c2_macro_candidate(star_point, tensor_border) {
        Ok(mapped) => {
            diagnostics.push(primary_fact(
                n,
                "positive_nonoverlap_c2_primary",
                mapped.detail(),
                Some(mapped.summary.min_jacobian_normal_cos),
            ));
            if mapped.is_acceptable() {
                return (diagnostics, Some(mapped.patches));
            }
            baseline = Some(mapped.summary);
        }
        Err(e) => diagnostics.push(primary_decline(n, "positive_nonoverlap_c2_primary", e)),
    }
    match try_support_domain_primary_candidate(star_point, tensor_border, baseline.as_ref()) {
        Ok(mapped) => {
            diagnostics.push(primary_fact(
                n,
                "positive_nonoverlap_c2_support_domain_primary",
                mapped.detail(),
                Some(mapped.summary.min_jacobian_normal_cos),
            ));
            if mapped.is_acceptable() {
                return (diagnostics, mapped.patches);
            }
        }
        Err(e) => diagnostics.push(primary_decline(
            n,
            "positive_nonoverlap_c2_support_domain_primary",
            e,
        )),
    }
    (diagnostics, None)
}

pub(super) struct Kp20PrimaryMappedCandidate {
    pub(super) patches: Vec<GregoryPatch>,
    pub(super) summary: Kp20SamplerFoldSummary,
    pub(super) validation: Kp20Bi6ValidationReport,
    domain: Kp20PieceMapDomainReport,
    stitch: Kp20PieceMapStitchReport,
}

impl Kp20PrimaryMappedCandidate {
    pub(super) fn is_acceptable(&self) -> bool {
        validation_is_clean(self.validation) && self.summary.fold_count == 0
    }

    pub(super) fn detail(&self) -> String {
        format!(
            "map_family=positive_nonoverlap_c2;status=sampled;folds={};piece_folds={};{};validation:{};inside_triangle={};non_positive_orientation_maps={};c2_edge={:.3e};origin_strip_edge={:.3e};origin_strip_reversed={};min_center_sep={:.3e}",
            self.summary.fold_count,
            piece_counts_label(self.summary.piece_fold_counts),
            fold_ownership_label(&self.summary),
            validation_label(self.validation),
            self.domain.stays_inside_triangle(),
            self.domain.non_positive_orientation_maps,
            self.stitch.max_internal_c2_edge_mismatch,
            self.stitch.max_origin_strip_edge_mismatch,
            self.stitch.origin_strip_reversed,
            self.stitch.min_piece_center_separation,
        )
    }
}

struct Kp20SupportDomainPrimaryCandidate {
    patches: Option<Vec<GregoryPatch>>,
    summary: Kp20SamplerFoldSummary,
    validation: Kp20Bi6ValidationReport,
    domain: Kp20PieceMapDomainReport,
    stitch: Kp20PieceMapStitchReport,
    params: Kp20PrimarySupportDomainParams,
    evaluated: usize,
    accepted: usize,
    rejected: usize,
}

impl Kp20SupportDomainPrimaryCandidate {
    fn is_acceptable(&self) -> bool {
        self.patches.is_some()
    }

    fn detail(&self) -> String {
        format!(
            "map_family=positive_nonoverlap_c2;fold_check={};status=sampled;folds={};piece_folds={};{};validation:{};inside_triangle={};non_positive_orientation_maps={};c2_edge={:.3e};origin_strip_edge={:.3e};origin_strip_reversed={};min_center_sep={:.3e};support_domain_evaluated={};support_domain_accepted={};support_domain_rejected={};support_domain_target_piece={:?};support_domain_origin=({:.3e},{:.3e});support_domain_span=({:.3e},{:.3e});support_domain_skew={:.3e};decline_reason={}",
            if self.is_acceptable() {
                "exact"
            } else {
                "changed_piece_jacobian"
            },
            self.summary.fold_count,
            piece_counts_label(self.summary.piece_fold_counts),
            fold_ownership_label(&self.summary),
            validation_label(self.validation),
            self.domain.stays_inside_triangle(),
            self.domain.non_positive_orientation_maps,
            self.stitch.max_internal_c2_edge_mismatch,
            self.stitch.max_origin_strip_edge_mismatch,
            self.stitch.origin_strip_reversed,
            self.stitch.min_piece_center_separation,
            self.evaluated,
            self.accepted,
            self.rejected,
            TARGET_PIECE,
            self.params.origin[0],
            self.params.origin[1],
            self.params.span_u,
            self.params.span_v,
            self.params.skew,
            if self.is_acceptable() {
                "none"
            } else {
                "support_domain_primary_fold_remaining"
            },
        )
    }
}

pub(super) fn try_positive_nonoverlap_c2_macro_candidate(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<Kp20PrimaryMappedCandidate, String> {
    mapped_candidate_from_maps(
        star_point,
        tensor_border,
        kp20_positive_nonoverlap_c2_sector_piece_maps(),
    )
}

fn try_support_domain_primary_candidate(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    baseline: Option<&Kp20SamplerFoldSummary>,
) -> Result<Kp20SupportDomainPrimaryCandidate, String> {
    let Some(baseline) = baseline else {
        return Err(
            "positive_nonoverlap_c2 support-domain primary missing baseline summary".to_string(),
        );
    };
    let mut best: Option<(
        Kp20Bi6MacroCapInput,
        Kp20PieceMapDomainReport,
        Kp20PieceMapStitchReport,
        (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
        Kp20PrimarySupportDomainParams,
    )> = None;
    let mut evaluated = 0;
    for params in Kp20PrimarySupportDomainParams::family() {
        evaluated += 1;
        let sampled = sampled_candidate_from_maps(
            star_point,
            tensor_border,
            support_domain_piece_maps(params),
        )?;
        let changed = changed_target_pieces(tensor_border.sectors.len());
        let fast = fast_candidate_summary_from_changed_patches(
            tensor_border,
            sampled.cap.sampled_macro_patches.clone(),
            baseline,
            &changed,
            baseline.fold_count,
        )?;
        if validation_is_clean(fast.1) && fast.0.fold_count == 0 {
            let candidate = mapped_candidate_from_sampled(
                tensor_border,
                sampled.cap,
                sampled.domain,
                sampled.stitch,
            )?;
            if !candidate.is_acceptable() {
                continue;
            }
            return Ok(Kp20SupportDomainPrimaryCandidate {
                patches: Some(candidate.patches),
                summary: candidate.summary,
                validation: candidate.validation,
                domain: candidate.domain,
                stitch: candidate.stitch,
                params,
                evaluated,
                accepted: 1,
                rejected: evaluated - 1,
            });
        }
        let best_summary = best
            .as_ref()
            .map(|(_, _, _, summary, _)| (summary.0.clone(), summary.1));
        if is_better_summary(&fast, best_summary.as_ref()) {
            best = Some((sampled.cap, sampled.domain, sampled.stitch, fast, params));
        }
    }
    let Some((cap, domain, stitch, _, params)) = best else {
        return Err("positive_nonoverlap_c2 support-domain primary had no candidates".to_string());
    };
    let (summary, validation) = fast_candidate_summary_from_changed_patches(
        tensor_border,
        cap.sampled_macro_patches,
        baseline,
        &changed_target_pieces(tensor_border.sectors.len()),
        baseline.fold_count,
    )?;
    Ok(Kp20SupportDomainPrimaryCandidate {
        patches: None,
        summary,
        validation,
        domain,
        stitch,
        params,
        evaluated,
        accepted: 0,
        rejected: evaluated,
    })
}

fn mapped_candidate_from_maps(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    maps: Kp20SectorPieceMaps,
) -> Result<Kp20PrimaryMappedCandidate, String> {
    let sampled = sampled_candidate_from_maps(star_point, tensor_border, maps)?;
    mapped_candidate_from_sampled(tensor_border, sampled.cap, sampled.domain, sampled.stitch)
}

struct Kp20SampledMappedCandidate {
    cap: Kp20Bi6MacroCapInput,
    domain: Kp20PieceMapDomainReport,
    stitch: Kp20PieceMapStitchReport,
}

fn sampled_candidate_from_maps(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    maps: Kp20SectorPieceMaps,
) -> Result<Kp20SampledMappedCandidate, String> {
    let domain = kp20_piece_map_domain_report(maps);
    if !domain.stays_inside_triangle() || domain.non_positive_orientation_maps > 0 {
        return Err(format!(
            "positive_nonoverlap_c2 map invalid: inside_triangle={} non_positive_orientation_maps={}",
            domain.stays_inside_triangle(),
            domain.non_positive_orientation_maps,
        ));
    }
    let stitch = kp20_piece_map_stitch_report(maps);
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let sampled_macro_patches = guide
        .sectors
        .iter()
        .map(|sector| sample_kp20_macro_patch_from_piece_maps(sector, maps))
        .collect::<Vec<_>>();
    let cap = Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches,
    };
    Ok(Kp20SampledMappedCandidate {
        cap,
        domain,
        stitch,
    })
}

fn mapped_candidate_from_sampled(
    tensor_border: &TensorBorderBi3,
    cap: Kp20Bi6MacroCapInput,
    domain: Kp20PieceMapDomainReport,
    stitch: Kp20PieceMapStitchReport,
) -> Result<Kp20PrimaryMappedCandidate, String> {
    let cap = assemble_kp20_bi6_macro_cap(&cap).map_err(|e| e.to_string())?;
    let (summary, validation) =
        candidate_summary_from_cap(tensor_border, cap.macro_patches.clone())?;
    Ok(Kp20PrimaryMappedCandidate {
        patches: gregory_patches_from_macro_patches(&cap.macro_patches),
        summary,
        validation,
        domain,
        stitch,
    })
}

fn changed_target_pieces(n: usize) -> Vec<[bool; 4]> {
    let mut changed = vec![[false; 4]; n];
    for pieces in &mut changed {
        pieces[TARGET_PIECE.index()] = true;
    }
    changed
}

#[derive(Debug, Clone, Copy)]
struct Kp20PrimarySupportDomainParams {
    origin: [f64; 2],
    span_u: f64,
    span_v: f64,
    skew: f64,
}

impl Kp20PrimarySupportDomainParams {
    const fn family() -> [Self; 8] {
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
}

fn support_domain_piece_maps(params: Kp20PrimarySupportDomainParams) -> Kp20SectorPieceMaps {
    let mut maps = kp20_positive_nonoverlap_c2_sector_piece_maps();
    for (piece, corner_maps) in &mut maps {
        if *piece == TARGET_PIECE {
            *corner_maps = support_domain_corner_maps(params);
        }
    }
    maps
}

fn support_domain_corner_maps(params: Kp20PrimarySupportDomainParams) -> Kp20PieceCornerMaps {
    square_maps(
        params.origin,
        [params.span_u, params.skew],
        [params.skew, params.span_v],
    )
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

fn primary_fact(
    n: usize,
    reason_code: &str,
    detail: String,
    residual: Option<f64>,
) -> ApexConstructionDiagnostic {
    super::super::kp20_polar_diagnostics::dispatch_fact(
        reason_code,
        n,
        Some(detail),
        residual,
        Some(super::super::kp_polar_bi6_promote::SELF_INT_TOL),
    )
}

fn primary_decline(n: usize, reason_code: &str, detail: String) -> ApexConstructionDiagnostic {
    super::super::kp20_polar_diagnostics::dispatch_declined(
        reason_code,
        n,
        Some(detail),
        None,
        Some(super::super::kp_polar_bi6_promote::SELF_INT_TOL),
    )
}

pub(crate) fn try_kp20_polar_macro_candidate(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<Vec<GregoryPatch>, String> {
    let cap = assemble_kp20_bi6_cap(&Kp20Bi6CapInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let report = measure_kp20_bi6_macro_cap(&cap, tensor_border).map_err(|e| e.to_string())?;
    if !report.passes(1e-8) {
        return Err(format!("validation residuals exceed tolerance: {report:?}"));
    }

    Ok(cap
        .patches_bi6
        .into_iter()
        .map(elevate_bi6_grid_to_bi9)
        .map(|grid| GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid })))
        .collect())
}
