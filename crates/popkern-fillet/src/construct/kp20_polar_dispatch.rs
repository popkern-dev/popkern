//! Gated KP20 polar macro-cap candidate construction for fillet dispatch.

use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20Bi6ValidationReport, Kp20PieceId, PolarBoundaryInput, Side, Spoke, SpokeEnd,
    TensorBorder, extract_bi3_tensor_border_from_polar_input, kp20_default_piece_map_domain_report,
};

use super::kp_polar_bi6_promote::SELF_INT_TOL;
use super::kp_polar_bi6_promote::extract_fillet_radius;
use super::kp20_polar_comparison::sampler_comparison_diagnostic;
use super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, summarize_candidate_intersections, summarize_intersections,
};
use super::kp20_polar_support_lift::Kp20SupportNormalLiftReport;
use super::kp20_polar_support_reparameterization::{
    Kp20CoupledDomainReport, Kp20SplitDomainReport, Kp20SupportDomainRefineReport,
    Kp20SupportDomainReport, Kp20SupportReparameterizationReport, Kp20SupportTransitionReport,
};
use super::three_quad_apply::extract_direct_boundary_input;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::ApexConstructionDiagnostic;
use coupled_lift::try_nonoverlap_coupled_lift_candidate;
use intersections::acceptance_screen;
use labels::{nonoverlap_coupled_lift_detail, piece_counts_label, validation_label};
use limited_lift::try_limited_nonoverlap_lift;
use origin_inner_row::{
    Kp20OriginInnerRowReport, Kp20PairedInnerRowReport, Kp20SupportInnerRowReport,
};
use primary_candidate::evaluate_positive_primary_candidates;
pub(crate) use primary_candidate::try_kp20_polar_macro_candidate;
use report::Kp20PolarDispatchReport;
use timing::Kp20DispatchStageTiming;

#[path = "kp20_polar_coupled_lift.rs"]
mod coupled_lift;
#[path = "kp20_polar_inner_field.rs"]
mod inner_field;
#[path = "kp20_polar_dispatch_intersections.rs"]
mod intersections;
#[path = "kp20_polar_dispatch_labels.rs"]
mod labels;
#[path = "kp20_polar_limited_lift.rs"]
mod limited_lift;
#[path = "kp20_polar_origin_inner_row.rs"]
mod origin_inner_row;
#[path = "kp20_polar_primary_candidate.rs"]
mod primary_candidate;
#[path = "kp20_polar_dispatch_report.rs"]
mod report;
#[path = "kp20_polar_support_domain.rs"]
mod support_domain;
#[path = "kp20_polar_support_transition.rs"]
mod support_transition;
#[path = "kp20_polar_dispatch_timing.rs"]
mod timing;

pub(super) fn try_kp20_polar_dispatch_with_diagnostics(
    sub_patches: &[GregoryPatch],
    fillet_surfaces: &[FilletSurfaceResult],
) -> Kp20PolarDispatchReport {
    let n = sub_patches.len();
    let mut diagnostics = Vec::new();
    let mut two_stripe = Vec::with_capacity(n);
    for (idx, patch) in sub_patches.iter().enumerate() {
        let GregoryPatch::RectangularTwoStripe(ts) = patch else {
            let detail = format!("sub-patch {idx}");
            record_decline(
                &mut diagnostics,
                "non_two_stripe_input",
                n,
                Some(detail),
                None,
            );
            return Kp20PolarDispatchReport::declined(diagnostics);
        };
        two_stripe.push((**ts).clone());
    }

    let spokes: Vec<Spoke> = (0..n)
        .map(|k| Spoke {
            a: SpokeEnd {
                patch_idx: k,
                side: Side::UOne,
                reversed: false,
            },
            b: SpokeEnd {
                patch_idx: (k + 1) % n,
                side: Side::VZero,
                reversed: true,
            },
        })
        .collect();
    let Some(mut direct) = extract_direct_boundary_input(&two_stripe, &spokes, None) else {
        record_decline(
            &mut diagnostics,
            "direct_boundary_unavailable",
            n,
            None,
            None,
        );
        return Kp20PolarDispatchReport::declined(diagnostics);
    };
    if fillet_surfaces.len() >= n {
        for (k, spoke) in direct.spokes.iter_mut().enumerate() {
            if let Some(radius) = extract_fillet_radius(&fillet_surfaces[k]) {
                spoke.radius_at_c = radius;
                spoke.radius_at_v = radius;
            }
        }
    }
    let polar = match PolarBoundaryInput::from_direct_boundary_input(&direct) {
        Ok(polar) => polar,
        Err(e) => {
            record_decline(&mut diagnostics, "polar_boundary_input", n, Some(e), None);
            return Kp20PolarDispatchReport::declined(diagnostics);
        }
    };
    let tensor_border = match extract_bi3_tensor_border_from_polar_input(&polar) {
        Ok(TensorBorder::Bi3(tensor_border)) => tensor_border,
        Ok(TensorBorder::Bi5(_)) => {
            record_decline(
                &mut diagnostics,
                "tensor_border_kind",
                n,
                Some("expected bi-3".into()),
                None,
            );
            return Kp20PolarDispatchReport::declined(diagnostics);
        }
        Err(e) => {
            record_decline(&mut diagnostics, "tensor_border_extract", n, Some(e), None);
            return Kp20PolarDispatchReport::declined(diagnostics);
        }
    };
    let (mut primary_diagnostics, primary_patches) =
        evaluate_positive_primary_candidates(n, polar.central_vertex, &tensor_border);
    diagnostics.append(&mut primary_diagnostics);
    if let Some(patches) = primary_patches {
        return Kp20PolarDispatchReport::accepted(patches, diagnostics);
    }
    let candidate = match try_kp20_polar_macro_candidate(polar.central_vertex, &tensor_border) {
        Ok(candidate) => candidate,
        Err(e) => {
            record_decline(&mut diagnostics, "candidate", n, Some(e), None);
            return Kp20PolarDispatchReport::declined(diagnostics);
        }
    };
    if candidate.len() != 4 * n {
        let detail = format!("expected {}, got {}", 4 * n, candidate.len());
        record_decline(
            &mut diagnostics,
            "output_length_mismatch",
            n,
            Some(detail),
            None,
        );
        return Kp20PolarDispatchReport::declined(diagnostics);
    }
    let map_report = kp20_default_piece_map_domain_report();
    if !map_report.stays_inside_triangle() {
        record_decline(
            &mut diagnostics,
            "piece_map_extrapolates",
            n,
            Some(format!(
                "outside_triangle_corners={};max_simplex_excess={:.3e};non_positive_orientation_maps={};{}",
                map_report.outside_triangle_corners,
                map_report.max_simplex_excess,
                map_report.non_positive_orientation_maps,
                map_report.offending_corners.join("|")
            )),
            Some(map_report.max_simplex_excess),
        );
    }
    if map_report.non_positive_orientation_maps > 0 {
        diagnostics.push(super::kp20_polar_diagnostics::dispatch_fact(
            "piece_map_orientation",
            n,
            Some(format!(
                "non_positive_orientation_maps={};{}",
                map_report.non_positive_orientation_maps,
                map_report.orientation_offenders.join("|")
            )),
            None,
            None,
        ));
    }
    let intersections = acceptance_screen(&candidate);
    if !intersections.is_empty() {
        let default_summary = summarize_intersections(&intersections);
        if n >= 6 {
            let mut bounded_summary = None;
            match try_limited_nonoverlap_lift(polar.central_vertex, &tensor_border) {
                Ok(lifted) => {
                    let lifted_intersections = acceptance_screen(&lifted.patches);
                    let lifted_summary =
                        summarize_candidate_intersections(&lifted.patches, &lifted_intersections);
                    diagnostics.push(super::kp20_polar_diagnostics::dispatch_fact(
                        "nonoverlap_limited_lift",
                        n,
                        Some(lifted.detail(&default_summary, &lifted_summary)),
                        Some(lifted_summary.min_jacobian_normal_cos),
                        Some(SELF_INT_TOL),
                    ));
                    if let Some(expanded_patches) = lifted.accepted_expanded_support_patches() {
                        let expanded_intersections = acceptance_screen(expanded_patches);
                        if expanded_intersections.is_empty() {
                            return Kp20PolarDispatchReport::accepted_expanded_support(
                                expanded_patches.to_vec(),
                                diagnostics,
                            );
                        }
                    }
                    if let Some(support_cell_patches) = lifted.accepted_support_cell_patches()
                        && support_cell_patches.len() % n == 0
                        && matches!(support_cell_patches.len() / n, 4 | 8 | 16)
                    {
                        if std::env::var_os("POPKERN_DETERMINISM_TRACE").is_some() {
                            let mut hash = 0xcbf29ce484222325u64;
                            for patch in support_cell_patches {
                                if let popkern_core::geometry::gregory::GregoryPatch::RectangularBiNonic(r) = patch {
                                    for row in &r.grid {
                                        for point in row {
                                            for value in [point.x, point.y, point.z] {
                                                hash ^= value.to_bits();
                                                hash = hash.wrapping_mul(0x100000001b3);
                                            }
                                        }
                                    }
                                }
                            }
                            eprintln!("dtrace: lift_support_cell_patches {hash:016x}");
                        }
                        return Kp20PolarDispatchReport::accepted_support_cell(
                            support_cell_patches.to_vec(),
                            lifted
                                .support_cell_floor_patches()
                                .map(|patches| patches.to_vec()),
                            diagnostics,
                        );
                    }
                    if lifted_intersections.is_empty() {
                        return Kp20PolarDispatchReport::accepted(lifted.patches, diagnostics);
                    }
                    let decline_reason = if lifted_summary.fold_count >= default_summary.fold_count
                        && lifted.report.optimizer_stats.accepted_candidates == 0
                    {
                        "nonoverlap_limited_lift_no_improving_candidate"
                    } else {
                        "nonoverlap_limited_lift_self_intersection"
                    };
                    record_candidate_decline(
                        &mut diagnostics,
                        decline_reason,
                        n,
                        &lifted_summary,
                        lifted.report.validation,
                    );
                    bounded_summary = Some(lifted_summary);
                }
                Err(error) => record_decline(
                    &mut diagnostics,
                    "nonoverlap_limited_lift",
                    n,
                    Some(error),
                    None,
                ),
            }
            record_lift_budget_decline(
                &mut diagnostics,
                n,
                &default_summary,
                bounded_summary.as_ref(),
            );
            return Kp20PolarDispatchReport::declined(diagnostics);
        }
        match try_nonoverlap_coupled_lift_candidate(polar.central_vertex, &tensor_border) {
            Ok(lifted) => {
                let lifted_intersections = acceptance_screen(&lifted.patches);
                let lifted_summary =
                    summarize_candidate_intersections(&lifted.patches, &lifted_intersections);
                diagnostics.push(nonoverlap_coupled_lift_fact(
                    n,
                    &default_summary,
                    &lifted_summary,
                    &lifted.report,
                    &lifted.origin_inner,
                    &lifted.paired_inner,
                    &lifted.support_inner,
                    &lifted.support_reparameterization,
                    &lifted.support_transition,
                    &lifted.support_domain,
                    &lifted.support_domain_refine,
                    &lifted.coupled_domain,
                    &lifted.split_domain,
                    lifted.timing,
                ));
                if lifted_intersections.is_empty() {
                    return Kp20PolarDispatchReport::accepted(lifted.patches, diagnostics);
                }
                record_candidate_decline(
                    &mut diagnostics,
                    "nonoverlap_coupled_lift_self_intersection",
                    n,
                    &lifted_summary,
                    lifted.report.validation,
                );
                if std::env::var_os("POPKERN_KP20_POLAR_COMPARE").is_some() {
                    diagnostics.push(sampler_comparison_diagnostic(
                        polar.central_vertex,
                        &tensor_border,
                        default_summary.clone(),
                    ));
                }
                for (idx, detail) in lifted_intersections {
                    record_self_intersection(&mut diagnostics, idx, &detail, n);
                }
                return Kp20PolarDispatchReport::declined(diagnostics);
            }
            Err(e) => record_decline(
                &mut diagnostics,
                "nonoverlap_coupled_lift",
                n,
                Some(e),
                None,
            ),
        }
        if map_report.non_positive_orientation_maps == 0 {
            diagnostics.push(super::kp20_polar_diagnostics::dispatch_fact(
                "piece_map_orientation",
                n,
                Some("non_positive_orientation_maps=0".to_string()),
                None,
                None,
            ));
        }
        for (idx, detail) in intersections {
            record_self_intersection(&mut diagnostics, idx, &detail, n);
        }
        return Kp20PolarDispatchReport::declined(diagnostics);
    }
    Kp20PolarDispatchReport::accepted(candidate, diagnostics)
}

#[allow(clippy::too_many_arguments)]
fn nonoverlap_coupled_lift_fact(
    n: usize,
    default: &Kp20SamplerFoldSummary,
    lifted: &Kp20SamplerFoldSummary,
    report: &Kp20SupportNormalLiftReport,
    origin_inner: &Kp20OriginInnerRowReport,
    paired_inner: &Kp20PairedInnerRowReport,
    support_inner: &Kp20SupportInnerRowReport,
    support_reparameterization: &Kp20SupportReparameterizationReport,
    support_transition: &Kp20SupportTransitionReport,
    support_domain: &Kp20SupportDomainReport,
    support_domain_refine: &Kp20SupportDomainRefineReport,
    coupled_domain: &Kp20CoupledDomainReport,
    split_domain: &Kp20SplitDomainReport,
    timing: Kp20DispatchStageTiming,
) -> ApexConstructionDiagnostic {
    super::kp20_polar_diagnostics::dispatch_fact(
        "nonoverlap_coupled_lift",
        n,
        Some(nonoverlap_coupled_lift_detail(
            default,
            lifted,
            report,
            origin_inner,
            paired_inner,
            support_inner,
            support_reparameterization,
            support_transition,
            support_domain,
            support_domain_refine,
            coupled_domain,
            split_domain,
            timing,
        )),
        Some(lifted.min_jacobian_normal_cos),
        Some(SELF_INT_TOL),
    )
}

fn record_candidate_decline(
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
    reason_code: &str,
    n_sides: usize,
    summary: &Kp20SamplerFoldSummary,
    validation: Kp20Bi6ValidationReport,
) {
    record_decline(
        diagnostics,
        reason_code,
        n_sides,
        Some(format!(
            "folds={};piece_folds={};{}",
            summary.fold_count,
            piece_counts_label(summary.piece_fold_counts),
            validation_label(validation),
        )),
        Some(summary.min_jacobian_normal_cos),
    );
}

fn record_self_intersection(
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
    idx: usize,
    detail: &popkern_core::geometry::gregory::GregorySelfIntersectionReport,
    n: usize,
) {
    let label = detail.diagnostic_label();
    let sector = idx / 4;
    let piece = Kp20PieceId::all()[idx % 4];
    record_decline(
        diagnostics,
        "self_intersection",
        n,
        Some(format!(
            "sector={sector};piece={piece:?};patch_index={idx};{label}"
        )),
        Some(detail.min_jacobian_normal_cos),
    );
}

fn record_decline(
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
    reason_code: &str,
    n_sides: usize,
    detail: Option<String>,
    residual: Option<f64>,
) {
    diagnostics.push(super::kp20_polar_diagnostics::dispatch_declined(
        reason_code,
        n_sides,
        detail,
        residual,
        Some(SELF_INT_TOL),
    ));
}

fn record_lift_budget_decline(
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
    n: usize,
    default_summary: &Kp20SamplerFoldSummary,
    bounded_summary: Option<&Kp20SamplerFoldSummary>,
) {
    let bounded_detail = bounded_summary
        .map(|summary| {
            format!(
                ";bounded_folds={};bounded_piece_folds={}",
                summary.fold_count,
                piece_counts_label(summary.piece_fold_counts),
            )
        })
        .unwrap_or_default();
    record_decline(
        diagnostics,
        "nonoverlap_coupled_lift_budget_exceeded",
        n,
        Some(format!(
            "folds={};piece_folds={};lift_solver_budget=limited_lift_only_for_valence_ge_6{}",
            default_summary.fold_count,
            piece_counts_label(default_summary.piece_fold_counts),
            bounded_detail,
        )),
        Some(default_summary.min_jacobian_normal_cos),
    );
}
#[cfg(test)]
#[path = "kp20_polar_dispatch_tests.rs"]
mod tests;
