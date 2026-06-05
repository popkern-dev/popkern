//! Gated KP20 polar macro-cap candidate construction for fillet dispatch.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20Bi6CapInput, Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport,
    Kp20GuideBuildInput, Kp20PieceId, PolarBoundaryInput, RectangularGregoryBiNonic, Side, Spoke,
    SpokeEnd, TensorBorder, TensorBorderBi3, assemble_kp20_bi6_cap, assemble_kp20_bi6_macro_cap,
    build_kp20_guide_ring, elevate_bi6_grid_to_bi9, extract_bi3_tensor_border_from_polar_input,
    kp20_default_piece_map_domain_report, kp20_nonoverlap_coupled_sector_piece_maps,
    measure_kp20_bi6_macro_cap, sample_kp20_macro_patch_from_piece_maps,
};

use super::kp_polar_bi6_promote::SELF_INT_TOL;
use super::kp_polar_bi6_promote::extract_fillet_radius;
use super::kp20_polar_comparison::sampler_comparison_diagnostic;
use super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, gregory_patches_from_macro_patches, summarize_candidate_intersections,
    summarize_intersections, validation_is_clean,
};
use super::kp20_polar_support_lift::{
    Kp20SupportNormalLiftReport, lift_candidate_from_macro_patches,
};
use super::three_quad_apply::extract_direct_boundary_input;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::ApexConstructionDiagnostic;
use labels::{fold_ownership_label, piece_counts_label, validation_label};
use origin_inner_row::{
    Kp20OriginInnerRowReport, Kp20PairedInnerRowReport, refine_origin_inner_interior_row,
    refine_paired_inner_interior_row,
};

#[path = "kp20_polar_dispatch_labels.rs"]
mod labels;
#[path = "kp20_polar_origin_inner_row.rs"]
mod origin_inner_row;

pub(super) struct Kp20PolarDispatchReport {
    pub patches: Option<Vec<GregoryPatch>>,
    pub diagnostics: Vec<ApexConstructionDiagnostic>,
}

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
            trace("non-two-stripe input", n, Some(detail.clone()));
            record_decline(
                &mut diagnostics,
                "non_two_stripe_input",
                n,
                Some(detail),
                None,
            );
            return Kp20PolarDispatchReport {
                patches: None,
                diagnostics,
            };
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
        trace("extract_direct_boundary_input", n, None);
        record_decline(
            &mut diagnostics,
            "direct_boundary_unavailable",
            n,
            None,
            None,
        );
        return Kp20PolarDispatchReport {
            patches: None,
            diagnostics,
        };
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
            trace("from_direct_boundary_input", n, Some(e.clone()));
            record_decline(&mut diagnostics, "polar_boundary_input", n, Some(e), None);
            return Kp20PolarDispatchReport {
                patches: None,
                diagnostics,
            };
        }
    };
    let tensor_border = match extract_bi3_tensor_border_from_polar_input(&polar) {
        Ok(TensorBorder::Bi3(tensor_border)) => tensor_border,
        Ok(TensorBorder::Bi5(_)) => {
            trace("tensor-border kind", n, Some("expected bi-3".into()));
            record_decline(
                &mut diagnostics,
                "tensor_border_kind",
                n,
                Some("expected bi-3".into()),
                None,
            );
            return Kp20PolarDispatchReport {
                patches: None,
                diagnostics,
            };
        }
        Err(e) => {
            trace("tensor-border extract", n, Some(e.clone()));
            record_decline(&mut diagnostics, "tensor_border_extract", n, Some(e), None);
            return Kp20PolarDispatchReport {
                patches: None,
                diagnostics,
            };
        }
    };
    let candidate = match try_kp20_polar_macro_candidate(polar.central_vertex, &tensor_border) {
        Ok(candidate) => candidate,
        Err(e) => {
            trace("candidate", n, Some(e.clone()));
            record_decline(&mut diagnostics, "candidate", n, Some(e), None);
            return Kp20PolarDispatchReport {
                patches: None,
                diagnostics,
            };
        }
    };
    if candidate.len() != 4 * n {
        let detail = format!("expected {}, got {}", 4 * n, candidate.len());
        trace("output length mismatch", n, Some(detail.clone()));
        record_decline(
            &mut diagnostics,
            "output_length_mismatch",
            n,
            Some(detail),
            None,
        );
        return Kp20PolarDispatchReport {
            patches: None,
            diagnostics,
        };
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
    let intersections = self_intersecting_patches(&candidate);
    if !intersections.is_empty() {
        let default_summary = summarize_intersections(&intersections);
        match try_nonoverlap_coupled_lift_candidate(polar.central_vertex, &tensor_border) {
            Ok(lifted) => {
                let lifted_intersections = self_intersecting_patches(&lifted.patches);
                let lifted_summary =
                    summarize_candidate_intersections(&lifted.patches, &lifted_intersections);
                diagnostics.push(nonoverlap_coupled_lift_fact(
                    n,
                    &default_summary,
                    &lifted_summary,
                    &lifted.report,
                    &lifted.origin_inner,
                    &lifted.paired_inner,
                ));
                if lifted_intersections.is_empty() {
                    trace(
                        "accepted nonoverlap_coupled_lift",
                        n,
                        Some("4N bi-9 patches".into()),
                    );
                    return Kp20PolarDispatchReport {
                        patches: Some(lifted.patches),
                        diagnostics,
                    };
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
                return Kp20PolarDispatchReport {
                    patches: None,
                    diagnostics,
                };
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
        return Kp20PolarDispatchReport {
            patches: None,
            diagnostics,
        };
    }
    trace("accepted", n, Some(format!("4N={} bi-9 patches", 4 * n)));
    Kp20PolarDispatchReport {
        patches: Some(candidate),
        diagnostics,
    }
}

struct Kp20LiftedDispatchCandidate {
    patches: Vec<GregoryPatch>,
    report: Kp20SupportNormalLiftReport,
    origin_inner: Kp20OriginInnerRowReport,
    paired_inner: Kp20PairedInnerRowReport,
}

fn try_nonoverlap_coupled_lift_candidate(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<Kp20LiftedDispatchCandidate, String> {
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let sampled_macro_patches = guide
        .sectors
        .iter()
        .map(|sector| {
            sample_kp20_macro_patch_from_piece_maps(
                sector,
                kp20_nonoverlap_coupled_sector_piece_maps(),
            )
        })
        .collect::<Vec<_>>();
    let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches,
    })
    .map_err(|e| e.to_string())?;
    let report = lift_candidate_from_macro_patches(star_point, tensor_border, cap.macro_patches)?;
    let (report, origin_inner) =
        refine_origin_inner_interior_row(star_point, tensor_border, report)?;
    let (report, paired_inner) =
        refine_paired_inner_interior_row(star_point, tensor_border, report)?;
    if !validation_is_clean(report.validation) {
        return Err(format!(
            "validation residuals exceed tolerance: {}",
            validation_label(report.validation)
        ));
    }
    let patches = gregory_patches_from_macro_patches(&report.macro_patches);
    Ok(Kp20LiftedDispatchCandidate {
        patches,
        report,
        origin_inner,
        paired_inner,
    })
}

fn nonoverlap_coupled_lift_fact(
    n: usize,
    default: &Kp20SamplerFoldSummary,
    lifted: &Kp20SamplerFoldSummary,
    report: &Kp20SupportNormalLiftReport,
    origin_inner: &Kp20OriginInnerRowReport,
    paired_inner: &Kp20PairedInnerRowReport,
) -> ApexConstructionDiagnostic {
    super::kp20_polar_diagnostics::dispatch_fact(
        "nonoverlap_coupled_lift",
        n,
        Some(format!(
            "status=sampled;default_folds={};folds={};piece_folds={};{};{};origin_inner_evaluated={};origin_inner_accepted={};origin_inner_rejected={};origin_inner_start_folds={};origin_inner_end_folds={};origin_inner_mode={};origin_inner_scale={:.3e};paired_inner_evaluated={};paired_inner_accepted={};paired_inner_rejected={};paired_inner_start_folds={};paired_inner_end_folds={};paired_inner_mode={};paired_inner_scale={:.3e};optimizer_active_fold_sites={};optimizer_row_variables={};optimizer_evaluated={};optimizer_accepted={};optimizer_rejected={};accepted_steps={}",
            default.fold_count,
            lifted.fold_count,
            piece_counts_label(lifted.piece_fold_counts),
            fold_ownership_label(lifted),
            validation_label(report.validation),
            origin_inner.evaluated,
            origin_inner.accepted,
            origin_inner.rejected,
            origin_inner.start_folds,
            origin_inner.end_folds,
            origin_inner.mode,
            origin_inner.scale,
            paired_inner.evaluated,
            paired_inner.accepted,
            paired_inner.rejected,
            paired_inner.start_folds,
            paired_inner.end_folds,
            paired_inner.mode,
            paired_inner.scale,
            report.optimizer_stats.active_fold_sites,
            report.optimizer_stats.row_variable_count,
            report.optimizer_stats.evaluated_candidates,
            report.optimizer_stats.accepted_candidates,
            report.optimizer_stats.rejected_candidates,
            report.accepted_steps.len(),
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

pub(super) fn try_kp20_polar_macro_candidate(
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

pub(super) fn self_intersecting_patches(
    candidate: &[GregoryPatch],
) -> Vec<(
    usize,
    popkern_core::geometry::gregory::GregorySelfIntersectionReport,
)> {
    candidate
        .iter()
        .enumerate()
        .filter_map(|(idx, gp)| {
            let report = gp.self_intersection_report(SELF_INT_TOL);
            report.intersects().then_some((idx, report))
        })
        .collect()
}

fn trace(stage: &str, n: usize, detail: Option<String>) {
    if std::env::var_os("POPKERN_KP20_POLAR_TRACE").is_some() {
        match detail {
            Some(d) => eprintln!("[kp20_polar N={n}] {stage}: {d}"),
            None => eprintln!("[kp20_polar N={n}] {stage}"),
        }
    }
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
    trace(
        "self-intersection",
        n,
        Some(format!(
            "sector {sector}; piece {piece:?}; patch index {idx}; {label}"
        )),
    );
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
