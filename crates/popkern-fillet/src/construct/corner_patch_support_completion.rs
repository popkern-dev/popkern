//! Result-level completion for Gregory cap/support seam ownership.

use popkern_core::topology::BRep;

use super::kp20_polar_support_cell_shared_trim_quality as shared_trim_quality;
use crate::types::{ApexConstructionDiagnostic, CornerPatchInfo, CornerSubPatchRole};

#[path = "kp20_polar_support_cell_boundary_reconcile.rs"]
mod support_cell_boundary_reconcile;
#[path = "kp20_polar_support_cell_completion_apply.rs"]
mod support_cell_completion_apply;
#[path = "kp20_polar_support_cell_pcurve_author.rs"]
mod support_cell_pcurve_author;

const BOUNDARY_SAFE_G0_PREFERENCE_RELATIVE_BAND: f64 = 5.0e-3;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SupportCellCompletionCandidate {
    PcurveOnly,
    BoundaryOnly,
    PcurveThenBoundary,
    PcurveThenBoundaryIfG0Safe,
    PcurveThenWorstBoundaryIfG0Safe,
    PcurveThenTop3BoundariesIfG0Safe,
    PcurveThenRankedBoundaryPrefixIfG0Safe,
    PcurveThenWorstBoundaryNeighborhoodIfG0Safe,
    PcurveThenWorstEndpointClusterIfG0Safe,
    PcurveThenRankedBoundaryPrefixThenWorstEndpointClusterIfG0Safe,
    PcurveThenRankedBoundaryEndpointClusterThenRowGroupIfG0Safe,
    PcurveThenRankedBoundaryEndpointClusterThenTwoRowGroupsIfG0Safe,
    PcurveThenRankedBoundaryEndpointRowGroupThenTangentIfG0Safe,
    PcurveThenWorstBoundaryRowGroupIfG0Safe,
}

impl SupportCellCompletionCandidate {
    const ALL: [Self; 14] = [
        Self::PcurveOnly,
        Self::BoundaryOnly,
        Self::PcurveThenBoundary,
        Self::PcurveThenBoundaryIfG0Safe,
        Self::PcurveThenWorstBoundaryIfG0Safe,
        Self::PcurveThenTop3BoundariesIfG0Safe,
        Self::PcurveThenRankedBoundaryPrefixIfG0Safe,
        Self::PcurveThenWorstBoundaryNeighborhoodIfG0Safe,
        Self::PcurveThenWorstEndpointClusterIfG0Safe,
        Self::PcurveThenRankedBoundaryPrefixThenWorstEndpointClusterIfG0Safe,
        Self::PcurveThenRankedBoundaryEndpointClusterThenRowGroupIfG0Safe,
        Self::PcurveThenRankedBoundaryEndpointClusterThenTwoRowGroupsIfG0Safe,
        Self::PcurveThenRankedBoundaryEndpointRowGroupThenTangentIfG0Safe,
        Self::PcurveThenWorstBoundaryRowGroupIfG0Safe,
    ];

    fn label(self) -> &'static str {
        match self {
            Self::PcurveOnly => "pcurve_only",
            Self::BoundaryOnly => "boundary_only",
            Self::PcurveThenBoundary => "pcurve_then_boundary",
            Self::PcurveThenBoundaryIfG0Safe => "pcurve_then_boundary_if_g0_safe",
            Self::PcurveThenWorstBoundaryIfG0Safe => "pcurve_then_worst_boundary_if_g0_safe",
            Self::PcurveThenTop3BoundariesIfG0Safe => "pcurve_then_top3_boundaries_if_g0_safe",
            Self::PcurveThenRankedBoundaryPrefixIfG0Safe => {
                "pcurve_then_ranked_boundary_prefix_if_g0_safe"
            }
            Self::PcurveThenWorstBoundaryNeighborhoodIfG0Safe => {
                "pcurve_then_worst_boundary_neighborhood_if_g0_safe"
            }
            Self::PcurveThenWorstEndpointClusterIfG0Safe => {
                "pcurve_then_worst_endpoint_cluster_if_g0_safe"
            }
            Self::PcurveThenRankedBoundaryPrefixThenWorstEndpointClusterIfG0Safe => {
                "pcurve_then_ranked_boundary_prefix_then_worst_endpoint_cluster_if_g0_safe"
            }
            Self::PcurveThenRankedBoundaryEndpointClusterThenRowGroupIfG0Safe => {
                "pcurve_then_ranked_boundary_endpoint_cluster_then_row_group_if_g0_safe"
            }
            Self::PcurveThenRankedBoundaryEndpointClusterThenTwoRowGroupsIfG0Safe => {
                "pcurve_then_ranked_boundary_endpoint_cluster_then_two_row_groups_if_g0_safe"
            }
            Self::PcurveThenRankedBoundaryEndpointRowGroupThenTangentIfG0Safe => {
                "pcurve_then_ranked_boundary_endpoint_row_group_then_tangent_if_g0_safe"
            }
            Self::PcurveThenWorstBoundaryRowGroupIfG0Safe => {
                "pcurve_then_worst_boundary_row_group_if_g0_safe"
            }
        }
    }

    fn is_boundary_safe(self) -> bool {
        matches!(
            self,
            Self::PcurveThenBoundaryIfG0Safe
                | Self::PcurveThenWorstBoundaryIfG0Safe
                | Self::PcurveThenTop3BoundariesIfG0Safe
                | Self::PcurveThenRankedBoundaryPrefixIfG0Safe
                | Self::PcurveThenWorstBoundaryNeighborhoodIfG0Safe
                | Self::PcurveThenWorstEndpointClusterIfG0Safe
                | Self::PcurveThenRankedBoundaryPrefixThenWorstEndpointClusterIfG0Safe
                | Self::PcurveThenRankedBoundaryEndpointClusterThenRowGroupIfG0Safe
                | Self::PcurveThenRankedBoundaryEndpointClusterThenTwoRowGroupsIfG0Safe
                | Self::PcurveThenRankedBoundaryEndpointRowGroupThenTangentIfG0Safe
                | Self::PcurveThenWorstBoundaryRowGroupIfG0Safe
        )
    }
}

#[derive(Clone)]
struct SupportCellCompletionSelection {
    brep: BRep,
    candidate: SupportCellCompletionCandidate,
    quality: shared_trim_quality::Kp20SupportCellSharedTrimQuality,
}

#[derive(Clone, Copy)]
struct SupportCellCompletionCandidateScore {
    candidate: SupportCellCompletionCandidate,
    status: &'static str,
    operation_status: &'static str,
    quality: shared_trim_quality::Kp20SupportCellSharedTrimQuality,
}

#[derive(Clone, Copy)]
pub(super) struct SupportCellCompletionApplication {
    operator_succeeded: bool,
    operation_status: &'static str,
}

pub(crate) fn complete_corner_patch_support_pcurves(
    brep: &mut BRep,
    corner_patches: &[CornerPatchInfo],
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
) {
    complete_support_cell_pcurves_transactionally(brep, corner_patches, diagnostics);
    let mut completed_cap_faces = Vec::new();
    for patch in corner_patches {
        if patch
            .sub_patch_provenance
            .is_some_and(|provenance| provenance.role == CornerSubPatchRole::SupportCell)
        {
            continue;
        }
        if completed_cap_faces.contains(&patch.face) {
            continue;
        }
        super::three_quad_support_pcurve_completion::complete_missing_cap_support_seams(
            brep, patch.face,
        );
        completed_cap_faces.push(patch.face);
    }
}

fn complete_support_cell_pcurves_transactionally(
    brep: &mut BRep,
    corner_patches: &[CornerPatchInfo],
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
) {
    let selected = shared_trim_quality::support_cell_selected_label(diagnostics);
    for group in shared_trim_quality::support_cell_face_groups(corner_patches) {
        let before =
            shared_trim_quality::measure_support_cell_shared_trim_quality(brep, &group.faces, 9);
        let (selection, scores) =
            select_support_cell_completion_candidate(brep, &group.faces, before);
        let (accepted, selected_candidate, after) = if let Some(selection) = selection {
            let candidate = selection.candidate.label();
            let quality = selection.quality;
            *brep = selection.brep;
            (true, candidate, quality)
        } else {
            (false, "none", before)
        };
        for score in scores {
            diagnostics.push(
                shared_trim_quality::support_cell_shared_trim_candidate_fact(
                    group.n_sides,
                    selected,
                    score.candidate.label(),
                    score.status,
                    score.operation_status,
                    before,
                    score.quality,
                ),
            );
        }
        diagnostics.push(shared_trim_quality::support_cell_shared_trim_repair_fact(
            group.n_sides,
            selected,
            selected_candidate,
            accepted,
            before,
            after,
        ));
    }
}

fn select_support_cell_completion_candidate(
    brep: &BRep,
    faces: &[popkern_core::topology::FaceHandle],
    before: shared_trim_quality::Kp20SupportCellSharedTrimQuality,
) -> (
    Option<SupportCellCompletionSelection>,
    Vec<SupportCellCompletionCandidateScore>,
) {
    let mut best = None;
    let mut scores = Vec::new();
    for candidate in SupportCellCompletionCandidate::ALL {
        let mut candidate_brep = brep.clone();
        let application =
            apply_support_cell_completion_candidate(&mut candidate_brep, faces, candidate);
        if !application.operator_succeeded {
            scores.push(SupportCellCompletionCandidateScore {
                candidate,
                status: "operator_failed",
                operation_status: application.operation_status,
                quality: before,
            });
            continue;
        }
        let quality = shared_trim_quality::measure_support_cell_shared_trim_quality(
            &candidate_brep,
            faces,
            9,
        );
        let status = if candidate_brep.validate().is_err() {
            "invalid_brep"
        } else if !quality.is_no_worse_than(before) {
            "quality_regressed"
        } else {
            "accepted_not_selected"
        };
        scores.push(SupportCellCompletionCandidateScore {
            candidate,
            status,
            operation_status: application.operation_status,
            quality,
        });
        if status != "accepted_not_selected" {
            continue;
        }
        let selection = SupportCellCompletionSelection {
            brep: candidate_brep,
            candidate,
            quality,
        };
        if support_cell_completion_selection_is_better(&selection, best.as_ref()) {
            best = Some(selection);
        }
    }
    if let Some(selection) = best.as_ref() {
        for score in &mut scores {
            if score.candidate == selection.candidate {
                score.status = "accepted_selected";
            }
        }
    }
    (best, scores)
}

fn apply_support_cell_completion_candidate(
    brep: &mut BRep,
    faces: &[popkern_core::topology::FaceHandle],
    candidate: SupportCellCompletionCandidate,
) -> SupportCellCompletionApplication {
    match candidate {
        SupportCellCompletionCandidate::PcurveOnly => support_cell_operation_result(
            support_cell_pcurve_author::author_assigned_support_cell_pcurves(brep, faces).is_ok(),
            "pcurve_authored",
            "pcurve_operator_failed",
        ),
        SupportCellCompletionCandidate::BoundaryOnly => support_cell_operation_result(
            support_cell_boundary_reconcile::reconcile_assigned_support_cell_full_boundaries(
                brep, faces,
            )
            .is_ok(),
            "boundary_reconciled",
            "boundary_operator_failed",
        ),
        SupportCellCompletionCandidate::PcurveThenBoundary => {
            if support_cell_pcurve_author::author_assigned_support_cell_pcurves(brep, faces)
                .is_err()
            {
                return SupportCellCompletionApplication {
                    operator_succeeded: false,
                    operation_status: "pcurve_operator_failed",
                };
            }
            support_cell_operation_result(
                support_cell_boundary_reconcile::reconcile_assigned_support_cell_full_boundaries(
                    brep, faces,
                )
                .is_ok(),
                "pcurve_and_boundary_reconciled",
                "boundary_operator_failed",
            )
        }
        SupportCellCompletionCandidate::PcurveThenBoundaryIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_support_cell_full_boundaries(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenWorstBoundaryIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_support_cell_worst_full_boundary(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenTop3BoundariesIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_support_cell_top_full_boundaries(
                    brep, faces, 3,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryPrefixIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_support_cell_ranked_boundary_prefix(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenWorstBoundaryNeighborhoodIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_support_cell_worst_boundary_neighborhood(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenWorstEndpointClusterIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_support_cell_worst_endpoint_cluster_with_status(
                        brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryPrefixThenWorstEndpointClusterIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointClusterThenRowGroupIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_row_group_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointClusterThenTwoRowGroupsIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_two_row_groups_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointRowGroupThenTangentIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_row_group_tangent_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenWorstBoundaryRowGroupIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_support_cell_worst_boundary_row_group(
                    brep, faces,
                )
            })
        }
    }
}

fn support_cell_operation_result(
    operator_succeeded: bool,
    success_status: &'static str,
    failure_status: &'static str,
) -> SupportCellCompletionApplication {
    SupportCellCompletionApplication {
        operator_succeeded,
        operation_status: if operator_succeeded {
            success_status
        } else {
            failure_status
        },
    }
}

fn support_cell_completion_selection_is_better(
    candidate: &SupportCellCompletionSelection,
    current: Option<&SupportCellCompletionSelection>,
) -> bool {
    let Some(current) = current else {
        return true;
    };
    let a = candidate.quality;
    let b = current.quality;
    match (a.pcurve_out_of_unit_count, a.degenerate_pcurve_count)
        .cmp(&(b.pcurve_out_of_unit_count, b.degenerate_pcurve_count))
    {
        std::cmp::Ordering::Less => true,
        std::cmp::Ordering::Greater => false,
        std::cmp::Ordering::Equal => g0_tie_break_prefers_candidate(candidate, current),
    }
}

fn g0_tie_break_prefers_candidate(
    candidate: &SupportCellCompletionSelection,
    current: &SupportCellCompletionSelection,
) -> bool {
    if boundary_safe_quality_is_equivalent(candidate, current) {
        return true;
    }
    if boundary_safe_quality_is_equivalent(current, candidate) {
        return false;
    }
    candidate.quality.worst_g0_gap + 1.0e-12 < current.quality.worst_g0_gap
}

fn boundary_safe_quality_is_equivalent(
    boundary_safe: &SupportCellCompletionSelection,
    other: &SupportCellCompletionSelection,
) -> bool {
    boundary_safe.candidate.is_boundary_safe()
        && !other.candidate.is_boundary_safe()
        && boundary_safe.quality.worst_g0_gap
            <= other.quality.worst_g0_gap * (1.0 + BOUNDARY_SAFE_G0_PREFERENCE_RELATIVE_BAND)
}
