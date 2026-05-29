use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::budget::ExactBoundaryMultiPatchBudget;
use super::report::{
    ExactBoundaryMultiPatchCandidateOutcome, ExactBoundaryMultiPatchDeclineReason,
    ExactBoundaryMultiPatchValidationOrder,
};
use super::strip::ExactBoundaryStripBoundaryDirection;
use super::{ExactBoundaryMultiPatchCandidateRequest, prepare_multipatch_candidate};
use crate::apex_construction::exact_boundary_cap::ribbon::ExactBoundaryRibbonSample;
use crate::apex_construction::exact_boundary_cap::sector::{
    ExactBoundaryCapSectorInput, ExactBoundaryOuterSeam, ExactBoundaryOuterSide,
};

#[test]
fn orchestrator_builds_commit_eligible_candidate_without_live_brep_mutation() {
    let sector = sector_input();
    let boundary = cubic_boundary();
    let budget = tight_budget();
    let attempt = prepare_multipatch_candidate(ExactBoundaryMultiPatchCandidateRequest {
        sector: &sector,
        selected_sides: &[ExactBoundaryOuterSide::FaceEdge],
        boundary: &boundary,
        boundary_direction: ExactBoundaryStripBoundaryDirection::Forward,
        support_samples: &support_samples(Vec3::new(0.0, 0.0, 1.0)),
        support_authority_boundary: None,
        pole: Point3::new(0.5, 0.7, 0.0),
        budget,
    });

    assert_eq!(
        attempt.report.outcome,
        ExactBoundaryMultiPatchCandidateOutcome::CommitEligible
    );
    assert!(attempt.report.diagnostics.is_empty());
    assert_eq!(attempt.report.performance.sector_count, 1);
    assert_eq!(attempt.report.performance.outer_strip_patch_count, 1);
    assert_eq!(attempt.report.performance.interface_count, 1);
    assert_eq!(attempt.report.performance.inner_patch_count, 1);
    assert_eq!(attempt.report.performance.strip_control_point_count, 15);
    assert_eq!(attempt.report.performance.inner_control_point_count, 15);
    assert!(attempt.report.performance.bounded_iteration_count >= 12);
    assert_eq!(attempt.report.performance.deterministic_skip_code, None);
    let candidate = attempt.candidate.expect("commit-eligible candidate");
    assert_eq!(candidate.inventory_report.actual_outer_strip_count, 1);
    assert!(candidate.strip_validation.support_strip_sample_count >= 32);
    assert!(candidate.inner_patch.interface_g1_validated);
    assert!(candidate.inner_validation.interface_g1_validated);
}

#[test]
fn orchestrator_declines_oversized_strip_before_surface_allocation() {
    let sector = sector_input();
    let boundary = oversized_boundary();
    let mut budget = tight_budget();
    budget.max_width_direction_controls = 4;
    let attempt = prepare_multipatch_candidate(ExactBoundaryMultiPatchCandidateRequest {
        sector: &sector,
        selected_sides: &[ExactBoundaryOuterSide::FaceEdge],
        boundary: &boundary,
        boundary_direction: ExactBoundaryStripBoundaryDirection::Forward,
        support_samples: &support_samples(Vec3::new(0.0, 0.0, 1.0)),
        support_authority_boundary: None,
        pole: Point3::new(0.5, 0.7, 0.0),
        budget,
    });

    assert_eq!(
        attempt.report.outcome,
        ExactBoundaryMultiPatchCandidateOutcome::Declined
    );
    assert_eq!(
        attempt.report.diagnostics[0].reason,
        ExactBoundaryMultiPatchDeclineReason::SupportStripBandUnsatisfied
    );
    assert_eq!(attempt.report.performance.strip_control_point_count, 0);
    assert_eq!(attempt.report.performance.inner_control_point_count, 0);
    assert_eq!(attempt.report.performance.bounded_iteration_count, 0);
}

#[test]
fn orchestrator_robustness_matrix_covers_accept_and_decline_classes() {
    let sector = sector_input();
    let accepted = [
        (
            "linear-forward",
            line_boundary(),
            ExactBoundaryStripBoundaryDirection::Forward,
        ),
        (
            "cubic-forward",
            cubic_boundary(),
            ExactBoundaryStripBoundaryDirection::Forward,
        ),
        (
            "cubic-reversed",
            cubic_boundary(),
            ExactBoundaryStripBoundaryDirection::Reversed,
        ),
    ];
    for (name, boundary, direction) in accepted {
        let attempt = prepare_multipatch_candidate(ExactBoundaryMultiPatchCandidateRequest {
            sector: &sector,
            selected_sides: &[ExactBoundaryOuterSide::FaceEdge],
            boundary: &boundary,
            boundary_direction: direction,
            support_samples: &support_samples(Vec3::new(0.0, 0.0, 1.0)),
            support_authority_boundary: None,
            pole: Point3::new(0.5, 0.7, 0.0),
            budget: tight_budget(),
        });
        assert_eq!(
            attempt.report.outcome,
            ExactBoundaryMultiPatchCandidateOutcome::CommitEligible,
            "{name}"
        );
        assert_eq!(
            attempt.report.performance.deterministic_skip_code, None,
            "{name}"
        );
    }

    let declined = [(
        "incompatible-support-normal",
        cubic_boundary(),
        ExactBoundaryStripBoundaryDirection::Forward,
        Vec3::new(1.0, 0.0, 0.0),
        ExactBoundaryMultiPatchDeclineReason::BoundaryTangentSupportIncompatible,
    )];
    for (name, boundary, direction, support_normal, reason) in declined {
        let attempt = prepare_multipatch_candidate(ExactBoundaryMultiPatchCandidateRequest {
            sector: &sector,
            selected_sides: &[ExactBoundaryOuterSide::FaceEdge],
            boundary: &boundary,
            boundary_direction: direction,
            support_samples: &support_samples(support_normal),
            support_authority_boundary: None,
            pole: Point3::new(0.5, 0.7, 0.0),
            budget: tight_budget(),
        });
        assert_eq!(
            attempt.report.outcome,
            ExactBoundaryMultiPatchCandidateOutcome::Declined,
            "{name}"
        );
        assert_eq!(attempt.report.diagnostics[0].reason, reason, "{name}");
        assert!(attempt.candidate.is_none(), "{name}");
    }
}

#[test]
fn orchestrator_reports_first_failed_gate_with_context() {
    let sector = sector_input();
    let boundary = cubic_boundary();
    let attempt = prepare_multipatch_candidate(ExactBoundaryMultiPatchCandidateRequest {
        sector: &sector,
        selected_sides: &[ExactBoundaryOuterSide::FaceEdge],
        boundary: &boundary,
        boundary_direction: ExactBoundaryStripBoundaryDirection::Forward,
        support_samples: &support_samples(Vec3::new(1.0, 0.0, 0.0)),
        support_authority_boundary: None,
        pole: Point3::new(0.5, 0.7, 0.0),
        budget: tight_budget(),
    });

    assert_eq!(
        attempt.report.outcome,
        ExactBoundaryMultiPatchCandidateOutcome::Declined
    );
    assert!(attempt.candidate.is_none());
    assert!(attempt.report.inventory.is_some());
    assert_eq!(
        attempt
            .report
            .performance
            .deterministic_skip_code
            .as_deref(),
        Some("exact_boundary_cap.multipatch.boundary_tangent_support_incompatible")
    );
    let diagnostic = &attempt.report.diagnostics[0];
    assert_eq!(
        diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::BoundaryTangentSupportIncompatible
    );
    assert_eq!(
        diagnostic.validation_order,
        ExactBoundaryMultiPatchValidationOrder::SupportStripBand
    );
    assert_eq!(
        diagnostic.support_half_edge,
        Some(sector.support_half_edges[1])
    );
    assert_eq!(diagnostic.support_face, Some(sector.bridge_support_face));
}

#[test]
fn orchestrator_declines_two_selected_outer_sides_before_construction() {
    let sector = sector_input();
    let boundary = cubic_boundary();
    let attempt = prepare_multipatch_candidate(ExactBoundaryMultiPatchCandidateRequest {
        sector: &sector,
        selected_sides: &[
            ExactBoundaryOuterSide::FilletEnd,
            ExactBoundaryOuterSide::FaceEdge,
        ],
        boundary: &boundary,
        boundary_direction: ExactBoundaryStripBoundaryDirection::Forward,
        support_samples: &support_samples(Vec3::new(0.0, 0.0, 1.0)),
        support_authority_boundary: None,
        pole: Point3::new(0.5, 0.7, 0.0),
        budget: tight_budget(),
    });

    assert_eq!(
        attempt.report.outcome,
        ExactBoundaryMultiPatchCandidateOutcome::Declined
    );
    assert!(attempt.candidate.is_none());
    assert_eq!(attempt.report.performance.outer_strip_patch_count, 0);
    let diagnostic = &attempt.report.diagnostics[0];
    assert_eq!(
        diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::TwoOuterSidesUnsupported
    );
    assert_eq!(
        diagnostic.validation_order,
        ExactBoundaryMultiPatchValidationOrder::OuterSideSelection
    );
}

fn tight_budget() -> ExactBoundaryMultiPatchBudget {
    let mut budget = ExactBoundaryMultiPatchBudget::default();
    budget.support_normal_deficit = 1.0e-8;
    budget
}

fn cubic_boundary() -> NurbsCurve {
    NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.33, 0.1, 0.0),
            Point3::new(0.66, 0.1, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0; 4],
    }
}

fn line_boundary() -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0; 2],
    }
}

fn oversized_boundary() -> NurbsCurve {
    NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.2, 0.08, 0.0),
            Point3::new(0.4, 0.12, 0.0),
            Point3::new(0.7, 0.08, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0; 5],
    }
}

fn support_samples(support_normal: Vec3) -> Vec<ExactBoundaryRibbonSample> {
    [0.0, 0.5, 1.0]
        .into_iter()
        .map(|parameter| ExactBoundaryRibbonSample {
            parameter,
            support_uv: (parameter, 0.0),
            boundary_point: Point3::new(parameter, 0.0, 0.0),
            support_point: Point3::new(parameter, 0.0, 0.0),
            boundary_tangent: Vec3::new(1.0, 0.0, 0.0),
            support_normal,
            cross_derivative: Vec3::new(0.0, 0.2, 0.0),
        })
        .collect()
}

fn sector_input() -> ExactBoundaryCapSectorInput {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube");
    let edges: Vec<EdgeHandle> = brep.edges.keys().collect();
    let half_edges: Vec<HalfEdgeHandle> = brep.half_edges.keys().collect();
    let faces: Vec<FaceHandle> = brep.faces.keys().collect();
    ExactBoundaryCapSectorInput {
        sector_index: 0,
        sector_name: "direct_apex.cap_sector.0".to_string(),
        fillet_index: 0,
        cap_face: faces[0],
        fillet_end_edge: edges[0],
        face_edge: edges[1],
        fillet_support_face: faces[1],
        bridge_support_face: faces[2],
        cap_half_edges: [half_edges[0], half_edges[2]],
        support_half_edges: [half_edges[1], half_edges[3]],
        outer_seams: [
            seam(
                ExactBoundaryOuterSide::FilletEnd,
                edges[0],
                half_edges[0],
                half_edges[1],
                faces[1],
            ),
            seam(
                ExactBoundaryOuterSide::FaceEdge,
                edges[1],
                half_edges[2],
                half_edges[3],
                faces[2],
            ),
        ],
    }
}

fn seam(
    side: ExactBoundaryOuterSide,
    edge: EdgeHandle,
    cap_half_edge: HalfEdgeHandle,
    support_half_edge: HalfEdgeHandle,
    support_face: FaceHandle,
) -> ExactBoundaryOuterSeam {
    ExactBoundaryOuterSeam {
        side,
        edge,
        cap_half_edge,
        support_half_edge,
        support_face,
        cap_pcurve_present: true,
        support_pcurve_present: true,
    }
}
