use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::preflight_multipatch_commit;
use crate::apex_construction::exact_boundary_cap::multipatch::budget::ExactBoundaryMultiPatchBudget;
use crate::apex_construction::exact_boundary_cap::multipatch::report::{
    ExactBoundaryMultiPatchCandidateOutcome, ExactBoundaryMultiPatchDeclineReason,
};
use crate::apex_construction::exact_boundary_cap::multipatch::strip::ExactBoundaryStripBoundaryDirection;
use crate::apex_construction::exact_boundary_cap::multipatch::{
    ExactBoundaryMultiPatchCandidate, ExactBoundaryMultiPatchCandidateRequest,
    prepare_multipatch_candidate,
};
use crate::apex_construction::exact_boundary_cap::ribbon::ExactBoundaryRibbonSample;
use crate::apex_construction::exact_boundary_cap::sector::{
    ExactBoundaryCapSectorInput, ExactBoundaryOuterSeam, ExactBoundaryOuterSide,
};

#[test]
fn commit_preflight_reports_candidate_topology_delta_without_mutation() {
    let sector = sector_input();
    let candidate = commit_eligible_candidate(&sector);

    let report = preflight_multipatch_commit(&sector, &candidate).expect("preflight passes");

    assert_eq!(
        report.logical_feature_name,
        "direct_apex.cap_sector.0.multipatch_cap"
    );
    assert_eq!(
        report.generated_face_names,
        vec![
            candidate.strip_inventory.patch_name.clone(),
            candidate.inner_patch.patch_name.clone()
        ]
    );
    assert_eq!(
        report.internal_interface_names,
        vec![candidate.interface.interface_name.clone()]
    );
    assert_eq!(
        report.reused_outer_edges,
        vec![candidate.strip_inventory.edge]
    );
    assert_eq!(report.pcurve_count, 3);
    assert!(report.generated_spoke_count >= 2);
}

#[test]
fn commit_preflight_declines_missing_surface_with_source_handles() {
    let sector = sector_input();
    let mut candidate = commit_eligible_candidate(&sector);
    candidate.inner_patch.surface.control_points.clear();

    let diagnostic =
        preflight_multipatch_commit(&sector, &candidate).expect_err("missing surface must decline");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::CommitPreflightMissingSurface
    );
    assert_eq!(
        diagnostic.code,
        "exact_boundary_cap.multipatch.commit_preflight_missing_surface"
    );
    assert_eq!(
        diagnostic.patch_name.as_deref(),
        Some(candidate.inner_patch.patch_name.as_str())
    );
    assert_eq!(diagnostic.edge, Some(candidate.strip_inventory.edge));
    assert_eq!(
        diagnostic.support_half_edge,
        Some(candidate.strip_inventory.support_half_edge)
    );
    assert_eq!(
        diagnostic.support_face,
        Some(candidate.strip_inventory.support_face)
    );
}

#[test]
fn commit_preflight_declines_missing_pcurve_or_interface_or_spokes() {
    let sector = sector_input();
    let mut missing_pcurve = commit_eligible_candidate(&sector);
    missing_pcurve.inner_patch.interface_pcurve.t_end =
        missing_pcurve.inner_patch.interface_pcurve.t_start;
    let pcurve_diagnostic = preflight_multipatch_commit(&sector, &missing_pcurve)
        .expect_err("missing pcurve must decline");
    assert_eq!(
        pcurve_diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::CommitPreflightMissingPcurve
    );

    let mut missing_interface = commit_eligible_candidate(&sector);
    missing_interface.interface.curve.control_points.clear();
    let interface_diagnostic = preflight_multipatch_commit(&sector, &missing_interface)
        .expect_err("missing interface must decline");
    assert_eq!(
        interface_diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::CommitPreflightMissingInterface
    );

    let mut missing_spokes = commit_eligible_candidate(&sector);
    for index in 1..missing_spokes.inner_patch.surface.num_u {
        let pole_index = index * missing_spokes.inner_patch.surface.num_v + 1;
        missing_spokes.inner_patch.surface.control_points[pole_index].x += index as f64;
    }
    let spoke_diagnostic = preflight_multipatch_commit(&sector, &missing_spokes)
        .expect_err("missing generated spokes must decline");
    assert_eq!(
        spoke_diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::CommitPreflightMissingGeneratedSpoke
    );
}

fn commit_eligible_candidate(
    sector: &ExactBoundaryCapSectorInput,
) -> ExactBoundaryMultiPatchCandidate {
    let boundary = cubic_boundary();
    let mut budget = ExactBoundaryMultiPatchBudget::default();
    budget.support_normal_deficit = 1.0e-8;
    let attempt = prepare_multipatch_candidate(ExactBoundaryMultiPatchCandidateRequest {
        sector,
        selected_sides: &[ExactBoundaryOuterSide::FaceEdge],
        boundary: &boundary,
        boundary_direction: ExactBoundaryStripBoundaryDirection::Forward,
        support_samples: &support_samples(),
        pole: Point3::new(0.5, 0.7, 0.0),
        budget,
    });
    assert_eq!(
        attempt.report.outcome,
        ExactBoundaryMultiPatchCandidateOutcome::CommitEligible
    );
    attempt.candidate.expect("commit-eligible candidate")
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

fn support_samples() -> Vec<ExactBoundaryRibbonSample> {
    [0.0, 0.5, 1.0]
        .into_iter()
        .map(|parameter| ExactBoundaryRibbonSample {
            parameter,
            support_uv: (parameter, 0.0),
            boundary_point: Point3::new(parameter, 0.0, 0.0),
            support_point: Point3::new(parameter, 0.0, 0.0),
            boundary_tangent: Vec3::new(1.0, 0.0, 0.0),
            support_normal: Vec3::new(0.0, 0.0, 1.0),
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
