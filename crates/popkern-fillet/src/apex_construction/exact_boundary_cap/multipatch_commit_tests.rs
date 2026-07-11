use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::{commit_multipatch_candidate_transactionally, preflight_multipatch_commit};
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
        let pole_index = index * missing_spokes.inner_patch.surface.num_v
            + missing_spokes.inner_patch.surface.num_v
            - 1;
        missing_spokes.inner_patch.surface.control_points[pole_index].x += index as f64;
    }
    let spoke_diagnostic = preflight_multipatch_commit(&sector, &missing_spokes)
        .expect_err("missing generated spokes must decline");
    assert_eq!(
        spoke_diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::CommitPreflightMissingGeneratedSpoke
    );
}

#[test]
fn transactional_commit_materializes_split_faces_and_validates_brep() {
    let (mut brep, sector, boundary, samples, pole) = live_sector_fixture();
    let initial_face_count = brep.faces.len();
    let initial_edge_count = brep.edges.len();
    let initial_vertex_count = brep.vertices.len();
    let support_pcurve_before = brep.half_edges[sector.support_half_edges[1]].pcurve.clone();
    let candidate = commit_eligible_candidate_from_boundary(&sector, &boundary, &samples, pole);

    let report = commit_multipatch_candidate_transactionally(&mut brep, &sector, &candidate)
        .expect("validated candidate should commit transactionally");

    assert_eq!(report.strip_face, sector.cap_face);
    assert!(brep.faces.contains_key(report.inner_face));
    assert_eq!(report.internal_interface_edges.len(), 1);
    assert_eq!(report.generated_connector_edges.len(), 2);
    assert_eq!(report.committed_pcurve_count, 3);
    assert_eq!(brep.faces.len(), initial_face_count + 1);
    assert_eq!(brep.edges.len(), initial_edge_count + 3);
    assert_eq!(brep.vertices.len(), initial_vertex_count + 2);
    assert!(brep.faces[report.strip_face].surface.is_some());
    assert!(brep.faces[report.inner_face].surface.is_some());
    assert!(brep.half_edges[sector.cap_half_edges[1]].pcurve.is_some());
    assert_eq!(
        brep.half_edges[sector.support_half_edges[1]]
            .pcurve
            .is_some(),
        support_pcurve_before.is_some()
    );
    brep.validate()
        .expect("committed multipatch split must leave valid BRep topology");
}

#[test]
fn transactional_commit_authors_support_authority_edge_curve() {
    let (mut brep, sector, boundary, samples, pole) = live_sector_fixture();
    let mut candidate = commit_eligible_candidate_from_boundary(&sector, &boundary, &samples, pole);
    let authority_boundary = quadratic_boundary_on_same_line(&boundary);
    candidate.support_authority_boundary = Some(authority_boundary.clone());

    commit_multipatch_candidate_transactionally(&mut brep, &sector, &candidate)
        .expect("validated candidate should commit transactionally");

    let edge_curve = brep.edges[sector.face_edge]
        .curve
        .as_ref()
        .expect("selected edge must keep 3-D curve authority");
    assert_eq!(edge_curve.t_start, authority_boundary.parameter_range().0);
    assert_eq!(edge_curve.t_end, authority_boundary.parameter_range().1);
    let CurveGeometry::Nurbs(committed_boundary) = &edge_curve.geometry else {
        panic!("support-authority boundary must be authored as selected edge NURBS");
    };
    assert_eq!(committed_boundary.degree, authority_boundary.degree);
    assert_eq!(
        committed_boundary.control_points,
        authority_boundary.control_points
    );
    brep.validate()
        .expect("support-authority edge curve must preserve valid topology");
}

#[test]
fn transactional_commit_rolls_back_when_source_face_is_stale() {
    let (mut brep, mut sector, boundary, samples, pole) = live_sector_fixture();
    let candidate = commit_eligible_candidate_from_boundary(&sector, &boundary, &samples, pole);
    let before_counts = topology_counts(&brep);
    sector.cap_face = FaceHandle::default();

    let diagnostic = commit_multipatch_candidate_transactionally(&mut brep, &sector, &candidate)
        .expect_err("stale live topology must decline");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::CommitMaterializationFailed
    );
    assert_eq!(topology_counts(&brep), before_counts);
}

#[test]
fn transactional_commit_rolls_back_post_validation_failure() {
    let (mut brep, sector, boundary, samples, pole) = live_sector_fixture();
    let mut candidate = commit_eligible_candidate_from_boundary(&sector, &boundary, &samples, pole);
    let before_counts = topology_counts(&brep);
    let source_surface_before = brep.faces[sector.cap_face].surface.clone();
    let collapsed = candidate.inner_patch.surface.control_points[0];
    for point in &mut candidate.inner_patch.surface.control_points {
        *point = collapsed;
    }

    let diagnostic = commit_multipatch_candidate_transactionally(&mut brep, &sector, &candidate)
        .expect_err("degenerate committed surface must fail post-validation");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::CommitValidationFailed
    );
    assert_eq!(topology_counts(&brep), before_counts);
    assert_eq!(
        brep.faces[sector.cap_face].surface.is_some(),
        source_surface_before.is_some()
    );
    brep.validate()
        .expect("rollback must restore the source BRep to validity");
}

fn commit_eligible_candidate(
    sector: &ExactBoundaryCapSectorInput,
) -> ExactBoundaryMultiPatchCandidate {
    let boundary = cubic_boundary();
    commit_eligible_candidate_from_boundary(
        sector,
        &boundary,
        &support_samples(),
        Point3::new(0.5, 0.7, 0.0),
    )
}

fn commit_eligible_candidate_from_boundary(
    sector: &ExactBoundaryCapSectorInput,
    boundary: &NurbsCurve,
    samples: &[ExactBoundaryRibbonSample],
    pole: Point3,
) -> ExactBoundaryMultiPatchCandidate {
    let budget = ExactBoundaryMultiPatchBudget {
        support_normal_deficit: 1.0e-8,
        ..Default::default()
    };
    let attempt = prepare_multipatch_candidate(ExactBoundaryMultiPatchCandidateRequest {
        sector,
        selected_sides: &[ExactBoundaryOuterSide::FaceEdge],
        boundary,
        boundary_direction: ExactBoundaryStripBoundaryDirection::Forward,
        support_samples: samples,
        support_authority_boundary: None,
        pole,
        budget,
    });
    assert_eq!(
        attempt.report.outcome,
        ExactBoundaryMultiPatchCandidateOutcome::CommitEligible
    );
    attempt.candidate.expect("commit-eligible candidate")
}

fn live_sector_fixture() -> (
    BRep,
    ExactBoundaryCapSectorInput,
    NurbsCurve,
    Vec<ExactBoundaryRibbonSample>,
    Point3,
) {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube");
    let cap_face = brep.faces.keys().next().expect("face");
    let loop_half_edges = loop_half_edges(&brep, cap_face);
    let selected = loop_half_edges[0];
    let companion = loop_half_edges[1];
    let selected_twin = brep.half_edges[selected].twin;
    let companion_twin = brep.half_edges[companion].twin;
    let selected_edge = brep.half_edges[selected].edge;
    let companion_edge = brep.half_edges[companion].edge;
    let selected_support_face = face_for_half_edge(&brep, selected_twin);
    let companion_support_face = face_for_half_edge(&brep, companion_twin);
    let start = brep.vertices[brep.half_edges[selected].origin].position;
    let end = brep.vertices[destination_vertex(&brep, selected)].position;
    let centroid = face_centroid(&brep, cap_face);
    let midpoint = start + (end - start) * 0.5;
    let tangent = (end - start).normalized().expect("non-degenerate edge");
    let inward = (centroid - midpoint)
        .normalized()
        .expect("selected edge has interior direction");
    let normal = tangent.cross(&inward).normalized().expect("face normal");
    let boundary = line_boundary(start, end);
    let samples = live_support_samples(start, end, normal, inward);
    let pole = midpoint + inward * 0.45;
    let sector = ExactBoundaryCapSectorInput {
        sector_index: 0,
        sector_name: "direct_apex.cap_sector.live".to_string(),
        fillet_index: 0,
        cap_face,
        fillet_end_edge: companion_edge,
        face_edge: selected_edge,
        fillet_support_face: companion_support_face,
        bridge_support_face: selected_support_face,
        cap_half_edges: [companion, selected],
        support_half_edges: [companion_twin, selected_twin],
        outer_seams: [
            seam(
                ExactBoundaryOuterSide::FilletEnd,
                companion_edge,
                companion,
                companion_twin,
                companion_support_face,
            ),
            seam(
                ExactBoundaryOuterSide::FaceEdge,
                selected_edge,
                selected,
                selected_twin,
                selected_support_face,
            ),
        ],
    };
    (brep, sector, boundary, samples, pole)
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

fn live_support_samples(
    start: Point3,
    end: Point3,
    support_normal: Vec3,
    inward: Vec3,
) -> Vec<ExactBoundaryRibbonSample> {
    [0.0, 0.5, 1.0]
        .into_iter()
        .map(|parameter| {
            let boundary_point = start + (end - start) * parameter;
            ExactBoundaryRibbonSample {
                parameter,
                support_uv: (parameter, 0.0),
                boundary_point,
                support_point: boundary_point,
                boundary_tangent: end - start,
                support_normal,
                cross_derivative: inward * 0.2,
            }
        })
        .collect()
}

fn line_boundary(start: Point3, end: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![start, end],
        weights: vec![1.0, 1.0],
    }
}

fn quadratic_boundary_on_same_line(boundary: &NurbsCurve) -> NurbsCurve {
    let (start, end) = boundary.parameter_range();
    let p0 = boundary.point_at(start);
    let p2 = boundary.point_at(end);
    let p1 = p0 + (p2 - p0) * 0.5;
    NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![p0, p1, p2],
        weights: vec![1.0, 1.0, 1.0],
    }
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

fn loop_half_edges(brep: &BRep, face: FaceHandle) -> Vec<HalfEdgeHandle> {
    let start = brep.loops[brep.faces[face].outer_loop].half_edge;
    let mut out = Vec::new();
    let mut cursor = start;
    loop {
        out.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            return out;
        }
    }
}

fn destination_vertex(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> popkern_core::topology::VertexHandle {
    brep.half_edges[brep.half_edges[half_edge].twin].origin
}

fn face_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> FaceHandle {
    brep.loops[brep.half_edges[half_edge].loop_].face
}

fn face_centroid(brep: &BRep, face: FaceHandle) -> Point3 {
    let points = loop_half_edges(brep, face)
        .into_iter()
        .map(|half_edge| brep.vertices[brep.half_edges[half_edge].origin].position)
        .collect::<Vec<_>>();
    let sum = points.iter().fold(Vec3::new(0.0, 0.0, 0.0), |acc, point| {
        acc + Vec3::new(point.x, point.y, point.z)
    });
    let scale = 1.0 / points.len() as f64;
    Point3::new(sum.x * scale, sum.y * scale, sum.z * scale)
}

fn topology_counts(brep: &BRep) -> (usize, usize, usize, usize, usize) {
    (
        brep.vertices.len(),
        brep.edges.len(),
        brep.half_edges.len(),
        brep.loops.len(),
        brep.faces.len(),
    )
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
