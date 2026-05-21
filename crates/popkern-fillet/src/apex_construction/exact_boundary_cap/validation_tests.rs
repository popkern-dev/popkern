use super::boundary::{ExactBoundaryCurveBudget, prepare_exact_boundary_curves};
use super::corner::{ExactBoundaryInnerSpokeBudget, prepare_exact_boundary_inner_spokes};
use super::pcurves::{PreparedExactBoundaryPcurvePair, PreparedExactBoundaryPcurves};
use super::ribbon::{ExactBoundaryRibbonBudget, prepare_exact_boundary_ribbons};
use super::sector::extract_exact_boundary_cap_sectors;
use super::solve::{ExactBoundarySolveBudget, solve_exact_boundary_surfaces};
use super::topology::materialize_exact_boundary_sector_topology;
use super::validation::{
    ExactBoundaryValidationBudget, ExactBoundaryValidationDeclineReason,
    validate_exact_boundary_candidate,
};
use crate::apex_construction::cap_kind::CapKind;
use crate::apex_construction::gregory_cap::GregoryCapSector;
use crate::apex_construction::plan::ApexPlan;
use crate::apex_construction::target_topology::{
    TargetCapFace, TargetEdge, TargetEdgeRole, TargetTopology,
};
use crate::apex_construction::telemetry::PlanTelemetry;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

#[test]
fn validation_accepts_materialized_exact_boundary_candidate() {
    let (mut brep, cap_face, edges) = hex_cap_fixture();
    let fixture = candidate_fixture(&mut brep, cap_face, &edges);

    let report = validate_exact_boundary_candidate(
        &brep,
        &fixture.topology,
        &fixture.topology.sectors,
        &fixture.spokes,
        &fixture.solved,
        &fixture.pcurves,
        ExactBoundaryValidationBudget::default(),
    )
    .expect("valid materialized candidate passes exact-boundary gates");

    assert_eq!(report.validated_sector_count, 3);
    assert_eq!(report.validated_pcurve_pair_count, 6);
    assert_eq!(report.validated_spoke_edge_count, 3);
}

#[test]
fn validation_declines_pair_residual_over_budget() {
    let (mut brep, cap_face, edges) = hex_cap_fixture();
    let mut fixture = candidate_fixture(&mut brep, cap_face, &edges);
    fixture.pcurves.pairs[0].max_pair_residual = 1.0;

    let diagnostic = validate_exact_boundary_candidate(
        &brep,
        &fixture.topology,
        &fixture.topology.sectors,
        &fixture.spokes,
        &fixture.solved,
        &fixture.pcurves,
        ExactBoundaryValidationBudget::default(),
    )
    .expect_err("oversized pair residual must decline before commit");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryValidationDeclineReason::PairBoundaryResidualExceeded
    );
    assert_eq!(
        diagnostic.code,
        "exact_boundary_cap.validation.pair_boundary_residual_exceeded"
    );
}

#[test]
fn validation_declines_foldover_budget_failure() {
    let (mut brep, cap_face, edges) = hex_cap_fixture();
    let fixture = candidate_fixture(&mut brep, cap_face, &edges);

    let diagnostic = validate_exact_boundary_candidate(
        &brep,
        &fixture.topology,
        &fixture.topology.sectors,
        &fixture.spokes,
        &fixture.solved,
        &fixture.pcurves,
        ExactBoundaryValidationBudget {
            min_normal_length: 1.0e9,
            ..ExactBoundaryValidationBudget::default()
        },
    )
    .expect_err("oversized normal budget must decline");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryValidationDeclineReason::FoldoverDetected
    );
    assert_eq!(
        diagnostic.code,
        "exact_boundary_cap.validation.foldover_detected"
    );
}

struct CandidateFixture {
    topology: super::topology::ExactBoundarySectorTopology,
    spokes: super::corner::PreparedExactBoundaryInnerSpokes,
    solved: super::solve::ExactBoundarySolvedSurfaces,
    pcurves: super::pcurves::PreparedExactBoundaryPcurves,
}

fn candidate_fixture(
    brep: &mut BRep,
    cap_face: FaceHandle,
    edges: &[EdgeHandle],
) -> CandidateFixture {
    let plan = trihedral_plan(brep, edges);
    let sectors =
        extract_exact_boundary_cap_sectors(brep, cap_face, &plan, &trihedral_sectors(edges))
            .unwrap();
    let boundaries =
        prepare_exact_boundary_curves(brep, &sectors, ExactBoundaryCurveBudget::default()).unwrap();
    let spokes = prepare_exact_boundary_inner_spokes(
        &sectors,
        &boundaries,
        ExactBoundaryInnerSpokeBudget::default(),
    )
    .unwrap();
    let ribbons = prepare_exact_boundary_ribbons(
        brep,
        &sectors,
        &boundaries,
        ExactBoundaryRibbonBudget::default(),
    )
    .unwrap();
    let solved =
        solve_exact_boundary_surfaces(&spokes, &ribbons, ExactBoundarySolveBudget::default())
            .unwrap();
    let topology =
        materialize_exact_boundary_sector_topology(brep, cap_face, &sectors, &spokes).unwrap();
    let pcurves = prepared_pcurves(&topology.sectors);
    CandidateFixture {
        topology,
        spokes,
        solved,
        pcurves,
    }
}

fn prepared_pcurves(
    sectors: &[super::sector::ExactBoundaryCapSectorInput],
) -> PreparedExactBoundaryPcurves {
    let mut pairs = Vec::with_capacity(sectors.len() * 2);
    for sector in sectors {
        pairs.push(pcurve_pair(
            sector,
            super::sector::ExactBoundaryOuterSide::FilletEnd,
        ));
        pairs.push(pcurve_pair(
            sector,
            super::sector::ExactBoundaryOuterSide::FaceEdge,
        ));
    }
    PreparedExactBoundaryPcurves { pairs }
}

fn pcurve_pair(
    sector: &super::sector::ExactBoundaryCapSectorInput,
    side: super::sector::ExactBoundaryOuterSide,
) -> PreparedExactBoundaryPcurvePair {
    PreparedExactBoundaryPcurvePair {
        sector_index: sector.sector_index,
        sector_name: sector.sector_name.clone(),
        side,
        edge: match side {
            super::sector::ExactBoundaryOuterSide::FilletEnd => sector.fillet_end_edge,
            super::sector::ExactBoundaryOuterSide::FaceEdge => sector.face_edge,
        },
        support_half_edge: match side {
            super::sector::ExactBoundaryOuterSide::FilletEnd => sector.support_half_edges[0],
            super::sector::ExactBoundaryOuterSide::FaceEdge => sector.support_half_edges[1],
        },
        support_face: match side {
            super::sector::ExactBoundaryOuterSide::FilletEnd => sector.fillet_support_face,
            super::sector::ExactBoundaryOuterSide::FaceEdge => sector.bridge_support_face,
        },
        cap_pcurve: Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0),
        support_pcurve: Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0),
        max_cap_residual: 0.0,
        max_support_residual: 0.0,
        max_pair_residual: 0.0,
    }
}

fn hex_cap_fixture() -> (BRep, FaceHandle, Vec<EdgeHandle>) {
    let brep = BRep::n_gonal_pyramid(6, 1.0, 0.8, Tolerance::default()).unwrap();
    let cap_face = brep
        .faces
        .keys()
        .find(|&face| loop_edges(&brep, face).len() == 6)
        .expect("hexagonal base face exists");
    let edges = loop_edges(&brep, cap_face);
    (brep, cap_face, edges)
}

fn loop_edges(brep: &BRep, face: FaceHandle) -> Vec<EdgeHandle> {
    let loop_handle = brep.faces[face].outer_loop;
    let start = brep.loops[loop_handle].half_edge;
    let mut cursor = start;
    let mut edges = Vec::new();
    for _ in 0..=brep.half_edges.len() {
        let half_edge = &brep.half_edges[cursor];
        edges.push(half_edge.edge);
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    edges
}

fn trihedral_sectors(edges: &[EdgeHandle]) -> Vec<GregoryCapSector> {
    (0..3)
        .map(|sector_index| GregoryCapSector {
            sector_index,
            fillet_index: sector_index,
            fillet_end_edge: edges[2 * sector_index],
            face_edge: edges[2 * sector_index + 1],
        })
        .collect()
}

fn trihedral_plan(brep: &BRep, edges: &[EdgeHandle]) -> ApexPlan {
    ApexPlan {
        apex_vertex: brep.vertices.keys().next().unwrap(),
        cap_kind: CapKind::Gregory {
            per_fillet_terminations: vec![
                [Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0)];
                3
            ],
            per_fillet_caps: vec![(Point3::new(0.0, 0.0, 0.0), 1.0); 3],
        },
        target_topology: TargetTopology {
            apex_vertices: Vec::new(),
            apex_edges: (0..6)
                .map(|index| TargetEdge {
                    origin_vertex_idx: index,
                    dest_vertex_idx: (index + 1) % 6,
                    curve: dummy_curve(),
                    role: if index % 2 == 0 {
                        TargetEdgeRole::CapVsCyl {
                            fillet_idx: index / 2,
                        }
                    } else {
                        TargetEdgeRole::CapVsPlane {
                            shared_face_idx: index / 2,
                        }
                    },
                })
                .collect(),
            cap_face: TargetCapFace::Gregory {
                outer_loop_edge_indices: (0..edges.len()).collect(),
            },
            face_boundary_updates: Vec::new(),
        },
        telemetry: PlanTelemetry::default(),
    }
}

fn dummy_curve() -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    }
}
