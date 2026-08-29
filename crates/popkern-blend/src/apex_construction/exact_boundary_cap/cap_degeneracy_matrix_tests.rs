use super::facade::{ExactBoundaryCapCandidateStatus, run_exact_boundary_cap_candidate};
use super::sector::{ExactBoundarySectorDeclineReason, extract_exact_boundary_cap_sectors};
use crate::apex_construction::cap_kind::CapKind;
use crate::apex_construction::gregory_cap::GregoryCapSector;
use crate::apex_construction::plan::ApexPlan;
use crate::apex_construction::target_topology::{
    TargetCapFace, TargetEdge, TargetEdgeRole, TargetTopology,
};
use crate::apex_construction::telemetry::PlanTelemetry;
use crate::types::{
    ContinuityLevel, ExactBoundaryCandidateOutcome, ExactBoundaryReportTransactionStatus,
};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

#[test]
fn candidate_degeneracy_matrix_produces_structured_reports() {
    for case in [
        CapFixtureCase::new("tiny_sector", 3, 1.0e-4, 0.8),
        CapFixtureCase::new("high_valence", 16, 1.0, 0.8),
        CapFixtureCase::new("nearly_tangent_supports", 6, 1.0, 1.0e-4),
    ] {
        let (mut brep, cap_face, edges) = n_sector_cap_fixture(case);
        let plan = n_sector_plan(&brep, &edges, case.sector_count);
        let sectors = n_sector_inventory(&edges, case.sector_count);

        let result = run_exact_boundary_cap_candidate(&mut brep, cap_face, &plan, &sectors)
            .unwrap_or_else(|err| {
                panic!("{} should produce a structured result: {err}", case.name)
            });

        assert_eq!(result.report.schema_version, 1, "{}", case.name);
        assert_eq!(
            result.report.sectors.len(),
            case.sector_count,
            "{}",
            case.name
        );
        assert_eq!(
            result.report.performance, result.performance,
            "{}",
            case.name
        );
        assert!(result.report.preview_available, "{}", case.name);
        assert_eq!(
            result.report.continuity.target.requested,
            ContinuityLevel::G1,
            "{}",
            case.name
        );
        assert!(
            result.performance.elapsed_micros.is_some_and(|us| us > 0),
            "{}",
            case.name
        );
        if result.status == ExactBoundaryCapCandidateStatus::Declined {
            assert_eq!(
                result.report.candidate_status,
                ExactBoundaryCandidateOutcome::Declined,
                "{}",
                case.name
            );
            assert!(
                result.performance.deterministic_skip_code.is_some(),
                "{} must carry a stable skip code when declined",
                case.name
            );
            assert!(!result.report.diagnostics.is_empty(), "{}", case.name);
        }
    }
}

#[test]
fn candidate_pcurve_adversaries_are_repaired_and_committed() {
    for adversary in [
        CandidateAdversary::MissingSupportPcurve,
        CandidateAdversary::StaleSupportPcurve,
    ] {
        let (mut brep, cap_face, edges) = n_sector_cap_fixture(CapFixtureCase::standard());
        let plan = n_sector_plan(&brep, &edges, 3);
        let sectors = n_sector_inventory(&edges, 3);
        let exact_sectors =
            extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &sectors).unwrap();
        let support_half_edge = exact_sectors[0].fillet_seam.support_half_edge;
        match adversary {
            CandidateAdversary::MissingSupportPcurve => {
                brep.half_edges[support_half_edge].pcurve = None;
            }
            CandidateAdversary::StaleSupportPcurve => {
                brep.half_edges[support_half_edge].pcurve =
                    Some(Pcurve::linear((100.0, 100.0), (101.0, 100.0), 0.0, 1.0));
            }
        }

        let result = run_exact_boundary_cap_candidate(&mut brep, cap_face, &plan, &sectors)
            .expect("pcurve adversary should return a structured result");

        assert_eq!(result.status, ExactBoundaryCapCandidateStatus::Accepted);
        assert_eq!(
            result.report.provenance.transaction_status,
            ExactBoundaryReportTransactionStatus::Committed
        );
        assert!(!result.report.provenance.topology_rollback_performed);
        assert_eq!(result.committed_faces.len(), 3);
        assert_eq!(result.report.committed_faces.len(), 3);
        assert_eq!(
            result.report.candidate_status,
            ExactBoundaryCandidateOutcome::Accepted
        );
        assert!(result.performance.deterministic_skip_code.is_none());
        assert!(
            result
                .report
                .diagnostics
                .iter()
                .all(|d| crate::types::diagnostic_code_is_informational(&d.code)),
            "{adversary:?}"
        );
        assert!(
            result
                .report
                .diagnostics
                .iter()
                .all(|diagnostic| diagnostic.schema_version == 1)
        );
    }
}

#[test]
fn topology_ownership_adversaries_decline_with_stable_reasons() {
    for adversary in [
        OwnershipAdversary::MissingEdgeCurve,
        OwnershipAdversary::MissingCapHalfEdge,
    ] {
        let (mut brep, cap_face, edges) = n_sector_cap_fixture(CapFixtureCase::standard());
        let plan = n_sector_plan(&brep, &edges, 3);
        let mut sectors = n_sector_inventory(&edges, 3);
        match adversary {
            OwnershipAdversary::MissingEdgeCurve => {
                brep.edges[edges[0]].curve = None;
            }
            OwnershipAdversary::MissingCapHalfEdge => {
                sectors[0].fillet_end_edge = first_non_cap_edge(&brep, &edges);
            }
        }

        let diagnostic =
            extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &sectors).unwrap_err();

        assert_eq!(diagnostic.reason, adversary.expected_reason());
        assert_eq!(diagnostic.sector_index, Some(0));
        assert_eq!(
            diagnostic.sector_name.as_deref(),
            Some("direct_apex.cap_sector.0")
        );
        assert_eq!(diagnostic.cap_face, Some(cap_face));
        assert!(diagnostic.edge.is_some());
    }
}

#[derive(Clone, Copy)]
struct CapFixtureCase {
    name: &'static str,
    sector_count: usize,
    radius: f64,
    height: f64,
}

impl CapFixtureCase {
    const fn new(name: &'static str, sector_count: usize, radius: f64, height: f64) -> Self {
        Self {
            name,
            sector_count,
            radius,
            height,
        }
    }

    const fn standard() -> Self {
        Self::new("standard", 3, 1.0, 0.8)
    }
}

#[derive(Debug, Clone, Copy)]
enum CandidateAdversary {
    MissingSupportPcurve,
    StaleSupportPcurve,
}

#[derive(Debug, Clone, Copy)]
enum OwnershipAdversary {
    MissingEdgeCurve,
    MissingCapHalfEdge,
}

impl OwnershipAdversary {
    const fn expected_reason(self) -> ExactBoundarySectorDeclineReason {
        match self {
            Self::MissingEdgeCurve => ExactBoundarySectorDeclineReason::MissingEdgeCurve,
            Self::MissingCapHalfEdge => ExactBoundarySectorDeclineReason::MissingCapHalfEdge,
        }
    }
}

fn n_sector_cap_fixture(case: CapFixtureCase) -> (BRep, FaceHandle, Vec<EdgeHandle>) {
    let edge_count = case.sector_count * 2;
    let brep =
        BRep::n_gonal_pyramid(edge_count, case.radius, case.height, Tolerance::default()).unwrap();
    let cap_face = brep
        .faces
        .keys()
        .find(|&face| loop_edges(&brep, face).len() == edge_count)
        .unwrap();
    let edges = loop_edges(&brep, cap_face);
    (brep, cap_face, edges)
}

fn first_non_cap_edge(brep: &BRep, cap_edges: &[EdgeHandle]) -> EdgeHandle {
    brep.edges
        .keys()
        .find(|edge| !cap_edges.contains(edge))
        .expect("pyramid fixture should contain side edges outside the cap loop")
}

fn loop_edges(brep: &BRep, face: FaceHandle) -> Vec<EdgeHandle> {
    let start = brep.loops[brep.faces[face].outer_loop].half_edge;
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

fn n_sector_inventory(edges: &[EdgeHandle], sector_count: usize) -> Vec<GregoryCapSector> {
    (0..sector_count)
        .map(|sector_index| GregoryCapSector {
            sector_index,
            fillet_index: sector_index,
            fillet_end_edge: edges[2 * sector_index],
            face_edge: Some(edges[2 * sector_index + 1]),
        })
        .collect()
}

fn n_sector_plan(brep: &BRep, edges: &[EdgeHandle], sector_count: usize) -> ApexPlan {
    ApexPlan {
        apex_vertex: brep.vertices.keys().next().unwrap(),
        cap_kind: CapKind::Gregory {
            per_fillet_terminations: vec![
                [Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0)];
                sector_count
            ],
            per_fillet_caps: vec![(Point3::new(0.0, 0.0, 0.0), 1.0); sector_count],
            corner_kinds: Vec::new(),
        },
        target_topology: TargetTopology {
            apex_vertices: Vec::new(),
            apex_edges: (0..edges.len())
                .map(|index| TargetEdge {
                    origin_vertex_idx: index,
                    dest_vertex_idx: (index + 1) % edges.len(),
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
