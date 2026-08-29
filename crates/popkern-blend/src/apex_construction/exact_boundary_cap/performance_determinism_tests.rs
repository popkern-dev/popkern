use super::facade::{ExactBoundaryCapCandidateStatus, run_exact_boundary_cap_candidate};
use super::ribbon::ExactBoundaryRibbonBudget;
use crate::apex_construction::cap_kind::CapKind;
use crate::apex_construction::gregory_cap::GregoryCapSector;
use crate::apex_construction::plan::ApexPlan;
use crate::apex_construction::target_topology::{
    TargetCapFace, TargetEdge, TargetEdgeRole, TargetTopology,
};
use crate::apex_construction::telemetry::PlanTelemetry;
use crate::types::ExactBoundaryCandidatePerformanceReport;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

#[test]
fn exact_boundary_candidate_reports_linear_scaling_counts() {
    let mut surface_control_points_per_sector = None;
    for sector_count in [3, 4, 8, 12] {
        let result = run_candidate(sector_count);
        let performance = &result.performance;

        assert_eq!(performance.sector_count, sector_count);
        assert_eq!(performance.prepared_curve_count, sector_count * 2);
        assert_eq!(performance.inner_spoke_count, sector_count * 2);
        assert_eq!(performance.support_ribbon_count, sector_count * 2);
        assert_eq!(
            performance.support_sample_count,
            performance.support_ribbon_count * ExactBoundaryRibbonBudget::default().sample_count
        );
        assert_eq!(performance.solved_surface_count, sector_count);
        assert_eq!(performance.surface_control_point_count % sector_count, 0);
        let per_sector = performance.surface_control_point_count / sector_count;
        if let Some(expected) = surface_control_points_per_sector {
            assert_eq!(per_sector, expected);
        } else {
            surface_control_points_per_sector = Some(per_sector);
        }
        assert!(per_sector >= 100);
        assert!(performance.pcurve_pair_count <= sector_count * 2);
        assert!(performance.validation_sector_count <= sector_count);
        assert!(performance.validation_pcurve_pair_count <= sector_count * 2);
        assert!(performance.validation_spoke_edge_count <= sector_count);
        assert!(performance.elapsed_micros.is_some_and(|us| us > 0));
        assert_eq!(result.report.performance, result.performance);

        if result.status == ExactBoundaryCapCandidateStatus::Declined {
            assert!(
                performance.deterministic_skip_code.is_some(),
                "declined candidate must expose a deterministic skip code"
            );
        }
    }
}

#[test]
fn exact_boundary_candidate_output_is_repeatable_except_elapsed_time() {
    let first = run_candidate(8);
    let second = run_candidate(8);

    assert_eq!(first.status, second.status);
    assert_eq!(
        without_elapsed(&first.performance),
        without_elapsed(&second.performance)
    );
    assert_eq!(first.report.sectors, second.report.sectors);
    assert_eq!(
        diagnostic_codes(&first.report.diagnostics),
        diagnostic_codes(&second.report.diagnostics)
    );
    assert_eq!(
        first.report.tolerance_evidence.deterministic_skip_code,
        second.report.tolerance_evidence.deterministic_skip_code
    );
    assert_eq!(first.report.continuity, second.report.continuity);
}

fn run_candidate(sector_count: usize) -> super::facade::ExactBoundaryCapCandidateResult {
    let (mut brep, cap_face, edges) = n_sector_cap_fixture(sector_count);
    let plan = n_sector_plan(&brep, &edges, sector_count);
    let sectors = n_sector_inventory(&edges, sector_count);
    run_exact_boundary_cap_candidate(&mut brep, cap_face, &plan, &sectors)
        .expect("fixture candidate should produce a structured result")
}

fn without_elapsed(
    performance: &ExactBoundaryCandidatePerformanceReport,
) -> ExactBoundaryCandidatePerformanceReport {
    let mut stable = performance.clone();
    stable.elapsed_micros = None;
    stable
}

fn diagnostic_codes(
    diagnostics: &[crate::types::ExactBoundaryReportDiagnostic],
) -> Vec<(String, String)> {
    diagnostics
        .iter()
        .map(|diagnostic| (diagnostic.stage.clone(), diagnostic.code.clone()))
        .collect()
}

fn n_sector_cap_fixture(sector_count: usize) -> (BRep, FaceHandle, Vec<EdgeHandle>) {
    let edge_count = sector_count * 2;
    let brep = BRep::n_gonal_pyramid(edge_count, 1.0, 0.8, Tolerance::default()).unwrap();
    let cap_face = brep
        .faces
        .keys()
        .find(|&face| loop_edges(&brep, face).len() == edge_count)
        .unwrap();
    let edges = loop_edges(&brep, cap_face);
    (brep, cap_face, edges)
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
