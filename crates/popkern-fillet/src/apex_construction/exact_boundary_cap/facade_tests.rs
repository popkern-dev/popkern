use super::facade::{ExactBoundaryCapCandidateStatus, build_exact_boundary_report};
use super::sector::extract_exact_boundary_cap_sectors;
use crate::apex_construction::cap_kind::CapKind;
use crate::apex_construction::gregory_cap::GregoryCapSector;
use crate::apex_construction::plan::ApexPlan;
use crate::apex_construction::target_topology::{
    TargetCapFace, TargetEdge, TargetEdgeRole, TargetTopology,
};
use crate::apex_construction::telemetry::PlanTelemetry;
use crate::types::{
    ApexConstructionDiagnostic, ExactBoundaryCandidateOutcome,
    ExactBoundaryCandidatePerformanceReport, ExactBoundaryCapDiagnostic,
    ExactBoundaryReportTransactionStatus,
};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

#[test]
fn report_sector_pcurve_flags_are_read_from_live_topology() {
    let (mut brep, cap_face, edges) = hex_cap_fixture();
    let plan = trihedral_plan(&brep, &edges);
    let sectors = trihedral_sectors(&edges);
    let mut exact_sectors =
        extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &sectors).unwrap();
    for sector in &mut exact_sectors {
        for seam in &mut sector.outer_seams {
            seam.cap_pcurve_present = false;
            seam.support_pcurve_present = false;
        }
        for half_edge in sector
            .cap_half_edges
            .into_iter()
            .chain(sector.support_half_edges)
        {
            brep.half_edges[half_edge].pcurve =
                Some(Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0));
        }
    }

    let report = build_exact_boundary_report(
        &brep,
        ExactBoundaryCapCandidateStatus::Accepted,
        &exact_sectors,
        &[cap_face],
        &ExactBoundaryCandidatePerformanceReport {
            sector_count: exact_sectors.len(),
            pcurve_pair_count: exact_sectors.len() * 2,
            ..Default::default()
        },
        &[],
    );

    assert_eq!(
        report.candidate_status,
        ExactBoundaryCandidateOutcome::Accepted
    );
    assert_eq!(
        report.provenance.transaction_status,
        ExactBoundaryReportTransactionStatus::Committed
    );
    assert!(!report.provenance.topology_rollback_performed);
    assert_eq!(report.sectors.len(), 3);
    for sector in report.sectors {
        assert!(
            sector
                .cap_pcurves_present
                .into_iter()
                .all(|present| present)
        );
        assert!(
            sector
                .support_pcurves_present
                .into_iter()
                .all(|present| present)
        );
    }
}

#[test]
fn declined_report_records_rollback_and_versioned_diagnostics() {
    let (brep, cap_face, edges) = hex_cap_fixture();
    let plan = trihedral_plan(&brep, &edges);
    let sectors = trihedral_sectors(&edges);
    let exact_sectors =
        extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &sectors).unwrap();
    let diagnostic = ApexConstructionDiagnostic::ExactBoundaryCap(ExactBoundaryCapDiagnostic {
        stage: "exact_boundary_cap.ribbon".to_string(),
        code: "exact_boundary_cap.ribbon.support_pcurve_stale".to_string(),
        sector_name: "direct_apex.cap_sector.0".to_string(),
        edge: Some(edges[0]),
        support_half_edge: Some(exact_sectors[0].support_half_edges[0]),
        support_face: Some(exact_sectors[0].fillet_support_face),
        residual: Some(1.0),
        shared_fact_ids: Vec::new(),
    });

    let report = build_exact_boundary_report(
        &brep,
        ExactBoundaryCapCandidateStatus::Declined,
        &exact_sectors,
        &[],
        &ExactBoundaryCandidatePerformanceReport {
            sector_count: exact_sectors.len(),
            deterministic_skip_code: Some(
                "exact_boundary_cap.ribbon.support_pcurve_stale".to_string(),
            ),
            ..Default::default()
        },
        &[diagnostic],
    );

    assert_eq!(
        report.provenance.transaction_status,
        ExactBoundaryReportTransactionStatus::RolledBack
    );
    assert!(report.provenance.topology_rollback_performed);
    assert_eq!(report.diagnostics[0].schema_version, 1);
    assert_eq!(report.diagnostics[0].tolerance, None);
    assert_eq!(
        report.tolerance_evidence.deterministic_skip_code.as_deref(),
        Some("exact_boundary_cap.ribbon.support_pcurve_stale")
    );
}

fn hex_cap_fixture() -> (BRep, FaceHandle, Vec<EdgeHandle>) {
    let brep = BRep::n_gonal_pyramid(6, 1.0, 0.8, Tolerance::default()).unwrap();
    let cap_face = brep
        .faces
        .keys()
        .find(|&face| loop_edges(&brep, face).len() == 6)
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
