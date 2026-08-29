use super::boundary::{ExactBoundaryCurveBudget, prepare_exact_boundary_curves};
use super::corner::{ExactBoundaryInnerSpokeBudget, prepare_exact_boundary_inner_spokes_with_pole};
use super::facade::{ExactBoundaryCapCandidateStatus, run_exact_boundary_cap_candidate};
use super::sector::{
    ExactBoundaryOuterSide, ExactBoundarySectorDeclineReason, extract_exact_boundary_cap_sectors,
};
use super::topology::materialize_exact_boundary_sector_topology;
use crate::apex_construction::cap_kind::CapKind;
use crate::apex_construction::gregory_cap::GregoryCapSector;
use crate::apex_construction::plan::ApexPlan;
use crate::apex_construction::target_topology::{
    TargetCapFace, TargetEdge, TargetEdgeRole, TargetTopology,
};
use crate::apex_construction::telemetry::PlanTelemetry;
use crate::types::ExactBoundaryCandidateOutcome;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

#[test]
fn extracts_four_sector_cap_with_stable_names_and_ownership() {
    let (brep, cap_face, edges) = n_sector_cap_fixture(4, 1.0);
    let plan = n_sector_plan(&brep, &edges, 4);
    let sectors = n_sector_inventory(&edges, 4);

    let exact_sectors =
        extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &sectors).unwrap();

    assert_eq!(exact_sectors.len(), 4);
    for (index, sector) in exact_sectors.iter().enumerate() {
        assert_eq!(sector.sector_index, index);
        assert_eq!(sector.fillet_index, index);
        assert_eq!(
            sector.sector_name,
            format!("direct_apex.cap_sector.{index}")
        );
        assert_eq!(sector.cap_face, cap_face);
        assert_ne!(Some(sector.fillet_end_edge), sector.face_edge);
        assert_eq!(sector.fillet_seam.side, ExactBoundaryOuterSide::FilletEnd);
        assert_eq!(
            sector.bridge_seam.as_ref().expect("bridged fixture").side,
            ExactBoundaryOuterSide::FaceEdge
        );
        assert_eq!(sector.fillet_seam.edge, edges[2 * index]);
        assert_eq!(
            sector.bridge_seam.as_ref().expect("bridged fixture").edge,
            edges[2 * index + 1]
        );
        assert_ne!(sector.fillet_support_face, cap_face);
        assert_ne!(sector.bridge_support_face, Some(cap_face));
        assert_ne!(
            sector.fillet_seam.cap_half_edge,
            sector
                .bridge_seam
                .as_ref()
                .expect("bridged fixture")
                .cap_half_edge
        );
        assert_ne!(
            sector.fillet_seam.support_half_edge,
            sector
                .bridge_seam
                .as_ref()
                .expect("bridged fixture")
                .support_half_edge
        );
    }
}

#[test]
fn four_sector_inner_spokes_and_topology_preserve_cyclic_ownership() {
    let (mut brep, cap_face, edges) = n_sector_cap_fixture(4, 1.0);
    let plan = n_sector_plan(&brep, &edges, 4);
    let sectors = n_sector_inventory(&edges, 4);
    let exact_sectors =
        extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &sectors).unwrap();
    let boundary_curves =
        prepare_exact_boundary_curves(&brep, &exact_sectors, ExactBoundaryCurveBudget::default())
            .unwrap();

    let spokes = prepare_exact_boundary_inner_spokes_with_pole(
        &exact_sectors,
        &boundary_curves,
        ExactBoundaryInnerSpokeBudget::default(),
        None,
    )
    .unwrap();

    assert_eq!(spokes.sectors.len(), 4);
    for index in 0..spokes.sectors.len() {
        let current = &spokes.sectors[index];
        let next = &spokes.sectors[(index + 1) % spokes.sectors.len()];
        assert_eq!(
            current.sector_name,
            format!("direct_apex.cap_sector.{index}")
        );
        assert!(
            current
                .face_inner_spoke
                .name
                .contains("face_anchor_to_pole")
        );
        assert!(
            current
                .fillet_inner_spoke
                .name
                .contains("fillet_anchor_to_pole")
        );
        assert!(
            current
                .face_inner_spoke
                .anchor
                .distance_to(&next.fillet_inner_spoke.anchor)
                <= 1.0e-7
        );
        assert_eq!(current.face_inner_spoke.pole, spokes.pole);
        assert_eq!(current.fillet_inner_spoke.pole, spokes.pole);
    }

    let topology =
        materialize_exact_boundary_sector_topology(&mut brep, cap_face, &exact_sectors, &spokes)
            .unwrap();

    assert_eq!(topology.sectors.len(), 4);
    assert_eq!(topology.spoke_edges.len(), 4);
    unique_faces(&topology.sectors);
    brep.validate().unwrap();
}

#[test]
fn four_sector_candidate_reports_linear_counts_and_structured_status() {
    let (mut brep, cap_face, edges) = n_sector_cap_fixture(4, 1.0);
    let plan = n_sector_plan(&brep, &edges, 4);
    let sectors = n_sector_inventory(&edges, 4);

    let result = run_exact_boundary_cap_candidate(&mut brep, cap_face, &plan, &sectors)
        .expect("four-sector candidate should produce a structured result");

    assert_eq!(result.sectors.len(), 4);
    assert_eq!(result.performance.sector_count, 4);
    assert_eq!(result.prepared_curve_count(), 8);
    assert_eq!(result.inner_spoke_count(), 8);
    assert_eq!(result.ribbon_count(), 8);
    assert!(result.performance.support_sample_count >= 40);
    assert!(result.performance.elapsed_micros.is_some_and(|us| us > 0));
    match result.status {
        ExactBoundaryCapCandidateStatus::Accepted => {
            assert_eq!(
                result.report.candidate_status,
                ExactBoundaryCandidateOutcome::Accepted
            );
            assert_eq!(result.committed_faces.len(), 4);
            assert_eq!(result.report.committed_faces.len(), 4);
            assert_eq!(result.performance.validation_sector_count, 4);
        }
        ExactBoundaryCapCandidateStatus::Declined => {
            assert_eq!(
                result.report.candidate_status,
                ExactBoundaryCandidateOutcome::Declined
            );
            assert!(!result.diagnostics.is_empty());
            assert!(!result.report.diagnostics.is_empty());
            assert!(result.performance.deterministic_skip_code.is_some());
            assert!(result.committed_faces.is_empty());
        }
    }
    assert_eq!(result.report.schema_version, 1);
    assert_eq!(result.report.sectors.len(), 4);
    assert_eq!(result.report.performance, result.performance);
    assert!(result.report.preview_available);
    assert!(
        result
            .report
            .tolerance_evidence
            .model_linear_tolerance
            .is_some()
    );
}

#[test]
fn small_angle_sector_candidate_keeps_bounded_linear_evidence() {
    let sector_count = 12;
    let (mut brep, cap_face, edges) = n_sector_cap_fixture(sector_count, 1.0);
    let plan = n_sector_plan(&brep, &edges, sector_count);
    let sectors = n_sector_inventory(&edges, sector_count);

    let result = run_exact_boundary_cap_candidate(&mut brep, cap_face, &plan, &sectors)
        .expect("small-angle N-sector candidate should produce a structured result");

    assert_eq!(result.performance.sector_count, sector_count);
    assert_eq!(result.prepared_curve_count(), sector_count * 2);
    assert_eq!(
        result.sectors[sector_count - 1].sector_name,
        "direct_apex.cap_sector.11"
    );
    assert!(result.performance.elapsed_micros.is_some_and(|us| us > 0));
    if result.status == ExactBoundaryCapCandidateStatus::Declined {
        assert!(!result.diagnostics.is_empty());
        assert!(result.performance.deterministic_skip_code.is_some());
    }
}

#[test]
fn incompatible_four_sector_inventory_declines_with_stable_diagnostic() {
    let (brep, cap_face, edges) = n_sector_cap_fixture(4, 1.0);
    let plan = n_sector_plan(&brep, &edges, 4);
    let mut sectors = n_sector_inventory(&edges, 4);
    sectors[2].face_edge = Some(sectors[2].fillet_end_edge);

    let diagnostic =
        extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &sectors).unwrap_err();

    assert_eq!(
        diagnostic.reason,
        ExactBoundarySectorDeclineReason::DuplicateOuterSeam
    );
    assert_eq!(
        diagnostic.code,
        "exact_boundary_cap.sector.duplicate_outer_seam"
    );
    assert_eq!(diagnostic.sector_index, Some(2));
    assert_eq!(
        diagnostic.sector_name,
        Some("direct_apex.cap_sector.2".to_string())
    );
    assert_eq!(diagnostic.cap_face, Some(cap_face));
    assert_eq!(diagnostic.edge, Some(edges[4]));
}

fn n_sector_cap_fixture(sector_count: usize, radius: f64) -> (BRep, FaceHandle, Vec<EdgeHandle>) {
    let edge_count = sector_count * 2;
    let brep = BRep::n_gonal_pyramid(edge_count, radius, 0.8, Tolerance::default()).unwrap();
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

fn unique_faces(sectors: &[super::sector::ExactBoundaryCapSectorInput]) {
    for index in 0..sectors.len() {
        for other in (index + 1)..sectors.len() {
            assert_ne!(sectors[index].cap_face, sectors[other].cap_face);
        }
    }
}
