use super::*;
use crate::geometry::curve::{CurveEval, CurveOnEdge};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{
    BRep, BrepRepairActionOutcome, BrepSewingOptions, BrepSewingProvenance, Edge, EdgeHandle, Face,
    FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, Shell, ShellHandle, Solid, Vertex,
    build_brep_repair_action_artifact,
};

fn pipeline_options(diagnose_loop_shell: bool) -> BrepRepairPipelineOptions {
    let clean_sewing = BrepSewingOptions {
        tolerance: 1e-9,
        samples: 3,
        require_surface_agreement: false,
        validate_output: false,
        allow_nonmanifold_edges: false,
    };
    BrepRepairPipelineOptions {
        classification: BrepSewingClassificationOptions {
            exact_tolerance: 1e-9,
            snap_tolerance: 1e-3,
            samples: 3,
            require_surface_agreement: false,
            allow_nonmanifold_edges: false,
        },
        snap_sewing: BrepSnapSewingOptions {
            max_vertex_motion: 1e-3,
            clean_sewing,
        },
        trim_repair: BrepTrimRepairOptions {
            endpoint_tolerance: 1e-6,
            sampled_tolerance: 1e-6,
            samples: 5,
            validate_output: true,
        },
        edge_curve_completion: BrepEdgeCurveCompletionOptions {
            tolerance: 1e-6,
            samples: 5,
            allow_vertex_segments_without_surface: true,
            require_all_completed: true,
            validate_output: true,
        },
        curve_trim_reconstruction: BrepCurveTrimReconstructionOptions {
            tolerance: 1e-6,
            samples: 5,
            replace_existing: false,
            replace_stale_existing: false,
            expand_bounded_surface_domains: false,
            require_all_reconstructed: false,
            validate_output: true,
        },
        edge_surface_compatibility: BrepEdgeSurfaceCompatibilityOptions {
            tolerance: 1e-6,
            max_anchor_motion: 0.0,
            samples: 7,
            require_all_compatible: true,
            validate_output: true,
        },
        degenerate_edges: BrepDegenerateEdgeOptions {
            tolerance: 1e-6,
            samples: 5,
        },
        require_all_pairs_repaired: true,
        complete_edge_curves: false,
        reconstruct_curve_trims: false,
        repair_edge_surfaces: false,
        repair_trims: true,
        diagnose_loop_shell,
        topology_output: BrepTopologyOutputPolicy::ClosedShell,
    }
}

fn open_sheet_pipeline_options() -> BrepRepairPipelineOptions {
    BrepRepairPipelineOptions {
        topology_output: BrepTopologyOutputPolicy::OpenSheet,
        ..pipeline_options(true)
    }
}

fn pair(first: EdgeHandle, second: EdgeHandle) -> BrepSewingEdgePair {
    BrepSewingEdgePair {
        first,
        second,
        provenance: BrepSewingProvenance {
            first: "pipeline left".to_string(),
            second: "pipeline right".to_string(),
        },
    }
}

#[test]
fn pipeline_repairs_clean_non_pns_pair_transactionally() {
    let (brep, left_edge, right_edge) = two_square_sheet();
    let result = repair_brep_boundary_edges_transactionally(
        &brep,
        &[pair(left_edge, right_edge)],
        pipeline_options(false),
    )
    .unwrap();

    assert_eq!(brep.boundary_edges.len(), 8, "input BRep is untouched");
    assert_eq!(result.report.candidate_pairs, 1);
    assert_eq!(result.report.policy_preflight.len(), 1);
    assert!(result.report.policy_preflight[0].is_accepted());
    assert_eq!(result.report.classification.unwrap().exact_clean_pairs, 1);
    assert_eq!(result.report.snap_sewing.unwrap().exact_pairs, 1);
    assert_eq!(result.report.actions.authored_actions, 1);
    assert_eq!(result.report.actions.declined_actions, 0);
    assert_eq!(
        result.report.actions.decisions[0].stage,
        BrepRepairPipelineStage::SnapSewing
    );
    assert!(
        !result.report.actions.decisions[0]
            .shared_fact_ids
            .is_empty(),
        "pipeline action should preserve snap/sewing shared facts"
    );
    assert_eq!(result.report.final_boundary_edges, Some(6));
    assert_eq!(result.brep.boundary_edges.len(), 6);
}

#[test]
fn pipeline_declines_nonmanifold_candidate_before_mutation() {
    let (mut brep, left_edge, right_edge) = two_square_sheet();
    let left_half_edge = brep.edges[left_edge].half_edge;
    brep.nonmanifold_radial
        .insert(left_half_edge, left_half_edge);

    let failure = repair_brep_boundary_edges_transactionally(
        &brep,
        &[pair(left_edge, right_edge)],
        pipeline_options(false),
    )
    .unwrap_err();

    assert_eq!(failure.stage, BrepRepairPipelineStage::PolicyPreflight);
    assert_eq!(
        failure.reason,
        BrepRepairPipelineFailureReason::PolicyDeclined
    );
    assert_eq!(failure.report.policy_preflight.len(), 1);
    assert!(
        failure.report.policy_preflight[0]
            .explanation
            .contains("radial non-manifold")
    );
    assert_eq!(brep.boundary_edges.len(), 8, "input BRep is untouched");
    assert_eq!(failure.report.actions.authored_actions, 0);
    assert_eq!(failure.report.actions.declined_actions, 0);
}

#[test]
fn pipeline_repairs_missing_trims_after_sewing_stage() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let half_edges = face_half_edges(&brep, face);
    for half_edge in &half_edges {
        brep.half_edges[*half_edge].pcurve = None;
    }

    let result =
        repair_brep_boundary_edges_transactionally(&brep, &[], pipeline_options(true)).unwrap();

    assert_eq!(result.report.candidate_pairs, 0);
    assert_eq!(result.report.snap_sewing.unwrap().classified_pairs, 0);
    let trim_report = result.report.trim_repair.unwrap();
    assert_eq!(trim_report.before.missing_pcurves, half_edges.len());
    assert!(
        result
            .report
            .actions
            .decisions
            .iter()
            .any(
                |decision| decision.stage == BrepRepairPipelineStage::TrimRepair
                    && !decision.shared_fact_ids.is_empty()
            ),
        "trim repair actions should preserve shared facts"
    );
    assert!(trim_report.after.is_clean());
    assert!(result.report.loop_shell.unwrap().has_no_errors());
    assert!(
        brep.half_edges[half_edges[0]].pcurve.is_none(),
        "input BRep is untouched"
    );
    assert!(result.brep.half_edges[half_edges[0]].pcurve.is_some());
}

#[test]
fn pipeline_reports_geometric_rejection_without_partial_output() {
    let (brep, left_edge, _) = two_square_sheet();
    let far_edge = brep
        .boundary_edges
        .iter()
        .copied()
        .find(|edge| *edge != left_edge && edge_segment_gap(&brep, left_edge, *edge) > 0.5)
        .unwrap();

    let failure = repair_brep_boundary_edges_transactionally(
        &brep,
        &[pair(left_edge, far_edge)],
        pipeline_options(false),
    )
    .unwrap_err();

    assert_eq!(failure.stage, BrepRepairPipelineStage::Classification);
    assert_eq!(
        failure.reason,
        BrepRepairPipelineFailureReason::ClassificationRejected
    );
    assert_eq!(brep.boundary_edges.len(), 8, "input BRep is untouched");
}

#[test]
fn pipeline_allows_open_sheet_boundary_diagnostics_under_open_sheet_policy() {
    let (brep, left_edge, right_edge) = two_square_sheet();
    let result = repair_brep_boundary_edges_transactionally(
        &brep,
        &[pair(left_edge, right_edge)],
        open_sheet_pipeline_options(),
    )
    .unwrap();

    let diagnostics = result.report.loop_shell.unwrap();
    assert!(diagnostics.boundary_edges > 0);
    assert!(diagnostics.has_no_errors());
    assert_eq!(result.report.final_boundary_edges, Some(6));
    assert_eq!(brep.boundary_edges.len(), 8, "input BRep is untouched");
}

#[test]
fn pipeline_rejects_same_open_sheet_when_closed_shell_required() {
    let (brep, left_edge, right_edge) = two_square_sheet();
    let failure = repair_brep_boundary_edges_transactionally(
        &brep,
        &[pair(left_edge, right_edge)],
        pipeline_options(true),
    )
    .unwrap_err();

    assert_eq!(failure.stage, BrepRepairPipelineStage::LoopShellDiagnostics);
    assert_eq!(
        failure.reason,
        BrepRepairPipelineFailureReason::LoopShellErrors
    );
    assert!(failure.report.loop_shell.unwrap().boundary_edges > 0);
    assert_eq!(brep.boundary_edges.len(), 8, "input BRep is untouched");
}

#[test]
fn pipeline_keeps_open_sheet_fragmentation_warning_observable() {
    let (brep, _, _) = two_square_sheet();
    let result =
        repair_brep_boundary_edges_transactionally(&brep, &[], open_sheet_pipeline_options())
            .unwrap();

    let diagnostics = result.report.loop_shell.unwrap();
    assert!(diagnostics.findings.iter().any(|finding| {
        finding.kind == BrepTopologyDiagnosticKind::FragmentedShell
            && finding.severity == BrepTopologyDiagnosticSeverity::Warning
    }));
    assert!(diagnostics.has_no_errors());
}

#[test]
fn pipeline_repairs_edge_surface_compatibility_after_trim_stages() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    corrupt_edge_curve(&mut brep, edge, Vec3::new(0.11, 0.07, 0.13));
    let mut options = pipeline_options(true);
    options.repair_edge_surfaces = true;
    options.repair_trims = false;

    let result = repair_brep_boundary_edges_transactionally(&brep, &[], options).unwrap();

    let report = result.report.edge_surface_compatibility.unwrap();
    assert_eq!(report.before.incompatible_half_edges, 2);
    assert_eq!(report.repaired_edges, 1);
    assert_eq!(report.after.incompatible_half_edges, 0);
    assert!(result.report.actions.decisions.iter().any(|decision| {
        decision.stage == BrepRepairPipelineStage::EdgeSurfaceCompatibility
            && decision.outcome == BrepRepairActionOutcome::Authored
    }));
    let artifact = build_brep_repair_action_artifact(&result.report.actions);
    assert!(artifact.actions.iter().any(|action| {
        action.stage == "edge_surface_compatibility" && action.outcome == "authored"
    }));
    assert!(
        brep.edges[edge]
            .curve
            .as_ref()
            .unwrap()
            .geometry
            .point_at(0.0)
            .distance_to(
                &result.brep.edges[edge]
                    .curve
                    .as_ref()
                    .unwrap()
                    .geometry
                    .point_at(0.0)
            )
            > 1e-2,
        "input BRep is untouched"
    );
}

#[test]
fn pipeline_reports_unresolved_edge_surface_compatibility_without_partial_output() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    corrupt_edge_curve(&mut brep, edge, Vec3::new(0.11, 0.07, 0.13));
    let half_edge = brep.edges[edge].half_edge;
    let face = brep.loops[brep.half_edges[half_edge].loop_].face;
    brep.faces[face].surface = None;
    let mut options = pipeline_options(false);
    options.repair_edge_surfaces = true;
    options.repair_trims = false;

    let failure = repair_brep_boundary_edges_transactionally(&brep, &[], options).unwrap_err();

    assert_eq!(
        failure.stage,
        BrepRepairPipelineStage::EdgeSurfaceCompatibility
    );
    assert!(matches!(
        failure.reason,
        BrepRepairPipelineFailureReason::EdgeSurfaceCompatibilityFailed(_)
    ));
    assert_eq!(
        failure
            .report
            .edge_surface_compatibility
            .unwrap()
            .repaired_edges,
        0
    );
    assert!(
        failure
            .report
            .actions
            .decisions
            .iter()
            .any(|decision| decision.outcome == BrepRepairActionOutcome::RolledBack)
    );
    let artifact = build_brep_repair_action_artifact(&failure.report.actions);
    assert!(artifact.actions.iter().any(|action| {
        action.outcome == "rolled_back"
            && action.reason.as_deref() == Some("unresolved_edge_surface_compatibility")
    }));
}

fn two_square_sheet() -> (BRep, EdgeHandle, EdgeHandle) {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: ShellHandle::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: FaceHandle::default(),
        solid,
    });
    brep.solids[solid].shell = shell;

    let left = add_square(
        &mut brep,
        shell,
        [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
    );
    let right = add_square(
        &mut brep,
        shell,
        [
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
    );
    brep.shells[shell].face = left.face;
    (brep, left.edges[1], right.edges[2])
}

struct SquareHandles {
    face: FaceHandle,
    edges: [EdgeHandle; 4],
}

fn add_square(brep: &mut BRep, shell: ShellHandle, points: [Point3; 4]) -> SquareHandles {
    let vertices = points.map(|position| {
        brep.vertices.insert(Vertex {
            position,
            outgoing: HalfEdgeHandle::default(),
            name: None,
            tolerance: None,
        })
    });
    let face = brep.faces.insert(Face {
        outer_loop: LoopHandle::default(),
        inner_loops: Vec::new(),
        shell,
        surface: None,
        name: None,
        tolerance: None,
        role: Default::default(),
    });
    let loop_ = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face,
    });
    brep.faces[face].outer_loop = loop_;

    let half_edges = vertices.map(|origin| {
        brep.half_edges.insert(HalfEdge {
            origin,
            twin: HalfEdgeHandle::default(),
            next: HalfEdgeHandle::default(),
            prev: HalfEdgeHandle::default(),
            edge: EdgeHandle::default(),
            loop_,
            pcurve: None,
        })
    });
    let edges = half_edges.map(|half_edge| {
        let edge = brep.edges.insert(Edge {
            half_edge,
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });
        brep.half_edges[half_edge].edge = edge;
        brep.boundary_edges.insert(edge);
        edge
    });

    for index in 0..4 {
        let prev = (index + 3) % 4;
        let next = (index + 1) % 4;
        brep.half_edges[half_edges[index]].prev = half_edges[prev];
        brep.half_edges[half_edges[index]].next = half_edges[next];
        brep.vertices[vertices[index]].outgoing = half_edges[index];
    }
    brep.loops[loop_].half_edge = half_edges[0];
    SquareHandles { face, edges }
}

fn face_half_edges(brep: &BRep, face: FaceHandle) -> Vec<HalfEdgeHandle> {
    let start = brep.loops[brep.faces[face].outer_loop].half_edge;
    let mut out = Vec::new();
    let mut cursor = start;
    loop {
        out.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            break;
        }
    }
    out
}

fn edge_segment_gap(brep: &BRep, first: EdgeHandle, second: EdgeHandle) -> f64 {
    let first = crate::topology::edge_segment(brep, first).unwrap();
    let second = crate::topology::edge_segment(brep, second).unwrap();
    first
        .0
        .distance_to(&second.1)
        .max(first.1.distance_to(&second.0))
}

fn corrupt_edge_curve(brep: &mut BRep, edge: EdgeHandle, offset: Vec3) {
    let curve = brep.edges[edge].curve.as_ref().unwrap();
    let start = curve.geometry.point_at(curve.t_start) + offset;
    let end = curve.geometry.point_at(curve.t_end) + offset;
    brep.edges[edge].curve = Some(
        CurveOnEdge::line_segment(start, end, 1e-12)
            .expect("offset keeps cube edge non-degenerate"),
    );
}
