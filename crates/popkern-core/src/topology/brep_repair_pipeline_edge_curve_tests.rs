use super::*;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, BrepRepairActionOutcome, BrepRepairTopologyHandle, BrepSewingOptions};

fn edge_curve_options() -> BrepRepairPipelineOptions {
    let tolerance = 1.0e-6;
    BrepRepairPipelineOptions {
        classification: BrepSewingClassificationOptions {
            exact_tolerance: tolerance,
            snap_tolerance: tolerance,
            samples: 3,
            require_surface_agreement: false,
            allow_nonmanifold_edges: false,
        },
        snap_sewing: BrepSnapSewingOptions {
            max_vertex_motion: 0.0,
            clean_sewing: BrepSewingOptions {
                tolerance,
                samples: 3,
                require_surface_agreement: false,
                validate_output: false,
                allow_nonmanifold_edges: false,
            },
        },
        trim_repair: BrepTrimRepairOptions {
            endpoint_tolerance: tolerance,
            sampled_tolerance: tolerance,
            samples: 5,
            validate_output: true,
        },
        edge_curve_completion: BrepEdgeCurveCompletionOptions {
            tolerance,
            samples: 5,
            allow_vertex_segments_without_surface: true,
            require_all_completed: true,
            validate_output: true,
        },
        curve_trim_reconstruction: BrepCurveTrimReconstructionOptions {
            tolerance,
            samples: 5,
            replace_existing: false,
            replace_stale_existing: false,
            expand_bounded_surface_domains: false,
            require_all_reconstructed: false,
            validate_output: true,
        },
        edge_surface_compatibility: BrepEdgeSurfaceCompatibilityOptions {
            tolerance,
            max_anchor_motion: 0.0,
            samples: 5,
            require_all_compatible: false,
            validate_output: true,
        },
        degenerate_edges: BrepDegenerateEdgeOptions {
            tolerance,
            samples: 5,
        },
        require_all_pairs_repaired: true,
        complete_edge_curves: true,
        reconstruct_curve_trims: false,
        repair_edge_surfaces: false,
        repair_trims: false,
        diagnose_loop_shell: true,
        topology_output: BrepTopologyOutputPolicy::ClosedShell,
    }
}

#[test]
fn pipeline_reports_edge_curve_completion_with_affected_edges() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    brep.edges[edge].curve = None;

    let result =
        repair_brep_boundary_edges_transactionally(&brep, &[], edge_curve_options()).unwrap();

    let completion = result.report.edge_curve_completion.unwrap();
    assert_eq!(completion.completed_edges, 1);
    assert!(result.report.actions.decisions.iter().any(|decision| {
        decision.stage == BrepRepairPipelineStage::EdgeCurveCompletion
            && decision.outcome == BrepRepairActionOutcome::Authored
            && decision.affected_handles == vec![BrepRepairTopologyHandle::Edge(edge)]
            && !decision.shared_fact_ids.is_empty()
    }));
}
