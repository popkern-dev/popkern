use super::*;
use crate::tolerance::Tolerance;
use crate::topology::{
    BRep, BrepRepairActionOutcome, BrepRepairTopologyHandle, BrepSewingOptions,
    analyze_brep_face_use_trim_authority, build_brep_repair_action_artifact,
};

fn face_use_options() -> BrepRepairPipelineOptions {
    let tolerance = 1e-6;
    BrepRepairPipelineOptions {
        classification: BrepSewingClassificationOptions {
            exact_tolerance: tolerance,
            snap_tolerance: tolerance,
            samples: 5,
            require_surface_agreement: false,
            allow_nonmanifold_edges: false,
        },
        snap_sewing: BrepSnapSewingOptions {
            max_vertex_motion: 0.0,
            clean_sewing: BrepSewingOptions {
                tolerance,
                samples: 5,
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
            require_all_completed: false,
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
        face_use_trim_authority: BrepFaceUseTrimAuthorityOptions {
            tolerance,
            samples: 5,
        },
        degenerate_edges: BrepDegenerateEdgeOptions {
            tolerance,
            samples: 5,
        },
        require_all_pairs_repaired: true,
        complete_edge_curves: false,
        reconstruct_curve_trims: false,
        repair_edge_surfaces: false,
        repair_face_use_trim_authority: true,
        repair_trims: false,
        diagnose_loop_shell: false,
        topology_output: BrepTopologyOutputPolicy::ClosedShell,
    }
}

#[test]
fn pipeline_materializes_face_use_trim_authority_transactionally() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let faces = brep.faces.keys().take(2).collect::<Vec<_>>();
    let reused_loop = brep.faces[faces[0]].outer_loop;
    brep.faces[faces[1]].outer_loop = reused_loop;
    brep.faces[faces[1]].surface = brep.faces[faces[0]].surface.clone();

    let before =
        analyze_brep_face_use_trim_authority(&brep, face_use_options().face_use_trim_authority);
    assert!(before.has_blocking_findings());

    let result = repair_brep_boundary_edges_transactionally(&brep, &[], face_use_options())
        .expect("same-surface reused loop should be repaired by the generic pipeline");

    let report = result.report.face_use_trim_authority.unwrap();
    assert_eq!(
        report.before.loop_ownership_mismatches,
        before.loop_ownership_mismatches
    );
    assert!(report.materialized_face_uses > 0);
    assert!(report.after.is_clean(), "{:#?}", report.after);
    assert!(
        analyze_brep_face_use_trim_authority(&brep, face_use_options().face_use_trim_authority)
            .has_blocking_findings(),
        "pipeline must repair a clone and leave input untouched"
    );
    assert_eq!(
        result.brep.face_use_pcurves.len(),
        report.materialized_face_uses
    );
    assert!(result.report.actions.decisions.iter().any(|decision| {
        decision.stage == BrepRepairPipelineStage::FaceUseTrimAuthority
            && decision.outcome == BrepRepairActionOutcome::Authored
            && !decision.shared_fact_ids.is_empty()
            && decision.affected_handles.iter().any(|handle| {
                matches!(
                    handle,
                    BrepRepairTopologyHandle::Face(face) if *face == faces[1]
                )
            })
    }));
    let artifact = build_brep_repair_action_artifact(&result.report.actions);
    assert!(artifact.actions.iter().any(|action| {
        action.stage == "face_use_trim_authority"
            && action.outcome == "authored"
            && !action.shared_fact_ids.is_empty()
    }));
}
