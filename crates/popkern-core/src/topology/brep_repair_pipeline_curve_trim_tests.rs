use super::*;
use crate::geometry::SurfaceGeometry;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, BrepSewingOptions, build_brep_repair_action_artifact};

fn options() -> BrepRepairPipelineOptions {
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
        degenerate_edges: BrepDegenerateEdgeOptions {
            tolerance,
            samples: 5,
        },
        require_all_pairs_repaired: true,
        complete_edge_curves: false,
        reconstruct_curve_trims: true,
        repair_edge_surfaces: false,
        repair_trims: false,
        diagnose_loop_shell: false,
        topology_output: BrepTopologyOutputPolicy::ClosedShell,
    }
}

#[test]
fn pipeline_runs_curve_backed_trim_reconstruction_before_edge_surface_repair() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let half_edge = brep.half_edges.keys().next().unwrap();
    bound_owning_face_to_existing_uv_box(&mut brep, half_edge);
    brep.half_edges[half_edge].pcurve = None;

    let result = repair_brep_boundary_edges_transactionally(&brep, &[], options()).unwrap();
    let report = result.report.curve_trim_reconstruction.unwrap();

    assert_eq!(report.reconstructed_half_edges, 1);
    assert!(result.brep.half_edges[half_edge].pcurve.is_some());
    assert!(brep.half_edges[half_edge].pcurve.is_none());
    let artifact = build_brep_repair_action_artifact(&result.report.actions);
    assert!(artifact.actions.iter().any(|action| {
        action.stage == "curve_trim_reconstruction"
            && action.outcome == "authored"
            && !action.shared_fact_ids.is_empty()
    }));
}

fn bound_owning_face_to_existing_uv_box(
    brep: &mut BRep,
    half_edge: crate::topology::HalfEdgeHandle,
) {
    let face = brep.loops[brep.half_edges[half_edge].loop_].face;
    let start = brep.loops[brep.faces[face].outer_loop].half_edge;
    let mut u_range = (f64::INFINITY, f64::NEG_INFINITY);
    let mut v_range = (f64::INFINITY, f64::NEG_INFINITY);
    let mut cursor = start;
    loop {
        if let Some(pcurve) = brep.half_edges[cursor].pcurve.as_ref() {
            for t in [pcurve.t_start, pcurve.t_end] {
                let uv = pcurve.point_at(t);
                u_range = (u_range.0.min(uv.x), u_range.1.max(uv.x));
                v_range = (v_range.0.min(uv.y), v_range.1.max(uv.y));
            }
        }
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            break;
        }
    }
    let surface = brep.faces[face].surface.as_ref().unwrap().geometry.clone();
    brep.faces[face].surface.as_mut().unwrap().geometry =
        SurfaceGeometry::Nurbs(surface.to_nurbs_in_range(u_range, v_range));
}
