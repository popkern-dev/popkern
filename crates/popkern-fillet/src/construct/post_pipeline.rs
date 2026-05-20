//! Late post-processing passes for applied blend topology.

use super::*;

pub(super) struct PostBlendPipelineInputs<'a> {
    pub(super) source_brep: &'a BRep,
    pub(super) result: &'a mut BRep,
    pub(super) edge_geometries: &'a [EdgeGeometry],
    pub(super) tolerance: &'a Tolerance,
    pub(super) late_wire_jobs: &'a [(crate::blend::BlendBoundary, Vec<FaceHandle>)],
    pub(super) support_seam_jobs: &'a [(crate::blend::BlendBoundary, Vec<FaceHandle>)],
    pub(super) direct_apex_built_faces: &'a [FaceHandle],
    pub(super) warnings: &'a mut Vec<(EdgeHandle, BlendWarning)>,
    pub(super) corner_patches: &'a mut Vec<crate::types::CornerPatchInfo>,
    pub(super) repair_report: &'a mut crate::types::BlendRepairReport,
}

pub(super) fn run_post_blend_pipeline(inputs: PostBlendPipelineInputs<'_>) {
    let brep = inputs.source_brep;
    let result = inputs.result;
    let edge_geometries = inputs.edge_geometries;
    let tolerance = inputs.tolerance;
    let late_wire_jobs = inputs.late_wire_jobs;
    let support_seam_jobs = inputs.support_seam_jobs;
    let direct_apex_built_faces = inputs.direct_apex_built_faces;
    let warnings = inputs.warnings;
    let corner_patches = inputs.corner_patches;
    let repair_report = inputs.repair_report;

    // Correct surface assignments — safety net.
    let fillet_surfaces: Vec<_> = edge_geometries
        .iter()
        .map(|eg| eg.surface.clone())
        .collect();
    surface_fix::correct_face_surfaces(result, brep, &fillet_surfaces);

    // Reconcile edge curves with their adjacent face surfaces.
    // multi_chord and corner-blend mef calls write straight-line
    // edge curves that were correct against the original planar
    // face but stale once surface_fix replaced the surface with a
    // cylinder. Replace each cylinder × plane edge curve with the
    // analytic intersection, trimmed to the existing endpoints.
    // Strict tolerance — only rewrites that hit existing vertex
    // positions to within `tolerance.linear()`. Off-cylinder
    // vertices (apex-only filleting case) keep their chord-line
    // edges; topology_repair handles those separately.
    edge_reconcile::reconcile_edge_curves(result);

    // Late-pipeline clamp: any half-edge pcurve on a bounded-UV
    // surface (Gregory/NURBS/sphere/torus) whose evaluation samples
    // outside the surface's parameter range gets its control points
    // clamped into the box. Catches stale pcurves attached after
    // Gregory cap surface assignment by other pipeline stages
    // (e.g. attach_fillet_pcurves authoring in the parent face's UV
    // space). Reduces tess-nm divergence on n04 from ~0.25 → smaller.
    // Optional global pcurve re-fit (`POPKERN_PCURVE_REFIT=1`):
    // re-project every half-edge's pcurve onto its current face
    // surface. Currently regresses tess-nm in `fillet_zoo`
    // (likely because the closest-point + degree-3 fit produces
    // noisier curves than the carefully-authored analytic pcurves
    // on plane / cylinder pairs); keep the pass for investigation
    // and selective use, but don't run it by default.
    if std::env::var("POPKERN_PCURVE_REFIT").is_ok() {
        pcurve_refit::refit_all_pcurves(result);
    }

    // Topology repair for off-cylinder vertices in cylinder face
    // outer loops (apex-only filleting case where the cylinder
    // ribbon's base-side ends at an unfilleted original vertex via
    // a chord pair). Each repair attempt runs in a snapshot and
    // must pass `BRep::validate`; failed attempts roll back without
    // affecting downstream operations.
    //
    // Default-on only for apex boundary contexts, where leaving
    // off-cylinder vertices in analytic cylinder loops poisons
    // Gregory boundary evidence and STEP pcurve fidelity. Simple
    // single-edge fillets keep their historical topology so the
    // current algorithmic unfillet path does not need to understand
    // repair sliver faces. `POPKERN_TOPO_REPAIR=1` remains the
    // broad diagnostic mode; `POPKERN_TOPO_REPAIR=0` disables even
    // the apex-context default. Different-adjacent-plane chord
    // chains remain separately gated inside `topology_repair`
    // because they require stronger multi-face ownership repair
    // before they can be commercial grade.
    let apex_boundary_context = !support_seam_jobs.is_empty()
        || !late_wire_jobs.is_empty()
        || !direct_apex_built_faces.is_empty()
        || corner_patches
            .iter()
            .any(|patch| matches!(patch.method, crate::types::CornerMethod::Gregory));
    if topology_repair_enabled(apex_boundary_context) {
        let repaired = topology_repair::repair_cylinder_face_off_surface_vertices(result);
        if std::env::var_os("POPKERN_TOPO_REPAIR_TRACE").is_some() {
            eprintln!("[topo_repair] default post-blend repaired {repaired} chord pair(s)");
        }
    }

    // §13.1(2) sliver absorption: Classes A/B/C/D. 12-edge cube test is Class D regression guard (26 faces preserved).
    //
    // Phase 6c-D Stage 4.5: protect substituted-boundary cap faces
    // from absorption.  When sphere-fit injection has fired, each
    // cap's outer loop is the canonical 3-edge per-fillet topology
    // that the injected `hex_boundaries[k]` describes; absorbing
    // adjacent Plane fragments would expand the loop across 2-3
    // fillets' boundaries and the Gregory cap surface (row 0 = ONE
    // FilletEnd) cannot represent that.  The protection set is
    // empty when no vertex took the substitution path (all other
    // fixtures behave identically to before).
    let mut protected_caps: Vec<popkern_core::topology::FaceHandle> = late_wire_jobs
        .iter()
        .flat_map(|(_, faces)| faces.iter().copied())
        .collect();
    // §6c-direct-apex dispatch.d — merge direct-apex-built cap faces into the protected
    // set so `sliver_absorb` / `cap_insert` / `late_wire_jobs`
    // don't re-author them.  When `use_direct_apex_construction` is
    // off this list is empty.
    protected_caps.extend(direct_apex_built_faces.iter().copied());
    sliver_absorb::run_and_route_warnings_with_protected(
        result,
        tolerance,
        edge_geometries,
        &protected_caps,
        warnings,
    );

    cap_insert::insert_caps(result, brep, edge_geometries, warnings); // §8.3 caps

    // Phase 6c-D Stage 4.4: re-run the substituted-boundary wire
    // step now that sliver_absorb and cap_insert have settled the
    // cap-face outer-loop topology. The wire originally runs
    // mid-pipeline (inside `three_quad_apply`) and only sees the
    // raw poke output's outer loop. After sliver_absorb merges
    // adjacent Plane fragments INTO each cap face, the cap face's
    // outer loop expands with the merged fragments' boundaries —
    // including new cap-vs-cylinder edges that the wire never
    // visited. Without this second pass, those edges retain their
    // original (rolling-ball footprint) curves and `cap_outer_g0`
    // stays at ~0.12 m on pyramid4_mixed under injection.
    //
    // The second wire pass uses the same `BlendBoundary` (saved by
    // `insert_vertex_blend`) and the same per-vertex sub-face
    // handles. Stage 4.3's cap-face-membership filter naturally
    // excludes cap-vs-cap edges, so only newly-absorbed cap-vs-non-
    // cap edges get re-wired; previously-wired edges get their
    // analytic curves re-written (idempotent).
    for (boundary, sub_faces) in late_wire_jobs {
        let wired =
            sphere_fit_edge_wire::wire_substituted_boundary_to_edges(result, boundary, sub_faces);
        if std::env::var_os("POPKERN_KP_TRACE").is_some() {
            eprintln!(
                "[late-wire] re-ran on {} sub_faces, wrote {wired} edges",
                sub_faces.len()
            );
        }
    }

    // Phase 6c-D Stage 4.7: lift each apex-region pyramid Plane
    // face to a curved NURBS surface that contains the substituted
    // FaceEdge bridge.  Without this, the bridge sits 2-4 cm off
    // the original plane, and the plane's pcurve refit can't close
    // the cap-vs-twin G0 gap (every plane UV maps to a point on
    // the plane, never to the off-plane bridge).
    //
    // On by default since #207 cleanup.  Set
    // `POPKERN_APEX_FACE_LIFT=0` to disable for debugging.
    if std::env::var("POPKERN_APEX_FACE_LIFT")
        .map(|v| v != "0")
        .unwrap_or(true)
    {
        for (_boundary, sub_faces) in late_wire_jobs {
            let lifted = apex_face_lift::lift_all_apex_plane_faces(result, sub_faces);
            if std::env::var_os("POPKERN_KP_TRACE").is_some() {
                eprintln!("[apex-face-lift] lifted {lifted} pyramid plane faces");
            }
        }
    }

    // Phase 6c-X.A.4: relocate apex vertices to match the substituted
    // boundary's termination points T.  After A.3 (full-range wire),
    // every cap-loop edge's curve agrees on T at each apex vertex;
    // this pass moves the BRep vertex to T.  After A.4, the
    // vertex-to-curve-endpoint coherence on cap-loop edges is exact;
    // BASE-EDGE curves still terminate at the OLD vertex position
    // and need A.5's curve-replacement pass.
    if std::env::var("POPKERN_APEX_PATH_A")
        .map(|v| v != "0")
        .unwrap_or(true)
    {
        for (boundary, sub_faces) in late_wire_jobs {
            let (n_moved, worst) =
                apex_vertex_coherence::relocate_apex_vertices(result, sub_faces, boundary);
            // Phase 6c-X.A.5: rebuild BASE-EDGE curves (cyl-vs-Plane
            // apex generators going AWAY from the apex through V) as
            // straight lines from the new V→T position to v_other.
            // Restores vertex↔curve coherence on the BASE-EDGE side.
            let n_repaired = apex_vertex_coherence::repair_base_edges(result, sub_faces);
            // Phase 6c-X.A.5b: realign cap-loop edges' curve direction
            // to match BRep convention (curve(t_start) =
            // primary_origin's vertex position).  Validate_watertight
            // demands this strict alignment; the wire's stored
            // direction was set by an old vertex-position heuristic
            // that no longer applies after A.4's relocations.
            let n_realigned =
                apex_vertex_coherence::realign_cap_loop_curve_directions(result, sub_faces);
            if std::env::var_os("POPKERN_KP_TRACE").is_some() {
                eprintln!(
                    "[apex-path-a] relocated {n_moved} vertices, worst V→T = {worst:.3e} m; \
                     repaired {n_repaired} base edges; realigned {n_realigned} cap-loop edges"
                );
            }
        }
    }

    // Phase 6c-C.2d — apex-region topology rebuild for sphere-fit
    // caps (env-gated).  Replaces the 6-vertex cap loop with a
    // 3-edge loop of cap-tangent arcs by running
    // plan → mef-execute → lens-collapse on each corner patch's
    // cap face.  Default off; opt-in via `POPKERN_APEX_TOPOLOGY_REBUILD=1`
    // until the validation matrix flips it on by default
    // (Phase 6c-F).
    if std::env::var("POPKERN_APEX_TOPOLOGY_REBUILD")
        .map(|v| v == "1")
        .unwrap_or(false)
    {
        let cap_faces: Vec<popkern_core::topology::FaceHandle> =
            corner_patches.iter().map(|p| p.face).collect();
        let _n = crate::vertex_blend_topology::rebuild::rebuild_all_apex_caps(result, &cap_faces);
    }

    // Phase 6c-A continuation §4.2 surgery — variable-radius cap-
    // region re-routing. Default off until validated.
    //
    // Source the cap-face list and the per-fillet sphere parameters
    // from `corner_patches` (populated during apply_blend) rather
    // than from the BRep's `Sphere` faces — this lets the surgery
    // run on Gregory cap faces too (Phase 6c-D Gregory cap rung)
    // without requiring the sphere-cap rung to be applied first.
    if std::env::var("POPKERN_VAR_RADIUS_REROUTE")
        .map(|v| v == "1")
        .unwrap_or(false)
    {
        let mut cap_faces: Vec<popkern_core::topology::FaceHandle> = Vec::new();
        let mut per_fillet_caps: Vec<(popkern_core::geometry::Point3, f64)> = Vec::new();
        for cp in corner_patches.iter() {
            if let Some(sphere) = cp.cap_sphere {
                cap_faces.push(cp.face);
                per_fillet_caps.push(sphere);
            }
        }
        if cap_faces.len() == per_fillet_caps.len() && !cap_faces.is_empty() {
            let r = crate::vertex_blend_topology::var_radius::reroute_neighbour_cylinder_edges(
                result,
                &cap_faces,
                &per_fillet_caps,
            );
            if std::env::var("POPKERN_VAR_RADIUS_TRACE").is_ok_and(|v| v == "1") {
                eprintln!("[var-radius-reroute] caps={} result={r:?}", cap_faces.len());
            }
        }
    }

    if support_seam_relocation_enabled() {
        for (boundary, sub_faces) in support_seam_jobs {
            let report = three_quad_apply::author_model_level_support_face_boundary_seams(
                result,
                sub_faces,
                boundary,
                &[],
            );
            if std::env::var_os("POPKERN_KP_TRACE").is_some() {
                eprintln!(
                    "[support-seam-relocation] sub_faces={} authored={} declined={} \
                     candidates={} max_move={:.6e} max_disagreement={:.6e} rollback={:?}",
                    sub_faces.len(),
                    report.authored,
                    report.declined,
                    report.candidates,
                    report.max_endpoint_motion,
                    report.max_sampled_disagreement,
                    report.rollback_reason
                );
                for decision in report.decisions.iter().take(3) {
                    eprintln!(
                        "[support-seam-relocation] decision he={:?} cap={:?} support={:?} \
                         side={:?} accepted={} reason={:?}",
                        decision.half_edge,
                        decision.cap_face,
                        decision.support_face,
                        decision.boundary_side,
                        decision.accepted,
                        decision.reason
                    );
                }
            }
            repair_report.support_seam_relocation.merge(report);
        }
    }

    blend_insert::reclassify_corner_methods(result, corner_patches);

    // Final pass: clamp any out-of-range pcurve control points on
    // bounded-UV surfaces (Gregory / NURBS / sphere / torus). Runs
    // *last* so it catches pcurves attached by every prior pipeline
    // stage (attach_fillet_pcurves, sliver_absorb, cap_insert, …).
    let n_clamped = gregory_pcurve_clamp::clamp_bounded_pcurves(result);
    if std::env::var("POPKERN_CLAMP_TRACE").is_ok() {
        eprintln!("  [gregory_pcurve_clamp] clamped {n_clamped} pcurves");
    }
}

fn support_seam_relocation_enabled() -> bool {
    std::env::var("POPKERN_SUPPORT_SEAM_RELOCATION")
        .map(|value| value != "0")
        .unwrap_or(true)
}

fn topology_repair_enabled(apex_boundary_context: bool) -> bool {
    match std::env::var("POPKERN_TOPO_REPAIR") {
        Ok(value) => value != "0",
        Err(_) => apex_boundary_context,
    }
}
