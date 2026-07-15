//! Late post-processing passes for applied blend topology.

use super::setback_vertex_topology::fan_quality_stage_probe as fanq;
use super::*;

#[path = "determinism_trace.rs"]
pub(crate) mod determinism_trace_mod;
pub(crate) use determinism_trace_mod::{
    determinism_face_trace, determinism_trace, determinism_trace_enabled, determinism_value_hash,
};

pub(super) struct PostBlendPipelineInputs<'a> {
    /// A setback vertex blend (KP20 or PnS cap family) landed in
    /// this result — the setback surgery owns its apex topology.
    pub setback_blend_present: bool,
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
    determinism_trace("pre_surface_fix", result);
    surface_fix::correct_face_surfaces(result, brep, &fillet_surfaces);
    determinism_trace("post_surface_fix", result);
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
    determinism_trace("post_edge_reconcile", result);
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
        determinism_trace("post_pcurve_refit", result);
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
    // The cylinder off-surface-vertex repair belongs to the LEGACY
    // apex machinery. Setback-built blends own their apex topology
    // (PK-097: on the G1 analytic family this repair consumed 3
    // lateral faces and collapsed the solid — corner patches with
    // CornerMethod::Gregory tripped the context even though the
    // setback surgery, not the legacy builder, produced them).
    // Recorded limitation: a mixed model (setback at one vertex,
    // legacy Gregory at another) skips the repair for both.
    let apex_boundary_context = !inputs.setback_blend_present
        && (!support_seam_jobs.is_empty()
            || !late_wire_jobs.is_empty()
            || !direct_apex_built_faces.is_empty()
            || corner_patches
                .iter()
                .any(|patch| matches!(patch.method, crate::types::CornerMethod::Gregory)));
    if topology_repair_enabled(apex_boundary_context) {
        let repaired = topology_repair::repair_cylinder_face_off_surface_vertices(result);
        determinism_trace("post_topology_repair", result);
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
    // don't re-author them. Exact-boundary caps also protect their
    // semantic bridge support faces: those plane fragments are part of
    // the accepted cap boundary contract, so absorbing them after
    // acceptance turns FaceEdge ownership into a stale cap-vs-cylinder
    // adjacency and breaks strict edge/surface validation.
    protected_caps.extend(direct_apex_built_faces.iter().copied());
    protected_caps.extend(exact_boundary_protected_support_faces(repair_report));
    determinism_trace("pre_sliver_absorb", result);
    // Sliver absorption unifies LEGACY cap-loop topology (sphere /
    // Gregory cap twins). Setback-built blends own their apex
    // topology and run their own band absorption inside the surgery
    // — on the G1 analytic family this pass consumed 3 lateral
    // faces and collapsed the solid (PK-097; the G2 family was
    // untouched only because the pass happened to no-op there).
    let pp_census = |tag: &str, b: &popkern_core::topology::BRep| {
        if std::env::var_os("PK_CARVE_TRACE").is_some() {
            let mut cyls = Vec::new();
            for (fh, f) in b.faces.iter() {
                if let Some(popkern_core::geometry::SurfaceGeometry::Cylinder { axis, .. }) =
                    f.surface.as_ref().map(|sof| &sof.geometry)
                {
                    let dir = if axis.z.abs() > 0.9 {
                        "V"
                    } else if axis.z.abs() < 0.1 {
                        "H"
                    } else {
                        "?"
                    };
                    cyls.push(format!("{fh:?}{dir}"));
                }
            }
            eprintln!("  [pp] {tag}: faces={} cyls={cyls:?}", b.faces.len());
        }
    };
    pp_census("pre-sliver-absorb", result);
    if !inputs.setback_blend_present {
        sliver_absorb::run_and_route_warnings_with_protected(
            result,
            tolerance,
            edge_geometries,
            &protected_caps,
            warnings,
        );
    }

    pp_census("post-sliver-absorb", result);

    determinism_trace("post_sliver_absorb", result);
    fanq("post_sliver_absorb", result, corner_patches);
    cap_insert::insert_caps(result, brep, edge_geometries, warnings); // §8.3 caps
    determinism_trace("post_cap_insert", result);
    if std::env::var_os("PK_CARVE_TRACE").is_some()
        && let Err(e) = result.validate_loop_closure()
    {
        eprintln!("  [closure@post_cap_insert] {e:?}");
    }

    // Route B corner reconstruction, LAST — after sliver absorption
    // AND §8.3 cap insertion (cap_insert re-anchors corner topology
    // from precomputed edge geometries; contracting before it left
    // resurrected anchor positions and 0.1131 open loops). The
    // absorbed corner fragments contribute the on-sphere contact
    // vertices to the cap loops (pre-absorb caps are pure anchor
    // triangles). Contract anchors into contacts, then rebuild the
    // cap boundary chords as exact sphere-centered arcs.
    let collapsed =
        super::corner_arc_boundaries::collapse_corner_anchor_chords(result, corner_patches);
    if collapsed > 0 && std::env::var_os("PK_CARVE_TRACE").is_some() {
        eprintln!("  [contact-collapse] contracted {collapsed} anchor chords");
    }
    let arcs = super::corner_arc_boundaries::arcify_corner_cap_boundaries(result, corner_patches);
    if arcs > 0 && std::env::var_os("PK_CARVE_TRACE").is_some() {
        eprintln!("  [arcify] rewrote {arcs} corner-cap boundary chords");
    }
    determinism_trace("post_corner_arcify", result);
    // Strip END reconstruction (the corner passes' analogue at
    // planar end faces): contract off-cylinder end corners, author
    // exact cylinder ∩ plane runout arcs, remove leftover spur
    // slits. Gated on blend context — corner patches or a
    // multi-fillet result. SINGLE plain fillets keep their legacy
    // end shape for two measured contracts: the PK-097 in-place
    // modify machinery cannot reproduce topology edits
    // (open_in_place_matches_remove_recreate_geometry failed 10 vs
    // 12 vertices under unconditional gating), and cap-kind
    // invariance on closed solids (a cap REQUEST is not a cap
    // INSERTION — gating on requested caps broke
    // cap_insertion_is_a_no_op_on_closed_cube_topology, 17 vs 15).
    // Teaching in-place modify to re-author strip ends is that
    // machinery's follow-up; see strip_end_arcs module docs.
    if !corner_patches.is_empty() || edge_geometries.len() > 1 {
        let (end_contracted, end_arcs, end_spurs) =
            super::strip_end_arcs::reconstruct_strip_end_boundaries(result, &fillet_surfaces);
        if (end_contracted + end_arcs + end_spurs) > 0
            && std::env::var_os("PK_CARVE_TRACE").is_some()
        {
            eprintln!(
                "  [strip-end] {end_contracted} corners contracted, {end_arcs} runout arcs, {end_spurs} spurs removed"
            );
        }
        determinism_trace("post_strip_end_arcs", result);
    }
    if std::env::var_os("PK_CARVE_TRACE").is_some()
        && let Err(e) = result.validate_loop_closure()
    {
        eprintln!("  [closure@post_corner_arcify] {e:?}");
        if let popkern_core::error::ValidationError::OpenLoop { loop_, .. } = e
            && let Some(l) = result.loops.get(loop_)
        {
            use popkern_core::geometry::curve::CurveEval;
            let f = l.face;
            let kind = result.faces[f]
                .surface
                .as_ref()
                .map(|sof| match &sof.geometry {
                    popkern_core::geometry::SurfaceGeometry::Plane { .. } => "Plane",
                    popkern_core::geometry::SurfaceGeometry::Cylinder { .. } => "Cyl",
                    popkern_core::geometry::SurfaceGeometry::Sphere { .. } => "Sphere",
                    _ => "?",
                })
                .unwrap_or("None");
            eprintln!("    face {f:?} ({kind})");
            let start = l.half_edge;
            let mut he = start;
            for _ in 0..result.half_edges.len() + 1 {
                let d = &result.half_edges[he];
                if let Some(c) = result.edges[d.edge].curve.as_ref() {
                    let ceval: &dyn CurveEval = &c.geometry;
                    let p0 = ceval.point_at(c.t_start);
                    let p1 = ceval.point_at(c.t_end);
                    let k = match &c.geometry {
                        popkern_core::geometry::CurveGeometry::Line { .. } => "Line",
                        popkern_core::geometry::CurveGeometry::Arc { .. } => "Arc",
                        popkern_core::geometry::CurveGeometry::Nurbs(_) => "Nurbs",
                        _ => "?",
                    };
                    eprintln!(
                        "    e {:?} {k} ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                        d.edge, p0.x, p0.y, p0.z, p1.x, p1.y, p1.z
                    );
                }
                he = d.next;
                if he == start {
                    break;
                }
            }
        }
    }

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
        determinism_trace("post_wire_jobs", result);
        if std::env::var_os("POPKERN_KP_TRACE").is_some() {
            eprintln!(
                "[late-wire] re-ran on {} sub_faces, wrote {wired} edges",
                sub_faces.len()
            );
        }
    }

    // Phase 6c-D Stage 4.7: lift legacy apex-region pyramid Plane
    // faces to curved NURBS surfaces where late-wire jobs own the
    // substituted bridge. Exact-boundary sectors keep their original
    // production support faces; measured pcurve repair handles their
    // export-readiness without broad support-face lifting.
    //
    // On by default since #207 cleanup.  Set
    // `POPKERN_APEX_FACE_LIFT=0` to disable for debugging.
    if std::env::var("POPKERN_APEX_FACE_LIFT")
        .map(|v| v != "0")
        .unwrap_or(true)
    {
        for (_boundary, sub_faces) in late_wire_jobs {
            let lifted = apex_face_lift::lift_all_apex_plane_faces(result, sub_faces);
            determinism_trace("post_apex_lift", result);
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
            determinism_trace("post_relocate", result);
            // Phase 6c-X.A.5: rebuild BASE-EDGE curves (cyl-vs-Plane
            // apex generators going AWAY from the apex through V) as
            // straight lines from the new V→T position to v_other.
            // Restores vertex↔curve coherence on the BASE-EDGE side.
            let n_repaired = apex_vertex_coherence::repair_base_edges(result, sub_faces);
            determinism_trace("post_repair_base", result);
            // Phase 6c-X.A.5b: realign cap-loop edges' curve direction
            // to match BRep convention (curve(t_start) =
            // primary_origin's vertex position).  Validate_watertight
            // demands this strict alignment; the wire's stored
            // direction was set by an old vertex-position heuristic
            // that no longer applies after A.4's relocations.
            let n_realigned =
                apex_vertex_coherence::realign_cap_loop_curve_directions(result, sub_faces);
            determinism_trace("post_realign", result);
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
        determinism_trace("post_rebuild_caps", result);
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
            determinism_trace("post_var_radius", result);
            if std::env::var("POPKERN_VAR_RADIUS_TRACE").is_ok_and(|v| v == "1") {
                eprintln!("[var-radius-reroute] caps={} result={r:?}", cap_faces.len());
            }
        }
    }

    let exact_boundary_owns_apex_seams = repair_report.apex_routes.iter().any(|route| {
        route.exact_boundary_report.candidate_status
            == crate::types::ExactBoundaryCandidateOutcome::Accepted
            && !route.exact_boundary_report.committed_faces.is_empty()
    });
    if support_seam_relocation_enabled() && !exact_boundary_owns_apex_seams {
        determinism_trace("pre_seams", result);
        for (boundary, sub_faces) in support_seam_jobs {
            let report = three_quad_apply::author_model_level_support_face_boundary_seams(
                result,
                sub_faces,
                boundary,
                &[],
            );
            determinism_trace("post_seams", result);
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

    determinism_trace("pre_reclassify", result);
    blend_insert::reclassify_corner_methods(result, corner_patches);

    // Final pass: clamp any out-of-range pcurve control points on
    // bounded-UV surfaces (Gregory / NURBS / sphere / torus). Runs
    // *last* so it catches pcurves attached by every prior pipeline
    // stage (attach_fillet_pcurves, sliver_absorb, cap_insert, …).
    determinism_trace("pre_pcurve_clamp", result);
    if std::env::var_os("PK_CARVE_TRACE").is_some()
        && let Err(e) = result.validate_loop_closure()
    {
        eprintln!("  [closure@pre_pcurve_clamp] {e:?}");
    }
    let n_clamped = gregory_pcurve_clamp::clamp_bounded_pcurves(result);
    if std::env::var("POPKERN_CLAMP_TRACE").is_ok() {
        eprintln!("  [gregory_pcurve_clamp] clamped {n_clamped} pcurves");
    }

    // Restore semantic curve storage after late repair passes that may fit
    // straight trim spans through a NURBS path. Preserve authored trim triples:
    // exact-boundary edges may be nearly straight while still requiring their
    // fitted 3-D curve to remain coupled to support-side pcurves.
    result.normalize_semantic_line_edge_curves_preserving_pcurves(tolerance.linear().max(1.0e-8));

    determinism_trace("pre_export_repair", result);
    if std::env::var_os("PK_CARVE_TRACE").is_some()
        && let Err(e) = result.validate_loop_closure()
    {
        eprintln!("  [closure@pre_export_repair] {e:?}");
    }
    exact_boundary_export_repair::repair_exact_boundary_export_blockers(result, repair_report);
    fanq("post_export_repair", result, corner_patches);
    repair_report.edge_curve_endpoint_coherence_repairs += result.separate_nonmanifold_edge_uses();
    repair_report.edge_curve_endpoint_coherence_repairs +=
        result.repair_edge_curve_endpoint_coherence(tolerance.linear().max(1.0e-8));

    determinism_trace("pre_completion", result);
    if std::env::var_os("PK_CARVE_TRACE").is_some()
        && let Err(e) = result.validate_loop_closure()
    {
        eprintln!("  [closure@pre_completion] {e:?}");
    }
    corner_patch_support_completion::complete_corner_patch_support_pcurves(
        result,
        corner_patches,
        &mut repair_report.apex_diagnostics,
    );
    determinism_trace("post_completion", result);
    if std::env::var_os("PK_CARVE_TRACE").is_some()
        && let Err(e) = result.validate_loop_closure()
    {
        eprintln!("  [closure@post_completion] {e:?}");
    }
    determinism_face_trace("post_completion", result);
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

fn exact_boundary_protected_support_faces(
    repair_report: &crate::types::BlendRepairReport,
) -> Vec<popkern_core::topology::FaceHandle> {
    repair_report
        .apex_routes
        .iter()
        .filter(|route| {
            route.exact_boundary_report.candidate_status
                == crate::types::ExactBoundaryCandidateOutcome::Accepted
        })
        .flat_map(|route| {
            route
                .exact_boundary_report
                .sectors
                .iter()
                .flat_map(|sector| [sector.fillet_support_face, sector.bridge_support_face])
        })
        .collect()
}
