//! Topology construction — apply fillet/chamfer geometry to the BRep. 3-phase single-pass (A: perp splits, B: multi-chord face splits, C: kef merge + surface assign) then surface-fix, sliver-absorb (§13.1(2)), and §8.3 cap insertion post-passes.

mod a2_dispatch;
mod apex_face_lift;
mod apex_vertex_coherence;
pub(crate) mod blend_insert;
pub(crate) mod blend_topo;
mod bridge_warnings;
pub(crate) mod cap_insert;
pub(crate) mod closed_fillet;
pub(crate) mod edge_reconcile;
mod gregory_pcurve_clamp;

/// Re-export of the late-pipeline clamp pass for use by `fillet_edges`
/// after `fix_stale_pcurves` has run.
pub fn clamp_bounded_pcurves_post_phase5(brep: &mut popkern_core::topology::BRep) -> usize {
    gregory_pcurve_clamp::clamp_bounded_pcurves(brep)
}
#[cfg(test)]
mod apex_termination_spike;
pub(crate) mod helpers;
mod kp13_bi5_promote;
mod kp_bi5_macro_apply;
mod kp_binonic_optimizer;
mod kp_biquintic_optimizer;
mod kp_lm_solver;
pub(crate) mod kp_polar_bi6_apply;
pub(crate) mod kp_polar_bi6_promote;
mod kp_polar_deep_trace;
mod kp_polar_lm_optimizer;
mod kp_twist_optimizer;
pub(crate) mod merge;
mod multi_chord;
pub(crate) mod pcurve_attach;
mod pcurve_refit;
mod phase_split;
pub(crate) mod poke_4n;
pub(crate) mod poke_hex;
mod promotion_ladder;
mod sliver_absorb;
mod sliver_absorb_classifiers;
mod sphere_fit_edge_wire;
pub(crate) mod surface_fix;
pub(crate) mod three_quad_apply;
mod three_quad_apply_edge_triple;
pub(crate) mod three_quad_promotions;
pub(crate) mod topology_repair;

use crate::analysis::EdgeAnalysis;
use crate::error::{BlendFailure, BlendWarning};
use crate::types::VertexMode;
use phase_split::{FaceTrim, TrimSplitResult};
use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};
use std::collections::HashMap;

/// Geometry computed for one fillet/chamfer edge.
#[derive(Debug)]
pub struct EdgeGeometry {
    /// The edge analysis.
    pub analysis: EdgeAnalysis,
    /// Fillet or chamfer surface.
    pub surface: SurfaceGeometry,
    /// Trim curve on face 1 (3D).
    pub trim_3d_1: NurbsCurve,
    /// Trim curve on face 2 (3D).
    pub trim_3d_2: NurbsCurve,
    /// Pcurve on face 1.
    pub pcurve_1: NurbsCurve,
    /// Pcurve on face 2.
    pub pcurve_2: NurbsCurve,
    /// Average fillet/chamfer radius (kept for legacy callers).
    pub effective_radius: f64,
    /// Radius at the edge's start vertex (constant ⇒ `effective_radius`; variable ⇒ `start`).
    pub radius_at_start: f64,
    /// Radius at the edge's end vertex.
    pub radius_at_end: f64,
    /// §18 telemetry: walker-path Newton iterations.
    pub iterations_used: Option<usize>,
    /// §18 telemetry: line-search rejections.
    pub step_rejections_used: Option<usize>,
    /// §8.3 start-endpoint cap (`None` when `cap_kind == Open`).
    pub start_cap: Option<crate::cap::CapGeometry>,
    /// §8.3 end-endpoint cap (`None` when `cap_kind == Open`).
    pub end_cap: Option<crate::cap::CapGeometry>,
    /// §7 continuity downgrade record — `Some((achieved, reason))` when the engine fell
    /// back; caller emits [`BlendWarning::ContinuityDowngraded`](crate::error::BlendWarning::ContinuityDowngraded).
    pub continuity_downgrade: Option<(crate::types::ContinuityLevel, String)>,
}

/// Result of applying blend geometry, including warnings.
pub struct ApplyBlendResult {
    /// The modified BRep.
    pub brep: BRep,
    /// Warnings generated during vertex mode handling.
    pub warnings: Vec<(EdgeHandle, BlendWarning)>,
    /// §5.3 per-corner patch reports — one entry per vertex where a blend patch was inserted.
    pub corner_patches: Vec<crate::types::CornerPatchInfo>,
}

/// Apply computed fillet/chamfer geometry to a BRep. 3-phase: A perp splits, B multi-chord face splits, C kef merge + surface assign.
#[allow(clippy::too_many_arguments)]
pub fn apply_blend_geometry(
    brep: &BRep,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
    vertex_mode: &VertexMode,
    use_rational_gregory: bool,
    use_deg5_gregory: bool,
    use_biquintic_full_gregory: bool,
    use_biseptic_gregory: bool,
    use_two_stripe_gregory: bool,
    use_kp_binonic_gregory: bool,
    use_kp_polar_gregory: bool,
    use_a2_apex_construction: bool,
) -> Result<ApplyBlendResult, Vec<(EdgeHandle, BlendFailure)>> {
    let mut result = brep.clone();
    let mut failures: Vec<(EdgeHandle, BlendFailure)> = Vec::new();
    let mut warnings: Vec<(EdgeHandle, BlendWarning)> = Vec::new();
    let mut corner_patches: Vec<crate::types::CornerPatchInfo> = Vec::new();
    // Phase 6c-D Stage 4.4: per-vertex (boundary, sub_faces) for
    // a post-`sliver_absorb` wire pass. Populated by
    // `insert_vertex_blend` when sphere-fit substitution fired.
    let mut late_wire_jobs: Vec<(
        crate::blend::BlendBoundary,
        Vec<popkern_core::topology::FaceHandle>,
    )> = Vec::new();

    // ── Phase A: ALL perpendicular edge splits ─ only modifies edges (inserts vertices), not faces.
    struct PerpSplitData {
        face_1_anchors: (VertexHandle, VertexHandle),
        face_2_anchors: (VertexHandle, VertexHandle),
    }
    let mut perp_data: Vec<Option<PerpSplitData>> = Vec::new();
    let mut closed_handled: std::collections::HashSet<usize> = std::collections::HashSet::new();
    // Round 5r+: snapshot before Phase A mutates result; anchor positions come from here.
    let phase_a_snapshot = result.clone();
    // Per-original-edge descendants after Phase A, for Phase C consumed-face lookup.
    let mut descendants: HashMap<EdgeHandle, Vec<EdgeHandle>> = HashMap::new();
    for eg in edge_geometries {
        descendants.insert(eg.analysis.edge, vec![eg.analysis.edge]);
    }

    for (idx, eg) in edge_geometries.iter().enumerate() {
        // Closed fillet edges (e.g. cylindrical hole rims) don't
        // have perpendicular boundary edges to split. Dispatch to
        // the closed-fillet surgery path and mark the edge as
        // already handled — skip Phase B/C for it.
        if closed_fillet::is_closed_fillet_edge(&result, eg.analysis.edge) {
            match closed_fillet::apply_closed_fillet(&mut result, eg) {
                Ok(()) => {
                    closed_handled.insert(idx);
                    perp_data.push(None);
                    continue;
                }
                Err(failure) => {
                    failures.push((eg.analysis.edge, failure));
                    perp_data.push(None);
                    continue;
                }
            }
        }

        match helpers::split_perpendicular_edges_by_radius(
            &phase_a_snapshot,
            &mut result,
            eg.analysis.face_1,
            eg.analysis.face_2,
            eg.analysis.edge,
            eg.radius_at_start,
            eg.radius_at_end,
            tolerance,
        ) {
            Ok(info) => {
                // Propagate every split that happened in this
                // iteration to the descendants map: any fillet
                // edge whose current descendants list contained
                // `old_edge` now tracks `edge_a` and `edge_b`
                // instead.
                for (old_edge, edge_a, edge_b) in &info.edge_splits {
                    for desc_list in descendants.values_mut() {
                        if let Some(pos) = desc_list.iter().position(|e| e == old_edge) {
                            desc_list.splice(pos..=pos, [*edge_a, *edge_b]);
                        }
                    }
                }
                perp_data.push(Some(PerpSplitData {
                    face_1_anchors: info.face_1_anchors,
                    face_2_anchors: info.face_2_anchors,
                }))
            }
            Err(e) => {
                let loc = crate::geometry_compute::edge_midpoint_by_handle(brep, eg.analysis.edge);
                failures.push((
                    eg.analysis.edge,
                    BlendFailure::InvalidTopology {
                        reason: e,
                        location: loc,
                    },
                ));
                perp_data.push(None);
            }
        }
    }

    // Collect trim data from Phase A results (BTreeMap for deterministic iteration).
    let mut face_trims: std::collections::BTreeMap<FaceHandle, Vec<FaceTrim>> =
        std::collections::BTreeMap::new();

    for (idx, eg) in edge_geometries.iter().enumerate() {
        let Some(ref pd) = perp_data[idx] else {
            continue;
        };

        face_trims
            .entry(eg.analysis.face_1)
            .or_default()
            .push(FaceTrim {
                edge_idx: idx,
                side: 1,
                anchors: pd.face_1_anchors,
            });
        face_trims
            .entry(eg.analysis.face_2)
            .or_default()
            .push(FaceTrim {
                edge_idx: idx,
                side: 2,
                anchors: pd.face_2_anchors,
            });
    }

    // Sort trims within each face by edge_idx for deterministic ordering.
    for trims in face_trims.values_mut() {
        trims.sort_by_key(|t| (t.edge_idx, t.side));
    }

    // ── Phase B: split faces by trim curves via `multi_chord`. ───
    let mut trim_results: Vec<TrimSplitResult> = Vec::new();
    for (original_face, trims) in &face_trims {
        let fillet_edges_by_trim: Vec<EdgeHandle> = trims
            .iter()
            .map(|t| edge_geometries[t.edge_idx].analysis.edge)
            .collect();
        let outcome = multi_chord::apply_multi_chord_for_face(
            &mut result,
            *original_face,
            trims,
            &fillet_edges_by_trim,
            &descendants,
            tolerance,
        );
        // Fold multi-chord's sub-edge splits (from inserting
        // crossing-point vertices on earlier trims' sub-chords)
        // into the global descendants map. Without this, Phase C
        // can't find the fillet edge on a face whose owning
        // sub-chord edge has been split.
        for (old_edge, edge_a, edge_b) in &outcome.extra_edge_splits {
            for desc_list in descendants.values_mut() {
                if let Some(pos) = desc_list.iter().position(|e| e == old_edge) {
                    desc_list.splice(pos..=pos, [*edge_a, *edge_b]);
                }
            }
        }
        trim_results.extend(outcome.trim_results);
        failures.extend(outcome.failures);
    }

    // ── Re-identify consumed faces from descendants (longest-sub-edge). ──
    for tr in &mut trim_results {
        let original_edge = edge_geometries[tr.edge_idx].analysis.edge;
        let fillet_descendants: &[EdgeHandle] = descendants
            .get(&original_edge)
            .map(|v| v.as_slice())
            .unwrap_or(std::slice::from_ref(&original_edge));

        // The trim's original face is face_1 or face_2 of the
        // fillet edge depending on `tr.side`.
        let eg = &edge_geometries[tr.edge_idx];
        let desired_parent = if tr.side == 1 {
            eg.analysis.face_1
        } else {
            eg.analysis.face_2
        };

        // Find the longest descendant and, of its two adjacent
        // faces, pick the one that's a descendant of
        // `desired_parent`.
        let mut longest: Option<(f64, FaceHandle, FaceHandle)> = None;
        for &d in fillet_descendants {
            let Some(edata) = result.edges.get(d) else {
                continue;
            };
            let he_h = edata.half_edge;
            let Some(he) = result.half_edges.get(he_h) else {
                continue;
            };
            let Some(twin) = result.half_edges.get(he.twin) else {
                continue;
            };
            let p0 = result.vertices[he.origin].position;
            let p1 = result.vertices[twin.origin].position;
            let len = p0.distance_to(&p1);
            let f_a = result.loops[he.loop_].face;
            let f_b = result.loops[twin.loop_].face;
            if longest.map(|(l, _, _)| len > l).unwrap_or(true) {
                longest = Some((len, f_a, f_b));
            }
        }
        if let Some((_, f_a, f_b)) = longest {
            // Prefer whichever of f_a, f_b descends from `desired_parent`.
            let pick = multi_chord::choose_side_face(&result, brep, desired_parent, f_a, f_b);
            tr.consumed = pick;
        }
    }

    // ── Phase C: merge consumed face pairs via kef ────────────────
    let mut edge_trims: std::collections::BTreeMap<usize, Vec<&TrimSplitResult>> =
        std::collections::BTreeMap::new();
    for tr in &trim_results {
        edge_trims.entry(tr.edge_idx).or_default().push(tr);
    }
    let mut merged_consumed: std::collections::HashSet<FaceHandle> =
        std::collections::HashSet::new();

    for (idx, trims) in &edge_trims {
        let eg = &edge_geometries[*idx];
        let fillet_edge = eg.analysis.edge;

        let side_1 = trims.iter().find(|t| t.side == 1);
        let side_2 = trims.iter().find(|t| t.side == 2);

        let (Some(s1), Some(s2)) = (side_1, side_2) else {
            // Missing one side — consumed face(s) cleaned up below.
            continue;
        };

        let fillet_descendants: &[EdgeHandle] = descendants
            .get(&fillet_edge)
            .map(|v| v.as_slice())
            .unwrap_or(std::slice::from_ref(&fillet_edge));
        match helpers::merge_consumed_faces_any(
            &mut result,
            s1.consumed,
            s2.consumed,
            fillet_descendants,
            &eg.surface,
        ) {
            Ok(fillet_face) => {
                merged_consumed.insert(s1.consumed);
                merged_consumed.insert(s2.consumed);

                // Use pcurves computed during Phase B (before kef).
                let params = pcurve_attach::FilletPcurveParams {
                    fillet_face,
                    fillet_surface: &eg.surface,
                    trim_edges_1: &s1.trim_edges,
                    trim_edges_2: &s2.trim_edges,
                    eg_pcurve_1: &s1.pcurve,
                    eg_pcurve_2: &s2.pcurve,
                    trim_3d_1: &s1.trim_3d,
                    trim_3d_2: &s2.trim_3d,
                };
                pcurve_attach::attach_fillet_pcurves(&mut result, &params);
            }
            Err(e) => {
                let loc = crate::geometry_compute::edge_midpoint_by_handle(brep, fillet_edge);
                failures.push((
                    fillet_edge,
                    BlendFailure::InvalidTopology {
                        reason: e,
                        location: loc,
                    },
                ));
            }
        }
    }

    // ── Phase C cleanup: undo orphaned consumed faces ─────────────
    // If a fillet edge only had one side split, its consumed face is
    // dangling. Undo the split by kef-ing the trim edge, but only if
    // the resulting merged face stays small (≤ 8 vertices). This
    // prevents creating oversized faces that no surface can fit.
    for tr in &trim_results {
        if merged_consumed.contains(&tr.consumed) {
            continue;
        }
        if !result.faces.contains_key(tr.consumed) {
            continue;
        }
        for &trim_edge in &tr.trim_edges {
            if !result.edges.contains_key(trim_edge) {
                continue;
            }
            // Check vertex count of both faces before merging.
            let he = result.edges[trim_edge].half_edge;
            let he_face = result.loops[result.half_edges[he].loop_].face;
            let twin = result.half_edges[he].twin;
            let twin_face = result.loops[result.half_edges[twin].loop_].face;
            if he_face == twin_face {
                continue;
            }
            let verts_a = result.face_outer_loop_positions(he_face).len();
            let verts_b = result.face_outer_loop_positions(twin_face).len();
            // After kef, the merged face has approximately verts_a + verts_b - 2 vertices.
            if verts_a + verts_b > 10 {
                continue; // Would create an oversized face.
            }
            let _ = result.kef(he);
        }
    }

    if !failures.is_empty() && failures.len() == edge_geometries.len() {
        return Err(failures);
    }

    // Vertex blend handling. `BTreeMap` keyed on VertexHandle for
    // deterministic insertion order across runs.
    if *vertex_mode == VertexMode::Smooth && edge_geometries.len() > 1 {
        let mut vertex_edges: std::collections::BTreeMap<VertexHandle, Vec<usize>> =
            std::collections::BTreeMap::new();
        for (idx, eg) in edge_geometries.iter().enumerate() {
            vertex_edges
                .entry(eg.analysis.vertex_start)
                .or_default()
                .push(idx);
            vertex_edges
                .entry(eg.analysis.vertex_end)
                .or_default()
                .push(idx);
        }
        // §5.3 item #5 — classify on the ORIGINAL `brep` (pre-Phase
        // A-through-C). Post-modification, original edge handles are
        // stale and `classify_vertex_blend` silently defaults to
        // `SameSideConvex`, masking `OppositeSide` saddles.
        type VC = crate::vertex_classifier::VertexBlendClassification;
        let precomputed: std::collections::BTreeMap<VertexHandle, VC> = vertex_edges
            .iter()
            .map(|(v, eis)| {
                let edges: Vec<_> = eis
                    .iter()
                    .filter_map(|&i| edge_geometries.get(i).map(|g| g.analysis.edge))
                    .collect();
                let c =
                    crate::vertex_classifier::classify_vertex_blend(brep, *v, &edges, tolerance)
                        .classification;
                (*v, c)
            })
            .collect();
        for (vertex, edge_indices) in &vertex_edges {
            if edge_indices.len() >= 2 {
                let pre = precomputed
                    .get(vertex)
                    .copied()
                    .unwrap_or(crate::vertex_classifier::VertexBlendClassification::SameSideConvex);

                // §6c-A2.3.b — when use_a2_apex_construction is set,
                // try the A2 planner + builder for this apex.  On
                // success, skip the legacy insert_vertex_blend path
                // entirely.  On any A2 error (plan or build), fall
                // back to legacy with a typed
                // BlendWarning::A2FellBackToLegacy recording the
                // reason.
                //
                // A2's failure modes:
                //   * plan_apex_resolution returns ApexPlanError
                //     (e.g. SphereCylinderNoIntersection on extreme
                //     radii, SingularCapCentre on nearly-coplanar
                //     faces).
                //   * construct_apex returns ApexBuildError
                //     (currently the dominant case during A2.2.b/c →
                //     A2.2.d/e/f transition: Sphere/Gregory cap
                //     construction is not yet implemented and the
                //     builder errors with PlanMismatch by design).
                //
                // The result of a successful A2 build is currently
                // not bridged into corner_patches — that's A2.2.d/e/f
                // work.  Until then, A2 success on non-Open plans is
                // unreachable (construct_apex errors on Sphere /
                // Gregory paths), so this dispatch hits the fallback
                // path on every cube/pyramid fixture today.
                let a2_outcome: Option<String> = if use_a2_apex_construction {
                    a2_dispatch::try_a2_apex_construction(
                        &mut result,
                        *vertex,
                        edge_indices,
                        edge_geometries,
                        &descendants,
                        tolerance,
                    )
                } else {
                    None // flag not set — legacy path runs unchanged
                };
                if let Some(reason) = a2_outcome {
                    let warn_edge = edge_geometries[edge_indices[0]].analysis.edge;
                    warnings.push((
                        warn_edge,
                        BlendWarning::A2FellBackToLegacy {
                            vertex: *vertex,
                            reason,
                        },
                    ));
                    // Fall through to legacy below.
                }

                match blend_insert::insert_vertex_blend(
                    &mut result,
                    *vertex,
                    edge_indices,
                    edge_geometries,
                    tolerance,
                    use_rational_gregory,
                    use_deg5_gregory,
                    use_biquintic_full_gregory,
                    use_biseptic_gregory,
                    use_two_stripe_gregory,
                    use_kp_binonic_gregory,
                    use_kp_polar_gregory,
                    Some(pre),
                ) {
                    Ok((infos, boundary_opt)) => {
                        for emit in [
                            bridge_warnings::emit_bridge_fallback_warnings,
                            bridge_warnings::emit_tangent_corner_warning,
                            bridge_warnings::emit_overlap_remediation_warning,
                            bridge_warnings::emit_corner_continuity_relaxed_warning,
                        ] {
                            emit(
                                &infos,
                                *vertex,
                                edge_indices,
                                edge_geometries,
                                &mut warnings,
                            );
                        }
                        // Phase 6c-D Stage 4.4: capture the blend
                        // boundary + per-vertex sub-face handles so a
                        // post-`sliver_absorb` wire pass can
                        // re-write substituted boundary curves to
                        // any cap-loop edges that were *added* by
                        // the absorb step. The wire originally runs
                        // mid-pipeline (inside three_quad_apply)
                        // when the cap face's outer loop is the
                        // raw poke output; sliver_absorb expands
                        // that loop with absorbed Plane fragments,
                        // leaving 8/12 cap-vs-cylinder edges on
                        // pyramid4_mixed unwired.
                        if let Some(b) = boundary_opt
                            && b.sphere_fit_substituted
                        {
                            let face_handles: Vec<_> = infos.iter().map(|cp| cp.face).collect();
                            late_wire_jobs.push((b, face_handles));
                        }
                        corner_patches.extend(infos);
                    }
                    Err(_) => {
                        let edge_idx = edge_indices[0];
                        if edge_idx < edge_geometries.len() {
                            warnings.push((
                                edge_geometries[edge_idx].analysis.edge,
                                BlendWarning::VertexModeFallback {
                                    requested: VertexMode::Smooth,
                                    actual: VertexMode::Mitred,
                                    vertex: Some(*vertex),
                                },
                            ));
                        }
                    }
                }
            }
        }
    }

    // Correct surface assignments — safety net.
    let fillet_surfaces: Vec<_> = edge_geometries
        .iter()
        .map(|eg| eg.surface.clone())
        .collect();
    surface_fix::correct_face_surfaces(&mut result, brep, &fillet_surfaces);

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
    edge_reconcile::reconcile_edge_curves(&mut result);

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
        pcurve_refit::refit_all_pcurves(&mut result);
    }

    // Topology repair for off-cylinder vertices in cylinder face
    // outer loops (apex-only filleting case where the cylinder
    // ribbon's base-side ends at an unfilleted original vertex via
    // a chord pair). Each repair attempt runs in a snapshot —
    // failed attempts roll back without affecting downstream
    // operations.
    //
    // Currently OPT-IN via `POPKERN_TOPO_REPAIR=1` because the
    // surgery splits the cylinder face into a (cylinder, sliver)
    // pair and the sliver face co-resides with the adjacent base
    // plane — valid topology for OCCT-strict slicers, but the
    // `unfillet` path doesn't yet know how to recognise sliver
    // faces and would fail to round-trip.
    if std::env::var("POPKERN_TOPO_REPAIR").is_ok() {
        let _ = topology_repair::repair_cylinder_face_off_surface_vertices(&mut result);
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
    let protected_caps: Vec<popkern_core::topology::FaceHandle> = late_wire_jobs
        .iter()
        .flat_map(|(_, faces)| faces.iter().copied())
        .collect();
    sliver_absorb::run_and_route_warnings_with_protected(
        &mut result,
        tolerance,
        edge_geometries,
        &protected_caps,
        &mut warnings,
    );

    cap_insert::insert_caps(&mut result, brep, edge_geometries, &mut warnings); // §8.3 caps

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
    for (boundary, sub_faces) in &late_wire_jobs {
        let wired = sphere_fit_edge_wire::wire_substituted_boundary_to_edges(
            &mut result,
            boundary,
            sub_faces,
        );
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
    // Default off until validated on more fixtures — opt-in via
    // `POPKERN_APEX_FACE_LIFT=1`.
    if std::env::var("POPKERN_APEX_FACE_LIFT").is_ok_and(|v| v == "1") {
        for (_boundary, sub_faces) in &late_wire_jobs {
            let lifted = apex_face_lift::lift_all_apex_plane_faces(&mut result, sub_faces);
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
    if std::env::var("POPKERN_APEX_PATH_A").is_ok_and(|v| v == "1") {
        for (_boundary, sub_faces) in &late_wire_jobs {
            let (n_moved, worst) =
                apex_vertex_coherence::relocate_apex_vertices(&mut result, sub_faces);
            // Phase 6c-X.A.5: rebuild BASE-EDGE curves (cyl-vs-Plane
            // apex generators going AWAY from the apex through V) as
            // straight lines from the new V→T position to v_other.
            // Restores vertex↔curve coherence on the BASE-EDGE side.
            let n_repaired = apex_vertex_coherence::repair_base_edges(&mut result, sub_faces);
            // Phase 6c-X.A.5b: realign cap-loop edges' curve direction
            // to match BRep convention (curve(t_start) =
            // primary_origin's vertex position).  Validate_watertight
            // demands this strict alignment; the wire's stored
            // direction was set by an old vertex-position heuristic
            // that no longer applies after A.4's relocations.
            let n_realigned =
                apex_vertex_coherence::realign_cap_loop_curve_directions(&mut result, sub_faces);
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
        let _n =
            crate::vertex_blend_topology::rebuild::rebuild_all_apex_caps(&mut result, &cap_faces);
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
                &mut result,
                &cap_faces,
                &per_fillet_caps,
            );
            if std::env::var("POPKERN_VAR_RADIUS_TRACE").is_ok_and(|v| v == "1") {
                eprintln!("[var-radius-reroute] caps={} result={r:?}", cap_faces.len());
            }
        }
    }

    blend_insert::reclassify_corner_methods(&result, &mut corner_patches);

    // Final pass: clamp any out-of-range pcurve control points on
    // bounded-UV surfaces (Gregory / NURBS / sphere / torus). Runs
    // *last* so it catches pcurves attached by every prior pipeline
    // stage (attach_fillet_pcurves, sliver_absorb, cap_insert, …).
    let n_clamped = gregory_pcurve_clamp::clamp_bounded_pcurves(&mut result);
    if std::env::var("POPKERN_CLAMP_TRACE").is_ok() {
        eprintln!("  [gregory_pcurve_clamp] clamped {n_clamped} pcurves");
    }

    Ok(ApplyBlendResult {
        brep: result,
        warnings,
        corner_patches,
    })
}
