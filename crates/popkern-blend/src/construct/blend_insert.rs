//! Vertex blend face insertion into BRep topology.

use popkern_core::geometry::surface::SurfaceOnFace;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::blend;
use crate::blend::BlendBoundary;
use crate::construct::EdgeGeometry;
use crate::construct::blend_topo::{
    build_fillet_surface_list, collect_edges_through_vertex, collect_faces_around_vertex,
    repair_vertex_outgoing,
};
use crate::construct::poke_hex::face_has_2n_outer_half_edges;
use crate::construct::three_quad_apply::apply_three_quad_decomposition_to_blend;

#[path = "blend_insert_diagnostics.rs"]
mod blend_insert_diagnostics;
#[path = "blend_insert_pns_acceptance.rs"]
mod blend_insert_pns_acceptance;
#[path = "blend_insert_ring_rung.rs"]
mod blend_insert_ring_rung;
#[path = "blend_insert_setback_rung.rs"]
mod blend_insert_setback_rung;
#[path = "blend_insert_support_cell_acceptance.rs"]
mod blend_insert_support_cell_acceptance;
#[path = "blend_insert_surface.rs"]
mod blend_insert_surface;
#[path = "blend_insert_three_quad_acceptance.rs"]
mod blend_insert_three_quad_acceptance;
#[path = "blend_insert_topology.rs"]
mod blend_insert_topology;
use crate::construct::direct_apex_promotion::cylinder_axes_mutually_orthogonal;
use crate::types::{
    ApexConstructionDiagnostic, ContinuityLevel, CornerMethod, CornerPatchInfo, VertexMode,
};
use crate::vertex_classifier::VertexBlendClassification;
use blend_insert_diagnostics::{
    pns_dispatch_declined_diagnostic, rational_gregory_dispatch_declined_diagnostic,
    rational_gregory_dispatch_precondition_reason,
};
use blend_insert_surface::{assign_blend_surface, corner_radii_uniform, sphere_owns_corner};
pub(crate) use blend_insert_surface::{reclassify_corner_methods, setback_min_n_for};
use blend_insert_three_quad_acceptance::record_fallback_cap_support_continuity;
pub(crate) use blend_insert_topology::{kev_v_bridge, remove_orphaned_vertex};
use blend_insert_topology::{normalize_underresolved_face_loop_to_n, split_corner_triangles};

type InsertVertexBlendOutput = (
    Vec<CornerPatchInfo>,
    Option<BlendBoundary>,
    Vec<ApexConstructionDiagnostic>,
);

/// Insert a smooth G1 vertex blend face at a multi-fillet vertex.
///
/// For each face touching the vertex, splits off the corner triangle with `mef()`,
/// merges the triangles into one N-gon by removing edges through the vertex with
/// `kef()`, and assigns the blend surface geometry to the merged face.
///
/// # Errors
///
/// Returns a descriptive string if any topology operation fails.
#[allow(clippy::too_many_arguments)]
pub(crate) fn insert_vertex_blend(
    brep: &mut BRep,
    vertex: VertexHandle,
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
    requested_continuity: ContinuityLevel,
    use_rational_gregory: bool,
    use_deg5_gregory: bool,
    use_biquintic_full_gregory: bool,
    use_biseptic_gregory: bool,
    use_two_stripe_gregory: bool,
    use_kp_binonic_gregory: bool,
    use_kp_polar_gregory: bool,
    precomputed_classification: Option<VertexBlendClassification>,
    cap_presolve: Option<&crate::construct::setback_cap_presolve::SetbackCapPresolveEntry>,
    corner_cap: crate::types::CornerCapPreference,
    explicit_subpatch_family: bool,
) -> Result<InsertVertexBlendOutput, String> {
    // §5.3 item #5 — callers pre-compute classification against the pre-modification BRep; post-Phase-C original edge handles are stale. Fallback for tests that call us on an untouched BRep.
    let classification = precomputed_classification.unwrap_or_else(|| {
        let filleted_edges: Vec<_> = edge_indices
            .iter()
            .filter_map(|&i| edge_geometries.get(i).map(|g| g.analysis.edge))
            .collect();
        crate::vertex_classifier::classify_vertex_blend(brep, vertex, &filleted_edges, tolerance)
            .classification
    });

    // PK-089: N >= 5 corners are all-or-nothing — an accepted PnS or
    // N-patch decomposition, or a full rollback to the mitred corner.
    // A lone undecomposed N-gon Gregory misrepresents the tiling and
    // shipped unmeasured; snapshot so partial routes can restore.
    let mut n5_entry_snapshot = (edge_indices.len() >= 5).then(|| brep.clone());

    // Repair the vertex's outgoing half-edge if it has become stale
    // due to prior topology operations (edge splitting, face splitting).
    repair_vertex_outgoing(brep, vertex)?;

    if edge_geometries.is_empty() {
        return Err(format!("vertex {vertex:?} has no fillet geometries"));
    }
    let face_info = collect_faces_around_vertex(brep, vertex)?;
    if face_info.is_empty() {
        return Err(format!("no faces found around vertex {vertex:?}"));
    }

    // §5.3 item #5 — `OppositeSide` flows through Gregory with a
    // saddle-specific central point. `Mixed` stays on polar-NURBS.
    let fillet_surfaces_pre = build_fillet_surface_list(edge_indices, edge_geometries);
    // Setback rung (setback_vertex_topology docs; PK-225/226; order PK-603; gate PK-597);
    // declines are transactional. Seams: SETBACK_DISABLE / with_setback_forced_to_decline.
    let radii_uniform_pre = corner_radii_uniform(edge_indices, edge_geometries, tolerance);
    let ball_owned = sphere_owns_corner(radii_uniform_pre, &fillet_surfaces_pre, classification);
    let setback_min_n = setback_min_n_for(
        edge_indices,
        edge_geometries,
        &fillet_surfaces_pre,
        classification,
        tolerance,
        explicit_subpatch_family,
    );
    let setback_attempted =
        crate::construct::setback_vertex_topology::setback_topology_commit_enabled()
            && edge_indices.len() >= setback_min_n;
    let (setback_served, setback_diagnostics) = blend_insert_setback_rung::try_setback_rung(
        brep,
        vertex,
        edge_indices,
        &fillet_surfaces_pre,
        tolerance,
        requested_continuity,
        classification,
        setback_min_n,
        cap_presolve,
        corner_cap,
        blend_insert_surface::corner_min_radius(edge_indices, edge_geometries),
    );
    if let Some(outcome) = setback_served {
        return outcome;
    }
    // Reaching here with `setback_attempted` set means the surgery
    // ran and declined transactionally (a served setback returned
    // above) — the fact the boundary-ring gate consumes.
    let setback_declined = setback_attempted;

    let local_indices_pre: Vec<usize> = (0..fillet_surfaces_pre.len()).collect();
    let n_sides = edge_indices.len();

    // Cap-kind dispatch: when the apex's incident fillets all share
    // a uniform radius within tolerance, suppress sphere-fit
    // substitution for this corner. The post-loop apex topology
    // rebuild then collapses the rolling-ball 6-vertex cap loop to a
    // 3-vertex cap-tangent-arc loop directly. For non-uniform
    // corners, substitution and rebuild stay active.
    let auto_disable_substitution = radii_uniform_pre;
    let boundary_options = crate::blend::BlendBoundaryOptions {
        disable_sphere_fit_substitution: if auto_disable_substitution {
            Some(true)
        } else {
            None
        },
        // PK-234: saddle bridges run through consumed support bands
        // without the rail-crossing station pushback (convex ring
        // corners get their pushed boundary inside the ring rung).
        ring_station_pushback: matches!(classification, VertexBlendClassification::OppositeSide),
        ring_sense_from_solid: false,
    };
    let rational_gregory_compatible = matches!(
        classification,
        VertexBlendClassification::SameSideConvex
            | VertexBlendClassification::SameSideConcave
            | VertexBlendClassification::OppositeSide,
    );
    // Compute the n-gon BlendBoundary whenever the rational-Gregory
    // dispatch or the PnS path could consume it — both share the
    // same boundary structure downstream.
    let mut n_gon_boundary: Option<BlendBoundary> = if n_sides >= 3 && rational_gregory_compatible {
        match blend::compute_blend_boundary_with_options(
            vertex,
            &local_indices_pre,
            &fillet_surfaces_pre,
            brep,
            tolerance,
            boundary_options,
        ) {
            Ok(b) if b.segments.len() == 2 * n_sides => Some(b),
            _ => None,
        }
    } else {
        None
    };

    // Boundary-ring rung (PK-234 saddles, PK-239 convex extension) —
    // gate semantics, the convex pushed-boundary contract, and the
    // transactional carve live in `blend_insert_ring_rung`. On a
    // served ring `n_gon_boundary` is substituted with the carved
    // (pushed) boundary; on a decline dispatch below is bit-identical
    // to the mitred path.
    let sphere_claim_false =
        ball_owned && !cylinder_axes_mutually_orthogonal(edge_indices, edge_geometries);
    let (boundary_ring_eligible, ring_face, ring_decline) =
        blend_insert_ring_rung::boundary_ring_rung(
            brep,
            vertex,
            classification,
            use_rational_gregory,
            setback_declined,
            sphere_claim_false,
            &mut n_gon_boundary,
            boundary_options,
            &local_indices_pre,
            &fillet_surfaces_pre,
            n_sides,
            tolerance,
        );

    let blend_face = match ring_face {
        Some(face) => face,
        None => {
            // For each face, cut a chord between the face's two
            // fillet-trim-endpoint anchors. The V-containing sub-face is
            // the "corner region" (may be a triangle, quad, or pentagon
            // depending on how many intermediate cube-sub-split vertices
            // sit between V and the anchors).
            let _corner_faces =
                split_corner_triangles(brep, vertex, &face_info, edge_geometries, tolerance)?;

            // Kef every V-edge separating two distinct faces;
            // last V-edge may become a same-face "bridge" handled by
            // the kev pass's `kev` after surface assignment.
            kef_cross_face_v_edges(brep, vertex)?
        }
    };

    // Compute the blend boundary + assign the surface. The
    // ring carve consumed the vertex (its region merge removed it),
    // so that path reuses the pre-carve boundary it was built from —
    // recomputation is impossible AND unnecessary (the ring IS the
    // boundary by construction).
    let fillet_surfaces = build_fillet_surface_list(edge_indices, edge_geometries);
    let local_indices: Vec<usize> = (0..fillet_surfaces.len()).collect();
    let blend_boundary = if ring_face.is_some() {
        n_gon_boundary
            .clone()
            .ok_or_else(|| "saddle ring path requires the n-gon boundary".to_string())?
    } else {
        blend::compute_blend_boundary_with_options(
            vertex,
            &local_indices,
            &fillet_surfaces,
            brep,
            tolerance,
            boundary_options,
        )
        .map_err(|e| format!("blend boundary computation failed: {e}"))?
    };
    let (method, achieved_continuity) = assign_blend_surface(
        brep,
        blend_face,
        &blend_boundary,
        &fillet_surfaces,
        classification,
        tolerance,
    )?;

    // Kev the remaining V-bridge if any. `kev` kills V
    // and the bridge edge simultaneously, un-pinching the merged
    // loop into a manifold polygon. This is R8 C6's crucial
    // topology fix — the old code left the bridge as an internal
    // dangling edge because its `kef` guard rejected same-loop
    // half-edge pairs. The saddle ring carve already consumed the
    // vertex in its region merge, so that path skips both steps.
    if ring_face.is_none() {
        kev_v_bridge(brep, vertex)?;
        remove_orphaned_vertex(brep, vertex);
    }

    // Sphere-fit substitution (N-path arm): when sphere-fit substitution
    // succeeded for `n_gon_boundary` and we're about to take the
    // simpler N-edge bridge construction (faces_2n=false — the
    // uniform-radius collapse case), write the substituted
    // boundary curves into blend_face's N outer-loop edges.
    // Without this, cube_uniform / pyramid_uniform never see the
    // substituted curves: the rational-Gregory branch below where
    // the previously-shipped 6c-D wiring lives is unreachable for
    // them.
    if let Some(boundary) = n_gon_boundary.as_ref()
        && boundary.sphere_fit_substituted
        && !face_has_2n_outer_half_edges(brep, blend_face, n_sides)
    {
        crate::construct::sphere_fit_edge_wire::wire_substituted_boundary_to_edges(
            brep,
            boundary,
            &[blend_face],
        );
    }

    // §5.3 item #6 — propagate pre-merge bridge-fallback names.
    let dispatch_boundary = n_gon_boundary.as_ref().unwrap_or(&blend_boundary);
    let bridge_fallbacks = dispatch_boundary.bridge_fallbacks.clone();
    let base_patch = CornerPatchInfo {
        vertex,
        n_sides: edge_indices.len(),
        method,
        achieved_continuity,
        classification,
        face: blend_face,
        bridge_fallbacks: bridge_fallbacks.clone(),
        overlap_remediated: false,
        relaxed_continuity_coplanarity: None,
        cap_sphere: None,
        sub_patch_provenance: None,
        setback_cap: false,
        exact_boundary_cap: false,
    };
    let mut diagnostics = setback_diagnostics;
    // The ring rung records its outcome either way — a decline means
    // the mitred fallback is serving (PK-236 deviation classes).
    diagnostics.extend(blend_insert_ring_rung::boundary_ring_diagnostic(
        boundary_ring_eligible,
        n_sides,
        ring_face,
        ring_decline.as_ref(),
    ));

    // Record WHY the PnS route was not entered. Until this existed, a
    // corner that failed the boundary preconditions was
    // indistinguishable from one where PnS was never enabled — both
    // produced zero diagnostics — which is precisely what hid the
    // n=3 / n=4 dispatch gap.
    if n_sides >= crate::blend::pns_promote::PNS_MIN_VALENCE {
        let segments = n_gon_boundary.as_ref().map(|b| b.segments.len());
        let precondition_failure = match segments {
            None => Some("no_n_gon_boundary".to_string()),
            Some(len) if len != 2 * n_sides => {
                Some(format!("segment_count={len}, expected={}", 2 * n_sides))
            }
            Some(_) => None,
        };
        if let Some(detail) = precondition_failure {
            diagnostics.push(pns_dispatch_declined_diagnostic(
                n_sides,
                "precondition_failed",
                Some(detail),
            ));
        }
    }

    if n_sides >= crate::blend::pns_promote::PNS_MIN_VALENCE
        && let Some(boundary) = n_gon_boundary.as_ref()
        && boundary.segments.len() == 2 * n_sides
        && let Some(patches) = blend_insert_pns_acceptance::try_pns_dispatch_with_acceptance(
            blend_insert_pns_acceptance::PnsDispatchAcceptanceInput {
                brep,
                blend_face,
                boundary,
                fillet_surfaces: &fillet_surfaces_pre,
                vertex,
                classification,
                n_sides,
                diagnostics: &mut diagnostics,
            },
        )
    {
        return Ok((patches, n_gon_boundary, Vec::new()));
    }

    if use_rational_gregory
        && use_kp_polar_gregory
        && rational_gregory_compatible
        && n_sides >= 5
        && dispatch_boundary.segments.len() == 2 * n_sides
        && !crate::construct::poke_hex::face_has_2n_outer_half_edges_with(
            brep,
            blend_face,
            n_sides,
            dispatch_boundary.sphere_fit_substituted,
        )
        && let Err(reason) = normalize_underresolved_face_loop_to_n(brep, blend_face, n_sides)
    {
        diagnostics.push(rational_gregory_dispatch_declined_diagnostic(
            n_sides,
            "normalization_failed",
            Some(reason),
        ));
    }

    // PK-222: N = 3 corners whose the cross-face kef pass collapsed the arc+bridge
    // pairs carry the physical N-half-edge loop — demanding the 2N
    // form starved the rational dispatch (polar-NURBS heuristic at
    // g0 ~ 0.2). The N-form arm is new dispatch reach with no fixture
    // history, so its commits must pass the same cap/support
    // acceptance the N >= 5 routes answer to; 2N keeps its historical
    // unconditional acceptance.
    let accept_n3_he = n_sides == 3 && !ball_owned;
    let n3_n_form_route = accept_n3_he && !face_has_2n_outer_half_edges(brep, blend_face, n_sides);
    let rational_dispatch_preconditions_met = dispatch_boundary.segments.len() == 2 * n_sides
        && crate::construct::poke_hex::face_has_2n_outer_half_edges_with_n3(
            brep,
            blend_face,
            n_sides,
            dispatch_boundary.sphere_fit_substituted,
            accept_n3_he,
        );
    // Record WHY dispatch didn't apply for ANY rational route (was
    // polar-gated: pyramid4_mixed binonic bailed with zero diagnostics
    // while its fallback cap shipped at g0 = 1.34).
    if use_rational_gregory && rational_gregory_compatible && !rational_dispatch_preconditions_met {
        diagnostics.push(rational_gregory_dispatch_declined_diagnostic(
            n_sides,
            "precondition_failed",
            Some(rational_gregory_dispatch_precondition_reason(
                brep,
                blend_face,
                n_sides,
                dispatch_boundary,
                accept_n3_he,
            )),
        ));
    }

    if use_rational_gregory && rational_gregory_compatible && rational_dispatch_preconditions_met {
        use crate::construct::three_quad_apply::ThreeQuadFailure;
        let three_quad_snapshot = brep.clone();
        match apply_three_quad_decomposition_to_blend(
            brep,
            blend_face,
            dispatch_boundary,
            &fillet_surfaces_pre,
            vertex,
            edge_indices.len(),
            classification,
            use_deg5_gregory,
            use_biquintic_full_gregory,
            use_biseptic_gregory,
            use_two_stripe_gregory,
            use_kp_binonic_gregory,
            use_kp_polar_gregory,
        ) {
            Ok(result) => {
                match blend_insert_three_quad_acceptance::accept_three_quad_result(
                    brep,
                    result,
                    n_sides,
                    n3_n_form_route,
                ) {
                    blend_insert_three_quad_acceptance::ThreeQuadAcceptance::Accepted(result) => {
                        diagnostics.extend(result.diagnostics);
                        return Ok((result.patches, Some(dispatch_boundary.clone()), diagnostics));
                    }
                    blend_insert_three_quad_acceptance::ThreeQuadAcceptance::Declined(
                        decline_diagnostics,
                    ) => {
                        diagnostics.extend(decline_diagnostics);
                        *brep = three_quad_snapshot;
                    }
                }
            }
            Err(ThreeQuadFailure::Overlap) => {
                *brep = three_quad_snapshot;
                diagnostics.push(rational_gregory_dispatch_declined_diagnostic(
                    n_sides,
                    "self_intersection_overlap",
                    None,
                ));
                if let Some(snapshot) = n5_entry_snapshot.take() {
                    *brep = snapshot;
                    return Ok((Vec::new(), None, diagnostics));
                }
                let mut r = base_patch.clone();
                r.overlap_remediated = true;
                return Ok((vec![r], Some(dispatch_boundary.clone()), diagnostics));
            }
            Err(ThreeQuadFailure::Other(reason)) => {
                *brep = three_quad_snapshot;
                diagnostics.push(rational_gregory_dispatch_declined_diagnostic(
                    n_sides,
                    "construction_error",
                    Some(reason),
                ));
            }
        }
    }

    if let Some(snapshot) = n5_entry_snapshot.take() {
        // PK-089: no route produced an accepted N-patch tiling — restore
        // the corner (mitred) instead of shipping the partial base patch.
        // The dispatch-decline diagnostics above record why; the caller
        // records the vertex-mode fallback warning.
        *brep = snapshot;
        diagnostics.push(rational_gregory_dispatch_declined_diagnostic(
            n_sides,
            "n5_all_or_nothing_rollback",
            None,
        ));
        return Ok((Vec::new(), None, diagnostics));
    }
    record_fallback_cap_support_continuity(brep, base_patch.face, &mut diagnostics);
    Ok((vec![base_patch], n_gon_boundary, diagnostics))
}

/// Kef every V-edge that separates two distinct faces.
///
/// Iterates until no cross-face V-edges remain. For trihedral
/// corners this merges (N-1) face-pairs out of N; the last V-edge
/// becomes a same-face bridge which `kef` cannot remove (its
/// precondition requires distinct loops). The bridge is handled
/// by `kev_v_bridge` in the kev pass, after the surface has been
/// assigned to the merged face.
///
/// Returns the surviving merged-face handle.
///
/// `pub(crate)` because the §6c-DirectApex.2.e apex builder
/// (`crate::apex_construction::construct`) calls this after its
/// own mef sequence to merge the corner sub-faces.
pub(crate) fn kef_cross_face_v_edges(
    brep: &mut BRep,
    vertex: VertexHandle,
) -> Result<FaceHandle, String> {
    let mut blend_face: Option<FaceHandle> = None;
    let mut safety = brep.edges.len() + 16;

    loop {
        if safety == 0 {
            return Err("kef_cross_face_v_edges: iteration cap hit".to_string());
        }
        safety -= 1;

        if !brep.vertices.contains_key(vertex) {
            break;
        }

        // `kef` invalidates V's `outgoing` field when the removed
        // half-edge *was* V's outgoing. Re-pick a valid outgoing
        // before every incidence walk.
        repair_vertex_outgoing(brep, vertex)?;

        let edges_through_v = collect_edges_through_vertex(brep, vertex)?;
        if edges_through_v.is_empty() {
            break;
        }

        let mut kef_target: Option<HalfEdgeHandle> = None;
        for edge in &edges_through_v {
            if !brep.edges.contains_key(*edge) {
                continue;
            }
            let he = brep.edges[*edge].half_edge;
            let he_loop = brep.half_edges[he].loop_;
            let twin_loop = brep.half_edges[brep.half_edges[he].twin].loop_;
            if he_loop != twin_loop {
                kef_target = Some(he);
                break;
            }
        }

        let Some(he) = kef_target else {
            // Only same-face bridges remain — the cross-face kef pass is done.
            // Record the pinched face as the surviving blend face
            // (it's the loop that currently holds V).
            if blend_face.is_none() {
                let outgoing = brep.vertices[vertex].outgoing;
                if brep.half_edges.contains_key(outgoing) {
                    let loop_h = brep.half_edges[outgoing].loop_;
                    blend_face = Some(brep.loops[loop_h].face);
                }
            }
            break;
        };
        let he_face = brep.loops[brep.half_edges[he].loop_].face;
        if popkern_core::trace::trace_enabled("corner-carve") {
            let twin_face = brep.loops[brep.half_edges[brep.half_edges[he].twin].loop_].face;
            let count = |fh: FaceHandle| {
                let start = brep.loops[brep.faces[fh].outer_loop].half_edge;
                let mut he2 = start;
                let mut n = 0;
                for _ in 0..256 {
                    n += 1;
                    he2 = brep.half_edges[he2].next;
                    if he2 == start {
                        break;
                    }
                }
                n
            };
            eprintln!(
                "  [vblend] kef at {vertex:?}: merging {he_face:?} (loop {}) with {twin_face:?} (loop {})",
                count(he_face),
                count(twin_face)
            );
        }
        brep.kef(he)
            .map_err(|e| format!("kef failed through vertex {vertex:?}: {e}"))?;
        blend_face = Some(he_face);
    }

    blend_face.ok_or_else(|| {
        "kef_cross_face_v_edges: no merge performed — blend face not created".to_string()
    })
}
