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
#[path = "blend_insert_support_cell_acceptance.rs"]
mod blend_insert_support_cell_acceptance;
#[path = "blend_insert_surface.rs"]
mod blend_insert_surface;
#[path = "blend_insert_three_quad_acceptance.rs"]
mod blend_insert_three_quad_acceptance;
#[path = "blend_insert_topology.rs"]
mod blend_insert_topology;
use crate::types::{
    ApexConstructionDiagnostic, ContinuityLevel, CornerMethod, CornerPatchInfo, VertexMode,
};
use crate::vertex_classifier::VertexBlendClassification;
use blend_insert_diagnostics::{
    rational_gregory_dispatch_declined_diagnostic, rational_gregory_dispatch_precondition_reason,
};
pub(crate) use blend_insert_surface::reclassify_corner_methods;
use blend_insert_surface::{assign_blend_surface, corner_radii_uniform};
use blend_insert_three_quad_acceptance::record_fallback_cap_support_continuity;
use blend_insert_topology::{normalize_underresolved_face_loop_to_n, split_corner_triangles};

type InsertVertexBlendOutput = (
    Vec<CornerPatchInfo>,
    Option<BlendBoundary>,
    Vec<ApexConstructionDiagnostic>,
);

/// Insert a smooth G1 vertex blend face at a multi-fillet vertex.
///
/// For each face touching the vertex, splits off the corner triangle using
/// `mef()`, then merges all corner triangles into one N-gon by removing
/// edges through the vertex using `kef()`. Finally assigns the blend surface
/// geometry to the merged face.
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
    // PK-094 stage-2 WIP arm: setback boundary surgery (gated; see
    // setback_vertex_topology module docs). On success the blend face
    // is the setback n-gon; the cap face fan is the recorded
    // follow-up, so for now the outcome is only exercised by
    // integration tests and the arm falls through on any error.
    if super::setback_vertex_topology::setback_topology_commit_enabled() && edge_indices.len() >= 5
    {
        // PK-094 commit mode: the setback vertex blend replaces the
        // legacy corner machinery entirely for eligible apexes. On
        // any surgery error the snapshot restore inside the surgery
        // has already run; fall through to the legacy path with a
        // recorded reason.
        let commit_snapshot = brep.clone();
        // The KP20 cap is the default for BOTH continuity levels: it
        // is the proven watertight tiling and a G2 cap satisfies G1.
        // The PnS bi-4 family stays reachable as a diagnostic /
        // comparison path (POPKERN_G1_PNS_CAP=1); it integrates at
        // the legacy PnS quality class (5e-2 seams — its patches
        // overlap by construction, see PK-097).
        let family = match requested_continuity {
            ContinuityLevel::G2 => super::setback_vertex_topology::CapFamily::Kp20G2,
            _ if std::env::var_os("POPKERN_G1_PNS_CAP").is_some() => {
                super::setback_vertex_topology::CapFamily::PnsG1
            }
            _ => super::setback_vertex_topology::CapFamily::Kp20G2,
        };
        match super::setback_vertex_topology::try_setback_boundary_surgery(
            brep,
            vertex,
            &fillet_surfaces_pre,
            tolerance,
            family,
        ) {
            Ok(outcome) => {
                if let Err(e) = brep.validate() {
                    *brep = commit_snapshot;
                    *super::setback_vertex_topology::SETBACK_TOPOLOGY_LAST
                        .lock()
                        .unwrap_or_else(|p| p.into_inner()) =
                        Some(format!("setback commit failed validation: {e}"));
                    return Err(format!("setback vertex blend failed validation: {e}"));
                }
                let infos: Vec<CornerPatchInfo> = outcome
                    .cap_faces
                    .iter()
                    .map(|&face| CornerPatchInfo {
                        vertex,
                        n_sides: edge_indices.len(),
                        method: match family {
                            super::setback_vertex_topology::CapFamily::Kp20G2 => {
                                CornerMethod::Gregory
                            }
                            super::setback_vertex_topology::CapFamily::PnsG1 => {
                                CornerMethod::PolyhedralNetSpline
                            }
                        },
                        // The blend's continuity is capped by the
                        // strips: a KP20 (G2-interior) cap on G1
                        // sweeps yields a G1 blend and must say so.
                        achieved_continuity: match (family, requested_continuity) {
                            (
                                super::setback_vertex_topology::CapFamily::Kp20G2,
                                ContinuityLevel::G2,
                            ) => crate::types::ContinuityLevel::G2,
                            _ => crate::types::ContinuityLevel::G1,
                        },
                        classification,
                        face,
                        bridge_fallbacks: Vec::new(),
                        overlap_remediated: false,
                        relaxed_continuity_coplanarity: None,
                        cap_sphere: None,
                        sub_patch_provenance: None,
                        setback_cap: true,
                    })
                    .collect();
                return Ok((infos, None, Vec::new()));
            }
            Err(e) => {
                // Fall through to the legacy path.
                *super::setback_vertex_topology::SETBACK_TOPOLOGY_LAST
                    .lock()
                    .unwrap_or_else(|p| p.into_inner()) =
                    Some(format!("setback commit declined: {e}"));
            }
        }
    }

    let local_indices_pre: Vec<usize> = (0..fillet_surfaces_pre.len()).collect();
    let n_sides = edge_indices.len();

    // Cap-kind dispatch: when apex topology rebuild is enabled and
    // the apex's incident fillets all share a uniform radius within
    // tolerance, suppress sphere-fit substitution for this corner.
    // The post-loop rebuild pass then collapses the rolling-ball
    // 6-vertex cap loop to a 3-vertex cap-tangent-arc loop directly.
    // For non-uniform corners, substitution and rebuild stay active.
    let apex_topology_rebuild_active = std::env::var("POPKERN_APEX_PATH_A")
        .map(|v| v != "0")
        .unwrap_or(true);
    let radii_uniform = corner_radii_uniform(edge_indices, edge_geometries, tolerance);
    let auto_disable_substitution = apex_topology_rebuild_active && radii_uniform;
    let boundary_options = crate::blend::BlendBoundaryOptions {
        disable_sphere_fit_substitution: if auto_disable_substitution {
            Some(true)
        } else {
            None
        },
    };
    let rational_gregory_compatible = matches!(
        classification,
        VertexBlendClassification::SameSideConvex
            | VertexBlendClassification::SameSideConcave
            | VertexBlendClassification::OppositeSide,
    );
    // Compute the n-gon BlendBoundary when either the rational-
    // Gregory dispatch OR the env-gated PnS path is enabled — both
    // consume the same boundary structure downstream.
    let pns_wants_boundary = crate::blend::pns_promote::pns_env_enabled();
    let n_gon_boundary: Option<BlendBoundary> = if (use_rational_gregory || pns_wants_boundary)
        && n_sides >= 3
        && rational_gregory_compatible
    {
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

    // Phase 1: For each face, cut a chord between the face's two
    // fillet-trim-endpoint anchors. The V-containing sub-face is
    // the "corner region" (may be a triangle, quad, or pentagon
    // depending on how many intermediate cube-sub-split vertices
    // sit between V and the anchors).
    let _corner_faces =
        split_corner_triangles(brep, vertex, &face_info, edge_geometries, tolerance)?;

    // Phase 2a: kef every V-edge separating two distinct faces;
    // last V-edge may become a same-face "bridge" handled by
    // Phase 2b's `kev` after surface assignment.
    let blend_face = kef_cross_face_v_edges(brep, vertex)?;

    // Phase 3: Compute the blend boundary + assign the surface.
    let fillet_surfaces = build_fillet_surface_list(edge_indices, edge_geometries);
    let local_indices: Vec<usize> = (0..fillet_surfaces.len()).collect();
    let blend_boundary = blend::compute_blend_boundary_with_options(
        vertex,
        &local_indices,
        &fillet_surfaces,
        brep,
        tolerance,
        boundary_options,
    )
    .map_err(|e| format!("blend boundary computation failed: {e}"))?;
    let (method, achieved_continuity) = assign_blend_surface(
        brep,
        blend_face,
        &blend_boundary,
        &fillet_surfaces,
        classification,
        tolerance,
    )?;

    // Phase 2b: kev the remaining V-bridge if any. `kev` kills V
    // and the bridge edge simultaneously, un-pinching the merged
    // loop into a manifold polygon. This is R8 C6's crucial
    // topology fix — the old code left the bridge as an internal
    // dangling edge because its `kef` guard rejected same-loop
    // half-edge pairs.
    kev_v_bridge(brep, vertex)?;

    remove_orphaned_vertex(brep, vertex);

    // Phase 6c-D (N-path arm): when sphere-fit substitution
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
    };
    let mut diagnostics = Vec::new();

    if crate::blend::pns_promote::pns_env_enabled()
        && n_sides >= 5
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

    let rational_dispatch_preconditions_met = dispatch_boundary.segments.len() == 2 * n_sides
        && crate::construct::poke_hex::face_has_2n_outer_half_edges_with(
            brep,
            blend_face,
            n_sides,
            dispatch_boundary.sphere_fit_substituted,
        );
    if use_rational_gregory
        && use_kp_polar_gregory
        && rational_gregory_compatible
        && n_sides >= 5
        && !rational_dispatch_preconditions_met
    {
        diagnostics.push(rational_gregory_dispatch_declined_diagnostic(
            n_sides,
            "precondition_failed",
            Some(rational_gregory_dispatch_precondition_reason(
                brep,
                blend_face,
                n_sides,
                dispatch_boundary,
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
                    brep, result, n_sides,
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
                let mut r = base_patch.clone();
                r.overlap_remediated = true;
                diagnostics.push(rational_gregory_dispatch_declined_diagnostic(
                    n_sides,
                    "self_intersection_overlap",
                    None,
                ));
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

    record_fallback_cap_support_continuity(brep, base_patch.face, &mut diagnostics);
    Ok((vec![base_patch], n_gon_boundary, diagnostics))
}

/// Phase 2a: kef every V-edge that separates two distinct faces.
///
/// Iterates until no cross-face V-edges remain. For trihedral
/// corners this merges (N-1) face-pairs out of N; the last V-edge
/// becomes a same-face bridge which `kef` cannot remove (its
/// precondition requires distinct loops). The bridge is handled
/// by `kev_v_bridge` in Phase 2b, after the surface has been
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
            // Only same-face bridges remain — Phase 2a is done.
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
        if std::env::var_os("PK_CARVE_TRACE").is_some() {
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

/// Phase 2b: kev the remaining V-bridge edge (if any) to un-pinch
/// the merged face and remove V in one Euler step.
///
/// At entry, V may have zero or one incident edges. Zero means
/// Phase 2a completed cleanly (all V-edges were cross-face kefs);
/// one means there's a same-loop bridge whose both half-edges are
/// on the merged face. `kev` requires the edge's destination
/// vertex to be valence-1, so we pick the half-edge that ends at
/// V and call `kev(that)` to kill V plus the bridge together.
///
/// `pub(crate)` for §6c-DirectApex.2.e — same rationale as
/// [`kef_cross_face_v_edges`].
pub(crate) fn kev_v_bridge(brep: &mut BRep, vertex: VertexHandle) -> Result<(), String> {
    if !brep.vertices.contains_key(vertex) {
        return Ok(());
    }

    let outgoings: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, h)| h.origin == vertex)
        .map(|(h, _)| h)
        .collect();

    if outgoings.is_empty() {
        return Ok(());
    }
    if outgoings.len() != 1 {
        return Err(format!(
            "kev_v_bridge: V has valence {} (expected 0 or 1 for single-bridge teardown)",
            outgoings.len()
        ));
    }

    let he_out = outgoings[0];
    let he_in = brep.half_edges[he_out].twin;
    // kev kills the destination vertex of its argument; we want to
    // kill V, so pass the half-edge whose destination is V.
    debug_assert_eq!(brep.half_edges[brep.half_edges[he_in].twin].origin, vertex);
    let same_loop = brep.half_edges[he_in].loop_ == brep.half_edges[he_out].loop_;
    let is_only_edge = brep.half_edges[he_in].next == he_out
        && brep.half_edges[he_out].next == he_in
        && brep.half_edges[he_in].prev == he_out
        && brep.half_edges[he_out].prev == he_in;
    match brep.kev(he_in) {
        Ok(_) => Ok(()),
        Err(_) if same_loop && !is_only_edge => {
            brep.dissolve_same_loop_bridge(he_in)
                .map_err(|e| format!("same-loop bridge dissolve failed on V-bridge: {e}"))?;
            remove_orphaned_vertex(brep, vertex);
            Ok(())
        }
        Err(kev_error) => Err(format!("kev failed on V-bridge: {kev_error}")),
    }
}

/// Phase 3: Remove vertex V if it is no longer connected to any edges.
///
/// `pub(crate)` for §6c-DirectApex.2.e — same rationale as
/// [`kef_cross_face_v_edges`].
pub(crate) fn remove_orphaned_vertex(brep: &mut BRep, vertex: VertexHandle) {
    let v_has_edges = brep.half_edges.values().any(|he| he.origin == vertex);
    if !v_has_edges {
        brep.vertices.remove(vertex);
    }
}
