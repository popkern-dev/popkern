//! Vertex blend face insertion into BRep topology.
//!
//! At a vertex where multiple fillets meet, carves a blend face from the
//! surrounding faces using Euler operators (`mef` / `kef`). The blend face
//! replaces the small triangular region around the vertex.

use popkern_core::geometry::surface::SurfaceOnFace;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::blend;
use crate::blend::BlendBoundary;
use crate::construct::EdgeGeometry;
use crate::construct::blend_topo::{
    build_fillet_surface_list, collect_edges_through_vertex, collect_faces_around_vertex,
    find_face_trim_anchors, find_halfedge_from_vertex_in_face, identify_corner_face,
    repair_vertex_outgoing,
};
use crate::construct::poke_hex::face_has_2n_outer_half_edges;
use crate::construct::three_quad_apply::apply_three_quad_decomposition_to_blend;

#[path = "blend_insert_surface.rs"]
mod blend_insert_surface;
use crate::types::{ContinuityLevel, CornerMethod, CornerPatchInfo, VertexMode};
use crate::vertex_classifier::VertexBlendClassification;
pub(crate) use blend_insert_surface::reclassify_corner_methods;
use blend_insert_surface::{assign_blend_surface, corner_radii_uniform};

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
    use_rational_gregory: bool,
    use_deg5_gregory: bool,
    use_biquintic_full_gregory: bool,
    use_biseptic_gregory: bool,
    use_two_stripe_gregory: bool,
    use_kp_binonic_gregory: bool,
    use_kp_polar_gregory: bool,
    precomputed_classification: Option<VertexBlendClassification>,
) -> Result<(Vec<CornerPatchInfo>, Option<BlendBoundary>), String> {
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

    // Phase 3 orphan (idempotent safety net for no-bridge cases).
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
    let bridge_fallbacks = n_gon_boundary
        .as_ref()
        .map(|b| b.bridge_fallbacks.clone())
        .unwrap_or_default();
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
    };

    // PnS rung (env-gated). At N=3 the rational-Gregory ladder
    // delivers 40–85% lower boundary-κ deviation than PnS bi-4 on
    // asymmetric-radius fixtures (R8 C8b cross-fixture sweep), so
    // PnS only fires for N ≥ 5 where the Gregory ladder isn't
    // competitive. N=4 is rejected internally (σ singularity).
    if crate::blend::pns_promote::pns_env_enabled()
        && n_sides >= 5
        && let Some(boundary) = n_gon_boundary.as_ref()
        && boundary.segments.len() == 2 * n_sides
        && let Some(patches) = crate::blend::pns_promote::try_pns_dispatch(
            brep,
            blend_face,
            boundary,
            &fillet_surfaces_pre,
            vertex,
            classification,
            boundary.bridge_fallbacks.clone(),
        )
    {
        return Ok((patches, n_gon_boundary));
    }

    // §5.3 item #4 — N ≥ 3 rational-Gregory dispatch.
    if use_rational_gregory
        && rational_gregory_compatible
        && let Some(boundary) = n_gon_boundary.as_ref()
        && boundary.segments.len() == 2 * n_sides
        && crate::construct::poke_hex::face_has_2n_outer_half_edges_with(
            brep,
            blend_face,
            n_sides,
            boundary.sphere_fit_substituted,
        )
    {
        use crate::construct::three_quad_apply::ThreeQuadFailure;
        match apply_three_quad_decomposition_to_blend(
            brep,
            blend_face,
            boundary,
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
            Ok(patches) if patches.len() == n_sides => return Ok((patches, n_gon_boundary)),
            Err(ThreeQuadFailure::Overlap) => {
                let mut r = base_patch.clone();
                r.overlap_remediated = true;
                return Ok((vec![r], n_gon_boundary));
            }
            _ => {}
        }
    }

    Ok((vec![base_patch], n_gon_boundary))
}

/// Phase 1: Split off corner regions from each face around V using
/// the per-face fillet trim endpoints as mef anchors.
fn split_corner_triangles(
    brep: &mut BRep,
    vertex: VertexHandle,
    face_info: &[(FaceHandle, HalfEdgeHandle)],
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Result<Vec<FaceHandle>, String> {
    let mut corner_faces: Vec<FaceHandle> = Vec::new();

    for (face, _he_arriving_at_v) in face_info {
        let anchors = find_face_trim_anchors(brep, *face, vertex, edge_geometries, tolerance)?;
        if anchors.len() < 2 {
            // Falling back to the legacy nearest-to-V anchors keeps
            // low-N or atypical configurations working during the R8
            // transition. Chamfer/mitred paths and non-trihedral
            // corners hit this.
            return Err(format!(
                "face {face:?} has {} trim-endpoint anchors around vertex {vertex:?}, expected ≥2",
                anchors.len()
            ));
        }
        // For trihedral corners we expect exactly 2. If more appear
        // (degenerate), just take the first two — the V-containing
        // sub-face still gets cut off.
        let (anchor_a, anchor_b) = (anchors[0], anchors[1]);

        let he_at_a = find_halfedge_from_vertex_in_face(brep, *face, anchor_a)?;
        let he_at_b = find_halfedge_from_vertex_in_face(brep, *face, anchor_b)?;

        // Verify both half-edges are in the same loop before calling mef.
        let loop_a = brep.half_edges[he_at_a].loop_;
        let loop_b = brep.half_edges[he_at_b].loop_;
        if loop_a != loop_b {
            return Err(format!(
                "anchor half-edges on face {face:?} are in different loops"
            ));
        }

        let mef_result = brep
            .mef(he_at_a, he_at_b)
            .map_err(|e| format!("mef failed on face {face:?}: {e}"))?;

        let corner = identify_corner_face(brep, *face, mef_result.face, vertex)?;
        corner_faces.push(corner);
    }

    Ok(corner_faces)
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
/// `pub(crate)` because the §6c-A2.2.e apex builder
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
/// `pub(crate)` for §6c-A2.2.e — same rationale as
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
/// `pub(crate)` for §6c-A2.2.e — same rationale as
/// [`kef_cross_face_v_edges`].
pub(crate) fn remove_orphaned_vertex(brep: &mut BRep, vertex: VertexHandle) {
    let v_has_edges = brep.half_edges.values().any(|he| he.origin == vertex);
    if !v_has_edges {
        brep.vertices.remove(vertex);
    }
}
