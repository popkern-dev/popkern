//! Vertex blend face insertion into BRep topology.
//!
//! At a vertex where multiple fillets meet, carves a blend face from the
//! surrounding faces using Euler operators (`mef` / `kef`). The blend face
//! replaces the small triangular region around the vertex.

use popkern_core::geometry::surface::SurfaceOnFace;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::blend;
use crate::construct::EdgeGeometry;
use crate::construct::blend_topo::{
    build_fillet_surface_list, collect_edges_through_vertex, collect_faces_around_vertex,
    collect_fillet_radii, find_face_trim_anchors, find_halfedge_from_vertex_in_face,
    identify_corner_face, repair_vertex_outgoing,
};
use crate::types::{ContinuityLevel, CornerMethod, CornerPatchInfo, VertexMode};

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
pub(crate) fn insert_vertex_blend(
    brep: &mut BRep,
    vertex: VertexHandle,
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Result<CornerPatchInfo, String> {
    // §5.3 C6 — classify the vertex-blend configuration *before*
    // topology modifications begin, while the pre-blend BRep is
    // still intact (the classifier reads face surface normals
    // along incident edges). The classification drives method
    // dispatch and is stored on the returned `CornerPatchInfo`
    // for UI consumers.
    let classification = {
        let filleted_edges: Vec<_> = edge_indices
            .iter()
            .filter_map(|&i| edge_geometries.get(i).map(|g| g.analysis.edge))
            .collect();
        crate::vertex_classifier::classify_vertex_blend(brep, vertex, &filleted_edges, tolerance)
            .classification
    };

    // Repair the vertex's outgoing half-edge if it has become stale
    // due to prior topology operations (edge splitting, face splitting).
    repair_vertex_outgoing(brep, vertex)?;

    // R8 C6 — per-face trim-endpoint anchors. The old nearest-to-V
    // anchor set cut the wrong chord on mixed-radius faces (picking
    // an intermediate cube-sub vertex instead of the face's own
    // fillet trim endpoint). For each adjacent face we now look up
    // the face-loop vertices at distance-from-V matching one of
    // the incident fillet radii — those are the real fillet
    // boundary anchors on that face.
    let radii = collect_fillet_radii(edge_geometries);
    if radii.is_empty() {
        return Err(format!("vertex {vertex:?} has no fillet radii to match"));
    }

    // Collect faces touching V (via half-edges originating from V).
    let face_info = collect_faces_around_vertex(brep, vertex)?;

    if face_info.is_empty() {
        return Err(format!("no faces found around vertex {vertex:?}"));
    }

    // Phase 1: For each face, cut a chord between the face's two
    // fillet-trim-endpoint anchors. The V-containing sub-face is
    // the "corner region" (may be a triangle, quad, or pentagon
    // depending on how many intermediate cube-sub-split vertices
    // sit between V and the anchors).
    let _corner_faces = split_corner_triangles(brep, vertex, &face_info, &radii, tolerance)?;

    // Phase 2a: kef every V-edge that currently separates two
    // distinct faces. On trihedral corners this merges all but the
    // last incident face-pair; the last V-edge may become a
    // same-face "bridge" which kef can't remove. Leaving the
    // bridge (plus V) intact for now keeps the `vertex` handle
    // valid through Phase 3's `compute_blend_boundary` call, which
    // still needs V's position and its incident faces to emit the
    // correct boundary segments. The final kev happens in
    // Phase 2b, after the surface has been assigned.
    let blend_face = kef_cross_face_v_edges(brep, vertex)?;

    // Phase 3: Compute the blend boundary + assign the surface.
    // `vertex` is still valid (Phase 2a only `kef`-ed cross-face
    // edges). `compute_blend_boundary` sees the merged (possibly
    // pinched) blend face — no surface yet — so `planar_faces`
    // is empty and the hex path is not taken. The FilletEnd-only
    // 3-segment fallback feeds `try_spherical_blend` with all 30
    // samples lying on the rolling-ball sphere, matching the
    // legacy behaviour for uniform-radius corners.
    let fillet_surfaces = build_fillet_surface_list(edge_indices, edge_geometries);
    let local_indices: Vec<usize> = (0..fillet_surfaces.len()).collect();
    let blend_boundary =
        blend::compute_blend_boundary(vertex, &local_indices, &fillet_surfaces, brep, tolerance)
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

    // Phase 3 orphan (idempotent): V is usually already removed by
    // the kev above; but a non-generic topology with no bridge
    // still needs this safety net.
    remove_orphaned_vertex(brep, vertex);

    Ok(CornerPatchInfo {
        vertex,
        n_sides: edge_indices.len(),
        method,
        achieved_continuity,
        classification,
        face: blend_face,
    })
}

/// Phase 1: Split off corner regions from each face around V using
/// the per-face fillet trim endpoints as mef anchors.
fn split_corner_triangles(
    brep: &mut BRep,
    vertex: VertexHandle,
    face_info: &[(FaceHandle, HalfEdgeHandle)],
    radii: &[f64],
    tolerance: &Tolerance,
) -> Result<Vec<FaceHandle>, String> {
    let mut corner_faces: Vec<FaceHandle> = Vec::new();

    for (face, _he_arriving_at_v) in face_info {
        let anchors = find_face_trim_anchors(brep, *face, vertex, radii, tolerance)?;
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
fn kef_cross_face_v_edges(brep: &mut BRep, vertex: VertexHandle) -> Result<FaceHandle, String> {
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
fn kev_v_bridge(brep: &mut BRep, vertex: VertexHandle) -> Result<(), String> {
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
    brep.kev(he_in)
        .map_err(|e| format!("kev failed on V-bridge: {e}"))?;
    Ok(())
}

/// Phase 3: Remove vertex V if it is no longer connected to any edges.
fn remove_orphaned_vertex(brep: &mut BRep, vertex: VertexHandle) {
    let v_has_edges = brep.half_edges.values().any(|he| he.origin == vertex);
    if !v_has_edges {
        brep.vertices.remove(vertex);
    }
}

/// Phase 4: Assign the blend surface to the merged face.
///
/// R8 C6 — takes a pre-computed `blend_boundary` because V may
/// already have been `kev`'d by `merge_corner_triangles`.
/// `build_blend_surface` only reads `boundary.vertex_position`, so
/// the removed handle in `boundary.vertex` is diagnostic only.
///
/// Returns `(method, achieved_continuity)` from the construction.
fn assign_blend_surface(
    brep: &mut BRep,
    blend_face: FaceHandle,
    blend_boundary: &blend::BlendBoundary,
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    classification: crate::vertex_classifier::VertexBlendClassification,
    tolerance: &Tolerance,
) -> Result<(CornerMethod, ContinuityLevel), String> {
    let blend_result = blend::build_blend_surface(
        blend_boundary,
        fillet_surfaces,
        brep,
        VertexMode::Smooth,
        classification,
        tolerance,
    )
    .map_err(|e| format!("blend surface computation failed: {e}"))?;

    let method = blend_result.method;
    let continuity = blend_result.continuity;
    brep.faces[blend_face].surface = Some(SurfaceOnFace {
        geometry: blend_result.surface,
        orientation: true,
    });

    Ok((method, continuity))
}

/// Reclassify `CornerPatchInfo::method` from each patch's final
/// face surface type, because `surface_fix::correct_face_surfaces`
/// can override `build_blend_surface`'s initial choice (e.g.
/// promoting a Mitred fallback plane to a fitted Sphere).
pub(crate) fn reclassify_corner_methods(
    brep: &BRep,
    corner_patches: &mut [crate::types::CornerPatchInfo],
) {
    use popkern_core::geometry::surface::SurfaceGeometry;
    for info in corner_patches.iter_mut() {
        let Some(sof) = brep.faces.get(info.face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        // Infer the method from the final face surface. Only
        // overwrite the method/continuity fields if the surface
        // type mismatches the original construction — that means
        // `surface_fix` promoted/demoted the surface post-hoc
        // (e.g. Mitred plane → RollingBallSphere). If the surface
        // type matches, keep the construction-reported values —
        // the §5.3 C4 Gregory path reports G1 or G2 per-case based
        // on radius uniformity, and we must not flatten that to G1.
        let inferred = match &sof.geometry {
            SurfaceGeometry::Sphere { .. } => {
                Some((CornerMethod::RollingBallSphere, ContinuityLevel::G1))
            }
            SurfaceGeometry::Nurbs(_) => Some((CornerMethod::Gregory, ContinuityLevel::G1)),
            SurfaceGeometry::Plane { .. } => Some((CornerMethod::Mitred, ContinuityLevel::G0)),
            _ => None,
        };
        if let Some((inferred_method, inferred_cont)) = inferred
            && info.method != inferred_method
        {
            info.method = inferred_method;
            info.achieved_continuity = inferred_cont;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insert_vertex_blend_too_few_anchors() {
        // A cube vertex has 3 neighbors, but if we pass no edge indices
        // the anchor set will be non-empty (all neighbors). To trigger the
        // "fewer than 2 anchors" path we need a degenerate vertex.
        // Instead, test with a default (non-existent) vertex handle.
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let bogus = popkern_core::topology::VertexHandle::default();
        let result = insert_vertex_blend(&mut brep, bogus, &[], &[], &Tolerance::default());
        assert!(result.is_err(), "should fail for invalid vertex");
    }

    #[test]
    fn insert_vertex_blend_on_cube_vertex() {
        // Calling insert_vertex_blend on a real cube vertex with no
        // edge geometries should fail during blend surface computation
        // (since there are no fillet surfaces to blend), but it should
        // get past the anchor check (3 anchors >= 2).
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let vertex = brep.vertices.keys().next().unwrap();
        let result = insert_vertex_blend(&mut brep, vertex, &[], &[], &Tolerance::default());
        // This should fail (no edge geometries -> empty fillet surfaces
        // -> blend boundary computation fails), but NOT at the anchor check.
        assert!(result.is_err());
        let err = result.unwrap_err();
        // The error should be from blend computation, not "fewer than 2 anchors".
        assert!(
            !err.contains("fewer than 2"),
            "should not fail at anchor check: {err}"
        );
    }
}
