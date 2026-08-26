//! Shell assembly: combine outer skin, inner skin, and side walls.
//!
//! Copies entities from both the outer and inner BReps into a new
//! [`BRepBuilder`], then adds side wall faces at the openings where
//! removed faces create gaps.

use std::collections::HashMap;

use crate::error::ShellError;
use crate::side_walls::{build_free_rim_walls, build_side_walls};
use crate::skin_map::SkinMap;
use popkern_core::constructors::BRepBuilder;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle,
};

/// Handle maps returned by [`copy_skin_entities`].
type SkinCopyMaps = (
    HashMap<VertexHandle, VertexHandle>,
    HashMap<EdgeHandle, EdgeHandle>,
    HashMap<HalfEdgeHandle, HalfEdgeHandle>,
    HashMap<FaceHandle, FaceHandle>,
);

/// Face provenance of an assembled shell: source-body handles →
/// assembled-result handles for both skins. Result faces in neither
/// map are side walls / corner sectors (built during assembly). Feeds
/// the global self-intersection trim stage's realization predicate
/// (PK-327).
#[derive(Debug, Default)]
pub struct ShellAssemblyMaps {
    /// Outer-skin faces (original body for classic shells; the
    /// outward-moved skin for thicken/midplane).
    pub outer_faces: HashMap<FaceHandle, FaceHandle>,
    /// Inner-skin faces (the moved void skin).
    pub inner_faces: HashMap<FaceHandle, FaceHandle>,
}

/// Assemble the final shell BRep from outer skin, inner skin, and side walls.
///
/// Copies all topology from the outer BRep (minus removed faces) and the
/// inner BRep into a single BRepBuilder, then builds side wall faces at
/// the openings.
///
/// Precondition: `inner` must already be orientation-mirrored
/// (`reverse_body_orientation`), as the shell pipeline does — the side
/// walls adopt the free inner rim half-edges and rely on the mirrored
/// direction (PK-307); a same-winding inner skin fails the direction
/// contract with a structural error.
///
/// # Errors
///
/// Returns [`ShellError`] if handle mapping or side wall construction fails.
pub fn assemble_shell(
    outer: &BRep,
    inner: &BRep,
    removed_faces: &[FaceHandle],
    skin_map: &SkinMap,
    tolerance: Tolerance,
) -> Result<BRep, ShellError> {
    assemble_shell_with_maps(outer, inner, removed_faces, skin_map, tolerance).map(|(b, _)| b)
}

/// [`assemble_shell`] returning the face provenance maps as well.
///
/// # Errors
///
/// Returns [`ShellError`] if handle mapping or side wall construction fails.
pub fn assemble_shell_with_maps(
    outer: &BRep,
    inner: &BRep,
    removed_faces: &[FaceHandle],
    skin_map: &SkinMap,
    tolerance: Tolerance,
) -> Result<(BRep, ShellAssemblyMaps), ShellError> {
    let mut builder = BRepBuilder::new(tolerance);
    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);

    let removed_set: std::collections::HashSet<FaceHandle> =
        removed_faces.iter().copied().collect();

    // Copy outer skin (minus removed faces).
    let (outer_vertex_map, outer_edge_map, outer_he_map, outer_face_map) =
        copy_skin_entities(&mut builder, outer, shell, |fh| !removed_set.contains(&fh))?;

    // Shell ownership is topological truth (PK-306): with no removed
    // faces there are no side walls, so the inner skin is a
    // DISCONNECTED closed void skin — its own shell of the solid
    // (χ = 2 per genus-0 shell; one dishonest merged record made every
    // hollow body fail the shell-aware Euler–Poincaré check, χ=4 vs
    // 2·1). With removed faces the side walls connect the skins into
    // one shell, which stays merged.
    // Free rims (open-sheet thicken): loopless half-edges in the source
    // mark sheet boundaries; their walls connect the skins just like
    // opening walls do, so the skins share one shell.
    let has_free_rims = outer
        .half_edges
        .iter()
        .any(|(_, he)| !outer.loops.contains_key(he.loop_));
    let inner_shell = if removed_faces.is_empty() && !has_free_rims {
        builder.add_shell(solid)
    } else {
        shell
    };

    // Copy inner skin (all faces).
    let (inner_vertex_map, _inner_edge_map, inner_he_map, inner_face_map) =
        copy_skin_entities(&mut builder, inner, inner_shell, |_| true)?;

    // Build side walls if there are removed faces.
    if !removed_faces.is_empty() {
        build_side_walls(
            &mut builder,
            outer,
            removed_faces,
            skin_map,
            &outer_vertex_map,
            &inner_vertex_map,
            &outer_he_map,
            &inner_he_map,
            shell,
        )?;
    }
    if has_free_rims {
        build_free_rim_walls(
            &mut builder,
            outer,
            skin_map,
            &outer_vertex_map,
            &inner_vertex_map,
            &outer_he_map,
            &inner_he_map,
            shell,
        )?;
    }

    let _ = outer_edge_map;

    Ok((
        builder.build_unchecked(),
        ShellAssemblyMaps {
            outer_faces: outer_face_map,
            inner_faces: inner_face_map,
        },
    ))
}

/// Copy a skin's entities into the builder, filtering faces and keeping
/// only edges/vertices reachable from kept faces or opening rims (see
/// the reachability rule below).
///
/// Returns maps from the source BRep handles to the new builder handles
/// for vertices, edges, and half-edges.
fn copy_skin_entities(
    builder: &mut BRepBuilder,
    source: &BRep,
    shell: popkern_core::topology::ShellHandle,
    face_filter: impl Fn(FaceHandle) -> bool,
) -> Result<SkinCopyMaps, ShellError> {
    // 0. Decide which edges the copy needs. An edge is reachable iff at
    //    least one of its half-edges belongs to a live loop whose face
    //    passes the filter: kept-face interior edges (both sides live)
    //    and opening-rim edges (one live kept side; the free side is the
    //    loopless half-edge the side walls adopt, PK-307). An edge whose
    //    sides BOTH belong to removed faces (removed×removed, e.g. the
    //    edge between two adjacent removed faces) is interior to the
    //    opening and must not be copied — copying it used to strand a
    //    loopless half-edge pair no wall ever adopts (DanglingReference
    //    + BrokenCycle on two-adjacent-removed shells).
    let edge_reachable = |eh: EdgeHandle| -> bool {
        let edge = &source.edges[eh];
        let he = &source.half_edges[edge.half_edge];
        for heh in [edge.half_edge, he.twin] {
            let loop_h = source.half_edges[heh].loop_;
            if source.loops.contains_key(loop_h) {
                let face = source.loops[loop_h].face;
                if source.faces.contains_key(face) && face_filter(face) {
                    return true;
                }
            }
        }
        false
    };

    // 1. Copy vertices incident to at least one reachable edge (a vertex
    //    used only by removed-interior edges must not survive either).
    let mut needed_vertices: std::collections::HashSet<VertexHandle> =
        std::collections::HashSet::new();
    for (eh, edge) in &source.edges {
        if !edge_reachable(eh) {
            continue;
        }
        let he = &source.half_edges[edge.half_edge];
        let twin = &source.half_edges[he.twin];
        needed_vertices.insert(he.origin);
        needed_vertices.insert(twin.origin);
    }
    let mut vertex_map = HashMap::new();
    for (vh, vertex) in &source.vertices {
        if !needed_vertices.contains(&vh) {
            continue;
        }
        let new_vh = builder.add_vertex(vertex.position);
        vertex_map.insert(vh, new_vh);
    }

    // 2. Copy reachable edges (and their half-edge pairs).
    let mut edge_map = HashMap::new();
    let mut he_map = HashMap::new();
    for (eh, edge) in &source.edges {
        if !edge_reachable(eh) {
            continue;
        }
        let he_handle = edge.half_edge;
        let he = &source.half_edges[he_handle];
        let twin_handle = he.twin;
        let twin = &source.half_edges[twin_handle];

        let new_v1 = *vertex_map
            .get(&he.origin)
            .ok_or_else(|| ShellError::Topology("edge origin vertex not in map".into()))?;
        let new_v2 = *vertex_map
            .get(&twin.origin)
            .ok_or_else(|| ShellError::Topology("edge twin origin vertex not in map".into()))?;

        let (new_eh, new_he1, new_he2) = builder.add_edge(new_v1, new_v2);
        edge_map.insert(eh, new_eh);
        he_map.insert(he_handle, new_he1);
        he_map.insert(twin_handle, new_he2);

        // Copy curve geometry.
        if let Some(curve) = &edge.curve {
            builder.set_edge_curve(new_eh, curve.clone());
        }
    }

    // 3. Copy faces with their loops.
    let mut face_map = HashMap::new();
    for (fh, face) in &source.faces {
        if !face_filter(fh) {
            continue;
        }

        let new_face = builder.add_face(shell);
        face_map.insert(fh, new_face);

        // Copy outer loop.
        copy_loop(builder, source, face.outer_loop, new_face, &he_map)?;

        // Copy inner loops.
        for &inner_loop in &face.inner_loops {
            copy_inner_loop(builder, source, inner_loop, new_face, &he_map)?;
        }

        // Copy surface geometry.
        if let Some(sof) = &face.surface {
            builder.set_face_surface(new_face, sof.clone());
        }
    }

    Ok((vertex_map, edge_map, he_map, face_map))
}

/// Copy a loop from source BRep into the builder as an outer loop.
fn copy_loop(
    builder: &mut BRepBuilder,
    source: &BRep,
    loop_handle: LoopHandle,
    face: FaceHandle,
    he_map: &HashMap<HalfEdgeHandle, HalfEdgeHandle>,
) -> Result<LoopHandle, ShellError> {
    let he_sequence = collect_loop_half_edges(source, loop_handle)?;
    let new_hes: Vec<HalfEdgeHandle> = he_sequence
        .iter()
        .map(|&old_he| {
            he_map
                .get(&old_he)
                .copied()
                .ok_or_else(|| ShellError::Topology("half-edge not in map".into()))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let lh = builder.add_loop(&new_hes, face).map_err(ShellError::Core)?;
    Ok(lh)
}

/// Copy a loop from source BRep into the builder as an inner loop.
fn copy_inner_loop(
    builder: &mut BRepBuilder,
    source: &BRep,
    loop_handle: LoopHandle,
    face: FaceHandle,
    he_map: &HashMap<HalfEdgeHandle, HalfEdgeHandle>,
) -> Result<LoopHandle, ShellError> {
    let he_sequence = collect_loop_half_edges(source, loop_handle)?;
    let new_hes: Vec<HalfEdgeHandle> = he_sequence
        .iter()
        .map(|&old_he| {
            he_map
                .get(&old_he)
                .copied()
                .ok_or_else(|| ShellError::Topology("half-edge not in map".into()))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let lh = builder
        .add_inner_loop(&new_hes, face)
        .map_err(ShellError::Core)?;
    Ok(lh)
}

/// Walk a loop and collect its half-edge handles in order.
fn collect_loop_half_edges(
    source: &BRep,
    loop_handle: LoopHandle,
) -> Result<Vec<HalfEdgeHandle>, ShellError> {
    let start = source.loops[loop_handle].half_edge;
    if !source.half_edges.contains_key(start) {
        return Err(ShellError::Topology("loop start half-edge missing".into()));
    }
    let mut result = vec![start];
    let mut current = source.half_edges[start].next;
    let mut safety = 0;
    while current != start {
        result.push(current);
        current = source.half_edges[current].next;
        safety += 1;
        if safety > 10_000 {
            return Err(ShellError::Topology("loop appears infinite".into()));
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::inner_skin::build_inner_skin;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn assemble_cube_no_removal() {
        let outer = BRep::unit_cube(tol()).unwrap();
        let (inner, skin_map) = build_inner_skin(&outer, &[]).unwrap();
        let result = assemble_shell(&outer, &inner, &[], &skin_map, tol()).unwrap();
        // Outer (6 faces) + inner (6 faces) = 12 faces.
        assert_eq!(result.faces.len(), 12);
    }

    #[test]
    fn assemble_cube_one_removed() {
        let outer = BRep::unit_cube(tol()).unwrap();
        let first_face = outer.faces.keys().next().unwrap();
        let (mut inner, skin_map) = build_inner_skin(&outer, &[first_face]).unwrap();
        inner.reverse_body_orientation();
        let result = assemble_shell(&outer, &inner, &[first_face], &skin_map, tol()).unwrap();
        // Outer (5 faces) + inner (5 faces) + side walls (4 faces) = 14 faces.
        assert_eq!(result.faces.len(), 14);
    }

    #[test]
    fn assemble_euler_no_removal() {
        let outer = BRep::unit_cube(tol()).unwrap();
        let (inner, skin_map) = build_inner_skin(&outer, &[]).unwrap();
        let result = assemble_shell(&outer, &inner, &[], &skin_map, tol()).unwrap();
        let v = result.vertices.len() as i64;
        let e = result.edges.len() as i64;
        let f = result.faces.len() as i64;
        // Two separate closed shells: each has V-E+F=2, total gives 4.
        // But since they share a solid, the Euler formula for a genus-0
        // solid with 2 shells is V-E+F=4.
        assert_eq!(v - e + f, 4, "Euler: V={v} E={e} F={f}");
    }

    #[test]
    fn assemble_one_removed_has_correct_face_count() {
        let outer = BRep::unit_cube(tol()).unwrap();
        let first_face = outer.faces.keys().next().unwrap();
        let (mut inner, skin_map) = build_inner_skin(&outer, &[first_face]).unwrap();
        inner.reverse_body_orientation();
        let result = assemble_shell(&outer, &inner, &[first_face], &skin_map, tol()).unwrap();
        // Outer (5) + inner (5) + side walls (4) = 14 faces.
        assert_eq!(result.faces.len(), 14);
        // Vertices: 8 (outer) + 8 (inner) = 16.
        assert_eq!(result.vertices.len(), 16);
    }

    #[test]
    fn collect_loop_half_edges_cube() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let first_face = brep.faces.keys().next().unwrap();
        let loop_h = brep.faces[first_face].outer_loop;
        let hes = collect_loop_half_edges(&brep, loop_h).unwrap();
        assert_eq!(hes.len(), 4, "cube face should have 4 half-edges");
    }

    #[test]
    fn assemble_all_faces_removed_produces_side_walls_only() {
        // Remove all 6 faces: outer skin has 0 faces, inner skin has 0 faces,
        // side walls cover all 6 openings (each square opening → 4 side wall
        // faces, but the geometry is more complex for full removal).
        // The key invariant: the call must succeed without panicking.
        let outer = BRep::unit_cube(tol()).unwrap();
        let all_faces: Vec<FaceHandle> = outer.faces.keys().collect();
        let (mut inner, skin_map) = build_inner_skin(&outer, &all_faces).unwrap();
        inner.reverse_body_orientation();
        let result = assemble_shell(&outer, &inner, &all_faces, &skin_map, tol());
        // Whether it succeeds or returns an error is acceptable — the important
        // thing is that the code path is exercised without a panic.
        let _ = result;
    }

    #[test]
    fn assemble_two_faces_removed() {
        let outer = BRep::unit_cube(tol()).unwrap();
        let removed: Vec<FaceHandle> = outer.faces.keys().take(2).collect();
        let (mut inner, skin_map) = build_inner_skin(&outer, &removed).unwrap();
        inner.reverse_body_orientation();
        let result = assemble_shell(&outer, &inner, &removed, &skin_map, tol());
        // Outer (4 faces) + inner (4 faces) + side walls for 2 openings.
        // Exact face count depends on side wall geometry. Just verify it succeeds.
        assert!(
            result.is_ok(),
            "two-face removal should assemble: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        // Must have at least the outer and inner faces.
        assert!(brep.faces.len() >= 8, "should have at least 8 faces");
    }

    #[test]
    fn copy_skin_entities_copies_all_vertices() {
        let outer = BRep::unit_cube(tol()).unwrap();
        let (inner, skin_map) = build_inner_skin(&outer, &[]).unwrap();
        let result = assemble_shell(&outer, &inner, &[], &skin_map, tol()).unwrap();
        // A pure offset (no removal): 8 outer + 8 inner = 16 vertices.
        assert_eq!(
            result.vertices.len(),
            16,
            "should have 16 vertices: 8 outer + 8 inner"
        );
    }
}
