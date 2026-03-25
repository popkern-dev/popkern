//! Shell assembly: combine outer skin, inner skin, and side walls.
//!
//! Copies entities from both the outer and inner BReps into a new
//! [`BRepBuilder`], then adds side wall faces at the openings where
//! removed faces create gaps.

use std::collections::HashMap;

use crate::error::ShellError;
use crate::side_walls::build_side_walls;
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
);

/// Assemble the final shell BRep from outer skin, inner skin, and side walls.
///
/// Copies all topology from the outer BRep (minus removed faces) and the
/// inner BRep into a single BRepBuilder, then builds side wall faces at
/// the openings.
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
    let mut builder = BRepBuilder::new(tolerance);
    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);

    let removed_set: std::collections::HashSet<FaceHandle> =
        removed_faces.iter().copied().collect();

    // Copy outer skin (minus removed faces).
    let (outer_vertex_map, outer_edge_map, outer_he_map) =
        copy_skin_entities(&mut builder, outer, shell, |fh| !removed_set.contains(&fh))?;

    // Copy inner skin (all faces).
    let (inner_vertex_map, _inner_edge_map, _inner_he_map) =
        copy_skin_entities(&mut builder, inner, shell, |_| true)?;

    // Build side walls if there are removed faces.
    if !removed_faces.is_empty() {
        build_side_walls(
            &mut builder,
            outer,
            removed_faces,
            skin_map,
            &outer_vertex_map,
            &inner_vertex_map,
            shell,
            &Tolerance::default(),
        )?;
    }

    let _ = outer_edge_map;
    let _ = outer_he_map;

    Ok(builder.build_unchecked())
}

/// Copy all entities from a BRep into the builder, filtering faces.
///
/// Returns maps from the source BRep handles to the new builder handles
/// for vertices, edges, and half-edges.
fn copy_skin_entities(
    builder: &mut BRepBuilder,
    source: &BRep,
    shell: popkern_core::topology::ShellHandle,
    face_filter: impl Fn(FaceHandle) -> bool,
) -> Result<SkinCopyMaps, ShellError> {
    // 1. Copy vertices.
    let mut vertex_map = HashMap::new();
    for (vh, vertex) in &source.vertices {
        let new_vh = builder.add_vertex(vertex.position);
        vertex_map.insert(vh, new_vh);
    }

    // 2. Copy edges (and their half-edge pairs).
    let mut edge_map = HashMap::new();
    let mut he_map = HashMap::new();
    for (eh, edge) in &source.edges {
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
    for (fh, face) in &source.faces {
        if !face_filter(fh) {
            continue;
        }

        let new_face = builder.add_face(shell);

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

    Ok((vertex_map, edge_map, he_map))
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
        let (inner, skin_map) = build_inner_skin(&outer, &[first_face]).unwrap();
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
        let (inner, skin_map) = build_inner_skin(&outer, &[first_face]).unwrap();
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
}
