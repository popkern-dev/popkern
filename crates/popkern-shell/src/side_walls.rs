//! Side wall construction at shell openings.
//!
//! When faces are removed from a shell, the gap between the outer and
//! inner skins must be closed with side wall faces. This module
//! identifies boundary edges (edges where a removed face meets a kept
//! face) and builds quadrilateral side wall faces connecting the two
//! skins.
//!
//! Topology contract (PK-307): a wall ADOPTS existing rim topology
//! instead of fabricating parallel edges. Each boundary edge already
//! carries two free half-edges in the assembled builder — the removed
//! face's side on the outer skin and its image on the inner skin
//! (whose face was likewise never copied) — and the wall's loop claims
//! exactly those two, plus ONE vertical edge per boundary vertex
//! shared with the adjacent wall. Before PK-307 each wall added four
//! fresh edges, which left the rim half-edges loopless (structural
//! `BrokenCycle`), duplicated every vertical, and attached the walls
//! to the skins at vertices only — the open-shell output was never
//! edge-manifold.

use std::collections::{HashMap, HashSet};

#[path = "side_walls_support.rs"]
mod side_walls_support;
use side_walls_support::{build_corner_sectors, set_line_curve, set_wall_surface};

use crate::error::ShellError;
use crate::skin_map::SkinMap;
use popkern_core::constructors::BRepBuilder;
use popkern_core::topology::{
    BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, ShellHandle, VertexHandle,
};

/// Result of side wall construction.
#[derive(Debug, Clone)]
pub struct SideWallResult {
    /// New faces created as side walls.
    pub faces: Vec<FaceHandle>,
}

/// One vertical (thickness-direction) edge at a boundary vertex pair,
/// shared by the walls adjacent along the opening's rim ring. Keyed by
/// (builder outer vertex, builder inner vertex): rim-corner arc splits
/// (PK-326/PK-340) give one outer vertex TWO distinct inner partners.
pub(super) struct VerticalEdge {
    /// Half-edge running outer → inner.
    pub(super) down: HalfEdgeHandle,
    /// Half-edge running inner → outer.
    pub(super) up: HalfEdgeHandle,
}

/// Build side wall faces connecting outer and inner skins at shell openings.
///
/// For each removed face, walks its outer loop edges. An edge whose
/// opposite face is NOT in the removed set is a boundary edge; each
/// yields one wall quad whose loop is
/// `outer rim (removed side) → vertical down → inner rim (free side) →
/// vertical up`.
///
/// # Errors
///
/// Returns [`ShellError::Topology`] if handle lookups fail or a rim
/// half-edge does not run the direction the topology promises (the
/// error names the half-edge and the expected/actual vertices).
#[allow(clippy::too_many_arguments)]
pub fn build_side_walls(
    builder: &mut BRepBuilder,
    outer: &BRep,
    removed_faces: &[FaceHandle],
    skin_map: &SkinMap,
    outer_vertex_map: &HashMap<VertexHandle, VertexHandle>,
    inner_vertex_map: &HashMap<VertexHandle, VertexHandle>,
    outer_he_map: &HashMap<HalfEdgeHandle, HalfEdgeHandle>,
    inner_he_map: &HashMap<HalfEdgeHandle, HalfEdgeHandle>,
    shell_handle: ShellHandle,
) -> Result<SideWallResult, ShellError> {
    let removed_set: HashSet<FaceHandle> = removed_faces.iter().copied().collect();
    let mut faces = Vec::new();
    let mut verticals: HashMap<(VertexHandle, VertexHandle), VerticalEdge> = HashMap::new();
    let mut vertical_faces: HashMap<(VertexHandle, VertexHandle), FaceHandle> = HashMap::new();
    // Each boundary edge yields exactly one wall (it is only ever seen
    // from its single removed side); the set guards adjacent removed
    // faces re-listing a shared rim.
    let mut processed_edges: HashSet<EdgeHandle> = HashSet::new();

    for &removed_fh in removed_faces {
        for (outer_eh, rim_he) in find_boundary_edges(outer, removed_fh, &removed_set)? {
            if !processed_edges.insert(outer_eh) {
                continue;
            }
            let face = build_one_side_wall(
                builder,
                outer,
                rim_he,
                skin_map,
                outer_vertex_map,
                inner_vertex_map,
                outer_he_map,
                inner_he_map,
                &mut verticals,
                &mut vertical_faces,
                shell_handle,
            )?;
            if let Some(face) = face {
                faces.push(face);
            }
        }
    }

    // PK-340: rim-corner arcs on the INNER skin have no original
    // preimage, so the edge iteration above never visits them. They
    // are the remaining loopless half-edges whose two endpoints both
    // carry verticals under ONE shared outer vertex — close each with
    // a planar sector [arc, vertical up, vertical down] lying in the
    // removed face's surface (adopted from an adjacent wall).
    faces.extend(build_corner_sectors(
        builder,
        &verticals,
        &vertical_faces,
        shell_handle,
    )?);

    Ok(SideWallResult { faces })
}

/// Build wall faces closing the FREE rims of an open sheet's thicken
/// (design §6 open-sheet row, PK-324): every source edge with a
/// loopless half-edge is a rim; the wall adopts that loopless
/// half-edge's builder image on the outer skin and its (mirrored)
/// image on the inner skin, exactly like an opening wall — the free
/// rim IS an opening whose "removed face" never existed.
///
/// # Errors
///
/// Returns [`ShellError::Topology`] when a rim wall would be
/// non-planar (ruled free-rim walls are Stage-3 work) or a handle
/// contract is violated.
#[allow(clippy::too_many_arguments)]
pub fn build_free_rim_walls(
    builder: &mut BRepBuilder,
    outer: &BRep,
    skin_map: &SkinMap,
    outer_vertex_map: &HashMap<VertexHandle, VertexHandle>,
    inner_vertex_map: &HashMap<VertexHandle, VertexHandle>,
    outer_he_map: &HashMap<HalfEdgeHandle, HalfEdgeHandle>,
    inner_he_map: &HashMap<HalfEdgeHandle, HalfEdgeHandle>,
    shell_handle: ShellHandle,
) -> Result<SideWallResult, ShellError> {
    let mut faces = Vec::new();
    let mut verticals: HashMap<(VertexHandle, VertexHandle), VerticalEdge> = HashMap::new();
    let mut vertical_faces: HashMap<(VertexHandle, VertexHandle), FaceHandle> = HashMap::new();
    let mut rim_hes: Vec<HalfEdgeHandle> = Vec::new();
    for (eh, edge) in &outer.edges {
        let he = &outer.half_edges[edge.half_edge];
        let he_loopless = !outer.loops.contains_key(he.loop_);
        let twin_loopless = !outer.loops.contains_key(outer.half_edges[he.twin].loop_);
        match (he_loopless, twin_loopless) {
            (true, false) => rim_hes.push(edge.half_edge),
            (false, true) => rim_hes.push(he.twin),
            (true, true) => {
                return Err(ShellError::Topology(format!(
                    "edge {eh:?} has no live face on either side — not a sheet rim"
                )));
            }
            (false, false) => {}
        }
    }
    // Deterministic wall order (handle-sorted; the map iteration above
    // is insertion-ordered per slotmap but sort anyway for the PK-277
    // contract).
    rim_hes.sort();
    for rim_he in rim_hes {
        let face = build_one_side_wall(
            builder,
            outer,
            rim_he,
            skin_map,
            outer_vertex_map,
            inner_vertex_map,
            outer_he_map,
            inner_he_map,
            &mut verticals,
            &mut vertical_faces,
            shell_handle,
        )?;
        if let Some(face) = face {
            faces.push(face);
        }
    }
    Ok(SideWallResult { faces })
}

/// Find boundary edges of a removed face.
///
/// Returns `(EdgeHandle, HalfEdgeHandle)` pairs where the half-edge belongs
/// to the removed face and the twin's face is NOT removed.
fn find_boundary_edges(
    outer: &BRep,
    removed_fh: FaceHandle,
    removed_set: &HashSet<FaceHandle>,
) -> Result<Vec<(EdgeHandle, HalfEdgeHandle)>, ShellError> {
    let face = &outer.faces[removed_fh];
    // ALL loops: a removed face pierced by a bore carries the bore's
    // rim on an INNER loop — that rim is an opening boundary exactly
    // like the outer one, and skipping it left its half-edges loopless
    // (structural BrokenCycle, PK-346).
    let loop_handles: Vec<_> = std::iter::once(face.outer_loop)
        .chain(face.inner_loops.iter().copied())
        .collect();
    let mut result = Vec::new();
    for lh in loop_handles {
        let Some(loop_data) = outer.loops.get(lh) else {
            continue;
        };
        let start_he = loop_data.half_edge;
        if !outer.half_edges.contains_key(start_he) {
            continue;
        }
        let mut current = start_he;
        loop {
            let he = &outer.half_edges[current];
            let twin = &outer.half_edges[he.twin];
            let twin_loop = twin.loop_;

            if outer.loops.contains_key(twin_loop) {
                let twin_face = outer.loops[twin_loop].face;
                if !removed_set.contains(&twin_face) {
                    result.push((he.edge, current));
                }
            }

            current = he.next;
            if current == start_he {
                break;
            }
        }
    }

    Ok(result)
}

/// The shared vertical edge at one boundary vertex, created on first
/// use (with a straight-line curve) and reused by the adjacent wall.
fn vertical_at(
    builder: &mut BRepBuilder,
    verticals: &mut HashMap<(VertexHandle, VertexHandle), VerticalEdge>,
    builder_outer_v: VertexHandle,
    builder_inner_v: VertexHandle,
) -> (HalfEdgeHandle, HalfEdgeHandle) {
    if let Some(v) = verticals.get(&(builder_outer_v, builder_inner_v)) {
        return (v.down, v.up);
    }
    let p_out = builder.brep().vertices[builder_outer_v].position;
    let p_in = builder.brep().vertices[builder_inner_v].position;
    let (edge, down, up) = builder.add_edge(builder_outer_v, builder_inner_v);
    set_line_curve(builder, edge, p_out, p_in);
    verticals.insert(
        (builder_outer_v, builder_inner_v),
        VerticalEdge { down, up },
    );
    (down, up)
}

/// Expect `he`'s origin in the builder to be `expected`; a mismatch is
/// a structural error naming the site (never silently mis-wired).
fn expect_origin(
    builder: &BRepBuilder,
    he: HalfEdgeHandle,
    expected: VertexHandle,
    role: &str,
) -> Result<(), ShellError> {
    let actual = builder.brep().half_edges[he].origin;
    if actual != expected {
        return Err(ShellError::Topology(format!(
            "side wall {role} half-edge {he:?} runs from {actual:?}, expected origin \
             {expected:?} — rim adoption direction contract violated"
        )));
    }
    Ok(())
}

/// Build a single side wall quad face by adopting the free rim
/// half-edges and the shared verticals.
///
/// `rim_he` is the boundary half-edge on the ORIGINAL outer body,
/// belonging to the removed face and running `v0 → v1` in that face's
/// loop. Its builder image is loopless (the removed face was never
/// copied), as is the builder image of its inner-skin counterpart —
/// which runs `v1' → v0'` because the inner skin's orientation was
/// reversed. The wall loop is therefore
/// `v0 → v1 → v1' → v0'` using existing rim half-edges on both skins.
#[allow(clippy::too_many_arguments)]
fn build_one_side_wall(
    builder: &mut BRepBuilder,
    outer: &BRep,
    rim_he: HalfEdgeHandle,
    skin_map: &SkinMap,
    outer_vertex_map: &HashMap<VertexHandle, VertexHandle>,
    inner_vertex_map: &HashMap<VertexHandle, VertexHandle>,
    outer_he_map: &HashMap<HalfEdgeHandle, HalfEdgeHandle>,
    inner_he_map: &HashMap<HalfEdgeHandle, HalfEdgeHandle>,
    verticals: &mut HashMap<(VertexHandle, VertexHandle), VerticalEdge>,
    vertical_faces: &mut HashMap<(VertexHandle, VertexHandle), FaceHandle>,
    shell_handle: ShellHandle,
) -> Result<Option<FaceHandle>, ShellError> {
    let he = &outer.half_edges[rim_he];
    let twin = &outer.half_edges[he.twin];
    let (v0_orig, v1_orig) = (he.origin, twin.origin);

    let map_err = |what: &str| ShellError::Topology(format!("side wall: {what} not in map"));

    // Builder-space skin vertices.
    let b_v0 = *outer_vertex_map
        .get(&v0_orig)
        .ok_or_else(|| map_err("outer vertex v0"))?;
    let b_v1 = *outer_vertex_map
        .get(&v1_orig)
        .ok_or_else(|| map_err("outer vertex v1"))?;
    // Mapped inner vertices are the UNSPLIT expectation only: a rim
    // corner peeled by PK-326/PK-340 ends a rim edge at a new vertex
    // with no map entry — the quad path reads the ACTUAL endpoints of
    // the adopted inner half-edge instead.
    let b_iv0 = skin_map
        .vertices
        .get(&v0_orig)
        .and_then(|iv| inner_vertex_map.get(iv))
        .copied();
    let b_iv1 = skin_map
        .vertices
        .get(&v1_orig)
        .and_then(|iv| inner_vertex_map.get(iv))
        .copied();

    // The free rim half-edges to adopt.
    let b_outer_rim = *outer_he_map
        .get(&rim_he)
        .ok_or_else(|| map_err("outer rim half-edge"))?;
    let inner_rim_orig = match skin_map.half_edges.get(&rim_he) {
        Some(h) => *h,
        None if b_iv0.is_some() && b_iv0 == b_iv1 => {
            // Rim-corner arc (PK-326 boundary corner): no partner rim
            // half-edge exists — both arc endpoints resolve to the ONE
            // inner corner vertex and the wall is a planar sector
            // [arc, vertical down, vertical up].
            let b_iv0 = b_iv0.ok_or_else(|| map_err("inner vertex v0"))?;
            let b_iv1 = b_iv1.ok_or_else(|| map_err("inner vertex v1"))?;
            expect_origin(builder, b_outer_rim, b_v0, "outer rim arc")?;
            let (down_v1, _) = vertical_at(builder, verticals, b_v1, b_iv1);
            let (_, up_v0) = vertical_at(builder, verticals, b_v0, b_iv0);
            let p0 = builder.brep().vertices[b_v0].position;
            let p1 = builder.brep().vertices[b_v1].position;
            let pc = builder.brep().vertices[b_iv0].position;
            let face = builder.add_face(shell_handle);
            builder
                .add_loop(&[b_outer_rim, down_v1, up_v0], face)
                .map_err(ShellError::Core)?;
            set_wall_surface(builder, outer, rim_he, None, face, [p0, p1, pc, pc])?;
            return Ok(Some(face));
        }
        None if b_iv0.is_none() || b_iv1.is_none() => {
            // Outer-side rim-corner arc (midplane openings): the arc
            // rides in the removed face's loop with split endpoints
            // that have no partner-skin image. Its planar sector is
            // built by the post-pass once the adjacent quads have
            // created both verticals.
            return Ok(None);
        }
        None => {
            return Err(ShellError::Topology(format!(
                "side wall: skin rim half-edge not in map (rim_he {rim_he:?}, v0                  {v0_orig:?} mapped {b_iv0:?}, v1 {v1_orig:?} mapped {b_iv1:?})"
            )));
        }
    };
    let b_inner_rim = *inner_he_map
        .get(&inner_rim_orig)
        .ok_or_else(|| map_err("inner rim half-edge"))?;

    // Direction contract: outer rim runs v0 → v1; the inner image runs
    // v1' → v0' (the inner skin's winding was mirrored honestly). The
    // inner ENDPOINTS are read from the adopted half-edge itself
    // (PK-340): a rim-corner arc split re-homes an endpoint to a new
    // vertex with no original preimage, so the mapped vertices are
    // only the unsplit expectation — a mapped-but-reversed adoption
    // still errors.
    expect_origin(builder, b_outer_rim, b_v0, "outer rim")?;
    let actual_iv1 = builder.brep().half_edges[b_inner_rim].origin;
    let inner_twin = builder.brep().half_edges[b_inner_rim].twin;
    let actual_iv0 = builder.brep().half_edges[inner_twin].origin;
    if b_iv0.is_some() && b_iv0 != b_iv1 && Some(actual_iv1) == b_iv0 && Some(actual_iv0) == b_iv1 {
        return Err(ShellError::Topology(format!(
            "side wall inner rim half-edge {b_inner_rim:?} runs v0' → v1' — rim \
             adoption direction contract violated"
        )));
    }
    let (b_iv1, b_iv0) = (actual_iv1, actual_iv0);

    // Shared verticals: down (outer → inner) at v1, up (inner → outer)
    // at v0.
    let (down_v1, _) = vertical_at(builder, verticals, b_v1, b_iv1);
    let (_, up_v0) = vertical_at(builder, verticals, b_v0, b_iv0);

    // Positions for the wall plane, in loop order.
    let p0 = builder.brep().vertices[b_v0].position;
    let p1 = builder.brep().vertices[b_v1].position;
    let p2 = builder.brep().vertices[b_iv1].position;
    let p3 = builder.brep().vertices[b_iv0].position;

    let face = builder.add_face(shell_handle);
    builder
        .add_loop(&[b_outer_rim, down_v1, b_inner_rim, up_v0], face)
        .map_err(ShellError::Core)?;

    let removed_face = outer.loops.get(he.loop_).map(|l| l.face);
    set_wall_surface(builder, outer, rim_he, removed_face, face, [p0, p1, p2, p3])?;
    vertical_faces.insert((b_v1, b_iv1), face);
    vertical_faces.insert((b_v0, b_iv0), face);

    Ok(Some(face))
}

#[cfg(test)]
#[path = "side_walls_tests.rs"]
mod tests;
