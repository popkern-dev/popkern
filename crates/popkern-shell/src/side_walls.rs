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

use crate::error::ShellError;
use crate::skin_map::SkinMap;
use popkern_core::constructors::BRepBuilder;
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{
    BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, ShellHandle, VertexHandle,
};

/// Result of side wall construction.
#[derive(Debug, Clone)]
pub struct SideWallResult {
    /// New faces created as side walls.
    pub faces: Vec<FaceHandle>,
}

/// One vertical (thickness-direction) edge at a boundary vertex,
/// shared by the two walls adjacent along the opening's rim ring.
struct VerticalEdge {
    /// Half-edge running outer → inner.
    down: HalfEdgeHandle,
    /// Half-edge running inner → outer.
    up: HalfEdgeHandle,
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
    // One shared vertical edge per boundary vertex (builder outer handle).
    let mut verticals: HashMap<VertexHandle, VerticalEdge> = HashMap::new();
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
                shell_handle,
            )?;
            faces.push(face);
        }
    }

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
    let mut verticals: HashMap<VertexHandle, VerticalEdge> = HashMap::new();
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
            shell_handle,
        )?;
        faces.push(face);
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
    let start_he = outer.loops[face.outer_loop].half_edge;
    let mut result = Vec::new();
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

    Ok(result)
}

/// The shared vertical edge at one boundary vertex, created on first
/// use (with a straight-line curve) and reused by the adjacent wall.
fn vertical_at(
    builder: &mut BRepBuilder,
    verticals: &mut HashMap<VertexHandle, VerticalEdge>,
    builder_outer_v: VertexHandle,
    builder_inner_v: VertexHandle,
) -> (HalfEdgeHandle, HalfEdgeHandle) {
    if let Some(v) = verticals.get(&builder_outer_v) {
        return (v.down, v.up);
    }
    let p_out = builder.brep().vertices[builder_outer_v].position;
    let p_in = builder.brep().vertices[builder_inner_v].position;
    let (edge, down, up) = builder.add_edge(builder_outer_v, builder_inner_v);
    set_line_curve(builder, edge, p_out, p_in);
    verticals.insert(builder_outer_v, VerticalEdge { down, up });
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
    verticals: &mut HashMap<VertexHandle, VerticalEdge>,
    shell_handle: ShellHandle,
) -> Result<FaceHandle, ShellError> {
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
    let iv0_orig = *skin_map
        .vertices
        .get(&v0_orig)
        .ok_or_else(|| map_err("skin vertex v0"))?;
    let iv1_orig = *skin_map
        .vertices
        .get(&v1_orig)
        .ok_or_else(|| map_err("skin vertex v1"))?;
    let b_iv0 = *inner_vertex_map
        .get(&iv0_orig)
        .ok_or_else(|| map_err("inner vertex v0"))?;
    let b_iv1 = *inner_vertex_map
        .get(&iv1_orig)
        .ok_or_else(|| map_err("inner vertex v1"))?;

    // The free rim half-edges to adopt.
    let b_outer_rim = *outer_he_map
        .get(&rim_he)
        .ok_or_else(|| map_err("outer rim half-edge"))?;
    let inner_rim_orig = match skin_map.half_edges.get(&rim_he) {
        Some(h) => *h,
        None if b_iv0 == b_iv1 => {
            // Rim-corner arc (PK-326 boundary corner): no partner rim
            // half-edge exists — both arc endpoints resolve to the ONE
            // inner corner vertex and the wall is a planar sector
            // [arc, vertical down, vertical up].
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
            return Ok(face);
        }
        None => return Err(map_err("skin rim half-edge")),
    };
    let b_inner_rim = *inner_he_map
        .get(&inner_rim_orig)
        .ok_or_else(|| map_err("inner rim half-edge"))?;

    // Direction contract: outer rim runs v0 → v1; the inner image runs
    // v1' → v0' (the inner skin's winding was mirrored honestly).
    expect_origin(builder, b_outer_rim, b_v0, "outer rim")?;
    expect_origin(builder, b_inner_rim, b_iv1, "inner rim")?;

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

    Ok(face)
}

/// Author the wall face's surface.
///
/// With a removed face the wall lies IN that face's original surface
/// (see the block comment below). At a FREE rim (open-sheet thicken,
/// PK-324) there is no removed face; the wall is the ruled strip
/// between the two rim curves, authored as a verified plane — the
/// planar-rim case (sheet bodies with straight rims). A non-planar
/// free rim declines with the edge named: authoring general ruled
/// NURBS walls is recorded Stage-3 work (PK-325), and serving a wrong
/// plane there would be a lie.
fn set_wall_surface(
    builder: &mut BRepBuilder,
    outer: &BRep,
    rim_he: HalfEdgeHandle,
    removed_face: Option<FaceHandle>,
    face: FaceHandle,
    [p0, p1, p2, p3]: [Point3; 4],
) -> Result<(), ShellError> {
    // The wall lies IN the removed face's original surface: the outer
    // rim is that face's own boundary, and the inner rim was recomputed
    // as (kept offset × removed ORIGINAL surface) — flush by the rim
    // contract (PK-307). Inheriting the surface covers every rim shape;
    // the former quad-plane fit degenerated on CLOSED rims (v0 == v1
    // collapses the quad) and annulus walls served surface-less —
    // invisible to the structural gate, but tessellation skipped them
    // and every shelled round opening meshed as a hole.
    //
    // Orientation: the wall loop adopts the removed face's own rim
    // half-edge, so its traversal along the rim matches — but the strip
    // lies on the removed face's INTERIOR side of the rim for inward
    // shells (inset rims) and on its EXTERIOR side for outward ones
    // (outset rims), which mirrors the winding. The side is measured:
    // interior-left of the rim traversal under the flag-folded normal
    // vs the actual offset direction (outer rim → skin rim, in-plane).
    if let Some(sof) = removed_face
        .and_then(|fh| outer.faces.get(fh))
        .and_then(|f| f.surface.clone())
    {
        use popkern_core::geometry::surface::SurfaceEval;
        let strip_on_exterior_side = (|| -> Option<bool> {
            let pc = outer.half_edges.get(rim_he)?.pcurve.as_ref()?;
            let (t0, t1) = (pc.t_start, pc.t_end);
            let uv0 = pc.point_at(t0);
            let uv1 = pc.point_at(t0 + (t1 - t0) * 1e-3);
            let q0 = sof.geometry.point_at(uv0.x, uv0.y);
            let q1 = sof.geometry.point_at(uv1.x, uv1.y);
            let rim_dir = (q1 - q0).normalized()?;
            let folded_normal = sof.geometry.normal_at(uv0.x, uv0.y).normalized()?
                * if sof.orientation { 1.0 } else { -1.0 };
            let interior_left = folded_normal.cross(&rim_dir);
            let offset_dir = (p3 - p0).normalized()?;
            Some(interior_left.dot(&offset_dir) < 0.0)
        })()
        .unwrap_or(false);
        let mut wall_sof = sof;
        if strip_on_exterior_side {
            wall_sof.orientation = !wall_sof.orientation;
        }
        builder.set_face_surface(face, wall_sof);
        return Ok(());
    }
    let Some(surface) = plane_from_quad(p0, p1, p2, p3) else {
        return Err(ShellError::Topology(format!(
            "free-rim wall quad for half-edge {rim_he:?} is degenerate"
        )));
    };
    // Verification: the fourth corner and both rim curve midpoints
    // must lie ON the fitted plane, else the plane is a lie for this
    // rim (curved rims need ruled walls — Stage 3, PK-325).
    let SurfaceGeometry::Plane { origin, normal, .. } = &surface else {
        unreachable!("plane_from_quad returns Plane");
    };
    let plane_dist = |q: &Point3| ((*q - *origin).dot(normal)).abs();
    let lin = builder.brep().tolerance.linear();
    let mut witnesses = vec![p3];
    if let Some(edge) = outer.half_edges.get(rim_he).map(|he| he.edge)
        && let Some(curve) = outer.edges.get(edge).and_then(|e| e.curve.as_ref())
    {
        use popkern_core::geometry::curve::CurveEval;
        witnesses.push(
            curve
                .geometry
                .point_at(curve.t_start + (curve.t_end - curve.t_start) * 0.5),
        );
    }
    if let Some(bad) = witnesses.iter().find(|q| plane_dist(q) > lin) {
        return Err(ShellError::Topology(format!(
            "free-rim wall for half-edge {rim_he:?} is non-planar (witness off-plane by              {:.3e}); ruled free-rim walls are Stage-3 work (PK-325)",
            plane_dist(bad)
        )));
    }
    builder.set_face_surface(
        face,
        SurfaceOnFace {
            geometry: surface,
            orientation: true,
        },
    );

    Ok(())
}

/// Set a line curve on an edge.
fn set_line_curve(builder: &mut BRepBuilder, edge: EdgeHandle, start: Point3, end: Point3) {
    let direction = end - start;
    builder.set_edge_curve(
        edge,
        CurveOnEdge {
            geometry: CurveGeometry::Line {
                origin: start,
                direction,
            },
            t_start: 0.0,
            t_end: 1.0,
        },
    );
}

/// Compute a planar surface from 4 quad vertices.
///
/// Uses the first 3 vertices to define the plane. Returns `None` if the
/// points are collinear.
fn plane_from_quad(p0: Point3, p1: Point3, p2: Point3, _p3: Point3) -> Option<SurfaceGeometry> {
    let raw_u = p1 - p0;
    let raw_v = p2 - p0;
    let normal = raw_u.cross(&raw_v);
    let normal = normal.normalized()?;
    let u_axis = raw_u.normalized()?;
    let v_proj = raw_v - u_axis * raw_v.dot(&u_axis);
    let v_axis = v_proj.normalized()?;
    Some(SurfaceGeometry::Plane {
        origin: p0,
        normal,
        u_axis,
        v_axis,
    })
}

#[cfg(test)]
#[path = "side_walls_tests.rs"]
mod tests;
