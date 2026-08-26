//! Side-wall support: corner-sector construction for rim arcs
//! (PK-340) and wall surface authoring.

use std::collections::HashMap;

use crate::error::ShellError;
use popkern_core::constructors::BRepBuilder;
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{
    BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, ShellHandle, VertexHandle,
};

use super::VerticalEdge;

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
pub(super) fn set_wall_surface(
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
pub(super) fn set_line_curve(
    builder: &mut BRepBuilder,
    edge: EdgeHandle,
    start: Point3,
    end: Point3,
) {
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
pub(super) fn plane_from_quad(
    p0: Point3,
    p1: Point3,
    p2: Point3,
    _p3: Point3,
) -> Option<SurfaceGeometry> {
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

/// PK-340: close rim-corner arcs left on the inner skin after the
/// quad walls. An arc is a still-loopless builder half-edge whose two
/// endpoints both carry verticals under ONE shared outer vertex; its
/// sector [arc, vertical up, vertical down] is planar in the removed
/// face's surface, adopted from an adjacent wall.
pub(super) fn build_corner_sectors(
    builder: &mut BRepBuilder,
    verticals: &HashMap<(VertexHandle, VertexHandle), VerticalEdge>,
    vertical_faces: &HashMap<(VertexHandle, VertexHandle), FaceHandle>,
    shell_handle: ShellHandle,
) -> Result<Vec<FaceHandle>, ShellError> {
    use std::collections::HashMap as Map;
    let mut by_inner: Map<VertexHandle, Vec<VertexHandle>> = Map::new();
    let mut by_outer: Map<VertexHandle, Vec<VertexHandle>> = Map::new();
    for &(ov, iv) in verticals.keys() {
        by_inner.entry(iv).or_default().push(ov);
        by_outer.entry(ov).or_default().push(iv);
    }
    let mut candidates: Vec<HalfEdgeHandle> = builder
        .brep()
        .half_edges
        .iter()
        .filter(|(_, he)| !builder.brep().loops.contains_key(he.loop_))
        .map(|(h, _)| h)
        .collect();
    candidates.sort();
    let mut faces = Vec::new();
    for he in candidates {
        if builder
            .brep()
            .loops
            .contains_key(builder.brep().half_edges[he].loop_)
        {
            continue; // consumed by an earlier sector this pass
        }
        let x = builder.brep().half_edges[he].origin;
        let twin = builder.brep().half_edges[he].twin;
        if !builder.brep().half_edges.contains_key(twin) {
            continue;
        }
        let y = builder.brep().half_edges[twin].origin;
        // INNER-side arc (classic shell): endpoints are the inner ends
        // of two verticals under one outer vertex. Loop:
        // arc (x → y), up at y (y → ov), down at x (ov → x).
        let inner_side = (|| {
            let (ovs_x, ovs_y) = (by_inner.get(&x)?, by_inner.get(&y)?);
            let ov = *ovs_x.iter().find(|o| ovs_y.contains(o))?;
            let vx = verticals.get(&(ov, x))?;
            let vy = verticals.get(&(ov, y))?;
            Some(([he, vy.up, vx.down], (ov, x), (ov, y)))
        })();
        // OUTER-side arc (midplane outer skin): endpoints are the
        // outer ends of two verticals under one inner vertex. Loop:
        // arc (x → y), down at y (y → iv), up at x (iv → x).
        let outer_side = (|| {
            let (ivs_x, ivs_y) = (by_outer.get(&x)?, by_outer.get(&y)?);
            let iv = *ivs_x.iter().find(|i| ivs_y.contains(i))?;
            let vx = verticals.get(&(x, iv))?;
            let vy = verticals.get(&(y, iv))?;
            Some(([he, vy.down, vx.up], (x, iv), (y, iv)))
        })();
        let Some((loop_hes, key_a, key_b)) = inner_side.or(outer_side) else {
            continue;
        };
        let face = builder.add_face(shell_handle);
        builder
            .add_loop(&loop_hes, face)
            .map_err(ShellError::Core)?;
        let sof = vertical_faces
            .get(&key_a)
            .or_else(|| vertical_faces.get(&key_b))
            .and_then(|wf| builder.brep().faces.get(*wf))
            .and_then(|f| f.surface.clone())
            .ok_or_else(|| {
                ShellError::Topology(
                    "rim-corner sector has no adjacent wall surface to adopt".into(),
                )
            })?;
        builder.set_face_surface(face, sof);
        faces.push(face);
    }
    Ok(faces)
}
