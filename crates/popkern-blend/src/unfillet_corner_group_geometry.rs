//! Geometry and small topology helpers for the recordless
//! corner-group unfillet (see `unfillet_corner_group`).

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::error::BlendError;
use crate::unfillet_helpers::{PlaneData, intersect_line_plane, intersect_planes};

/// The two plane neighbours a cylinder strip is TANGENT to: normal
/// perpendicular to the fillet axis, plane at axis distance equal to
/// the fillet radius (relative band 10%). Far end faces have normals
/// along the axis and the cap is not a plane, so both are excluded
/// by construction.
pub(crate) fn cylinder_support_pair(
    brep: &BRep,
    strip: FaceHandle,
) -> Option<(FaceHandle, FaceHandle)> {
    let Some(SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    }) = surface_kind_of(brep, strip).cloned()
    else {
        return None;
    };
    let axis_n = axis.normalized()?;
    let mut found: Vec<FaceHandle> = Vec::new();
    for f in crate::identify::blend_neighbours(brep, strip) {
        let Some(SurfaceGeometry::Plane {
            origin: po, normal, ..
        }) = surface_kind_of(brep, f)
        else {
            continue;
        };
        if normal.dot(&axis_n).abs() > 0.1 {
            continue;
        }
        let axis_to_plane = ((origin - *po).dot(normal)).abs();
        if (axis_to_plane - radius).abs() <= radius * 0.1 && !found.contains(&f) {
            found.push(f);
        }
    }
    match found.as_slice() {
        [a, b] => Some((*a, *b)),
        _ => None,
    }
}

pub(crate) fn face_diameter(brep: &BRep, face: FaceHandle) -> f64 {
    let positions = brep.face_outer_loop_positions(face);
    let mut d = 0.0_f64;
    for (i, p) in positions.iter().enumerate() {
        for q in positions.iter().skip(i + 1) {
            d = d.max(p.distance_to(q));
        }
    }
    d
}

pub(crate) fn intersect_three_planes(planes: &[PlaneData]) -> Option<Point3> {
    let [a, b, c] = planes else { return None };
    let (o, d) = intersect_planes(a, b).ok()?;
    intersect_line_plane(&o, &d, c)
}

pub(crate) fn plane_deviation(g: &SurfaceGeometry, positions: &[Point3]) -> f64 {
    let SurfaceGeometry::Plane { origin, normal, .. } = g else {
        return f64::INFINITY;
    };
    positions
        .iter()
        .map(|p| ((*p - *origin).dot(normal)).abs())
        .fold(0.0, f64::max)
}

/// Distinct vertices on any loop of `face` within `gate` of `target`.
pub(crate) fn loop_vertices_near(
    brep: &BRep,
    face: FaceHandle,
    target: &Point3,
    gate: f64,
) -> Vec<VertexHandle> {
    let mut found: Vec<VertexHandle> = Vec::new();
    let Some(fd) = brep.faces.get(face) else {
        return found;
    };
    let mut loops = vec![fd.outer_loop];
    loops.extend(fd.inner_loops.iter().copied());
    for lh in loops {
        let Some(ld) = brep.loops.get(lh) else {
            continue;
        };
        let start = ld.half_edge;
        let mut cur = start;
        for _ in 0..brep.half_edges.len() + 1 {
            let Some(he) = brep.half_edges.get(cur) else {
                break;
            };
            let v = he.origin;
            if let Some(vd) = brep.vertices.get(v)
                && vd.position.distance_to(target) <= gate
                && !found.contains(&v)
            {
                found.push(v);
            }
            cur = he.next;
            if cur == start {
                break;
            }
        }
    }
    found
}

pub(crate) fn weld_cluster(
    out: &mut BRep,
    cluster: &[VertexHandle],
    target: Point3,
) -> Result<VertexHandle, BlendError> {
    let primary = cluster[0];
    out.vertices[primary].position = target;
    for &v in &cluster[1..] {
        out.vertices[v].position = target;
        out.weld_vertex_pair(primary, v)
            .map_err(BlendError::Topology)?;
    }
    Ok(primary)
}

/// A restored edge already exists between the two vertices when the
/// surviving slit connects them — then no mef is needed for that
/// sector boundary.
pub(crate) fn needs_mef(brep: &BRep, a: VertexHandle, b: VertexHandle) -> bool {
    !brep.edges.iter().any(|(_, ed)| {
        let Some(he) = brep.half_edges.get(ed.half_edge) else {
            return false;
        };
        let Some(tw) = brep.half_edges.get(he.twin) else {
            return false;
        };
        (he.origin == a && tw.origin == b) || (he.origin == b && tw.origin == a)
    })
}

/// Keep at most one corner-to-far spur before the sector split (PK-608).
///
/// The group merge leaves a same-face spur along a strip contact edge
/// per independent cycle of the group's face adjacency. The uniform
/// sphere corner (a 6-ring plus one cap) leaves one, from the corner
/// to one far vertex; the two sector mefs to the other far vertices
/// re-legitimise it as the third restored edge. A setback corner's cap
/// is several sector faces and leaves more: measured on the cube
/// 0.08/0.1/0.12 KP20 fan and the 0.05/0.2/0.3 fitted cap, spurs to two
/// far vertices (one of them doubled), so [`needs_mef`] skipped two
/// splits and the restore served a face folded over two support planes.
/// Every spur but the first is zero-area (`remove_slit_loop_edges`), so
/// dropping them changes no geometry; the dropped far vertices then
/// take their mefs.
pub(crate) fn drop_extra_corner_spurs(brep: &mut BRep, corner: VertexHandle, far: &[VertexHandle]) {
    let view: &BRep = brep;
    let spurs: Vec<EdgeHandle> = far
        .iter()
        .flat_map(|&f| {
            view.edges
                .iter()
                .filter(move |(_, ed)| {
                    let Some(he) = view.half_edges.get(ed.half_edge) else {
                        return false;
                    };
                    let Some(tw) = view.half_edges.get(he.twin) else {
                        return false;
                    };
                    let joins = (he.origin == corner && tw.origin == f)
                        || (he.origin == f && tw.origin == corner);
                    joins && he.loop_ == tw.loop_ && (he.next == he.twin || tw.next == ed.half_edge)
                })
                .map(|(eh, _)| eh)
                .collect::<Vec<_>>()
        })
        .collect();
    let extra: std::collections::BTreeSet<EdgeHandle> = spurs.into_iter().skip(1).collect();
    if !extra.is_empty() {
        crate::unfillet_surgery::remove_slit_loop_edges_where(brep, |eh| extra.contains(&eh));
    }
}

/// Pick a same-loop half-edge pair (origin `corner`, origin `far`)
/// for the sector mef. Tries every corner-origin candidate against
/// every far-origin candidate and returns the first same-loop pair.
pub(crate) fn pick_mef_pair(
    brep: &BRep,
    corner: VertexHandle,
    far: VertexHandle,
) -> Option<(HalfEdgeHandle, HalfEdgeHandle)> {
    let at: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, he)| he.origin == corner)
        .map(|(h, _)| h)
        .collect();
    let bt: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, he)| he.origin == far)
        .map(|(h, _)| h)
        .collect();
    for &ha in &at {
        for &hb in &bt {
            if brep.half_edges[ha].loop_ == brep.half_edges[hb].loop_ && ha != hb {
                return Some((ha, hb));
            }
        }
    }
    None
}

pub(crate) fn edge_between(brep: &BRep, a: VertexHandle, b: VertexHandle) -> Option<EdgeHandle> {
    brep.edges.iter().find_map(|(eh, ed)| {
        let he = brep.half_edges.get(ed.half_edge)?;
        let tw = brep.half_edges.get(he.twin)?;
        ((he.origin == a && tw.origin == b) || (he.origin == b && tw.origin == a)).then_some(eh)
    })
}

pub(crate) fn invalid(parameter: &str, value: String, reason: String) -> BlendError {
    BlendError::InvalidParameter {
        parameter: parameter.to_string(),
        value,
        reason,
    }
}

pub(crate) fn is_cap_surface(g: &SurfaceGeometry) -> bool {
    matches!(
        g,
        SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Nurbs(_) | SurfaceGeometry::Gregory(_)
    )
}

pub(crate) fn is_strip_surface(g: &SurfaceGeometry) -> bool {
    matches!(
        g,
        SurfaceGeometry::Cylinder { .. } | SurfaceGeometry::Torus { .. }
    )
}

pub(crate) fn surface_kind_of(brep: &BRep, f: FaceHandle) -> Option<&SurfaceGeometry> {
    brep.faces
        .get(f)
        .and_then(|d| d.surface.as_ref())
        .map(|s| &s.geometry)
}
