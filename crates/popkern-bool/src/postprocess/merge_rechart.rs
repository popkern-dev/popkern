//! Re-chart the pcurves a face inherits when a KEF merges two
//! co-surfacial faces (PK-806).

use popkern_core::topology::BRep;

/// The two faces an edge separates, with each face's half-edges and
/// surface, captured before a KEF merges them.
pub(crate) type MergeSide = (
    popkern_core::topology::FaceHandle,
    Vec<popkern_core::topology::HalfEdgeHandle>,
    Option<popkern_core::geometry::SurfaceGeometry>,
);

pub(crate) fn merge_sides(
    brep: &BRep,
    he: popkern_core::topology::HalfEdgeHandle,
) -> Vec<MergeSide> {
    let face_of = |h| {
        brep.half_edges
            .get(h)
            .and_then(|x| brep.loops.get(x.loop_))
            .map(|l| l.face)
    };
    let twin = brep.half_edges.get(he).and_then(|x| x.mate());
    [Some(he), twin]
        .into_iter()
        .flatten()
        .filter_map(face_of)
        .map(|face| {
            let f = &brep.faces[face];
            let hes = std::iter::once(f.outer_loop)
                .chain(f.inner_loops.iter().copied())
                .flat_map(|l| {
                    let start = brep.loops[l].half_edge;
                    let mut out = vec![start];
                    let mut cur = brep.half_edges[start].next;
                    while cur != start && out.len() <= brep.half_edges.len() {
                        out.push(cur);
                        cur = brep.half_edges[cur].next;
                    }
                    out
                })
                .collect();
            (face, hes, f.surface.as_ref().map(|s| s.geometry.clone()))
        })
        .collect()
}

/// After a KEF merged two co-surfacial faces, re-author the pcurves the
/// survivor inherited from `removed` when the two faces' surfaces are
/// not the same object: coplanar planes can carry different charts, and
/// the inherited pcurves were in the removed face's (PK-806: a post's
/// outer face merged into the base's side face kept post-chart pcurves,
/// the merged face's uv boundary covered only the base part, and the
/// next union clipped every record on the post part away). Identical
/// surfaces share their chart, so their pcurves stay exact.
pub(crate) fn rechart_inherited(
    brep: &mut BRep,
    sides: &[MergeSide],
    removed: popkern_core::topology::FaceHandle,
) {
    let (Some(gone), Some(kept)) = (
        sides.iter().find(|s| s.0 == removed),
        sides.iter().find(|s| s.0 != removed),
    ) else {
        return;
    };
    if format!("{:?}", gone.2) == format!("{:?}", kept.2) {
        return;
    }
    for &h in &gone.1 {
        if let Some(x) = brep.half_edges.get_mut(h) {
            x.pcurve = None;
        }
    }
    let _ = brep.recompute_face_pcurves(kept.0);
}

/// Whether two coincident surfaces share one chart, so faces on them
/// can merge into one face on either's surface. Planes always do (a
/// merged plane face's pcurves are re-derived, and a plane has no seam
/// or pole). A sphere, cylinder, cone or torus does when its axis is
/// the same: the chart's angle is measured from a direction its axis
/// fixes, and its poles and seam lie where the axis puts them. Two
/// faces of one sphere about different axes are the same locus in
/// different charts: merged, the face enclosed the survivor's pole and
/// crossed its seam with no seam edge, its pcurves 0.1 off its edges,
/// and the union of two half balls declined (PK-1138).
pub(crate) fn same_chart(
    a: &popkern_core::geometry::SurfaceGeometry,
    b: &popkern_core::geometry::SurfaceGeometry,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> bool {
    use popkern_core::geometry::SurfaceGeometry as S;
    let same_axis = |x: &popkern_core::geometry::Vec3, y: &popkern_core::geometry::Vec3| {
        x.dot(y) > 1.0 - tolerance.angular()
    };
    match (a, b) {
        (S::Plane { .. }, S::Plane { .. }) => true,
        (S::Sphere { axis: x, .. }, S::Sphere { axis: y, .. })
        | (S::Cylinder { axis: x, .. }, S::Cylinder { axis: y, .. })
        | (S::Cone { axis: x, .. }, S::Cone { axis: y, .. })
        | (S::Torus { axis: x, .. }, S::Torus { axis: y, .. }) => same_axis(x, y),
        _ => false,
    }
}

/// Stations per boundary edge at which [`faces_share_surface`] measures
/// a face's boundary against the other surface (the ends and three
/// interior points: a plane face's distance to another plane is affine
/// along a line, and a curved edge's deviation peaks between its ends).
const BOUNDARY_STATIONS: usize = 4;

/// Whether faces `fa` and `fb` lie on one surface over their extent, so
/// a merge may serve them as one face: their surfaces coincide
/// analytically and share a chart ([`same_chart`]), and every boundary
/// of each face lies within the linear tolerance of the other's
/// surface. The surface comparison alone is scale-blind: two planes
/// whose normals are within the angular tolerance (up to 0.8 degrees)
/// and one of whose origins lies on the other read coincident, and a
/// box united with a block whose top rises from the box's top edge by
/// 1e-3 to 0.1 over 10 served the two tops merged into one face, its
/// far edge 0.1 off its plane (PK-1148). Coincidence is decided by
/// distance over the faces, as plane x plane decides it (PK-1137).
pub(crate) fn faces_share_surface(
    brep: &BRep,
    fa: popkern_core::topology::FaceHandle,
    fb: popkern_core::topology::FaceHandle,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> bool {
    let surface = |f| {
        brep.faces
            .get(f)
            .and_then(|d| d.surface.as_ref())
            .map(|s| &s.geometry)
    };
    let (Some(sa), Some(sb)) = (surface(fa), surface(fb)) else {
        return false;
    };
    same_chart(sa, sb, tolerance) && faces_coincide(brep, (fa, sa), (fb, sb), tolerance)
}

/// [`faces_share_surface`] without the chart condition: the surfaces
/// coincide analytically and each face's boundary lies within the
/// linear tolerance of the other's surface.
pub(crate) fn faces_coincide(
    brep: &BRep,
    (fa, sa): (
        popkern_core::topology::FaceHandle,
        &popkern_core::geometry::SurfaceGeometry,
    ),
    (fb, sb): (
        popkern_core::topology::FaceHandle,
        &popkern_core::geometry::SurfaceGeometry,
    ),
    tolerance: &popkern_core::tolerance::Tolerance,
) -> bool {
    use popkern_core::geometry::intersection::coincident_analytic::are_surfaces_coincident_analytic;
    are_surfaces_coincident_analytic(sa, sb, tolerance)
        && boundary_on_surface(brep, fa, sb, tolerance)
        && boundary_on_surface(brep, fb, sa, tolerance)
}

/// Whether every edge bounding `face` lies within the linear tolerance
/// of `surface` at [`BOUNDARY_STATIONS`] stations (a vertex where the
/// edge has no curve). An edge whose stations project nowhere fails.
fn boundary_on_surface(
    brep: &BRep,
    face: popkern_core::topology::FaceHandle,
    surface: &popkern_core::geometry::SurfaceGeometry,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> bool {
    use popkern_core::geometry::curve::CurveEval as _;
    use popkern_core::geometry::surface::SurfaceEval as _;
    let Some(f) = brep.faces.get(face) else {
        return false;
    };
    let off = |p: &popkern_core::geometry::Point3| -> bool {
        surface
            .closest_point(p)
            .is_none_or(|(u, v)| surface.point_at(u, v).distance_to(p) > tolerance.linear())
    };
    for &lh in std::iter::once(&f.outer_loop).chain(f.inner_loops.iter()) {
        let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
            return false;
        };
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            let Some(h) = brep.half_edges.get(he) else {
                return false;
            };
            let curve = brep.edges.get(h.edge).and_then(|e| e.curve.as_ref());
            match curve {
                Some(c) => {
                    for i in 0..=BOUNDARY_STATIONS {
                        let t =
                            c.t_start + (c.t_end - c.t_start) * i as f64 / BOUNDARY_STATIONS as f64;
                        if off(&c.geometry.point_at(t)) {
                            return false;
                        }
                    }
                }
                None => {
                    if brep.vertices.get(h.origin).is_none_or(|v| off(&v.position)) {
                        return false;
                    }
                }
            }
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    true
}
