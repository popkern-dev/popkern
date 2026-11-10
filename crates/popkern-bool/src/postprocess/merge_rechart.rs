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
