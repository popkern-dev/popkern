use super::*;

/// Move a vertex and update the endpoint of any incident edge's
/// 3D curve that previously matched the old vertex position.
/// Preserves loop closure: each `he`'s origin is at the moved
/// vertex, and the corresponding curve endpoint is updated to
/// match.
pub(super) fn update_incident_edges_for_vertex_move(
    brep: &mut BRep,
    v: VertexHandle,
    _old_pos: popkern_core::geometry::Point3,
    new_pos: popkern_core::geometry::Point3,
) {
    use popkern_core::geometry::CurveGeometry;

    // Snapshot half-edges with origin = v.  For each, identify
    // whether it is the edge's "primary" half-edge (= edge.half_edge):
    // primary's origin matches CP[0] (curve start);
    // twin's origin matches CP[last] (curve end).
    //
    // This is more precise than position-based matching: it can
    // disambiguate the V vs V_new ends of a zero-length bridge
    // (where both CPs are at V_pos initially).
    let halves: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter_map(|(h, he)| (he.origin == v).then_some(h))
        .collect();

    for he_h in halves {
        let edge_h = brep.half_edges[he_h].edge;
        let primary = brep.edges.get(edge_h).map(|e| e.half_edge);
        let is_primary = primary == Some(he_h);
        let Some(edge) = brep.edges.get_mut(edge_h) else {
            continue;
        };
        let Some(coe) = edge.curve.as_mut() else {
            continue;
        };
        if let CurveGeometry::Nurbs(nc) = &mut coe.geometry
            && !nc.control_points.is_empty()
        {
            let n = nc.control_points.len();
            if is_primary {
                nc.control_points[0] = new_pos;
            } else {
                nc.control_points[n - 1] = new_pos;
            }
        }
    }

    if let Some(vd) = brep.vertices.get_mut(v) {
        vd.position = new_pos;
    }
}
