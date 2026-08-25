use std::collections::HashSet;

use popkern_core::constructors::BRepBuilder;
use popkern_core::topology::HalfEdgeHandle;

pub(super) fn normalize_open_shell_boundaries(builder: &mut BRepBuilder) {
    let used_half_edges: HashSet<HalfEdgeHandle> = builder
        .brep()
        .half_edges
        .iter()
        .filter_map(|(handle, half_edge)| (half_edge.loop_ != Default::default()).then_some(handle))
        .collect();

    let mut removals = Vec::new();
    for (edge, edge_data) in &builder.brep().edges {
        let first = edge_data.half_edge;
        let second = builder.brep().half_edges[first].twin;
        let first_used = used_half_edges.contains(&first);
        let second_used = used_half_edges.contains(&second);
        match (first_used, second_used) {
            (true, false) => removals.push((edge, first, second, false)),
            (false, true) => removals.push((edge, second, first, true)),
            _ => {}
        }
    }

    for (edge, surviving_half_edge, removed_half_edge, surviving_was_reverse) in removals {
        builder.brep_mut().edges[edge].half_edge = surviving_half_edge;
        if surviving_was_reverse && let Some(curve) = builder.brep_mut().edges[edge].curve.as_mut()
        {
            std::mem::swap(&mut curve.t_start, &mut curve.t_end);
        }
        builder.brep_mut().half_edges[surviving_half_edge].twin = HalfEdgeHandle::default();
        let removed_origin = builder.brep().half_edges[removed_half_edge].origin;
        if builder.brep().vertices[removed_origin].outgoing == removed_half_edge
            && let Some(replacement) = builder
                .brep()
                .half_edges
                .iter()
                .find(|(handle, half_edge)| {
                    *handle != removed_half_edge && half_edge.origin == removed_origin
                })
                .map(|(handle, _)| handle)
        {
            builder.brep_mut().vertices[removed_origin].outgoing = replacement;
        }
        builder
            .brep_mut()
            .remove_half_edge_reanchored(removed_half_edge);
        builder.brep_mut().boundary_edges.insert(edge);
    }
}
