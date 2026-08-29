use super::*;

pub(super) fn half_edge_is_degenerate(brep: &BRep, half_edge: HalfEdgeHandle) -> bool {
    let Some(start) = brep.half_edges.get(half_edge).map(|data| data.origin) else {
        return true;
    };
    let Some(end) = half_edge_end_vertex(brep, half_edge) else {
        return true;
    };
    let Some(start_point) = brep.vertices.get(start).map(|vertex| vertex.position) else {
        return true;
    };
    let Some(end_point) = brep.vertices.get(end).map(|vertex| vertex.position) else {
        return true;
    };
    start_point.distance_to(&end_point) <= brep.tolerance.linear().max(1.0e-9)
}

pub(super) fn vertices_are_coincident(brep: &BRep, start: VertexHandle, end: VertexHandle) -> bool {
    let Some(start_point) = brep.vertices.get(start).map(|vertex| vertex.position) else {
        return false;
    };
    let Some(end_point) = brep.vertices.get(end).map(|vertex| vertex.position) else {
        return false;
    };
    start_point.distance_to(&end_point) <= brep.tolerance.linear().max(1.0e-9)
}

pub(super) fn support_cycle_open_chain_diagnostic(
    brep: &BRep,
    half_edges: &[HalfEdgeHandle],
) -> String {
    let position = |vertex: VertexHandle| -> String {
        brep.vertices
            .get(vertex)
            .map(|data| {
                format!(
                    "({:.4},{:.4},{:.4})",
                    data.position.x, data.position.y, data.position.z
                )
            })
            .unwrap_or_else(|| "stale".to_string())
    };
    half_edges
        .iter()
        .map(|half_edge| {
            let start = brep
                .half_edges
                .get(*half_edge)
                .map(|data| data.origin)
                .unwrap_or_default();
            let end = support_split_half_edge_end_vertex(brep, *half_edge).unwrap_or_default();
            format!(
                "{half_edge:?}:{start:?}{}->{end:?}{}",
                position(start),
                position(end)
            )
        })
        .collect::<Vec<_>>()
        .join(", ")
}

pub(super) fn support_split_half_edge_end_vertex(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Option<VertexHandle> {
    half_edge_end_vertex(brep, half_edge)
        .or_else(|| boundary_curve_end_vertex(brep, half_edge))
        .or_else(|| {
            let next = brep.half_edges.get(half_edge)?.next;
            brep.half_edges
                .get(next)
                .map(|next_half_edge| next_half_edge.origin)
        })
}

fn boundary_curve_end_vertex(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<VertexHandle> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    if !brep.boundary_edges.contains(&half_edge_data.edge) {
        return None;
    }
    let curve = brep.edges.get(half_edge_data.edge)?.curve.as_ref()?;
    let start_point = brep.vertices.get(half_edge_data.origin)?.position;
    let p0 = curve.geometry.point_at(curve.t_start);
    let p1 = curve.geometry.point_at(curve.t_end);
    let end_point = if start_point.distance_to(&p0) <= start_point.distance_to(&p1) {
        p1
    } else {
        p0
    };
    let tol = brep.tolerance.linear().max(1.0e-9);
    brep.vertices
        .iter()
        .filter_map(|(vertex, data)| {
            let distance = data.position.distance_to(&end_point);
            (distance <= tol).then_some((vertex, distance))
        })
        .min_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(vertex, _)| vertex)
}

pub(super) fn closed_half_edge_cycles_from_unwired_set(
    brep: &BRep,
    half_edges: &[HalfEdgeHandle],
) -> Vec<Vec<HalfEdgeHandle>> {
    let mut by_start: std::collections::BTreeMap<VertexHandle, Vec<HalfEdgeHandle>> =
        std::collections::BTreeMap::new();
    for half_edge in half_edges {
        if let Some(start) = brep.half_edges.get(*half_edge).map(|data| data.origin) {
            by_start.entry(start).or_default().push(*half_edge);
        }
    }
    if by_start.values().any(|edges| edges.len() != 1) {
        return Vec::new();
    }
    let input = half_edges
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    let mut visited = std::collections::BTreeSet::new();
    let mut cycles = Vec::new();
    for start_half_edge in half_edges {
        if visited.contains(start_half_edge) {
            continue;
        }
        let mut cycle = Vec::new();
        let mut cursor = *start_half_edge;
        for _ in 0..=half_edges.len() {
            if !input.contains(&cursor) || !visited.insert(cursor) {
                return Vec::new();
            }
            cycle.push(cursor);
            let Some(end) = support_split_half_edge_end_vertex(brep, cursor) else {
                return Vec::new();
            };
            let Some(next) = support_split_next_half_edge_from_end(brep, &by_start, end) else {
                return Vec::new();
            };
            if next == *start_half_edge {
                cycles.push(cycle);
                break;
            }
            cursor = next;
        }
    }
    if visited.len() == half_edges.len() {
        cycles
    } else {
        Vec::new()
    }
}

fn support_split_next_half_edge_from_end(
    brep: &BRep,
    by_start: &std::collections::BTreeMap<VertexHandle, Vec<HalfEdgeHandle>>,
    end: VertexHandle,
) -> Option<HalfEdgeHandle> {
    if let Some(next) = by_start.get(&end).and_then(|edges| edges.first()).copied() {
        return Some(next);
    }
    let coincident = by_start
        .iter()
        .filter_map(|(start, edges)| {
            vertices_are_coincident(brep, end, *start).then_some(edges.first().copied())
        })
        .flatten()
        .collect::<Vec<_>>();
    if coincident.len() == 1 {
        Some(coincident[0])
    } else {
        None
    }
}
