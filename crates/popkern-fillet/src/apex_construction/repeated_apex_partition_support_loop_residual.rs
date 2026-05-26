use super::*;

pub(in crate::apex_construction) fn materialize_authored_support_boundary_residual_cycles(
    brep: &mut BRep,
    boundary_spans: &[HalfEdgeHandle],
    existing_residual_sides: &[HalfEdgeHandle],
) -> Result<(), ApexBuildError> {
    materialize_support_boundary_residual_cycles_candidate(
        brep,
        boundary_spans,
        existing_residual_sides,
    )
}

pub(super) fn materialize_support_boundary_residual_cycles(
    brep: &mut BRep,
    boundary_spans: &[HalfEdgeHandle],
) -> Result<(), ApexBuildError> {
    materialize_support_boundary_residual_cycles_candidate(brep, boundary_spans, &[])
}

fn materialize_support_boundary_residual_cycles_candidate(
    brep: &mut BRep,
    boundary_spans: &[HalfEdgeHandle],
    existing_residual_sides: &[HalfEdgeHandle],
) -> Result<(), ApexBuildError> {
    let boundary_spans = boundary_spans
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    let existing_residual_sides = existing_residual_sides
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    if boundary_spans.is_empty() && existing_residual_sides.is_empty() {
        return Ok(());
    }
    let rollback = brep.clone();
    let result = (|| {
        let mut residual_sides = existing_residual_sides
            .iter()
            .copied()
            .filter(|half_edge| {
                brep.half_edges.get(*half_edge).is_some_and(|data| {
                    data.loop_ == LoopHandle::default() && brep.half_edges.contains_key(data.twin)
                })
            })
            .collect::<Vec<_>>();
        let mut template_face = None;
        for half_edge in &residual_sides {
            template_face = template_face.or_else(|| {
                brep.half_edges
                    .get(*half_edge)
                    .and_then(|data| brep.half_edges.get(data.twin))
                    .and_then(|twin| brep.loops.get(twin.loop_))
                    .map(|loop_data| loop_data.face)
            });
        }
        for half_edge in boundary_spans {
            let half_edge_data = brep.half_edges.get(half_edge).cloned().ok_or_else(|| {
                ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition cannot close support boundary \
                         residual for stale half-edge {half_edge:?}"
                    ),
                }
            })?;
            if !brep.boundary_edges.contains(&half_edge_data.edge)
                || brep.half_edges.contains_key(half_edge_data.twin)
            {
                continue;
            }
            template_face = template_face.or_else(|| {
                brep.loops
                    .get(half_edge_data.loop_)
                    .map(|loop_data| loop_data.face)
            });
            let residual_origin =
                support_split_half_edge_end_vertex(brep, half_edge).ok_or_else(|| {
                    ApexBuildError::PlanMismatch {
                        reason: format!(
                            "planned repeated-apex partition cannot close support boundary \
                             residual for {half_edge:?}: missing endpoint"
                        ),
                    }
                })?;
            let residual = brep.half_edges.insert(HalfEdge {
                origin: residual_origin,
                twin: half_edge,
                next: HalfEdgeHandle::default(),
                prev: HalfEdgeHandle::default(),
                edge: half_edge_data.edge,
                loop_: LoopHandle::default(),
                pcurve: None,
            });
            brep.half_edges[half_edge].twin = residual;
            brep.boundary_edges.remove(&half_edge_data.edge);
            residual_sides.push(residual);
        }
        if residual_sides.is_empty() {
            return Ok(());
        }
        let cycles = residual_cycles_with_connector_closure(brep, &residual_sides);
        if cycles.is_empty() {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot close support boundary residual \
                     ownership: residual sides {:?} did not form closed cycles; {}",
                    residual_sides,
                    support_cycle_open_chain_diagnostic(brep, &residual_sides)
                ),
            });
        }
        let template_face = template_face.ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: "planned repeated-apex partition cannot close support boundary residual \
                     ownership without a template face"
                .to_string(),
        })?;
        for cycle in cycles {
            validate_half_edge_cycle(
                brep,
                &cycle,
                "support boundary residual cycle",
                template_face,
            )?;
            let residual_face = insert_face_with_outer_loop(brep, template_face, &cycle)?;
            brep.faces[residual_face].role = TopologyRole::ConstructionResidual;
            brep.faces[residual_face].surface = None;
        }
        Ok(())
    })();
    if result.is_err() {
        *brep = rollback;
    }
    result
}

fn closed_residual_cycles_with_shared_vertices(
    brep: &BRep,
    half_edges: &[HalfEdgeHandle],
) -> Vec<Vec<HalfEdgeHandle>> {
    let input = half_edges
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    if input.is_empty() {
        return Vec::new();
    }
    let mut by_start = std::collections::BTreeMap::<VertexHandle, Vec<HalfEdgeHandle>>::new();
    for half_edge in &input {
        let Some(start) = brep.half_edges.get(*half_edge).map(|data| data.origin) else {
            return Vec::new();
        };
        by_start.entry(start).or_default().push(*half_edge);
    }
    for half_edges in by_start.values_mut() {
        half_edges.sort_by_key(|half_edge| format!("{half_edge:?}"));
    }

    let mut visited = std::collections::BTreeSet::new();
    let mut cycles = Vec::new();
    for start_half_edge in &input {
        if visited.contains(start_half_edge) {
            continue;
        }
        let Some(start_vertex) = brep
            .half_edges
            .get(*start_half_edge)
            .map(|data| data.origin)
        else {
            return Vec::new();
        };
        let mut cursor = *start_half_edge;
        let mut cycle = Vec::new();
        for _ in 0..=input.len() {
            if !input.contains(&cursor) || !visited.insert(cursor) {
                return Vec::new();
            }
            cycle.push(cursor);
            let Some(end) = support_split_half_edge_end_vertex(brep, cursor) else {
                return Vec::new();
            };
            if end == start_vertex || vertices_are_coincident(brep, end, start_vertex) {
                cycles.push(cycle);
                break;
            }
            let Some(next) =
                next_unvisited_residual_side_from_vertex(brep, &by_start, &visited, end)
            else {
                return Vec::new();
            };
            cursor = next;
        }
    }
    if visited.len() == input.len() {
        cycles
    } else {
        Vec::new()
    }
}

fn residual_cycles_with_connector_closure(
    brep: &mut BRep,
    half_edges: &[HalfEdgeHandle],
) -> Vec<Vec<HalfEdgeHandle>> {
    let closed = closed_residual_cycles_with_shared_vertices(brep, half_edges);
    if !closed.is_empty() {
        return closed;
    }
    let input = half_edges
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    if input.is_empty() {
        return Vec::new();
    }
    let mut by_start = std::collections::BTreeMap::<VertexHandle, Vec<HalfEdgeHandle>>::new();
    for half_edge in &input {
        let Some(start) = brep.half_edges.get(*half_edge).map(|data| data.origin) else {
            return Vec::new();
        };
        by_start.entry(start).or_default().push(*half_edge);
    }
    for half_edges in by_start.values_mut() {
        half_edges.sort_by_key(|half_edge| format!("{half_edge:?}"));
    }

    let mut visited = std::collections::BTreeSet::new();
    let mut cycles = Vec::new();
    for start_half_edge in &input {
        if visited.contains(start_half_edge) {
            continue;
        }
        let Some(start_vertex) = brep
            .half_edges
            .get(*start_half_edge)
            .map(|data| data.origin)
        else {
            return Vec::new();
        };
        let mut cursor = *start_half_edge;
        let mut cycle = Vec::new();
        for _ in 0..=input.len() {
            if !input.contains(&cursor) || !visited.insert(cursor) {
                return Vec::new();
            }
            cycle.push(cursor);
            let Some(end) = support_split_half_edge_end_vertex(brep, cursor) else {
                return Vec::new();
            };
            if end == start_vertex || vertices_are_coincident(brep, end, start_vertex) {
                cycles.push(cycle);
                break;
            }
            if let Some(next) =
                next_unvisited_residual_side_from_vertex(brep, &by_start, &visited, end)
            {
                cursor = next;
                continue;
            }
            let (_, connector, _) = insert_unwired_edge_pair(brep, end, start_vertex);
            cycle.push(connector);
            cycles.push(cycle);
            break;
        }
    }
    if visited.len() == input.len() {
        cycles
    } else {
        Vec::new()
    }
}

fn next_unvisited_residual_side_from_vertex(
    brep: &BRep,
    by_start: &std::collections::BTreeMap<VertexHandle, Vec<HalfEdgeHandle>>,
    visited: &std::collections::BTreeSet<HalfEdgeHandle>,
    end: VertexHandle,
) -> Option<HalfEdgeHandle> {
    by_start
        .get(&end)
        .and_then(|edges| edges.iter().find(|edge| !visited.contains(edge)).copied())
        .or_else(|| {
            by_start.iter().find_map(|(start, edges)| {
                vertices_are_coincident(brep, end, *start)
                    .then(|| edges.iter().find(|edge| !visited.contains(edge)).copied())
                    .flatten()
            })
        })
}
