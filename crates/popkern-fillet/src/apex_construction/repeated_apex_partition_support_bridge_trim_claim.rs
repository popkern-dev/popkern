//! Claim support-continuation bridge boundary spans for cylinder trims.

use super::*;

pub(super) fn claim_side_support_bridge_boundary_for_cylinder_trim(
    brep: &mut BRep,
    bridge_loops: &[LoopHandle],
    start: VertexHandle,
    end: VertexHandle,
    curve: CurveOnEdge,
    label: &'static str,
    fillet_index: usize,
) -> Result<Option<(HalfEdgeHandle, FaceHandle)>, ApexBuildError> {
    let Some((bridge_face, support_half_edge)) = bridge_loops.iter().find_map(|loop_handle| {
        let loop_data = brep.loops.get(*loop_handle)?;
        collect_loop_half_edges(brep, *loop_handle)
            .ok()?
            .into_iter()
            .find(|half_edge| {
                let Some(data) = brep.half_edges.get(*half_edge) else {
                    return false;
                };
                data.origin == end
                    && support_split_half_edge_end_vertex(brep, *half_edge) == Some(start)
                    && brep.boundary_edges.contains(&data.edge)
                    && !brep.half_edges.contains_key(data.twin)
            })
            .map(|half_edge| (loop_data.face, half_edge))
    }) else {
        trace_side_support_bridge_claim_miss(brep, bridge_loops, start, end, label, fillet_index);
        return Ok(None);
    };
    let rollback = brep.clone();
    let result = (|| {
        let edge = brep
            .half_edges
            .get(support_half_edge)
            .map(|half_edge| half_edge.edge)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot claim \
                     stale side-support bridge boundary side {support_half_edge:?} for {label}"
                ),
            })?;
        if brep.edges[edge].curve.is_none() {
            brep.edges[edge].curve = Some(curve);
        }
        let half_edge = brep.half_edges.insert(HalfEdge {
            origin: start,
            twin: support_half_edge,
            next: HalfEdgeHandle::default(),
            prev: HalfEdgeHandle::default(),
            edge,
            loop_: LoopHandle::default(),
            pcurve: None,
        });
        brep.half_edges[support_half_edge].twin = half_edge;
        brep.boundary_edges.remove(&edge);
        validate_side_support_bridge_claim(
            brep,
            bridge_face,
            support_half_edge,
            half_edge,
            label,
            fillet_index,
        )?;
        if crate::a2_trace::enabled() {
            eprintln!(
                "[a2-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] claimed {label} \
                 {start:?}->{end:?} from side-support bridge loop on {bridge_face:?}; cylinder \
                 side {half_edge:?}, support side {support_half_edge:?}"
            );
        }
        Ok((half_edge, bridge_face))
    })();
    if result.is_err() {
        *brep = rollback;
    }
    result.map(Some)
}

fn trace_side_support_bridge_claim_miss(
    brep: &BRep,
    bridge_loops: &[LoopHandle],
    start: VertexHandle,
    end: VertexHandle,
    label: &'static str,
    fillet_index: usize,
) {
    if !crate::a2_trace::enabled() {
        return;
    }
    let reports = bridge_loops
        .iter()
        .map(|loop_handle| {
            let face = brep.loops.get(*loop_handle).map(|loop_data| loop_data.face);
            let spans = collect_loop_half_edges(brep, *loop_handle)
                .unwrap_or_default()
                .into_iter()
                .filter_map(|half_edge| {
                    let data = brep.half_edges.get(half_edge)?;
                    let span_end = support_split_half_edge_end_vertex(brep, half_edge)?;
                    let boundary = brep.boundary_edges.contains(&data.edge);
                    let twin_live = brep.half_edges.contains_key(data.twin);
                    Some(format!(
                        "{half_edge:?}:{:?}->{span_end:?} boundary={boundary} twin_live={twin_live}",
                        data.origin
                    ))
                })
                .collect::<Vec<_>>()
                .join(", ");
            format!("{loop_handle:?}/{face:?} [{spans}]")
        })
        .collect::<Vec<_>>()
        .join("; ");
    let path_diagnostic = bridge_boundary_path_for_cylinder_side(brep, bridge_loops, start, end)
        .map(|path| {
            let spans = path
                .iter()
                .map(|step| {
                    format!(
                        "{:?}/{:?}:{:?}->{:?}",
                        step.loop_, step.half_edge, step.start, step.end
                    )
                })
                .collect::<Vec<_>>()
                .join(", ");
            format!(
                " boundary_path_exists_but_requires_loop_regrouping=[{spans}]; this path cannot \
                 be claimed as one trim twin because its half-edges already belong to bridge \
                 loops"
            )
        })
        .unwrap_or_default();
    eprintln!(
        "[a2-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] found no side-support bridge \
         boundary side to claim for {label} {start:?}->{end:?}; bridge candidates {reports}{path_diagnostic}"
    );
}

#[derive(Clone, Copy)]
struct BridgeBoundaryPathStep {
    loop_: LoopHandle,
    half_edge: HalfEdgeHandle,
    start: VertexHandle,
    end: VertexHandle,
}

fn bridge_boundary_path_for_cylinder_side(
    brep: &BRep,
    bridge_loops: &[LoopHandle],
    start: VertexHandle,
    end: VertexHandle,
) -> Option<Vec<BridgeBoundaryPathStep>> {
    let steps = bridge_loops
        .iter()
        .flat_map(|loop_handle| {
            collect_loop_half_edges(brep, *loop_handle)
                .unwrap_or_default()
                .into_iter()
                .filter_map(|half_edge| {
                    let data = brep.half_edges.get(half_edge)?;
                    if !brep.boundary_edges.contains(&data.edge)
                        || brep.half_edges.contains_key(data.twin)
                    {
                        return None;
                    }
                    let span_end = support_split_half_edge_end_vertex(brep, half_edge)?;
                    Some(BridgeBoundaryPathStep {
                        loop_: *loop_handle,
                        half_edge,
                        start: data.origin,
                        end: span_end,
                    })
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    shortest_bridge_boundary_path(start, end, &steps)
}

fn shortest_bridge_boundary_path(
    start: VertexHandle,
    end: VertexHandle,
    steps: &[BridgeBoundaryPathStep],
) -> Option<Vec<BridgeBoundaryPathStep>> {
    let mut queue = std::collections::VecDeque::from([(start, Vec::new())]);
    let mut best_depth_by_vertex = std::collections::BTreeMap::from([(start, 0usize)]);
    while let Some((cursor, path)) = queue.pop_front() {
        if cursor == end && !path.is_empty() {
            return Some(path);
        }
        for step in steps.iter().filter(|step| step.start == cursor) {
            if path
                .iter()
                .any(|visited: &BridgeBoundaryPathStep| visited.half_edge == step.half_edge)
            {
                continue;
            }
            let next_depth = path.len() + 1;
            if best_depth_by_vertex
                .get(&step.end)
                .is_some_and(|best_depth| *best_depth <= next_depth)
            {
                continue;
            }
            best_depth_by_vertex.insert(step.end, next_depth);
            let mut next_path = path.clone();
            next_path.push(*step);
            queue.push_back((step.end, next_path));
        }
    }
    None
}

fn validate_side_support_bridge_claim(
    brep: &BRep,
    bridge_face: FaceHandle,
    support_half_edge: HalfEdgeHandle,
    cylinder_half_edge: HalfEdgeHandle,
    label: &'static str,
    fillet_index: usize,
) -> Result<(), ApexBuildError> {
    let support_owner = brep
        .half_edges
        .get(support_half_edge)
        .and_then(|half_edge| brep.loops.get(half_edge.loop_))
        .map(|loop_data| loop_data.face);
    if support_owner != Some(bridge_face) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] claimed {label} from \
                 side-support bridge half-edge {support_half_edge:?}, but its owner is \
                 {support_owner:?} instead of {bridge_face:?}"
            ),
        });
    }
    let Some(cylinder_data) = brep.half_edges.get(cylinder_half_edge) else {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] claimed stale cylinder \
                 half-edge {cylinder_half_edge:?} for {label}"
            ),
        });
    };
    if cylinder_data.twin != support_half_edge
        || brep.half_edges[support_half_edge].twin != cylinder_half_edge
        || brep.boundary_edges.contains(&cylinder_data.edge)
    {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] claimed {label} with \
                 invalid bridge ownership on {cylinder_half_edge:?}/{support_half_edge:?}"
            ),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use popkern_core::geometry::Point3;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::{Edge, Face, HalfEdge, Loop, Shell, Solid, Vertex};

    use super::*;

    #[test]
    fn bridge_boundary_path_detects_multi_loop_regrouping_need() {
        let mut brep = BRep::new(Tolerance::default());
        let solid = brep.solids.insert(Solid {
            shell: Default::default(),
        });
        let shell = brep.shells.insert(Shell {
            face: Default::default(),
            solid,
        });
        brep.solids[solid].shell = shell;
        let vertices = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ]
        .map(|position| {
            brep.vertices.insert(Vertex {
                position,
                outgoing: Default::default(),
                name: None,
                tolerance: None,
            })
        });
        let first_loop =
            add_boundary_loop(&mut brep, shell, &[vertices[0], vertices[1], vertices[3]]);
        let second_loop =
            add_boundary_loop(&mut brep, shell, &[vertices[1], vertices[2], vertices[3]]);

        let path = bridge_boundary_path_for_cylinder_side(
            &brep,
            &[first_loop, second_loop],
            vertices[0],
            vertices[2],
        )
        .expect("multi-loop boundary path is detected");

        assert_eq!(path.len(), 2);
        assert_eq!(path[0].start, vertices[0]);
        assert_eq!(path[0].end, vertices[1]);
        assert_eq!(path[1].start, vertices[1]);
        assert_eq!(path[1].end, vertices[2]);
    }

    fn add_boundary_loop(
        brep: &mut BRep,
        shell: popkern_core::topology::ShellHandle,
        vertices: &[VertexHandle],
    ) -> LoopHandle {
        let face = brep.faces.insert(Face {
            outer_loop: Default::default(),
            inner_loops: vec![],
            shell,
            surface: None,
            name: None,
            tolerance: None,
            role: Default::default(),
        });
        let loop_ = brep.loops.insert(Loop {
            half_edge: Default::default(),
            face,
        });
        brep.faces[face].outer_loop = loop_;
        let half_edges = vertices
            .iter()
            .copied()
            .map(|vertex| {
                let edge = brep.edges.insert(Edge {
                    half_edge: Default::default(),
                    curve: None,
                    name: None,
                    tolerance: None,
                    canonical: Default::default(),
                    role: Default::default(),
                });
                let half_edge = brep.half_edges.insert(HalfEdge {
                    origin: vertex,
                    twin: Default::default(),
                    next: Default::default(),
                    prev: Default::default(),
                    edge,
                    loop_,
                    pcurve: None,
                });
                brep.edges[edge].half_edge = half_edge;
                brep.boundary_edges.insert(edge);
                if brep.vertices[vertex].outgoing == HalfEdgeHandle::default() {
                    brep.vertices[vertex].outgoing = half_edge;
                }
                half_edge
            })
            .collect::<Vec<_>>();
        for (index, half_edge) in half_edges.iter().copied().enumerate() {
            brep.half_edges[half_edge].next = half_edges[(index + 1) % half_edges.len()];
            brep.half_edges[half_edge].prev =
                half_edges[(index + half_edges.len() - 1) % half_edges.len()];
        }
        brep.loops[loop_].half_edge = half_edges[0];
        loop_
    }
}
