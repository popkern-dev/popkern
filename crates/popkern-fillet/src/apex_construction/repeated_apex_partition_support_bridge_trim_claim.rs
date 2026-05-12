//! Claim support-continuation bridge boundary spans for cylinder trims.

use super::*;

pub(super) fn claim_side_support_bridge_boundary_for_cylinder_trim(
    brep: &mut BRep,
    request: SideSupportBridgeTrimClaimRequest<'_>,
) -> Result<Option<(HalfEdgeHandle, FaceHandle)>, ApexBuildError> {
    let SideSupportBridgeTrimClaimRequest {
        support_face,
        bridge_loops,
        start,
        end,
        curve,
        label,
        fillet_index,
    } = request;
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
        if let Some(result) = author_side_support_bridge_path_patch_for_cylinder_trim(
            brep,
            SideSupportBridgePathPatchRequest {
                support_face,
                bridge_loops,
                start,
                end,
                curve,
                label,
                fillet_index,
            },
        )? {
            return Ok(Some(result));
        }
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

pub(super) struct SideSupportBridgeTrimClaimRequest<'a> {
    pub(super) support_face: FaceHandle,
    pub(super) bridge_loops: &'a [LoopHandle],
    pub(super) start: VertexHandle,
    pub(super) end: VertexHandle,
    pub(super) curve: CurveOnEdge,
    pub(super) label: &'static str,
    pub(super) fillet_index: usize,
}

struct SideSupportBridgePathPatchRequest<'a> {
    support_face: FaceHandle,
    bridge_loops: &'a [LoopHandle],
    start: VertexHandle,
    end: VertexHandle,
    curve: CurveOnEdge,
    label: &'static str,
    fillet_index: usize,
}

fn author_side_support_bridge_path_patch_for_cylinder_trim(
    brep: &mut BRep,
    request: SideSupportBridgePathPatchRequest<'_>,
) -> Result<Option<(HalfEdgeHandle, FaceHandle)>, ApexBuildError> {
    let SideSupportBridgePathPatchRequest {
        support_face,
        bridge_loops,
        start,
        end,
        curve,
        label,
        fillet_index,
    } = request;
    let Some(path) = bridge_boundary_path_for_cylinder_side(brep, bridge_loops, start, end) else {
        return Ok(None);
    };
    if path.len() < 2 {
        return Ok(None);
    }
    let rollback = brep.clone();
    let before = BridgePathPatchCounts::capture(brep);
    let result = (|| {
        let (edge, cylinder_half_edge, support_half_edge) =
            insert_unwired_edge_pair(brep, start, end);
        brep.edges[edge].curve = Some(curve);
        let mut loop_half_edges = Vec::with_capacity(path.len() + 1);
        loop_half_edges.push(support_half_edge);
        for step in &path {
            loop_half_edges.push(author_bridge_path_patch_span(
                brep,
                *step,
                fillet_index,
                label,
            )?);
        }
        let support_patch_face = insert_face_with_outer_loop(brep, support_face, &loop_half_edges)?;
        brep.faces[support_patch_face].role = TopologyRole::ConstructionResidual;
        for half_edge in &loop_half_edges {
            let edge = brep.half_edges[*half_edge].edge;
            brep.edges[edge].role = TopologyRole::ConstructionResidual;
        }
        validate_bridge_path_patch_postconditions(
            brep,
            BridgePathPatchPostconditions {
                before,
                support_patch_face,
                bridge_loops,
                cylinder_half_edge,
                support_half_edge,
                loop_half_edges: &loop_half_edges,
                path: &path,
                label,
                fillet_index,
            },
        )?;
        if crate::a2_trace::enabled() {
            let path_label = path
                .iter()
                .map(|step| format!("{:?}->{:?}", step.start, step.end))
                .collect::<Vec<_>>()
                .join(", ");
            eprintln!(
                "[a2-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] authored {label} \
                 {start:?}->{end:?} as side-support bridge path patch {support_patch_face:?}; \
                 cylinder side {cylinder_half_edge:?}, support side {support_half_edge:?}, \
                 path [{path_label}]"
            );
        }
        Ok((cylinder_half_edge, support_patch_face))
    })();
    if result.is_err() {
        *brep = rollback;
    }
    result.map(Some)
}

#[derive(Clone, Copy)]
struct BridgePathPatchCounts {
    edges: usize,
    half_edges: usize,
    faces: usize,
    loops: usize,
}

impl BridgePathPatchCounts {
    fn capture(brep: &BRep) -> Self {
        Self {
            edges: brep.edges.len(),
            half_edges: brep.half_edges.len(),
            faces: brep.faces.len(),
            loops: brep.loops.len(),
        }
    }
}

fn author_bridge_path_patch_span(
    brep: &mut BRep,
    step: BridgeBoundaryPathStep,
    fillet_index: usize,
    label: &'static str,
) -> Result<HalfEdgeHandle, ApexBuildError> {
    let (edge, half_edge) = insert_boundary_half_edge(brep, step.start);
    brep.edges[edge].curve = bridge_path_patch_span_curve(brep, step, fillet_index, label)?;
    if support_split_half_edge_end_vertex(brep, half_edge) != Some(step.end) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] authored bridge-path \
                 patch span for {label} with wrong endpoint: expected {:?}->{:?}, got \
                 {half_edge:?}",
                step.start, step.end
            ),
        });
    }
    Ok(half_edge)
}

fn bridge_path_patch_span_curve(
    brep: &BRep,
    step: BridgeBoundaryPathStep,
    fillet_index: usize,
    label: &'static str,
) -> Result<Option<CurveOnEdge>, ApexBuildError> {
    if let Some(curve) = brep
        .half_edges
        .get(step.half_edge)
        .and_then(|half_edge| brep.edges.get(half_edge.edge))
        .and_then(|edge| edge.curve.clone())
    {
        return Ok(Some(curve));
    }
    let start = brep.vertices.get(step.start).ok_or_else(|| ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author bridge-path \
             patch span for {label}: stale start vertex {:?}",
            step.start
        ),
    })?;
    let end = brep.vertices.get(step.end).ok_or_else(|| ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author bridge-path \
             patch span for {label}: stale end vertex {:?}",
            step.end
        ),
    })?;
    CurveOnEdge::line_segment(start.position, end.position, brep.tolerance.linear())
        .map(Some)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                 degenerate bridge-path patch span {:?}->{:?} for {label}",
                step.start, step.end
            ),
        })
}

fn validate_bridge_path_patch_postconditions(
    brep: &BRep,
    request: BridgePathPatchPostconditions<'_>,
) -> Result<(), ApexBuildError> {
    let BridgePathPatchPostconditions {
        before,
        support_patch_face,
        bridge_loops,
        cylinder_half_edge,
        support_half_edge,
        loop_half_edges,
        path,
        label,
        fillet_index,
    } = request;
    let after = BridgePathPatchCounts::capture(brep);
    if after.edges != before.edges + path.len() + 1
        || after.half_edges != before.half_edges + path.len() + 2
        || after.faces != before.faces + 1
        || after.loops != before.loops + 1
    {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] bridge-path patch for \
                 {label} violated topology delta dV=0 dE={} dHE={} dF=1 dL=1",
                path.len() + 1,
                path.len() + 2
            ),
        });
    }
    let cylinder_data = brep.half_edges.get(cylinder_half_edge).ok_or_else(|| {
        ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] bridge-path patch for \
                 {label} produced stale cylinder side {cylinder_half_edge:?}"
            ),
        }
    })?;
    if cylinder_data.twin != support_half_edge {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] bridge-path patch for \
                 {label} did not twin cylinder side {cylinder_half_edge:?} to support side \
                 {support_half_edge:?}"
            ),
        });
    }
    let support_loop = brep
        .half_edges
        .get(support_half_edge)
        .map(|half_edge| half_edge.loop_)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] bridge-path patch for \
                 {label} produced stale support side {support_half_edge:?}"
            ),
        })?;
    if brep.loops.get(support_loop).map(|loop_data| loop_data.face) != Some(support_patch_face) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] bridge-path patch for \
                 {label} did not assign support side {support_half_edge:?} to support patch \
                 {support_patch_face:?}"
            ),
        });
    }
    validate_half_edge_cycle(
        brep,
        loop_half_edges,
        "side-support bridge path patch",
        support_patch_face,
    )?;
    for loop_handle in bridge_loops {
        if let Some(loop_data) = brep.loops.get(*loop_handle) {
            let bridge_cycle = collect_loop_half_edges(brep, *loop_handle).map_err(|err| {
                ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] bridge-path \
                         patch for {label} damaged bridge loop {loop_handle:?}: {err}"
                    ),
                }
            })?;
            if bridge_cycle.is_empty()
                || bridge_cycle
                    .iter()
                    .any(|half_edge| brep.half_edges[*half_edge].loop_ != *loop_handle)
            {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] bridge-path \
                         patch for {label} damaged source bridge loop {loop_handle:?} on \
                         {:?}",
                        loop_data.face
                    ),
                });
            }
        }
    }
    Ok(())
}

struct BridgePathPatchPostconditions<'a> {
    before: BridgePathPatchCounts,
    support_patch_face: FaceHandle,
    bridge_loops: &'a [LoopHandle],
    cylinder_half_edge: HalfEdgeHandle,
    support_half_edge: HalfEdgeHandle,
    loop_half_edges: &'a [HalfEdgeHandle],
    path: &'a [BridgeBoundaryPathStep],
    label: &'static str,
    fillet_index: usize,
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

    #[test]
    fn bridge_path_patch_authors_support_owned_trim_without_rewiring_source_loops() {
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
        let support_face =
            add_boundary_loop_face(&mut brep, shell, &[vertices[0], vertices[2], vertices[3]]);
        let first_loop =
            add_boundary_loop(&mut brep, shell, &[vertices[0], vertices[1], vertices[3]]);
        let second_loop =
            add_boundary_loop(&mut brep, shell, &[vertices[1], vertices[2], vertices[3]]);
        let first_loop_before = collect_loop_half_edges(&brep, first_loop).expect("first loop");
        let second_loop_before = collect_loop_half_edges(&brep, second_loop).expect("second loop");
        let curve = CurveOnEdge::line_segment(
            brep.vertices[vertices[0]].position,
            brep.vertices[vertices[2]].position,
            brep.tolerance.linear(),
        )
        .expect("nondegenerate trim curve");

        let (cylinder_half_edge, support_patch_face) =
            author_side_support_bridge_path_patch_for_cylinder_trim(
                &mut brep,
                SideSupportBridgePathPatchRequest {
                    support_face,
                    bridge_loops: &[first_loop, second_loop],
                    start: vertices[0],
                    end: vertices[2],
                    curve,
                    label: "test generator",
                    fillet_index: 0,
                },
            )
            .expect("patch authoring succeeds")
            .expect("path patch applies");

        let support_half_edge = brep.half_edges[cylinder_half_edge].twin;
        let support_loop = brep.half_edges[support_half_edge].loop_;
        assert_eq!(brep.loops[support_loop].face, support_patch_face);
        assert_eq!(
            half_edge_end_vertex(&brep, cylinder_half_edge),
            Some(vertices[2])
        );
        assert_eq!(
            support_split_half_edge_end_vertex(&brep, support_half_edge),
            Some(vertices[0])
        );
        assert_eq!(
            collect_loop_half_edges(&brep, first_loop).expect("first loop after"),
            first_loop_before
        );
        assert_eq!(
            collect_loop_half_edges(&brep, second_loop).expect("second loop after"),
            second_loop_before
        );
        validate_replacement_loop_support_side_ownership(
            &brep,
            FaceHandle::default(),
            0,
            &[cylinder_half_edge],
        )
        .expect("cylinder trim has a live non-cylinder support side");
    }

    fn add_boundary_loop(
        brep: &mut BRep,
        shell: popkern_core::topology::ShellHandle,
        vertices: &[VertexHandle],
    ) -> LoopHandle {
        let face = add_boundary_loop_face(brep, shell, vertices);
        brep.faces[face].outer_loop
    }

    fn add_boundary_loop_face(
        brep: &mut BRep,
        shell: popkern_core::topology::ShellHandle,
        vertices: &[VertexHandle],
    ) -> FaceHandle {
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
        face
    }
}
