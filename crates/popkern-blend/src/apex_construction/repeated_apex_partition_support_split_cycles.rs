use super::super::constructive_topology;
use super::*;

pub(super) struct SupportSplitCyclePlan {
    pub(super) cycles: Vec<Vec<HalfEdgeHandle>>,
    pub(super) authored_boundary_spans: Vec<HalfEdgeHandle>,
    missing_spans: Vec<(VertexHandle, VertexHandle)>,
    /// True when the trim set neither chained into closed cycles nor
    /// produced authorable two-sided spans (balanced endpoint
    /// multiplicities with degenerate or ambiguous chaining, e.g.
    /// coincident duplicate corner vertices on a shared-corner
    /// repeated apex). PK-216: this state used to read as complete —
    /// `is_complete()` checked only `missing_spans` — and the
    /// partition corrupt-accepted a half-wired shell (the glued
    /// duplicate apex membrane).
    degenerate_unclosed: bool,
}

impl SupportSplitCyclePlan {
    pub(super) fn is_complete(&self) -> bool {
        self.missing_spans.is_empty() && !self.degenerate_unclosed
    }

    pub(super) fn diagnostic(&self) -> String {
        if self.degenerate_unclosed {
            return "support trim sides neither chain into closed cycles nor leave authorable \
                    two-sided spans (balanced endpoint multiplicities with degenerate chaining \
                    — coincident duplicate corner vertices?)"
                .to_string();
        }
        if self.missing_spans.is_empty() {
            return "support split cycles are closed".to_string();
        }
        let spans = self
            .missing_spans
            .iter()
            .map(|(start, end)| format!("{start:?}->{end:?}"))
            .collect::<Vec<_>>()
            .join(", ");
        format!(
            "support split planner requires authored two-sided span(s) [{}] before it can create \
             support-owned loops",
            spans
        )
    }
}

pub(super) fn materialize_support_split_cycles(
    brep: &mut BRep,
    support_face: FaceHandle,
    half_edges: &[HalfEdgeHandle],
) -> Result<SupportSplitCyclePlan, ApexBuildError> {
    if let Some(cycle) =
        support_replacement_outer_loop_from_boundary_order(brep, support_face, half_edges)
        && half_edge_cycle_uses_unique_edges(&cycle)
        && validate_half_edge_cycle(
            brep,
            &cycle,
            "boundary-ordered support replacement loop",
            support_face,
        )
        .is_ok()
    {
        return Ok(SupportSplitCyclePlan {
            cycles: vec![cycle],
            authored_boundary_spans: Vec::new(),
            missing_spans: Vec::new(),
            degenerate_unclosed: false,
        });
    }
    let cycles = closed_half_edge_cycles_from_unwired_set(brep, half_edges);
    let covered = cycles.iter().map(Vec::len).sum::<usize>();
    if covered == half_edges.len() {
        return Ok(SupportSplitCyclePlan {
            cycles,
            authored_boundary_spans: Vec::new(),
            missing_spans: Vec::new(),
            degenerate_unclosed: false,
        });
    }
    // No coincidence filter here: a manufactured span between
    // coincident-but-distinct vertices is a degeneracy the caller
    // must refuse (the span authoring rejects it as "degenerate
    // endpoints"), not a request to silently drop. Filtering these
    // spans out was the PK-216 corrupt-accept: the plan read as
    // complete with zero cycles and the partition proceeded on a
    // half-wired shell.
    let missing_spans = missing_support_split_spans(brep, half_edges);
    if missing_spans.is_empty() {
        return Ok(SupportSplitCyclePlan {
            cycles: Vec::new(),
            authored_boundary_spans: Vec::new(),
            missing_spans: Vec::new(),
            degenerate_unclosed: true,
        });
    }
    let mut supplemented = half_edges.to_vec();
    let mut authored_boundary_spans = Vec::new();
    for (start, end) in &missing_spans {
        let half_edge = author_support_split_closing_span(brep, support_face, *start, *end)?;
        supplemented.push(half_edge);
        authored_boundary_spans.push(half_edge);
    }
    let cycles = closed_half_edge_cycles_from_unwired_set(brep, &supplemented);
    let covered = cycles.iter().map(Vec::len).sum::<usize>();
    if covered == supplemented.len() {
        return Ok(SupportSplitCyclePlan {
            cycles,
            authored_boundary_spans,
            missing_spans: Vec::new(),
            degenerate_unclosed: false,
        });
    }
    trace_support_split_incomplete_cycles(
        brep,
        support_face,
        &missing_spans,
        &supplemented,
        &cycles,
    );
    Ok(SupportSplitCyclePlan {
        cycles: Vec::new(),
        authored_boundary_spans,
        missing_spans,
        degenerate_unclosed: false,
    })
}

fn author_support_split_closing_span(
    brep: &mut BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Result<HalfEdgeHandle, ApexBuildError> {
    validate_support_split_span_geometry(brep, support_face, start, end)?;
    let rollback = brep.clone();
    let result = (|| {
        let start_point = brep
            .vertices
            .get(start)
            .map(|vertex| vertex.position)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot author support closing span on \
                         {support_face:?}: stale start vertex {start:?}"
                ),
            })?;
        let end_point = brep
            .vertices
            .get(end)
            .map(|vertex| vertex.position)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot author support closing span on \
                     {support_face:?}: stale end vertex {end:?}"
                ),
            })?;
        let curve = CurveOnEdge::line_segment(start_point, end_point, brep.tolerance.linear())
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot author degenerate support closing \
                     span {start:?}->{end:?} on {support_face:?}"
                ),
            })?;
        let before_edges = brep.edges.len();
        let before_half_edges = brep.half_edges.len();
        let before_faces = brep.faces.len();
        let before_loops = brep.loops.len();
        let (edge, half_edge) = insert_boundary_half_edge(brep, start);
        brep.edges[edge].curve = Some(curve);
        if brep.edges.len() != before_edges + 1
            || brep.half_edges.len() != before_half_edges + 1
            || brep.faces.len() != before_faces
            || brep.loops.len() != before_loops
        {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition support closing span {start:?}->{end:?} on \
                     {support_face:?} violated topology delta dV=0 dE=1 dHE=1 dF=0 dL=0"
                ),
            });
        }
        if !brep.boundary_edges.contains(&edge) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition support closing span {start:?}->{end:?} on \
                     {support_face:?} did not mark its unowned side as boundary"
                ),
            });
        }
        let actual_start = brep.half_edges[half_edge].origin;
        let actual_end = support_split_half_edge_end_vertex(brep, half_edge);
        if actual_start != start || actual_end != Some(end) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition support closing span on {support_face:?} \
                     has wrong orientation: expected {start:?}->{end:?}, got \
                     {actual_start:?}->{actual_end:?}"
                ),
            });
        }
        Ok(half_edge)
    })();
    if result.is_err() {
        *brep = rollback;
    }
    result
}

fn trace_support_split_incomplete_cycles(
    brep: &BRep,
    support_face: FaceHandle,
    missing_spans: &[(VertexHandle, VertexHandle)],
    supplemented: &[HalfEdgeHandle],
    cycles: &[Vec<HalfEdgeHandle>],
) {
    if !crate::direct_apex_trace::enabled() {
        return;
    }
    let spans = missing_spans
        .iter()
        .map(|(start, end)| format!("{start:?}->{end:?}"))
        .collect::<Vec<_>>()
        .join(", ");
    let supplemented_edges = support_cycle_open_chain_diagnostic(brep, supplemented);
    let cycle_sizes = cycles.iter().map(Vec::len).collect::<Vec<_>>();
    let current_loop = brep
        .faces
        .get(support_face)
        .and_then(|face| collect_loop_half_edges(brep, face.outer_loop).ok())
        .map(|loop_edges| support_cycle_open_chain_diagnostic(brep, &loop_edges))
        .unwrap_or_else(|| "unavailable".to_string());
    eprintln!(
        "[direct-apex-trace §14.9/support-split] support {support_face:?} supplemented spans [{spans}] \
         did not close all half-edges; supplemented=[{supplemented_edges}] cycle_sizes={cycle_sizes:?}; \
         current outer loop=[{current_loop}]"
    );
}

fn support_replacement_outer_loop_from_boundary_order(
    brep: &BRep,
    support_face: FaceHandle,
    half_edges: &[HalfEdgeHandle],
) -> Option<Vec<HalfEdgeHandle>> {
    if half_edges.is_empty() {
        return None;
    }
    let boundary_vertices = constructive_topology::collect_outer_loop_vertices(brep, support_face)?;
    let boundary_index = boundary_vertices
        .iter()
        .enumerate()
        .map(|(index, vertex)| (*vertex, index))
        .collect::<std::collections::BTreeMap<_, _>>();
    let mut ordered = half_edges
        .iter()
        .filter_map(|half_edge| {
            let start = brep.half_edges.get(*half_edge)?.origin;
            let end = half_edge_end_vertex(brep, *half_edge)?;
            let start_index = boundary_index.get(&start).copied()?;
            boundary_index.get(&end)?;
            Some((start_index, *half_edge))
        })
        .collect::<Vec<_>>();
    if ordered.len() != half_edges.len() {
        return None;
    }
    ordered.sort_by_key(|(start_index, half_edge)| (*start_index, *half_edge));
    let ordered_half_edges = ordered
        .iter()
        .map(|(_, half_edge)| *half_edge)
        .collect::<Vec<_>>();
    let mut cycle = Vec::new();
    for (index, half_edge) in ordered_half_edges.iter().enumerate() {
        cycle.push(*half_edge);
        let end = half_edge_end_vertex(brep, *half_edge)?;
        let next = ordered_half_edges[(index + 1) % ordered_half_edges.len()];
        let next_start = brep.half_edges.get(next)?.origin;
        if end == next_start {
            continue;
        }
        let path =
            support_outer_loop_forward_path_between_vertices(brep, support_face, end, next_start)?;
        cycle.extend(path);
    }
    Some(cycle)
}

fn half_edge_cycle_uses_unique_edges(cycle: &[HalfEdgeHandle]) -> bool {
    let unique = cycle
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    unique.len() == cycle.len()
}

fn support_outer_loop_forward_path_between_vertices(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    support_outer_loop_paths_between_vertices(brep, support_face, start, end)
        .into_iter()
        .next()
}

fn support_outer_loop_paths_between_vertices(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Vec<Vec<HalfEdgeHandle>> {
    let Some(face) = brep.faces.get(support_face) else {
        return Vec::new();
    };
    let Ok(loop_edges) = collect_loop_half_edges(brep, face.outer_loop) else {
        return Vec::new();
    };
    if loop_edges.iter().any(|half_edge| {
        brep.half_edges
            .get(*half_edge)
            .is_none_or(|half_edge_data| half_edge_data.loop_ != face.outer_loop)
    }) {
        return Vec::new();
    }
    let n = loop_edges.len();
    if n == 0 {
        return Vec::new();
    }
    let mut paths = Vec::new();
    for start_index in 0..n {
        if brep
            .half_edges
            .get(loop_edges[start_index])
            .map(|he| he.origin)
            != Some(start)
        {
            continue;
        }
        let mut path = Vec::new();
        for offset in 0..n {
            let index = (start_index + offset) % n;
            path.push(loop_edges[index]);
            let next_index = (index + 1) % n;
            if brep
                .half_edges
                .get(loop_edges[next_index])
                .map(|he| he.origin)
                == Some(end)
            {
                paths.push(path);
                break;
            }
        }
    }
    paths
}

fn missing_support_split_spans(
    brep: &BRep,
    half_edges: &[HalfEdgeHandle],
) -> Vec<(VertexHandle, VertexHandle)> {
    let mut starts = std::collections::BTreeMap::<VertexHandle, usize>::new();
    let mut ends = std::collections::BTreeMap::<VertexHandle, usize>::new();
    for half_edge in half_edges {
        let Some(start) = brep.half_edges.get(*half_edge).map(|data| data.origin) else {
            continue;
        };
        let Some(end) = half_edge_end_vertex(brep, *half_edge) else {
            continue;
        };
        *starts.entry(start).or_default() += 1;
        *ends.entry(end).or_default() += 1;
    }
    // `then(|| …)`, NOT `then_some(…)`: then_some evaluates its
    // argument eagerly, and `count - end_count` underflows (debug
    // panic) whenever the guard is false (measured: nearly-uniform
    // pyramid4 radii panicked the whole fillet op here).
    let mut open_starts = starts
        .iter()
        .filter_map(|(vertex, count)| {
            let end_count = ends.get(vertex).copied().unwrap_or_default();
            (count > &end_count).then(|| (*vertex, count - end_count))
        })
        .collect::<Vec<_>>();
    let mut open_ends = ends
        .iter()
        .filter_map(|(vertex, count)| {
            let start_count = starts.get(vertex).copied().unwrap_or_default();
            (count > &start_count).then(|| (*vertex, count - start_count))
        })
        .collect::<Vec<_>>();
    let mut spans = Vec::new();
    for (end_vertex, end_count) in open_ends.drain(..) {
        for _ in 0..end_count {
            let Some((start_vertex, start_count)) = open_starts.first_mut() else {
                return spans;
            };
            spans.push((end_vertex, *start_vertex));
            *start_count -= 1;
            if *start_count == 0 {
                open_starts.remove(0);
            }
        }
    }
    spans
}

fn validate_support_split_span_geometry(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Result<(), ApexBuildError> {
    let start_point = brep.vertices[start].position;
    let end_point = brep.vertices[end].position;
    let Some(SurfaceGeometry::Plane { origin, normal, .. }) = brep
        .faces
        .get(support_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author support split span {start:?}->{end:?}: \
                 support face {support_face:?} has no plane surface"
            ),
        });
    };
    let tol = brep.tolerance.linear().max(1.0e-9);
    let start_distance = (start_point - *origin).dot(normal).abs();
    let end_distance = (end_point - *origin).dot(normal).abs();
    if start_distance > tol || end_distance > tol {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author support split span {start:?}->{end:?}: \
                 endpoints are off support face {support_face:?} by {:.3e}/{:.3e}",
                start_distance, end_distance
            ),
        });
    }
    if start_point.distance_to(&end_point) <= tol {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author support split span {start:?}->{end:?}: \
                 degenerate endpoints"
            ),
        });
    }
    Ok(())
}

#[cfg(test)]
mod support_split_cycle_plan_honesty {
    use super::SupportSplitCyclePlan;

    /// PK-216 regression pin: a plan with zero cycles and zero
    /// missing spans (balanced-but-unchainable trim sides, the
    /// shared-corner degenerate state) must NOT read as complete.
    #[test]
    fn degenerate_unclosed_plan_is_not_complete() {
        let plan = SupportSplitCyclePlan {
            cycles: Vec::new(),
            authored_boundary_spans: Vec::new(),
            missing_spans: Vec::new(),
            degenerate_unclosed: true,
        };
        assert!(!plan.is_complete());
        assert!(plan.diagnostic().contains("degenerate chaining"));
    }

    #[test]
    fn closed_cycle_plan_is_complete() {
        let plan = SupportSplitCyclePlan {
            cycles: vec![Vec::new()],
            authored_boundary_spans: Vec::new(),
            missing_spans: Vec::new(),
            degenerate_unclosed: false,
        };
        assert!(plan.is_complete());
    }
}
