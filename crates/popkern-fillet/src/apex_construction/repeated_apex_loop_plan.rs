//! Classifier and typed plan for repeated-apex cap loops.

use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

pub(super) fn plan_repeated_apex_loop_normalization(
    brep: &BRep,
    face: FaceHandle,
    apex_vertex: VertexHandle,
    apex_vertex_handles: &[VertexHandle],
    n_fillets: usize,
) -> Option<RepeatedApexLoopPlan> {
    let face_data = brep.faces.get(face)?;
    let loop_handle = face_data.outer_loop;
    let loop_data = brep.loops.get(loop_handle)?;
    let start = loop_data.half_edge;
    let mut cursor = start;
    let mut vertices = Vec::new();
    let mut half_edges = Vec::new();
    let mut guard = 0usize;
    let max_iters = brep.half_edges.len() + 1;
    while let Some(half_edge) = brep.half_edges.get(cursor) {
        half_edges.push(cursor);
        vertices.push(half_edge.origin);
        cursor = half_edge.next;
        guard += 1;
        if cursor == start || guard >= max_iters {
            break;
        }
    }

    let apex_positions: Vec<usize> = vertices
        .iter()
        .enumerate()
        .filter_map(|(index, vertex)| (*vertex == apex_vertex).then_some(index))
        .collect();
    if apex_positions.len() < 2 || vertices.is_empty() {
        return None;
    }

    let mut arcs = Vec::with_capacity(apex_positions.len());
    for (arc_index, &start_pos) in apex_positions.iter().enumerate() {
        let end_pos = apex_positions[(arc_index + 1) % apex_positions.len()];
        let interior = arc_vertices_between(&vertices, start_pos, end_pos);
        let kind = classify_apex_arc(&interior, apex_vertex_handles, n_fillets);
        arcs.push(RepeatedApexArc {
            index: arc_index,
            start_pos,
            end_pos,
            interior,
            kind,
        });
    }

    // Welded corners are encoded by HANDLE IDENTITY: the authoring
    // layer aliases T_b[k] and T_a[k+1] to one vertex (PK-219).
    let welded_corners: Vec<bool> = (0..n_fillets)
        .map(|k| apex_vertex_handles[2 * k + 1] == apex_vertex_handles[2 * ((k + 1) % n_fillets)])
        .collect();

    let mut replacement_spans = Vec::with_capacity(n_fillets);
    for k in 0..n_fillets {
        let t_a = apex_vertex_handles[2 * k];
        let t_b = apex_vertex_handles[2 * k + 1];
        replacement_spans.push(FilletEndReplacementSpan {
            fillet_index: k,
            t_a_vertex: t_a,
            t_b_vertex: t_b,
            t_a_positions: vertex_positions(&vertices, t_a),
            t_b_positions: vertex_positions(&vertices, t_b),
        });
    }

    let normalized_boundary = if arcs
        .iter()
        .all(|arc| !matches!(arc.kind, RepeatedApexArcKind::Unresolved))
    {
        Some(plan_normalized_cap_boundary(
            &arcs,
            &replacement_spans,
            &welded_corners,
            n_fillets,
        )?)
    } else {
        None
    };

    Some(RepeatedApexLoopPlan {
        face,
        loop_handle,
        half_edges,
        vertices,
        apex_positions,
        arcs,
        replacement_spans,
        welded_corners,
        normalized_boundary,
    })
}

pub(super) struct RepeatedApexLoopPlan {
    pub(super) face: FaceHandle,
    pub(super) loop_handle: popkern_core::topology::LoopHandle,
    pub(super) half_edges: Vec<popkern_core::topology::HalfEdgeHandle>,
    pub(super) vertices: Vec<VertexHandle>,
    pub(super) apex_positions: Vec<usize>,
    pub(super) arcs: Vec<RepeatedApexArc>,
    pub(super) replacement_spans: Vec<FilletEndReplacementSpan>,
    /// `welded_corners[k]`: corner between fillet k and k+1 is a
    /// welded shared corner (aliased T handles, PK-219).
    pub(super) welded_corners: Vec<bool>,
    pub(super) normalized_boundary: Option<Vec<NormalizedCapBoundaryEdge>>,
}

impl RepeatedApexLoopPlan {
    pub(super) fn diagnostic_summary(&self, n_fillets: usize) -> String {
        let mut face_edge_arcs = Vec::new();
        let mut base_arcs = Vec::new();
        let mut unresolved_arcs = Vec::new();

        for arc in &self.arcs {
            match &arc.kind {
                RepeatedApexArcKind::FaceEdge { fillet_index } => {
                    face_edge_arcs.push(format!(
                        "FaceEdge[{fillet_index}] arc[{}] V@{}->V@{}: {}",
                        arc.index,
                        arc.start_pos,
                        arc.end_pos,
                        format_vertices(&arc.interior)
                    ));
                }
                RepeatedApexArcKind::SharedCorner { corner_index } => {
                    face_edge_arcs.push(format!(
                        "SharedCorner[{corner_index}] arc[{}] V@{}->V@{}: {}",
                        arc.index,
                        arc.start_pos,
                        arc.end_pos,
                        format_vertices(&arc.interior)
                    ));
                }
                RepeatedApexArcKind::Base => {
                    base_arcs.push(format!(
                        "arc[{}] V@{}->V@{}",
                        arc.index, arc.start_pos, arc.end_pos
                    ));
                }
                RepeatedApexArcKind::Unresolved => {
                    unresolved_arcs.push(format!(
                        "arc[{}] V@{}->V@{}: {}",
                        arc.index,
                        arc.start_pos,
                        arc.end_pos,
                        format_vertices(&arc.interior)
                    ));
                }
            }
        }

        let replacement_spans: Vec<String> = self
            .replacement_spans
            .iter()
            .map(|span| {
                format!(
                    "FilletEnd[{}] replacement {:?}@{:?}->{:?}@{:?}",
                    span.fillet_index,
                    span.t_a_vertex,
                    span.t_a_positions,
                    span.t_b_vertex,
                    span.t_b_positions
                )
            })
            .collect();

        let normalized_boundary = self
            .normalized_boundary
            .as_ref()
            .map(|edges| {
                edges
                    .iter()
                    .map(|edge| match edge {
                        NormalizedCapBoundaryEdge::FilletEnd { fillet_index } => {
                            format!("FilletEnd[{fillet_index}]")
                        }
                        NormalizedCapBoundaryEdge::FaceEdge { fillet_index } => {
                            format!("FaceEdge[{fillet_index}]")
                        }
                    })
                    .collect::<Vec<_>>()
                    .join(" -> ")
            })
            .unwrap_or_else(|| "<unplanned>".to_string());

        format!(
            "normalization required on {:?}/{:?}: candidate_ok={}, {} half-edges, {} vertices, {} apex occurrences; face-edge arcs=[{}]; base arcs=[{}]; unresolved arcs=[{}]; replacement spans=[{}]; normalized cap boundary=[{}]",
            self.face,
            self.loop_handle,
            self.is_commercial_grade_candidate(n_fillets),
            self.half_edges.len(),
            self.vertices.len(),
            self.apex_positions.len(),
            face_edge_arcs.join("; "),
            base_arcs.join("; "),
            unresolved_arcs.join("; "),
            replacement_spans.join("; "),
            normalized_boundary
        )
    }

    pub(super) fn is_commercial_grade_candidate(&self, n_fillets: usize) -> bool {
        let welded = self.welded_corners.iter().filter(|w| **w).count();
        self.apex_positions.len() == n_fillets + 1
            && self
                .arcs
                .iter()
                .filter(|arc| matches!(arc.kind, RepeatedApexArcKind::FaceEdge { .. }))
                .count()
                == n_fillets - welded
            && self
                .arcs
                .iter()
                .filter(|arc| matches!(arc.kind, RepeatedApexArcKind::SharedCorner { .. }))
                .count()
                == welded
            && self
                .arcs
                .iter()
                .filter(|arc| matches!(arc.kind, RepeatedApexArcKind::Base))
                .count()
                == 1
            && self
                .arcs
                .iter()
                .all(|arc| !matches!(arc.kind, RepeatedApexArcKind::Unresolved))
            && self
                .replacement_spans
                .iter()
                .all(|span| span.t_a_positions.len() == 1 && span.t_b_positions.len() == 1)
            && self.normalized_boundary.as_ref().is_some_and(|boundary| {
                boundary.len() == 2 * n_fillets - self.welded_corners.iter().filter(|w| **w).count()
            })
    }
}

pub(super) struct RepeatedApexArc {
    pub(super) index: usize,
    pub(super) start_pos: usize,
    pub(super) end_pos: usize,
    pub(super) interior: Vec<VertexHandle>,
    pub(super) kind: RepeatedApexArcKind,
}

pub(super) struct FilletEndReplacementSpan {
    fillet_index: usize,
    t_a_vertex: VertexHandle,
    t_b_vertex: VertexHandle,
    t_a_positions: Vec<usize>,
    t_b_positions: Vec<usize>,
}

pub(super) enum NormalizedCapBoundaryEdge {
    FilletEnd { fillet_index: usize },
    FaceEdge { fillet_index: usize },
}

#[derive(Clone, Copy, Debug)]
pub(super) enum RepeatedApexArcKind {
    FaceEdge {
        fillet_index: usize,
    },
    /// Welded shared corner (PK-219): the arc's single interior
    /// vertex is simultaneously T_b[corner_index] and
    /// T_a[corner_index + 1].
    SharedCorner {
        corner_index: usize,
    },
    Base,
    Unresolved,
}

fn plan_normalized_cap_boundary(
    arcs: &[RepeatedApexArc],
    replacement_spans: &[FilletEndReplacementSpan],
    welded_corners: &[bool],
    n_fillets: usize,
) -> Option<Vec<NormalizedCapBoundaryEdge>> {
    if replacement_spans.len() != n_fillets {
        return None;
    }
    // Each BRIDGED corner needs exactly one FaceEdge arc; each welded
    // corner exactly one SharedCorner arc.
    let mut face_edge_seen = vec![false; n_fillets];
    let mut shared_corner_seen = vec![false; n_fillets];
    for arc in arcs {
        match arc.kind {
            RepeatedApexArcKind::FaceEdge { fillet_index } => {
                let seen = face_edge_seen.get_mut(fillet_index)?;
                if *seen || welded_corners[fillet_index] {
                    return None;
                }
                *seen = true;
            }
            RepeatedApexArcKind::SharedCorner { corner_index } => {
                let seen = shared_corner_seen.get_mut(corner_index)?;
                if *seen || !welded_corners[corner_index] {
                    return None;
                }
                *seen = true;
            }
            RepeatedApexArcKind::Base | RepeatedApexArcKind::Unresolved => {}
        }
    }
    for k in 0..n_fillets {
        if welded_corners[k] != shared_corner_seen[k] || welded_corners[k] == face_edge_seen[k] {
            return None;
        }
    }

    let mut boundary = Vec::with_capacity(2 * n_fillets);
    for (k, welded) in welded_corners.iter().enumerate().take(n_fillets) {
        boundary.push(NormalizedCapBoundaryEdge::FilletEnd { fillet_index: k });
        if !welded {
            boundary.push(NormalizedCapBoundaryEdge::FaceEdge { fillet_index: k });
        }
    }
    Some(boundary)
}

fn classify_apex_arc(
    interior: &[VertexHandle],
    apex_vertex_handles: &[VertexHandle],
    n_fillets: usize,
) -> RepeatedApexArcKind {
    let planned: Vec<Vec<(usize, bool)>> = interior
        .iter()
        .map(|&vertex| planned_apex_vertex_roles(vertex, apex_vertex_handles, n_fillets))
        .filter(|roles| !roles.is_empty())
        .collect();
    if planned.is_empty() {
        return RepeatedApexArcKind::Base;
    }
    if let [roles] = planned.as_slice() {
        // One planned interior vertex: a welded shared corner iff the
        // SAME handle holds both T_b[k] and T_a[k+1] roles (PK-219).
        for &(k, is_a) in roles {
            if !is_a && roles.contains(&((k + 1) % n_fillets, true)) {
                return RepeatedApexArcKind::SharedCorner { corner_index: k };
            }
        }
    }
    if let [first_roles, second_roles] = planned.as_slice() {
        for &(first_k, first_is_a) in first_roles {
            let wanted = ((first_k + 1) % n_fillets, true);
            if !first_is_a && second_roles.contains(&wanted) {
                return RepeatedApexArcKind::FaceEdge {
                    fillet_index: first_k,
                };
            }
        }
    }
    RepeatedApexArcKind::Unresolved
}

/// ALL planned roles a loop vertex matches — a welded corner handle
/// occupies both its T_b slot and the next sector's T_a slot.
fn planned_apex_vertex_roles(
    vertex: VertexHandle,
    apex_vertex_handles: &[VertexHandle],
    n_fillets: usize,
) -> Vec<(usize, bool)> {
    let mut roles = Vec::new();
    for k in 0..n_fillets {
        if apex_vertex_handles[2 * k] == vertex {
            roles.push((k, true));
        }
        if apex_vertex_handles[2 * k + 1] == vertex {
            roles.push((k, false));
        }
    }
    roles
}

fn arc_vertices_between(
    vertices: &[VertexHandle],
    start_pos: usize,
    end_pos: usize,
) -> Vec<VertexHandle> {
    let mut out = Vec::new();
    if vertices.is_empty() {
        return out;
    }
    let mut pos = (start_pos + 1) % vertices.len();
    while pos != end_pos {
        out.push(vertices[pos]);
        pos = (pos + 1) % vertices.len();
    }
    out
}

fn vertex_positions(vertices: &[VertexHandle], needle: VertexHandle) -> Vec<usize> {
    vertices
        .iter()
        .enumerate()
        .filter_map(|(index, vertex)| (*vertex == needle).then_some(index))
        .collect()
}

fn format_vertices(vertices: &[VertexHandle]) -> String {
    if vertices.is_empty() {
        return "<empty>".to_string();
    }
    vertices
        .iter()
        .map(|vertex| format!("{vertex:?}"))
        .collect::<Vec<_>>()
        .join(" -> ")
}
