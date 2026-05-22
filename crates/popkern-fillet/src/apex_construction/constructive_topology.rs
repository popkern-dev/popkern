//! Loop walking and topology validation helpers for constructive apex building.

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use super::error::ApexBuildError;

pub(super) fn find_outgoing_half_edge(brep: &BRep, vertex: VertexHandle) -> Option<HalfEdgeHandle> {
    for (he, data) in &brep.half_edges {
        if data.origin == vertex {
            return Some(he);
        }
    }
    None
}

pub(super) fn find_two_origins_on_face(
    brep: &BRep,
    face: FaceHandle,
    vertex_a: VertexHandle,
    vertex_b: VertexHandle,
) -> Option<(HalfEdgeHandle, HalfEdgeHandle)> {
    let face_data = brep.faces.get(face)?;
    let mut loop_handles = Vec::with_capacity(1 + face_data.inner_loops.len());
    loop_handles.push(face_data.outer_loop);
    loop_handles.extend(face_data.inner_loops.iter().copied());

    for loop_handle in loop_handles {
        let start = brep.loops.get(loop_handle)?.half_edge;
        let mut cursor = start;
        let cap = brep.half_edges.len() + 1;
        let mut iters = 0;
        let mut he_a: Option<HalfEdgeHandle> = None;
        let mut he_b: Option<HalfEdgeHandle> = None;
        loop {
            let curr = brep.half_edges.get(cursor)?;
            if he_a.is_none() && curr.origin == vertex_a {
                he_a = Some(cursor);
            }
            if he_b.is_none() && curr.origin == vertex_b {
                he_b = Some(cursor);
            }
            if he_a.is_some() && he_b.is_some() {
                return Some((he_a?, he_b?));
            }
            cursor = curr.next;
            iters += 1;
            if cursor == start || iters > cap {
                break;
            }
        }
    }
    None
}

pub(super) fn find_two_origins_on_all_faces(
    brep: &BRep,
    vertex_a: VertexHandle,
    vertex_b: VertexHandle,
) -> Vec<(FaceHandle, HalfEdgeHandle, HalfEdgeHandle)> {
    let mut out = Vec::new();
    for (face, _) in &brep.faces {
        if let Some((he_a, he_b)) = find_two_origins_on_face(brep, face, vertex_a, vertex_b) {
            out.push((face, he_a, he_b));
        }
    }
    out
}

pub(super) fn loop_segment_contains_origin(
    brep: &BRep,
    start: HalfEdgeHandle,
    end: HalfEdgeHandle,
    vertex: VertexHandle,
) -> bool {
    let mut cursor = start;
    for _ in 0..(brep.half_edges.len() + 1) {
        if cursor == end {
            return false;
        }
        let Some(he) = brep.half_edges.get(cursor) else {
            return false;
        };
        if he.origin == vertex {
            return true;
        }
        cursor = he.next;
    }
    false
}

pub(super) fn collect_outer_loop_vertices(
    brep: &BRep,
    face: FaceHandle,
) -> Option<Vec<VertexHandle>> {
    let loop_handle = brep.faces.get(face)?.outer_loop;
    let start = brep.loops.get(loop_handle)?.half_edge;
    let mut cursor = start;
    let mut vertices = Vec::new();
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        vertices.push(half_edge.origin);
        cursor = half_edge.next;
        if cursor == start {
            return Some(vertices);
        }
    }
    None
}

pub(super) fn cyclic_vertices_match(actual: &[VertexHandle], expected: &[VertexHandle]) -> bool {
    if actual.len() != expected.len() {
        return false;
    }
    if actual.is_empty() {
        return true;
    }
    (0..actual.len())
        .any(|offset| (0..actual.len()).all(|i| actual[(offset + i) % actual.len()] == expected[i]))
}

fn cyclic_vertices_match_with_coalesced_expected(
    brep: &BRep,
    actual: &[VertexHandle],
    expected: &[VertexHandle],
) -> bool {
    if actual.len() != expected.len() {
        return false;
    }
    if actual.is_empty() {
        return true;
    }
    (0..actual.len()).any(|offset| {
        (0..actual.len()).all(|i| {
            let expected = expected_representative_after_coalescence(brep, expected, i);
            actual[(offset + i) % actual.len()] == expected
        })
    })
}

fn expected_representative_after_coalescence(
    brep: &BRep,
    expected: &[VertexHandle],
    index: usize,
) -> VertexHandle {
    let vertex = expected[index];
    if brep.vertices.contains_key(vertex) {
        return vertex;
    }
    for step in 1..=expected.len() {
        let candidate = expected[(index + step) % expected.len()];
        if brep.vertices.contains_key(candidate) {
            return candidate;
        }
    }
    vertex
}

pub(super) fn assert_normalized_cap_boundary(
    brep: &BRep,
    cap_face: FaceHandle,
    apex_vertex_handles: &[VertexHandle],
    n_fillets: usize,
) -> Result<(), ApexBuildError> {
    let vertices = collect_outer_loop_vertices(brep, cap_face).ok_or_else(|| {
        ApexBuildError::PlanMismatch {
            reason: format!(
                "repeated-apex normalization: cap face {cap_face:?} has no walkable outer loop"
            ),
        }
    })?;

    let mut expected = Vec::with_capacity(2 * n_fillets);
    for k in 0..n_fillets {
        expected.push(apex_vertex_handles[2 * k]);
        expected.push(apex_vertex_handles[2 * k + 1]);
    }

    if vertices.len() == expected.len()
        && (cyclic_vertices_match(&vertices, &expected)
            || cyclic_vertices_match_with_coalesced_expected(brep, &vertices, &expected))
    {
        return Ok(());
    }

    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "repeated-apex normalization produced unexpected cap boundary: got {:?}, expected cyclic {:?}",
            vertices, expected
        ),
    })
}

pub(super) fn validate_repeated_apex_normalization_structure(
    brep: &BRep,
) -> Result<(), ApexBuildError> {
    let result = brep
        .validate_euler_poincare()
        .and_then(|_| brep.validate_half_edge_pairing())
        .and_then(|_| brep.validate_half_edge_cycles())
        .and_then(|_| brep.validate_vertex_consistency())
        .and_then(|_| brep.validate_manifold_edges())
        .and_then(|_| brep.validate_vertex_valence());
    if let Err(error) = &result
        && crate::direct_apex_trace::enabled()
    {
        eprintln!(
            "[direct-apex-trace §14.9/normalization-validation] failure {error:?}; {}",
            topology_count_diagnostic(brep)
        );
        trace_repeated_apex_loop_walks(brep);
    }
    result.map_err(|e| ApexBuildError::PlanMismatch {
        reason: format!(
            "repeated-apex normalization produced invalid topology: {e:?}; {}",
            topology_count_diagnostic(brep)
        ),
    })
}

fn trace_repeated_apex_loop_walks(brep: &BRep) {
    for (loop_handle, loop_data) in &brep.loops {
        let mut cursor = loop_data.half_edge;
        let mut entries = Vec::new();
        for _ in 0..=brep.half_edges.len().min(40) {
            let Some(half_edge) = brep.half_edges.get(cursor) else {
                entries.push(format!("{cursor:?}:missing"));
                break;
            };
            let role = brep
                .edges
                .get(half_edge.edge)
                .map(|edge| format!("{:?}", edge.role))
                .unwrap_or_else(|| "missing-edge".to_string());
            entries.push(format!(
                "{cursor:?}(loop={:?} edge={:?} role={} next={:?} prev={:?})",
                half_edge.loop_, half_edge.edge, role, half_edge.next, half_edge.prev
            ));
            cursor = half_edge.next;
            if cursor == loop_data.half_edge {
                break;
            }
        }
        eprintln!(
            "[direct-apex-trace §14.9/normalization-validation] loop {loop_handle:?} face={:?} walk=[{}]",
            loop_data.face,
            entries.join(" -> ")
        );
    }
}

fn topology_count_diagnostic(brep: &BRep) -> String {
    let v = brep.vertices.len() as i64;
    let e = brep.edges.len() as i64;
    let f = brep.faces.len() as i64;
    let l = brep.loops.len() as i64;
    let s = brep.shells.len() as i64;
    let genus_zero_lhs = v - e + f - (l - f);
    let genus_zero_rhs = 2 * s;
    let any_genus_lhs = v - e + 2 * f - l;
    let two_genus = 2 * s - any_genus_lhs;
    format!(
        "topology counts V={v} E={e} F={f} L={l} S={s}; genus0 lhs={genus_zero_lhs} rhs={genus_zero_rhs}; any-genus lhs={any_genus_lhs} two_genus={two_genus}"
    )
}

pub(super) fn production_faces_containing_vertex(
    brep: &BRep,
    vertex: VertexHandle,
) -> Vec<FaceHandle> {
    brep.faces
        .iter()
        .filter_map(|(face_handle, face)| {
            face.surface.as_ref()?;
            collect_outer_loop_vertices(brep, face_handle)
                .is_some_and(|vertices| vertices.contains(&vertex))
                .then_some(face_handle)
        })
        .collect()
}
