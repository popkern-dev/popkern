//! Diagnostics for cap-boundary rewrite plans.

use super::*;

pub(super) struct CapBoundaryRewriteDiagnosticInput<'a> {
    pub(super) required: bool,
    pub(super) cap_face: FaceHandle,
    pub(super) cap_loop: LoopHandle,
    pub(super) fillet_end_count: usize,
    pub(super) face_edge_count: usize,
    pub(super) current_cap_cycle: &'a [HalfEdgeHandle],
    pub(super) reauthored_cap_cycle: &'a [HalfEdgeHandle],
    pub(super) adjacent_loop_rewrites: Vec<String>,
}

pub(super) fn cap_boundary_rewrite_diagnostic_label(
    input: CapBoundaryRewriteDiagnosticInput<'_>,
) -> String {
    if input.required {
        format!(
            "reauthor the whole cap boundary on {:?}/{:?} across {} FilletEnd seam(s) and {} \
             FaceEdge seam(s) so every cap-side half-edge remains opposite its support/cylinder-side \
             twin; current_cap_cycle={:?}; reauthored_cap_cycle={:?}; adjacent_loop_rewrites=[{}]",
            input.cap_face,
            input.cap_loop,
            input.fillet_end_count,
            input.face_edge_count,
            input.current_cap_cycle,
            input.reauthored_cap_cycle,
            input.adjacent_loop_rewrites.join("; ")
        )
    } else {
        format!(
            "reuse current cap boundary on {:?}/{:?} with {} FilletEnd seam(s) and {} FaceEdge \
             seam(s); cap_cycle={:?}",
            input.cap_face,
            input.cap_loop,
            input.fillet_end_count,
            input.face_edge_count,
            input.current_cap_cycle
        )
    }
}
