use super::super::partition_cap_boundary_rewrite_diagnostics::{
    CapBoundaryRewriteDiagnosticInput, cap_boundary_rewrite_diagnostic_label,
};
use super::{AdjacentLoopRewriteRequirement, CapBoundaryRewritePlan};

impl CapBoundaryRewritePlan {
    pub(in crate::apex_construction) fn diagnostic_label(&self) -> String {
        cap_boundary_rewrite_diagnostic_label(CapBoundaryRewriteDiagnosticInput {
            required: self.required,
            cap_face: self.cap_face,
            cap_loop: self.cap_loop,
            fillet_end_count: self.fillet_end_count,
            face_edge_count: self.face_edge_count,
            current_cap_cycle: &self.current_cap_cycle,
            reauthored_cap_cycle: &self.reauthored_cap_cycle,
            adjacent_loop_rewrites: self
                .adjacent_loop_rewrites
                .iter()
                .map(AdjacentLoopRewriteRequirement::diagnostic_label)
                .collect(),
        })
    }
}

impl AdjacentLoopRewriteRequirement {
    pub(in crate::apex_construction) fn diagnostic_label(&self) -> String {
        format!(
            "stolen_cap_side={:?}, replacement_neighbor_side={:?}, neighbor_loop={:?}, \
             neighbor_face={:?}, rewrite_kind={}, current_neighbor_cycle={:?}, \
             reauthored_neighbor_cycle={:?}: {}",
            self.stolen_cap_side,
            self.replacement_neighbor_side,
            self.neighbor_loop,
            self.neighbor_face,
            self.rewrite_kind.diagnostic_label(),
            self.current_neighbor_cycle,
            self.reauthored_neighbor_cycle,
            self.rewrite_kind.reason()
        )
    }
}
