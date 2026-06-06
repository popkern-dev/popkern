//! Structured timing facts for the apply-blend construction pipeline.

use std::time::Instant;

use crate::types::{ApexConstructionDiagnostic, ExactBoundaryCapDiagnostic};

#[derive(Debug, Clone, Copy, Default)]
pub(super) struct ApplyBlendPipelineTiming {
    pub(super) perpendicular_split_ms: u128,
    pub(super) face_trim_ms: u128,
    pub(super) vertex_blend_ms: u128,
    pub(super) post_pipeline_ms: u128,
}

impl ApplyBlendPipelineTiming {
    fn total_ms(self) -> u128 {
        self.perpendicular_split_ms
            + self.face_trim_ms
            + self.vertex_blend_ms
            + self.post_pipeline_ms
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn diagnostic(
        self,
        edge_count: usize,
        corner_patch_count: usize,
        warning_count: usize,
        direct_apex_built_face_count: usize,
        late_wire_job_count: usize,
        support_seam_job_count: usize,
    ) -> ApexConstructionDiagnostic {
        ApexConstructionDiagnostic::ExactBoundaryCap(ExactBoundaryCapDiagnostic {
            stage: "fillet.apply_blend.pipeline.timing".to_string(),
            code: "fillet.apply_blend.pipeline.fact".to_string(),
            detail_code: Some(format!(
                "timing_total_ms={};timing_perpendicular_split_ms={};timing_face_trim_ms={};timing_vertex_blend_ms={};timing_post_pipeline_ms={};edge_count={};corner_patch_count={};warning_count={};direct_apex_built_face_count={};late_wire_job_count={};support_seam_job_count={}",
                self.total_ms(),
                self.perpendicular_split_ms,
                self.face_trim_ms,
                self.vertex_blend_ms,
                self.post_pipeline_ms,
                edge_count,
                corner_patch_count,
                warning_count,
                direct_apex_built_face_count,
                late_wire_job_count,
                support_seam_job_count,
            )),
            sector_name: "fillet_apply_blend_pipeline".to_string(),
            edge: None,
            support_half_edge: None,
            support_face: None,
            residual: Some(self.total_ms() as f64),
            outer_side: None,
            sample_parameter: None,
            tolerance: None,
            shared_fact_ids: Vec::new(),
        })
    }
}

pub(super) fn timed<T>(f: impl FnOnce() -> T) -> (T, u128) {
    let start = Instant::now();
    let value = f();
    (value, start.elapsed().as_millis())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diagnostic_exposes_stable_pipeline_timing_fields() {
        let timing = ApplyBlendPipelineTiming {
            perpendicular_split_ms: 1,
            face_trim_ms: 2,
            vertex_blend_ms: 3,
            post_pipeline_ms: 4,
        };

        let ApexConstructionDiagnostic::ExactBoundaryCap(cap) =
            timing.diagnostic(5, 6, 7, 8, 9, 10)
        else {
            panic!("pipeline timing should emit an exact-boundary diagnostic fact");
        };

        assert_eq!(cap.stage, "fillet.apply_blend.pipeline.timing");
        assert_eq!(cap.code, "fillet.apply_blend.pipeline.fact");
        assert_eq!(cap.sector_name, "fillet_apply_blend_pipeline");
        assert_eq!(cap.residual, Some(10.0));
        let detail = cap.detail_code.expect("timing detail label");
        assert!(detail.contains("timing_total_ms=10"));
        assert!(detail.contains("timing_vertex_blend_ms=3"));
        assert!(detail.contains("support_seam_job_count=10"));
    }
}
