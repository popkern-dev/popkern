//! Env-gated (the fan-quality trace topic) pipeline-stage
//! observability for committed cap fans: worst live-seam g0/g1 of
//! the corner-patch faces. This tracer localized the pentagon5_blunt
//! post-pipeline pcurve clobber to `corner_patch_support_completion`
//! (PK-087). Split from `setback_vertex_topology` for the file-size
//! budget.

use popkern_core::topology::{BRep, FaceHandle};

use super::setback_fan_seam_refresh;

/// Print one pipeline-stage fan-quality row when the `fan-quality`
/// trace topic is enabled.
pub(crate) fn fan_quality_stage_probe(
    label: &str,
    brep: &BRep,
    corner_patches: &[crate::types::CornerPatchInfo],
) {
    if !popkern_core::trace::trace_enabled("fan-quality") {
        return;
    }
    let faces: Vec<FaceHandle> = corner_patches.iter().map(|cp| cp.face).collect();
    let (g0, g1, g2) = setback_fan_seam_refresh::worst_fan_seam_quality(brep, &faces);
    eprintln!("  [fan-quality] {label}: g0={g0:.3e} g1={g1:.3e} g2={g2:.3e}");
}
