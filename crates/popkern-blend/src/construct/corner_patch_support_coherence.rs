//! Stage-B1 coherence-validation diagnostic emission for corner-patch
//! support completion (nested module of
//! `corner_patch_support_completion` for the file-size budget).

use popkern_core::topology::BRep;

use super::surfaced_cap_face_set;
use crate::types::{ApexConstructionDiagnostic, CornerPatchInfo};

/// Stage B1 of `docs/decisions/2026-06-21-edge-curve-face-boundary-
/// coherence-invariants.md`: at the end of corner-patch support
/// completion, scope `validate_coherence` to the cap faces and emit
/// the result as an `ApexConstructionDiagnostic::CoherenceValidation`.
/// Stage B1 is observational — the report is recorded, no behavior
/// gating yet.
pub(super) fn emit_coherence_validation_diagnostic(
    brep: &BRep,
    corner_patches: &[CornerPatchInfo],
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
) {
    let cap_faces = surfaced_cap_face_set(brep, corner_patches);
    let report = brep.validate_coherence(&popkern_core::topology::CoherenceValidationScope::Faces(
        cap_faces,
    ));
    diagnostics.push(ApexConstructionDiagnostic::CoherenceValidation(
        crate::types::CoherenceValidationDiagnostic {
            stage: "corner_patch_support_completion".to_string(),
            report,
        },
    ));
}
