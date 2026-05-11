//! Exposure of twin-side blockers for cylinder bridge resources.

use crate::error::CoreError;
use crate::topology::BRep;

use super::apply_cylinder_resource_remap::CylinderBridgeResourceRemap;
use super::apply_exposure_candidates::{
    ExposureConstructionCandidate, ExposureConstructionFamily, ExposureConstructionFeasibility,
    ExposureRollbackBoundary, ExposureSemanticOwnership,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct CylinderResourceBlockerExposureResult {
    spans: usize,
    half_edges: usize,
    exposed_spans: usize,
    exposed_half_edges: usize,
}

impl CylinderResourceBlockerExposureResult {
    pub(super) fn candidate_count(self) -> usize {
        1
    }

    pub(super) fn candidate_gate_complete(self) -> bool {
        self.candidate().mutation_gate_complete()
    }

    pub(super) fn exposed_spans(self) -> usize {
        self.exposed_spans
    }

    pub(super) fn diagnostic_label(self) -> String {
        let candidate = self.candidate();
        format!(
            "cylinder_resource_blocker_exposure=spans={} half_edges={} exposed_spans={} exposed_half_edges={} gate_complete:{} {}",
            self.spans,
            self.half_edges,
            self.exposed_spans,
            self.exposed_half_edges,
            self.candidate_gate_complete(),
            candidate.diagnostic_label()
        )
    }

    fn candidate(self) -> ExposureConstructionCandidate {
        ExposureConstructionCandidate {
            family: ExposureConstructionFamily::CylinderResourceBlockerRetire,
            feasibility: if self.spans == 0 {
                ExposureConstructionFeasibility::NotApplicable
            } else {
                ExposureConstructionFeasibility::ExecutableByExistingOperator
            },
            reason: if self.spans == 0 {
                "cylinder-resource-blocker-exposure-not-required"
            } else {
                "cylinder-resource-blocker-boundary-fragment-retire-ready"
            },
            spans: self.spans,
            half_edges: self.half_edges,
            expected_edges: self.exposed_spans as i64,
            expected_faces: 0,
            expected_loops: 0,
            semantic_ownership: if self.spans == 0 {
                ExposureSemanticOwnership::None
            } else {
                ExposureSemanticOwnership::CylinderResourceBoundarySideExposed
            },
            rollback_boundary: if self.spans == 0 {
                ExposureRollbackBoundary::NoMutation
            } else {
                ExposureRollbackBoundary::CandidateBrepBeforeCommit
            },
        }
    }
}

pub(super) fn expose_cylinder_resource_blockers(
    brep: &mut BRep,
    remap: &CylinderBridgeResourceRemap,
) -> Result<CylinderResourceBlockerExposureResult, CoreError> {
    let spans = remap.blocking_spans();
    let mut exposed_spans = 0usize;
    let mut exposed_half_edges = 0usize;
    for span in spans {
        if !span.contiguous {
            return Err(CoreError::PreconditionFailed(format!(
                "cylinder resource blocker span on {:?}/{:?} is not contiguous",
                span.owner_face, span.owner_loop
            )));
        }
        let Some((&first, &last)) = span.half_edges.first().zip(span.half_edges.last()) else {
            continue;
        };
        let after_last = brep
            .half_edges
            .get(last)
            .ok_or_else(|| {
                CoreError::InvalidHandle(format!(
                    "cylinder resource blocker span has stale half-edge {last:?}"
                ))
            })?
            .next;
        brep.retire_face_boundary_fragment(span.owner_face, first, after_last)
            .map_err(|err| {
                CoreError::PreconditionFailed(format!(
                    "cylinder resource blocker exposure rejected for {:?}/{:?} {:?}->{:?}: {err}",
                    span.owner_face, span.owner_loop, first, after_last
                ))
            })?;
        exposed_spans += 1;
        exposed_half_edges += span.half_edges.len();
    }
    Ok(CylinderResourceBlockerExposureResult {
        spans: spans.len(),
        half_edges: spans.iter().map(|span| span.half_edges.len()).sum(),
        exposed_spans,
        exposed_half_edges,
    })
}
