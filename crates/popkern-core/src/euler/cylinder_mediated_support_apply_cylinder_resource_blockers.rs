//! Exposure of twin-side blockers for cylinder bridge resources.

use crate::error::CoreError;
use crate::topology::BRep;

use super::apply_cylinder_resource_remap::CylinderBridgeResourceRemap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct CylinderResourceBlockerExposureResult {
    spans: usize,
    half_edges: usize,
    exposed_spans: usize,
    exposed_half_edges: usize,
}

impl CylinderResourceBlockerExposureResult {
    pub(super) fn diagnostic_label(self) -> String {
        format!(
            "cylinder_resource_blocker_exposure=spans={} half_edges={} exposed_spans={} exposed_half_edges={}",
            self.spans, self.half_edges, self.exposed_spans, self.exposed_half_edges
        )
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
