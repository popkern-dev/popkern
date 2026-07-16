//! Owner-side span derivation for cylinder-mediated support continuation.

use crate::error::CoreError;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

use super::{
    CylinderMediatedOwnerRepairSpan, CylinderMediatedSupportContinuationSegmentAuthoringPlan,
    apply_exposure_candidates::{
        ExposureConstructionCandidate, ExposureConstructionFamily, ExposureConstructionFeasibility,
        ExposureRollbackBoundary, ExposureSemanticOwnership,
    },
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct OwnerRepairExposureResult {
    pub(super) exposed_spans: usize,
    exposed_half_edges: usize,
}

impl OwnerRepairExposureResult {
    pub(super) fn candidate_count(self) -> usize {
        1
    }

    pub(super) fn candidate_gate_complete(self) -> bool {
        self.candidate().mutation_gate_complete()
    }

    pub(super) fn diagnostic_label(self) -> String {
        let candidate = self.candidate();
        format!(
            "owner_repair_candidate=exposed_spans={} exposed_half_edges={} gate_complete:{} {}",
            self.exposed_spans,
            self.exposed_half_edges,
            self.candidate_gate_complete(),
            candidate.diagnostic_label()
        )
    }

    fn candidate(self) -> ExposureConstructionCandidate {
        ExposureConstructionCandidate {
            family: ExposureConstructionFamily::OwnerBoundaryFragmentRetire,
            feasibility: if self.exposed_spans == 0 {
                ExposureConstructionFeasibility::NotApplicable
            } else {
                ExposureConstructionFeasibility::ExecutableByExistingOperator
            },
            reason: if self.exposed_spans == 0 {
                "owner-repair-exposure-not-required"
            } else {
                "owner-repair-boundary-fragment-retire-ready"
            },
            spans: self.exposed_spans,
            half_edges: self.exposed_half_edges,
            expected_edges: self.exposed_spans as i64,
            expected_faces: 0,
            expected_loops: 0,
            semantic_ownership: if self.exposed_spans == 0 {
                ExposureSemanticOwnership::None
            } else {
                ExposureSemanticOwnership::OutsideBoundarySideExposed
            },
            rollback_boundary: if self.exposed_spans == 0 {
                ExposureRollbackBoundary::NoMutation
            } else {
                ExposureRollbackBoundary::CandidateBrepBeforeCommit
            },
        }
    }
}

pub(super) fn owner_repair_spans(
    brep: &BRep,
    segment: &CylinderMediatedSupportContinuationSegmentAuthoringPlan,
) -> Result<Vec<CylinderMediatedOwnerRepairSpan>, CoreError> {
    let mut spans = Vec::new();
    let mut current: Option<OwnerRepairSpanBuilder> = None;
    for &outside_half_edge in &segment.half_edges {
        if !segment_requires_owner_repair(segment, outside_half_edge) {
            push_owner_repair_span(brep, &mut spans, current.take())?;
            continue;
        }
        let twin = segment_twin_half_edge(brep, outside_half_edge)?;
        let owner_loop = segment_owner_loop(brep, twin)?;
        let owner_face = brep.loops[owner_loop].face;
        match current.as_mut() {
            Some(builder) if builder.owner_loop == owner_loop => {
                builder.outside_half_edges.push(outside_half_edge);
                builder.owner_half_edges.push(twin);
            }
            Some(_) => {
                push_owner_repair_span(brep, &mut spans, current.take())?;
                current = Some(OwnerRepairSpanBuilder::new(
                    owner_face,
                    owner_loop,
                    outside_half_edge,
                    twin,
                ));
            }
            None => {
                current = Some(OwnerRepairSpanBuilder::new(
                    owner_face,
                    owner_loop,
                    outside_half_edge,
                    twin,
                ));
            }
        }
    }
    push_owner_repair_span(brep, &mut spans, current)?;
    Ok(spans)
}

pub(super) fn expose_owner_repair_spans(
    brep: &mut BRep,
    authoring: &[CylinderMediatedSupportContinuationSegmentAuthoringPlan],
) -> Result<OwnerRepairExposureResult, CoreError> {
    let mut exposed_spans = 0usize;
    let mut exposed_half_edges = 0usize;
    if std::env::var_os("PK_CMS_TRACE").is_some() {
        for segment in authoring {
            eprintln!(
                "[cms] segment loop={:?} idx={} hes={:?} owned(cyl/out/sup/ctx)={:?}/{:?}/{:?}/{:?} boundary={:?} req={:?}",
                segment.outside_loop,
                segment.segment_index,
                segment.half_edges,
                segment.cylinder_owned_half_edges,
                segment.outside_owned_half_edges,
                segment.support_owned_half_edges,
                segment.context_owned_half_edges,
                segment.boundary_bridge_half_edges,
                segment.requirement,
            );
            for span in &segment.owner_repair_spans {
                eprintln!(
                    "[cms]   owner-span face={:?} loop={:?} outside={:?} owner={:?}",
                    span.owner_face,
                    span.owner_loop,
                    span.outside_half_edges,
                    span.owner_half_edges
                );
            }
        }
    }
    for span in authoring
        .iter()
        .flat_map(|segment| segment.owner_repair_spans.iter())
    {
        if span.owner_half_edges.is_empty() {
            continue;
        }
        // A slit edge interior to the outside loop (both half-edges on
        // the same loop) produces MIRRORED self-owner spans — (a; b)
        // and (b; a). Retiring the first already exposes both sides;
        // the mirror must be recognized as done, not re-retired into
        // the fresh boundary ("interval already contains a boundary
        // edge", measured on the PK-092 paired cross-owner batch whose
        // cap loop carries an apex slit).
        let already_exposed = span.outside_half_edges.iter().all(|&outside_half_edge| {
            match brep.half_edges.get(outside_half_edge) {
                None => true,
                Some(data) => {
                    brep.boundary_edges.contains(&data.edge)
                        && !brep.half_edges.contains_key(data.twin)
                }
            }
        });
        if already_exposed {
            exposed_spans += 1;
            exposed_half_edges += span.outside_half_edges.len();
            continue;
        }
        // Owner-side traversal runs OPPOSITE to the outside span
        // (twin half-edges reverse direction), so a multi-edge owned
        // run's fragment interval starts at the twin of the LAST
        // outside half-edge. Single-half-edge spans masked this until
        // the first multi-edge run (PK-092 paired cross-owner spans:
        // the wrong-way walk swept the rest of the owner loop and hit
        // "interval already contains a boundary edge"). Orientation
        // is read from the loop itself rather than assumed.
        let (first, last) = if span.owner_half_edges.len() >= 2 {
            let head = span.owner_half_edges[0];
            let second = span.owner_half_edges[1];
            let head_next = brep.half_edges.get(head).map(|data| data.next);
            if head_next == Some(second) {
                (head, span.owner_half_edges[span.owner_half_edges.len() - 1])
            } else {
                (span.owner_half_edges[span.owner_half_edges.len() - 1], head)
            }
        } else {
            (span.owner_half_edges[0], span.owner_half_edges[0])
        };
        let after_last = brep
            .half_edges
            .get(last)
            .ok_or_else(|| {
                CoreError::InvalidHandle(format!("owner repair span has stale half-edge {last:?}"))
            })?
            .next;
        brep.retire_face_boundary_fragment(span.owner_face, first, after_last)
            .map_err(|err| {
                CoreError::PreconditionFailed(format!(
                    "owner repair span exposure rejected for {:?}/{:?} {:?}->{:?}: {err}",
                    span.owner_face, span.owner_loop, first, after_last
                ))
            })?;
        for &outside_half_edge in &span.outside_half_edges {
            let data = brep.half_edges.get(outside_half_edge).ok_or_else(|| {
                CoreError::PreconditionFailed(format!(
                    "owner repair exposure removed outside half-edge {outside_half_edge:?}"
                ))
            })?;
            if !brep.boundary_edges.contains(&data.edge) || brep.half_edges.contains_key(data.twin)
            {
                return Err(CoreError::PreconditionFailed(format!(
                    "owner repair exposure did not expose outside half-edge {outside_half_edge:?}"
                )));
            }
        }
        exposed_spans += 1;
        exposed_half_edges += span.outside_half_edges.len();
    }
    Ok(OwnerRepairExposureResult {
        exposed_spans,
        exposed_half_edges,
    })
}

fn segment_requires_owner_repair(
    segment: &CylinderMediatedSupportContinuationSegmentAuthoringPlan,
    half_edge: HalfEdgeHandle,
) -> bool {
    segment.cylinder_owned_half_edges.contains(&half_edge)
        || segment.outside_owned_half_edges.contains(&half_edge)
        || segment.support_owned_half_edges.contains(&half_edge)
        || segment.context_owned_half_edges.contains(&half_edge)
}

fn segment_twin_half_edge(
    brep: &BRep,
    outside_half_edge: HalfEdgeHandle,
) -> Result<HalfEdgeHandle, CoreError> {
    let twin = brep
        .half_edges
        .get(outside_half_edge)
        .ok_or_else(|| {
            CoreError::InvalidHandle(format!(
                "owner repair span references stale outside half-edge {outside_half_edge:?}"
            ))
        })?
        .twin;
    if !brep.half_edges.contains_key(twin) {
        return Err(CoreError::PreconditionFailed(format!(
            "owner repair span outside half-edge {outside_half_edge:?} has no live twin"
        )));
    }
    Ok(twin)
}

fn segment_owner_loop(brep: &BRep, twin: HalfEdgeHandle) -> Result<LoopHandle, CoreError> {
    let owner_loop = brep.half_edges[twin].loop_;
    if !brep.loops.contains_key(owner_loop) {
        return Err(CoreError::PreconditionFailed(format!(
            "owner repair span twin half-edge {twin:?} has no live owner loop"
        )));
    }
    Ok(owner_loop)
}

fn push_owner_repair_span(
    brep: &BRep,
    spans: &mut Vec<CylinderMediatedOwnerRepairSpan>,
    builder: Option<OwnerRepairSpanBuilder>,
) -> Result<(), CoreError> {
    let Some(mut builder) = builder else {
        return Ok(());
    };
    builder.owner_half_edges.reverse();
    for pair in builder.owner_half_edges.windows(2) {
        if brep.half_edges[pair[0]].next != pair[1] {
            return Err(CoreError::PreconditionFailed(format!(
                "owner repair span on {:?}/{:?} is not contiguous at {:?}->{:?}",
                builder.owner_face, builder.owner_loop, pair[0], pair[1]
            )));
        }
    }
    spans.push(CylinderMediatedOwnerRepairSpan {
        owner_face: builder.owner_face,
        owner_loop: builder.owner_loop,
        outside_half_edges: builder.outside_half_edges,
        owner_half_edges: builder.owner_half_edges,
    });
    Ok(())
}

struct OwnerRepairSpanBuilder {
    owner_face: FaceHandle,
    owner_loop: LoopHandle,
    outside_half_edges: Vec<HalfEdgeHandle>,
    owner_half_edges: Vec<HalfEdgeHandle>,
}

impl OwnerRepairSpanBuilder {
    fn new(
        owner_face: FaceHandle,
        owner_loop: LoopHandle,
        outside_half_edge: HalfEdgeHandle,
        owner_half_edge: HalfEdgeHandle,
    ) -> Self {
        Self {
            owner_face,
            owner_loop,
            outside_half_edges: vec![outside_half_edge],
            owner_half_edges: vec![owner_half_edge],
        }
    }
}
