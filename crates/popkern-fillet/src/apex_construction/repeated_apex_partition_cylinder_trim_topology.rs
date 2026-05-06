//! Cylinder trim topology choices for repeated-apex partition rewrites.

use super::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum FarSideTrimStrategy {
    Forward,
    Reverse,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum CylinderLoopOrientationKind {
    CurrentSharedFilletEnd,
    ReauthoredSharedFilletEnd,
}

#[derive(Clone)]
pub(super) struct CylinderLoopOrientationCandidate {
    pub(super) kind: CylinderLoopOrientationKind,
    pub(super) far_side_strategy: FarSideTrimStrategy,
    pub(super) segments: Vec<OrientedCylinderTrimSegment>,
}

pub(super) struct CylinderTrimRewriteTopologyPlan {
    pub(super) fillet_index: usize,
    pub(super) cylinder_face: FaceHandle,
    pub(super) loop_candidate: CylinderLoopOrientationCandidate,
    pub(super) seam_rewrite: FilletEndSeamRewrite,
}

pub(super) struct PendingCylinderLoopReplacement {
    pub(super) fillet_index: usize,
    pub(super) cylinder_face: FaceHandle,
    pub(super) kind: CylinderLoopOrientationKind,
    pub(super) cycle: Vec<HalfEdgeHandle>,
}

pub(super) enum FilletEndSeamRewrite {
    ReuseCurrentTwinPair {
        edge: EdgeHandle,
        cap_half_edge: HalfEdgeHandle,
        current_cylinder_half_edge: HalfEdgeHandle,
    },
    ReauthorCapBoundaryAndCylinderSide {
        edge: EdgeHandle,
        cap_half_edge: HalfEdgeHandle,
        current_cylinder_half_edge: HalfEdgeHandle,
    },
}

#[derive(Clone)]
pub(super) struct OrientedCylinderTrimSegment {
    pub(super) label: &'static str,
    pub(super) start: VertexHandle,
    pub(super) end: VertexHandle,
}

impl CylinderTrimRewriteTopologyPlan {
    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "FilletEnd[{}] cylinder {:?}: {}, {}",
            self.fillet_index,
            self.cylinder_face,
            self.loop_candidate.diagnostic_label(),
            self.seam_rewrite.diagnostic_label()
        )
    }
}

impl FilletEndSeamRewrite {
    pub(super) fn diagnostic_label(&self) -> String {
        match self {
            FilletEndSeamRewrite::ReuseCurrentTwinPair {
                edge,
                cap_half_edge,
                current_cylinder_half_edge,
            } => format!(
                "reuse current FilletEnd cap/cylinder twin pair edge={edge:?} \
                 cap_he={cap_half_edge:?} current_cylinder_he={current_cylinder_half_edge:?}"
            ),
            FilletEndSeamRewrite::ReauthorCapBoundaryAndCylinderSide {
                edge,
                cap_half_edge,
                current_cylinder_half_edge,
            } => format!(
                "reauthor cap boundary and cylinder side for FilletEnd seam edge={edge:?} \
                 cap_he={cap_half_edge:?} current_cylinder_he={current_cylinder_half_edge:?}"
            ),
        }
    }
}

impl CylinderLoopOrientationCandidate {
    pub(super) fn requires_fillet_end_reorientation(&self) -> bool {
        self.kind == CylinderLoopOrientationKind::ReauthoredSharedFilletEnd
    }

    pub(super) fn validate_contiguous(&self, fillet_index: usize) -> Result<(), String> {
        for i in 0..self.segments.len() {
            let current = &self.segments[i];
            let next = &self.segments[(i + 1) % self.segments.len()];
            if current.end != next.start {
                return Err(format!(
                    "FilletEnd[{fillet_index}] segment {} ends at {:?}, but segment {} starts \
                     at {:?}",
                    current.label, current.end, next.label, next.start
                ));
            }
        }
        Ok(())
    }

    pub(super) fn diagnostic_label(&self) -> String {
        let segments = self
            .segments
            .iter()
            .map(OrientedCylinderTrimSegment::diagnostic_label)
            .collect::<Vec<_>>()
            .join(", ");
        format!(
            "{:?}/{:?} reorient_fillet_end={} segments=[{}]",
            self.kind,
            self.far_side_strategy,
            self.requires_fillet_end_reorientation(),
            segments
        )
    }
}

impl OrientedCylinderTrimSegment {
    pub(super) fn diagnostic_label(&self) -> String {
        format!("{} {:?}->{:?}", self.label, self.start, self.end)
    }
}
