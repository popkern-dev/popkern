//! Shared-trim ownership planning for KP20 support-cell boundaries.

use super::FullBoundaryTrim;
use popkern_core::topology::{EdgeHandle, FaceHandle, HalfEdgeHandle};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SupportCellTrimOwnership {
    FullBoundaryRow,
    BoundarySegment,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct OwnedSupportCellTrim {
    pub(super) trim: FullBoundaryTrim,
    pub(super) ownership: SupportCellTrimOwnership,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct SegmentOwnedBoundaryTrim {
    pub(super) trim: FullBoundaryTrim,
    pub(super) edge: EdgeHandle,
    pub(super) half_edge_a: HalfEdgeHandle,
    pub(super) face_a: FaceHandle,
    pub(super) interval_a: BoundaryInterval,
    pub(super) half_edge_b: HalfEdgeHandle,
    pub(super) face_b: FaceHandle,
    pub(super) interval_b: BoundaryInterval,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct BoundaryInterval {
    pub(super) start: f64,
    pub(super) end: f64,
}

impl SegmentOwnedBoundaryTrim {
    pub(super) fn has_topology_address(self) -> bool {
        self.edge != EdgeHandle::default()
            && self.half_edge_a != HalfEdgeHandle::default()
            && self.half_edge_b != HalfEdgeHandle::default()
            && self.face_a != FaceHandle::default()
            && self.face_b != FaceHandle::default()
            && self.interval_a.is_valid()
            && self.interval_b.is_valid()
    }

    pub(super) fn has_paired_half_edges(self) -> bool {
        self.half_edge_a != self.half_edge_b
    }

    pub(super) fn requires_edge_split(self) -> bool {
        false
    }
}

impl BoundaryInterval {
    fn is_valid(self) -> bool {
        self.start.is_finite()
            && self.end.is_finite()
            && (self.start - self.end).abs() > 1.0e-9
            && (0.0..=1.0).contains(&self.start)
            && (0.0..=1.0).contains(&self.end)
    }
}

pub(super) fn support_cell_trim_ownership_plan(
    trims: &[FullBoundaryTrim],
) -> Vec<OwnedSupportCellTrim> {
    trims
        .iter()
        .copied()
        .map(|trim| OwnedSupportCellTrim {
            trim,
            ownership: trim_ownership(trim),
        })
        .collect()
}

pub(super) fn full_row_owned_trims(trims: &[FullBoundaryTrim]) -> Vec<FullBoundaryTrim> {
    owned_trims(trims, SupportCellTrimOwnership::FullBoundaryRow)
}

pub(super) fn segment_owned_boundary_trims(
    trims: &[FullBoundaryTrim],
) -> Vec<SegmentOwnedBoundaryTrim> {
    support_cell_trim_ownership_plan(trims)
        .into_iter()
        .filter(|owned| owned.ownership == SupportCellTrimOwnership::BoundarySegment)
        .map(|owned| segment_owned_boundary_trim(owned.trim))
        .collect()
}

fn owned_trims(
    trims: &[FullBoundaryTrim],
    ownership: SupportCellTrimOwnership,
) -> Vec<FullBoundaryTrim> {
    support_cell_trim_ownership_plan(trims)
        .into_iter()
        .filter(|owned| owned.ownership == ownership)
        .map(|owned| owned.trim)
        .collect()
}

fn segment_owned_boundary_trim(trim: FullBoundaryTrim) -> SegmentOwnedBoundaryTrim {
    SegmentOwnedBoundaryTrim {
        trim,
        edge: trim.edge,
        half_edge_a: trim.half_edge_a,
        face_a: trim.face_a,
        interval_a: BoundaryInterval {
            start: trim.start_a,
            end: trim.end_a,
        },
        half_edge_b: trim.half_edge_b,
        face_b: trim.face_b,
        interval_b: BoundaryInterval {
            start: trim.start_b,
            end: trim.end_b,
        },
    }
}

fn trim_ownership(trim: FullBoundaryTrim) -> SupportCellTrimOwnership {
    if full_boundary_span(trim.start_a, trim.end_a) && full_boundary_span(trim.start_b, trim.end_b)
    {
        SupportCellTrimOwnership::FullBoundaryRow
    } else {
        SupportCellTrimOwnership::BoundarySegment
    }
}

fn full_boundary_span(start: f64, end: f64) -> bool {
    (close(start, 0.0) && close(end, 1.0)) || (close(start, 1.0) && close(end, 0.0))
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1.0e-9
}
