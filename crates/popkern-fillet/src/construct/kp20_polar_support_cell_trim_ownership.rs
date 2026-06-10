//! Shared-trim ownership planning for KP20 support-cell boundaries.

use super::FullBoundaryTrim;

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

pub(super) fn segment_owned_trims(trims: &[FullBoundaryTrim]) -> Vec<FullBoundaryTrim> {
    owned_trims(trims, SupportCellTrimOwnership::BoundarySegment)
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
