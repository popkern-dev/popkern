//! Canonical curve source choices for support-cell chain authoring.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ChainCurveSourcePreference {
    LeastSquaresFit,
    PreferSpanningOwnerRow,
}

impl ChainCurveSourcePreference {
    pub(super) fn prefers_spanning_owner_row(self) -> bool {
        matches!(self, Self::PreferSpanningOwnerRow)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ChainCurveSource {
    LeastSquaresFit,
    SpanningOwnerRow,
}
