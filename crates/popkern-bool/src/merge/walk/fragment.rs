//! A face fragment the walk produces (split from the walk's module
//! under the file budget).

use super::super::graph::FaceRef;
use super::FragmentChart;

/// A face fragment — loops forming a subset of an original face.
#[derive(Debug)]
pub(crate) struct Fragment {
    /// The face this fragment belongs to.
    pub face: FaceRef,
    /// Half-edge indices forming the outer loop (CCW).
    pub outer_loop: Vec<usize>,
    /// Half-edge indices forming inner loops (CW holes).
    pub inner_loops: Vec<Vec<usize>>,
    /// The chart walk's own view of this fragment, when the chart
    /// permutation walk produced it. `None` for 3D-Newell fallback
    /// walks, whose loops carry no chart handedness guarantee.
    pub chart: Option<FragmentChart>,
}

impl Fragment {
    /// A fragment the walk produced (no chart view yet, not identity).
    pub(super) fn walked(
        face: FaceRef,
        outer_loop: Vec<usize>,
        inner_loops: Vec<Vec<usize>>,
    ) -> Self {
        Self {
            face,
            outer_loop,
            inner_loops,
            chart: None,
        }
    }
}
