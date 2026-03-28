use serde::{Deserialize, Serialize};

slotmap::new_key_type! {
    /// Opaque handle to a [`Point`] inside a [`Sketch`](crate::Sketch).
    pub struct PointHandle;
}

slotmap::new_key_type! {
    /// Opaque handle to a [`Line`] inside a [`Sketch`](crate::Sketch).
    pub struct LineHandle;
}

slotmap::new_key_type! {
    /// Opaque handle to an [`Arc`] inside a [`Sketch`](crate::Sketch).
    pub struct ArcHandle;
}

/// A reference to any geometric element in the sketch.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ElementHandle {
    /// References a [`Point`].
    Point(PointHandle),
    /// References a [`Line`].
    Line(LineHandle),
    /// References an [`Arc`].
    Arc(ArcHandle),
}

impl From<PointHandle> for ElementHandle {
    fn from(h: PointHandle) -> Self {
        ElementHandle::Point(h)
    }
}

impl From<LineHandle> for ElementHandle {
    fn from(h: LineHandle) -> Self {
        ElementHandle::Line(h)
    }
}

impl From<ArcHandle> for ElementHandle {
    fn from(h: ArcHandle) -> Self {
        ElementHandle::Arc(h)
    }
}

/// A 2D point stored as two consecutive entries in the global parameter
/// vector: `params[param_start]` = x, `params[param_start + 1]` = y.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Point {
    /// Index of the x parameter in the sketch's parameter vector.
    pub param_start: usize,
    /// When `true` this element is reference geometry and is not exported
    /// as a model edge.
    pub construction: bool,
}

impl Point {
    /// Return the `[x_index, y_index]` pair into the global parameter vector.
    pub fn param_indices(&self) -> [usize; 2] {
        [self.param_start, self.param_start + 1]
    }
}

/// A straight line segment between two existing points.
///
/// Lines do not own any parameters themselves; all numeric values come from
/// the referenced [`Point`] elements.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Line {
    /// Handle to the start point.
    pub start: PointHandle,
    /// Handle to the end point.
    pub end: PointHandle,
    /// When `true` this element is reference geometry.
    pub construction: bool,
}

/// A circular arc stored as five consecutive entries in the global parameter
/// vector:
///
/// | offset | meaning      |
/// |--------|--------------|
/// | 0      | cx           |
/// | 1      | cy           |
/// | 2      | r            |
/// | 3      | start\_angle |
/// | 4      | end\_angle   |
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Arc {
    /// Index of the first (cx) parameter in the sketch's parameter vector.
    pub param_start: usize,
    /// When `true` this element is reference geometry.
    pub construction: bool,
}

impl Arc {
    /// Return the `[cx, cy, r, start_angle, end_angle]` indices into the
    /// global parameter vector.
    pub fn param_indices(&self) -> [usize; 5] {
        [
            self.param_start,
            self.param_start + 1,
            self.param_start + 2,
            self.param_start + 3,
            self.param_start + 4,
        ]
    }
}
