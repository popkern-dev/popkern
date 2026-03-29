use crate::rigid_set::RigidSetHandle;
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

slotmap::new_key_type! {
    /// Opaque handle to an [`Ellipse`] inside a [`Sketch`](crate::Sketch).
    pub struct EllipseHandle;
}

slotmap::new_key_type! {
    /// Opaque handle to a [`Parabola`] inside a [`Sketch`](crate::Sketch).
    pub struct ParabolaHandle;
}

slotmap::new_key_type! {
    /// Opaque handle to a [`Hyperbola`] inside a [`Sketch`](crate::Sketch).
    pub struct HyperbolaHandle;
}

slotmap::new_key_type! {
    /// Opaque handle to a [`Spline`] inside a [`Sketch`](crate::Sketch).
    pub struct SplineHandle;
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
    /// References an [`Ellipse`].
    Ellipse(EllipseHandle),
    /// References a [`Parabola`].
    Parabola(ParabolaHandle),
    /// References a [`Hyperbola`].
    Hyperbola(HyperbolaHandle),
    /// References a [`Spline`].
    Spline(SplineHandle),
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

impl From<EllipseHandle> for ElementHandle {
    fn from(h: EllipseHandle) -> Self {
        ElementHandle::Ellipse(h)
    }
}

impl From<ParabolaHandle> for ElementHandle {
    fn from(h: ParabolaHandle) -> Self {
        ElementHandle::Parabola(h)
    }
}

impl From<HyperbolaHandle> for ElementHandle {
    fn from(h: HyperbolaHandle) -> Self {
        ElementHandle::Hyperbola(h)
    }
}

impl From<SplineHandle> for ElementHandle {
    fn from(h: SplineHandle) -> Self {
        ElementHandle::Spline(h)
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
    /// Rigid set membership, if any.
    #[serde(default)]
    pub rigid_set: Option<RigidSetHandle>,
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
    /// Rigid set membership, if any.
    #[serde(default)]
    pub rigid_set: Option<RigidSetHandle>,
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

/// An ellipse stored as five consecutive entries in the global parameter
/// vector:
///
/// | offset | meaning    |
/// |--------|------------|
/// | 0      | cx         |
/// | 1      | cy         |
/// | 2      | semi\_a    |
/// | 3      | semi\_b    |
/// | 4      | rotation   |
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ellipse {
    /// Index of the first (cx) parameter in the sketch's parameter vector.
    pub param_start: usize,
    /// When `true` this element is reference geometry.
    pub construction: bool,
    /// Rigid set membership, if any.
    #[serde(default)]
    pub rigid_set: Option<RigidSetHandle>,
}

impl Ellipse {
    /// Return the `[cx, cy, semi_a, semi_b, rotation]` indices into the
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

/// A parabola stored as four consecutive entries in the global parameter
/// vector:
///
/// | offset | meaning      |
/// |--------|--------------|
/// | 0      | vx           |
/// | 1      | vy           |
/// | 2      | focal\_dist  |
/// | 3      | rotation     |
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Parabola {
    /// Index of the first (vx) parameter in the sketch's parameter vector.
    pub param_start: usize,
    /// When `true` this element is reference geometry.
    pub construction: bool,
    /// Rigid set membership, if any.
    #[serde(default)]
    pub rigid_set: Option<RigidSetHandle>,
}

impl Parabola {
    /// Return the `[vx, vy, focal_dist, rotation]` indices into the global
    /// parameter vector.
    pub fn param_indices(&self) -> [usize; 4] {
        [
            self.param_start,
            self.param_start + 1,
            self.param_start + 2,
            self.param_start + 3,
        ]
    }
}

/// A hyperbola stored as five consecutive entries in the global parameter
/// vector:
///
/// | offset | meaning    |
/// |--------|------------|
/// | 0      | cx         |
/// | 1      | cy         |
/// | 2      | semi\_a    |
/// | 3      | semi\_b    |
/// | 4      | rotation   |
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Hyperbola {
    /// Index of the first (cx) parameter in the sketch's parameter vector.
    pub param_start: usize,
    /// When `true` this element is reference geometry.
    pub construction: bool,
    /// Rigid set membership, if any.
    #[serde(default)]
    pub rigid_set: Option<RigidSetHandle>,
}

impl Hyperbola {
    /// Return the `[cx, cy, semi_a, semi_b, rotation]` indices into the
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

/// A NURBS spline whose geometry is fully determined by its control points.
///
/// The spline itself stores no parameters in the global parameter vector —
/// all numeric values come from the referenced [`Point`] elements in
/// `control_points`.  The `knots` and `weights` vectors are fixed structural
/// data that do not participate in the Newton–Raphson solver.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Spline {
    /// Ordered list of control point handles that define the spline shape.
    pub control_points: Vec<PointHandle>,
    /// Polynomial degree of each Bézier segment (≥ 1).
    pub degree: usize,
    /// Knot vector of length `control_points.len() + degree + 1`.
    pub knots: Vec<f64>,
    /// Per-control-point weights for rational (NURBS) evaluation.
    /// Must have the same length as `control_points`.
    pub weights: Vec<f64>,
    /// `true` if the spline is periodic (wraps around to form a closed loop).
    pub periodic: bool,
    /// When `true` this element is reference geometry and is not exported
    /// as a model edge.
    pub construction: bool,
}
