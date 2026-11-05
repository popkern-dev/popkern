//! Parametric primitives (PK-1083): a box, cylinder, cone or frustum,
//! sphere, regular prism, pyramid or its frustum, wedge or torus built
//! from its dimensions and placed by a frame,
//! [`Feature::Primitive`](crate::feature::Feature::Primitive).
//!
//! Every dimension is a parameter (lengths in the model's unit, angles
//! in degrees), so a primitive is driven like any feature. Replay
//! builds it in its local frame (axis +z, reference +x, origin at the
//! base centre or, for a box, a corner) and places it rigidly; faces
//! carry role names under the feature's operation (`replay_primitive`).

use popkern_core::naming::EntityName;
use serde::{Deserialize, Serialize};

use crate::parameter::ParameterId;

/// A primitive and where it stands.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PrimitiveDesc {
    /// The shape and its dimensions.
    pub shape: PrimitiveShape,
    /// The local frame in world coordinates.
    #[serde(default, skip_serializing_if = "PrimitiveFrame::is_world")]
    pub frame: PrimitiveFrame,
}

/// The primitive shapes.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum PrimitiveShape {
    /// A box `x` × `y` × `z`, a corner at the origin (or centred on it).
    Box {
        /// Extent along the local x.
        x: ParameterId,
        /// Extent along the local y.
        y: ParameterId,
        /// Extent along the local z (the axis).
        z: ParameterId,
        /// Centred on the origin instead of cornered there.
        #[serde(default, skip_serializing_if = "std::ops::Not::not")]
        centred: bool,
    },
    /// A cylinder on the axis, its base centred on the origin,
    /// optionally a sector of `angle` degrees from the reference.
    Cylinder {
        /// Radius.
        radius: ParameterId,
        /// Height along the axis.
        height: ParameterId,
        /// Sector angle in degrees from the reference direction
        /// (`None`: whole).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        angle: Option<ParameterId>,
    },
    /// A cone or frustum on the axis: base radius `r1` at the origin,
    /// top radius `r2` at `height` (either may be 0, not both),
    /// optionally a sector.
    Cone {
        /// Base radius.
        r1: ParameterId,
        /// Top radius.
        r2: ParameterId,
        /// Height along the axis.
        height: ParameterId,
        /// Sector angle in degrees from the reference direction
        /// (`None`: whole).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        angle: Option<ParameterId>,
    },
    /// A sphere centred on the origin, optionally only between two
    /// latitudes (degrees from the equator, −90 ≤ lo < hi ≤ 90) and
    /// only a sector.
    Sphere {
        /// Radius.
        radius: ParameterId,
        /// `[lo, hi]` latitude limits.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        latitudes: Option<[ParameterId; 2]>,
        /// Sector angle in degrees from the reference direction
        /// (`None`: whole).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        angle: Option<ParameterId>,
    },
    /// A regular prism: a `sides`-gon on the origin, corner 0 on the
    /// reference direction, standing `height` along the axis.
    Prism {
        /// Number of sides (at least 3).
        sides: u32,
        /// Centre to corner, or to the middle of a side when `flats`.
        radius: ParameterId,
        /// Height along the axis.
        height: ParameterId,
        /// `radius` is measured to the flats (the inscribed circle).
        #[serde(default, skip_serializing_if = "std::ops::Not::not")]
        flats: bool,
    },
    /// A regular pyramid or its frustum: base radius `r1` on the
    /// origin, top radius `r2` at `height` (either may be 0, an apex;
    /// not both), corner 0 of each on the reference direction.
    Pyramid {
        /// Number of sides (at least 3).
        sides: u32,
        /// Base radius (centre to corner, or to a flat when `flats`).
        r1: ParameterId,
        /// Top radius.
        r2: ParameterId,
        /// Height along the axis.
        height: ParameterId,
        /// The radii are measured to the flats.
        #[serde(default, skip_serializing_if = "std::ops::Not::not")]
        flats: bool,
    },
    /// A wedge: the box `x` × `y` × `z` on a corner at the origin
    /// whose top runs only `top_x` along x from the origin's side (0:
    /// a ramp, its section a right triangle).
    Wedge {
        /// Extent along the local x at the bottom.
        x: ParameterId,
        /// Extent along the local y.
        y: ParameterId,
        /// Height along the axis.
        z: ParameterId,
        /// Extent along x at the top.
        top_x: ParameterId,
    },
    /// A wire: straight edges through `points` (local frame, metres),
    /// the last joined back to the first when `closed`.
    WirePolyline {
        /// The corners, at least two.
        points: Vec<[f64; 3]>,
        /// Joined back to the first point.
        #[serde(default, skip_serializing_if = "std::ops::Not::not")]
        closed: bool,
    },
    /// A wire: the circular arc from `start` through `via` to `end`.
    WireArc {
        /// Start point.
        start: [f64; 3],
        /// A point on the arc between the ends.
        via: [f64; 3],
        /// End point.
        end: [f64; 3],
    },
    /// A wire: the circle of `radius` about the axis in the local xy
    /// plane, one closed edge starting on the reference direction.
    WireCircle {
        /// Radius.
        radius: ParameterId,
    },
    /// A wire: a helix of `radius` about the axis from the reference
    /// direction, rising `pitch` per turn for `turns`.
    WireHelix {
        /// Radius.
        radius: ParameterId,
        /// Rise per turn.
        pitch: ParameterId,
        /// Number of turns (fractional allowed).
        turns: ParameterId,
        /// Winding sense (`false`: right-hand).
        #[serde(default, skip_serializing_if = "std::ops::Not::not")]
        left_handed: bool,
    },
    /// A wire: the cubic (or `degree`) spline interpolating `points`
    /// (local frame, metres) at chord-length parameters, or, when
    /// `control`, the clamped uniform B-spline with `points` as its
    /// control polygon (four points of degree 3: a Bézier curve).
    WireSpline {
        /// The points it passes through (or its control points), at
        /// least `degree + 1`.
        points: Vec<[f64; 3]>,
        /// Polynomial degree.
        degree: u32,
        /// `points` are control points, not points on the curve.
        #[serde(default, skip_serializing_if = "std::ops::Not::not")]
        control: bool,
    },
    /// A sheet: the rectangle `x` × `y` in the local xy plane, a
    /// corner on the origin (or centred on it), facing the axis.
    SheetRect {
        /// Extent along the local x.
        x: ParameterId,
        /// Extent along the local y.
        y: ParameterId,
        /// Centred on the origin instead of cornered there.
        #[serde(default, skip_serializing_if = "std::ops::Not::not")]
        centred: bool,
    },
    /// A sheet: the disc of `radius` about the origin in the local xy
    /// plane, facing the axis.
    SheetDisc {
        /// Radius.
        radius: ParameterId,
    },
    /// A sheet: the polygon through `points` (local xy, metres),
    /// facing the axis.
    SheetPolygon {
        /// The corners, at least three.
        points: Vec<[f64; 2]>,
    },
    /// A sheet: the planar face a closed chain of wire edges bounds
    /// (in world coordinates; the frame is not used). The wire is
    /// referenced, not consumed.
    SheetFill {
        /// The bounding edges, in any order and direction.
        edges: Vec<EntityName>,
    },
    /// A torus about the axis, centred on the origin, optionally a
    /// sector of `angle` degrees from the reference direction.
    Torus {
        /// Distance from the axis to the tube's centre.
        major: ParameterId,
        /// Tube radius (less than `major`).
        minor: ParameterId,
        /// Sector angle in degrees (`None`: a whole ring).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        angle: Option<ParameterId>,
    },
}

/// A primitive's local frame: `origin`, its `axis` (local z) and a
/// `reference` direction (local x, made square to the axis).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PrimitiveFrame {
    /// The local origin.
    pub origin: [f64; 3],
    /// The local z.
    pub axis: [f64; 3],
    /// The local x (its part along `axis` is dropped).
    pub reference: [f64; 3],
}

impl Default for PrimitiveFrame {
    fn default() -> Self {
        PrimitiveFrame {
            origin: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
            reference: [1.0, 0.0, 0.0],
        }
    }
}

impl PrimitiveFrame {
    /// Whether this is the world frame (omitted on disk).
    pub fn is_world(&self) -> bool {
        *self == PrimitiveFrame::default()
    }
}

impl PrimitiveShape {
    /// The operation spelling (`box`, `cylinder`, `cone`, `sphere`,
    /// `prism`, `pyramid`, `wedge`, `torus`, and the wires `wire_line`,
    /// `wire_polyline`, `wire_arc`, `wire_circle`, `wire_helix`,
    /// `wire_spline`).
    pub fn spelling(&self) -> &'static str {
        match self {
            PrimitiveShape::Box { .. } => "box",
            PrimitiveShape::Cylinder { .. } => "cylinder",
            PrimitiveShape::Cone { .. } => "cone",
            PrimitiveShape::Sphere { .. } => "sphere",
            PrimitiveShape::Prism { .. } => "prism",
            PrimitiveShape::Pyramid { .. } => "pyramid",
            PrimitiveShape::Wedge { .. } => "wedge",
            PrimitiveShape::Torus { .. } => "torus",
            PrimitiveShape::WirePolyline { points, closed } => {
                if points.len() == 2 && !closed {
                    "wire_line"
                } else {
                    "wire_polyline"
                }
            }
            PrimitiveShape::WireArc { .. } => "wire_arc",
            PrimitiveShape::WireCircle { .. } => "wire_circle",
            PrimitiveShape::WireHelix { .. } => "wire_helix",
            PrimitiveShape::WireSpline { .. } => "wire_spline",
            PrimitiveShape::SheetRect { .. } => "sheet_rect",
            PrimitiveShape::SheetDisc { .. } => "sheet_disc",
            PrimitiveShape::SheetPolygon { .. } => "sheet_polygon",
            PrimitiveShape::SheetFill { .. } => "sheet_fill",
        }
    }

    /// Visit every entity name the shape references (a fill's edges).
    pub fn for_each_name_mut(&mut self, f: &mut dyn FnMut(&mut EntityName)) {
        if let PrimitiveShape::SheetFill { edges } = self {
            for name in edges.iter_mut() {
                f(name);
            }
        }
    }

    /// Whether the shape is a sheet (one face with free edges).
    pub fn is_sheet(&self) -> bool {
        matches!(
            self,
            PrimitiveShape::SheetRect { .. }
                | PrimitiveShape::SheetDisc { .. }
                | PrimitiveShape::SheetPolygon { .. }
                | PrimitiveShape::SheetFill { .. }
        )
    }

    /// Whether the shape is a wire (edges, no faces).
    pub fn is_wire(&self) -> bool {
        matches!(
            self,
            PrimitiveShape::WirePolyline { .. }
                | PrimitiveShape::WireArc { .. }
                | PrimitiveShape::WireCircle { .. }
                | PrimitiveShape::WireHelix { .. }
                | PrimitiveShape::WireSpline { .. }
        )
    }

    /// The sector angle, when the shape is only a sector of its
    /// revolution.
    pub fn angle(&self) -> Option<ParameterId> {
        match self {
            PrimitiveShape::Box { .. }
            | PrimitiveShape::Prism { .. }
            | PrimitiveShape::Pyramid { .. }
            | PrimitiveShape::Wedge { .. }
            | PrimitiveShape::WirePolyline { .. }
            | PrimitiveShape::WireArc { .. }
            | PrimitiveShape::WireCircle { .. }
            | PrimitiveShape::WireHelix { .. }
            | PrimitiveShape::WireSpline { .. }
            | PrimitiveShape::SheetRect { .. }
            | PrimitiveShape::SheetDisc { .. }
            | PrimitiveShape::SheetPolygon { .. }
            | PrimitiveShape::SheetFill { .. } => None,
            PrimitiveShape::Cylinder { angle, .. }
            | PrimitiveShape::Cone { angle, .. }
            | PrimitiveShape::Sphere { angle, .. }
            | PrimitiveShape::Torus { angle, .. } => *angle,
        }
    }

    /// Every dimension with its role.
    pub fn param_ids(&self) -> Vec<(&'static str, ParameterId)> {
        let mut v = self.dimensions();
        if let Some(a) = self.angle() {
            v.push(("angle", a));
        }
        v
    }

    fn dimensions(&self) -> Vec<(&'static str, ParameterId)> {
        match self {
            PrimitiveShape::Box { x, y, z, .. } => vec![("x", *x), ("y", *y), ("z", *z)],
            PrimitiveShape::Cylinder { radius, height, .. } => {
                vec![("radius", *radius), ("height", *height)]
            }
            PrimitiveShape::Cone { r1, r2, height, .. } => {
                vec![("r1", *r1), ("r2", *r2), ("height", *height)]
            }
            PrimitiveShape::Sphere {
                radius, latitudes, ..
            } => {
                let mut v = vec![("radius", *radius)];
                if let Some([lo, hi]) = latitudes {
                    v.push(("lat_lo", *lo));
                    v.push(("lat_hi", *hi));
                }
                v
            }
            PrimitiveShape::Prism { radius, height, .. } => {
                vec![("radius", *radius), ("height", *height)]
            }
            PrimitiveShape::Pyramid { r1, r2, height, .. } => {
                vec![("r1", *r1), ("r2", *r2), ("height", *height)]
            }
            PrimitiveShape::Wedge { x, y, z, top_x } => {
                vec![("x", *x), ("y", *y), ("z", *z), ("top_x", *top_x)]
            }
            PrimitiveShape::Torus { major, minor, .. } => {
                vec![("major", *major), ("minor", *minor)]
            }
            PrimitiveShape::WireCircle { radius } => vec![("radius", *radius)],
            PrimitiveShape::WireHelix {
                radius,
                pitch,
                turns,
                ..
            } => vec![("radius", *radius), ("pitch", *pitch), ("turns", *turns)],
            PrimitiveShape::SheetRect { x, y, .. } => vec![("x", *x), ("y", *y)],
            PrimitiveShape::SheetDisc { radius } => vec![("radius", *radius)],
            PrimitiveShape::WirePolyline { .. }
            | PrimitiveShape::WireArc { .. }
            | PrimitiveShape::WireSpline { .. }
            | PrimitiveShape::SheetPolygon { .. }
            | PrimitiveShape::SheetFill { .. } => Vec::new(),
        }
    }
}
