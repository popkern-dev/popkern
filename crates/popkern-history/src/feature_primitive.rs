//! Parametric primitives (PK-1083): a box, cylinder, cone or frustum,
//! sphere or torus built from its dimensions and placed by a frame,
//! [`Feature::Primitive`](crate::feature::Feature::Primitive).
//!
//! Every dimension is a parameter (lengths in the model's unit, angles
//! in degrees), so a primitive is driven like any feature. Replay
//! builds it in its local frame (axis +z, reference +x, origin at the
//! base centre or, for a box, a corner) and places it rigidly; faces
//! carry role names under the feature's operation (`replay_primitive`).

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
    /// A cylinder on the axis, its base centred on the origin.
    Cylinder {
        /// Radius.
        radius: ParameterId,
        /// Height along the axis.
        height: ParameterId,
    },
    /// A cone or frustum on the axis: base radius `r1` at the origin,
    /// top radius `r2` at `height` (either may be 0, not both).
    Cone {
        /// Base radius.
        r1: ParameterId,
        /// Top radius.
        r2: ParameterId,
        /// Height along the axis.
        height: ParameterId,
    },
    /// A sphere centred on the origin, optionally only between two
    /// latitudes (degrees from the equator, −90 ≤ lo < hi ≤ 90).
    Sphere {
        /// Radius.
        radius: ParameterId,
        /// `[lo, hi]` latitude limits.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        latitudes: Option<[ParameterId; 2]>,
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
    /// `torus`).
    pub fn spelling(&self) -> &'static str {
        match self {
            PrimitiveShape::Box { .. } => "box",
            PrimitiveShape::Cylinder { .. } => "cylinder",
            PrimitiveShape::Cone { .. } => "cone",
            PrimitiveShape::Sphere { .. } => "sphere",
            PrimitiveShape::Torus { .. } => "torus",
        }
    }

    /// Every dimension with its role.
    pub fn param_ids(&self) -> Vec<(&'static str, ParameterId)> {
        match self {
            PrimitiveShape::Box { x, y, z, .. } => vec![("x", *x), ("y", *y), ("z", *z)],
            PrimitiveShape::Cylinder { radius, height } => {
                vec![("radius", *radius), ("height", *height)]
            }
            PrimitiveShape::Cone { r1, r2, height } => {
                vec![("r1", *r1), ("r2", *r2), ("height", *height)]
            }
            PrimitiveShape::Sphere { radius, latitudes } => {
                let mut v = vec![("radius", *radius)];
                if let Some([lo, hi]) = latitudes {
                    v.push(("lat_lo", *lo));
                    v.push(("lat_hi", *hi));
                }
                v
            }
            PrimitiveShape::Torus {
                major,
                minor,
                angle,
            } => {
                let mut v = vec![("major", *major), ("minor", *minor)];
                if let Some(a) = angle {
                    v.push(("angle", *a));
                }
                v
            }
        }
    }
}
