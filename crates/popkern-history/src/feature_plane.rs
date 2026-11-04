//! Sketch-plane and direction descriptors of the feature tree (split
//! from `feature.rs` under the 500-line file budget).

use crate::parameter::ParameterId;
use popkern_core::naming::EntityName;
use serde::{Deserialize, Serialize};

/// Describes the plane a sketch is placed on.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SketchPlaneDesc {
    /// The global XY plane (Z = 0, normal = +Z).
    XY,
    /// The global XZ plane (Y = 0, normal = +Y).
    XZ,
    /// The global YZ plane (X = 0, normal = +X).
    YZ,
    /// A planar face in the model, referenced by stable name.
    OnFace(EntityName),
    /// A plane offset from a model face by the given distance (resolved via
    /// parameter table).
    OnFaceOffset(EntityName, ParameterId),
    /// An arbitrary plane defined by an origin point and a unit normal.
    Custom {
        /// Plane origin in world space.
        origin: [f64; 3],
        /// Plane unit normal in world space.
        normal: [f64; 3],
    },
    /// A base plane translated along its normal by a parameter
    /// (`xy_at(z)` and friends): keeps the base plane's u/v frame, so a
    /// profile at height z compiles exactly as on the base plane.
    BaseOffset {
        /// Which base plane.
        base: BasePlane,
        /// Signed offset along the base normal.
        distance: ParameterId,
    },
}

/// One of the three principal planes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BasePlane {
    /// The XY plane (normal +Z).
    XY,
    /// The XZ plane.
    XZ,
    /// The YZ plane.
    YZ,
}

/// A serialisable description of a 3-D direction vector.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Vec3Desc {
    /// A literal XYZ direction vector.
    Literal([f64; 3]),
    /// The outward normal of a named face.
    FaceNormal(EntityName),
    /// The tangent direction of a named edge at its midpoint.
    EdgeDirection(EntityName),
    /// The negation of another direction.
    Negated(Box<Vec3Desc>),
}
