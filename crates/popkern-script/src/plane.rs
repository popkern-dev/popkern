//! Sketch plane specifications (spec §1.4): the script-level form of
//! `SketchPlaneDesc`, with offsets and faces given as scalar bindings
//! and selections that the profile builder resolves at record time.

use popkern_history::feature::{BasePlane, SketchPlaneDesc};
use popkern_history::unit::Unit;

use crate::scalar::Scalar;
use crate::selection::FaceSelector;

/// Where a profile lives.
#[derive(Debug, Clone, PartialEq)]
pub enum Plane {
    /// The XY plane.
    XY,
    /// The XZ plane.
    XZ,
    /// The YZ plane.
    YZ,
    /// A base plane offset along its normal.
    Offset(BasePlane, Scalar),
    /// On a planar face of the active body.
    OnFace(FaceSelector),
    /// Offset from a planar face of the active body.
    OnFaceOffset(FaceSelector, Scalar),
    /// An arbitrary plane by origin and normal.
    Custom {
        /// A point on the plane.
        origin: [f64; 3],
        /// The plane normal.
        normal: [f64; 3],
    },
}

impl Plane {
    /// The XY plane at height `z`.
    pub fn xy_at(z: impl Into<Scalar>) -> Plane {
        Plane::Offset(BasePlane::XY, z.into())
    }
    /// The XZ plane at `y`.
    pub fn xz_at(y: impl Into<Scalar>) -> Plane {
        Plane::Offset(BasePlane::XZ, y.into())
    }
    /// The YZ plane at `x`.
    pub fn yz_at(x: impl Into<Scalar>) -> Plane {
        Plane::Offset(BasePlane::YZ, x.into())
    }
    /// On a face.
    pub fn on_face(faces: FaceSelector) -> Plane {
        Plane::OnFace(faces)
    }
    /// Offset from a face.
    pub fn on_face_offset(faces: FaceSelector, distance: impl Into<Scalar>) -> Plane {
        Plane::OnFaceOffset(faces, distance.into())
    }
    /// By origin and normal.
    pub fn custom(origin: [f64; 3], normal: [f64; 3]) -> Plane {
        Plane::Custom { origin, normal }
    }
    /// The descriptor for the base planes that need no resolution. A
    /// custom plane's origin is a position in the model's `unit`; the
    /// descriptor stores it in metres, as every bound length is.
    pub(crate) fn simple_desc(&self, unit: Unit) -> Option<SketchPlaneDesc> {
        match self {
            Plane::XY => Some(SketchPlaneDesc::XY),
            Plane::XZ => Some(SketchPlaneDesc::XZ),
            Plane::YZ => Some(SketchPlaneDesc::YZ),
            Plane::Custom { origin, normal } => Some(SketchPlaneDesc::Custom {
                origin: origin.map(|c| unit.to_si(c)),
                normal: *normal,
            }),
            _ => None,
        }
    }
}
