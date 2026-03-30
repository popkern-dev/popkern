//! The [`Feature`] enum and its supporting parameter types.
//!
//! Every variant stores the operation's parameters in a serialisable form.
//! Entity references use [`EntityName`] (stable names) rather than raw handles
//! so that references survive rollback and tree reordering.

use crate::node::BodyId;
use popkern_core::naming::EntityName;
use serde::{Deserialize, Serialize};

/// A single parametric modeling operation in the feature tree.
///
/// Each variant corresponds to one type of BRep operation.  The parameters
/// are stored in a serialisable form so the whole tree can be saved and
/// replayed from scratch.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Feature {
    /// A 2D sketch placed on a plane (or face), described by a JSON payload
    /// produced by `popkern-sketch`.
    Sketch {
        /// Description of the plane on which the sketch lives.
        plane: SketchPlaneDesc,
        /// JSON-encoded sketch data (geometry + constraints).
        sketch_data: String,
    },
    /// Extrude a face by sweeping it along a direction.
    Extrude {
        /// The face to extrude, identified by stable name.
        face: EntityName,
        /// Direction of extrusion.
        direction: Vec3Desc,
        /// Distance to extrude.
        distance: f64,
    },
    /// Round one or more edges with a constant radius.
    Fillet {
        /// The edges to fillet, identified by stable name.
        edges: Vec<EntityName>,
        /// Fillet radius.
        radius: f64,
    },
    /// Chamfer one or more edges with a constant distance.
    Chamfer {
        /// The edges to chamfer, identified by stable name.
        edges: Vec<EntityName>,
        /// Chamfer distance.
        distance: f64,
    },
    /// Boolean set operation combining the current body with a tool body.
    Boolean {
        /// Which boolean operation to perform.
        operation: BoolOpKind,
        /// The tool body to operate with.
        tool_body: BodyId,
    },
    /// Hollow out a solid by removing faces and offsetting the remaining surfaces.
    Shell {
        /// Faces to remove (the "opening"), identified by stable name.
        faces_to_remove: Vec<EntityName>,
        /// Wall thickness.
        thickness: f64,
        /// Direction to offset the shell.
        direction: ShellDirectionKind,
    },
    /// Apply a draft (taper) angle to faces relative to a neutral plane.
    Draft {
        /// Faces to draft, identified by stable name.
        faces: Vec<EntityName>,
        /// Draft angle in radians.
        angle: f64,
        /// Neutral plane face, identified by stable name.
        neutral_plane: EntityName,
    },
}

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
    /// A plane offset from a model face by the given distance.
    OnFaceOffset(EntityName, f64),
    /// An arbitrary plane defined by an origin point and a unit normal.
    Custom {
        /// Plane origin in world space.
        origin: [f64; 3],
        /// Plane unit normal in world space.
        normal: [f64; 3],
    },
}

/// A serialisable description of a 3-D direction vector.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Vec3Desc {
    /// A literal XYZ direction vector.
    Literal([f64; 3]),
    /// The outward normal of a named face.
    FaceNormal(EntityName),
    /// The negation of another direction.
    Negated(Box<Vec3Desc>),
}

/// Boolean set-operation kind.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum BoolOpKind {
    /// Unite (add) two bodies.
    Union,
    /// Retain only the intersection of two bodies.
    Intersection,
    /// Subtract the tool body from the target.
    Subtraction,
}

/// Shell offset direction.
///
/// This mirrors [`popkern_shell::types::ShellDirection`] but adds
/// `Serialize`/`Deserialize` so the feature tree can be persisted.
/// Convert to/from `ShellDirection` at the replay boundary.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum ShellDirectionKind {
    /// Offset inward (remove material from the inside).
    Inward,
    /// Offset outward (add material to the outside).
    Outward,
    /// Split the thickness equally: half inward, half outward.
    Midplane,
}
