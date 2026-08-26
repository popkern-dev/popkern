//! The [`Feature`] enum and its supporting parameter types.
//!
//! Every variant stores the operation's parameters in a serialisable form.
//! Entity references use [`EntityName`] (stable names) rather than raw handles
//! so that references survive rollback and tree reordering.

use crate::node::BodyId;
use crate::parameter::ParameterId;
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
        /// Direction of extrusion.  `None` means use the face normal.
        #[serde(default)]
        direction: Option<Vec3Desc>,
        /// Where and how far to extrude.
        end_condition: ExtrudeEnd,
        /// Whether to create a new body, boss, or cut.
        #[serde(default)]
        mode: OperationMode,
        /// Optional draft angle applied along the primary direction (radians,
        /// resolved via parameter table).
        #[serde(default)]
        draft_angle: Option<ParameterId>,
        /// Optional second direction for a symmetric or two-sided extrude.
        #[serde(default)]
        second_direction: Option<ExtrudeSecondDir>,
    },
    /// Round one or more edges with a constant radius.
    Fillet {
        /// The edges to fillet, identified by stable name.
        edges: Vec<EntityName>,
        /// Fillet radius (resolved via parameter table).
        radius: ParameterId,
    },
    /// Chamfer one or more edges with a constant distance.
    Chamfer {
        /// The edges to chamfer, identified by stable name.
        edges: Vec<EntityName>,
        /// Chamfer distance (resolved via parameter table).
        distance: ParameterId,
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
        /// Wall thickness (resolved via parameter table).
        thickness: ParameterId,
        /// Direction to offset the shell.
        direction: ShellDirectionKind,
        /// How gap edges are closed (PK-328; defaults keep the
        /// Minkowski-faithful rolling-ball behavior for trees saved
        /// before the knob existed).
        #[serde(default)]
        join_style: JoinStyleKind,
        /// How globally self-intersecting offset regions are handled
        /// (PK-328).
        #[serde(default)]
        trim_policy: TrimPolicyKind,
    },
    /// Thicken the body symmetrically by ±thickness/2 (PK-328).
    Thicken {
        /// Total wall thickness (resolved via parameter table).
        thickness: ParameterId,
        /// How gap edges are closed.
        #[serde(default)]
        join_style: JoinStyleKind,
        /// How globally self-intersecting offset regions are handled.
        #[serde(default)]
        trim_policy: TrimPolicyKind,
    },
    /// Offset every face of the solid: positive grows, negative
    /// shrinks (PK-328).
    OffsetSolid {
        /// Signed offset distance (resolved via parameter table).
        distance: ParameterId,
        /// How gap edges are closed.
        #[serde(default)]
        join_style: JoinStyleKind,
        /// How globally self-intersecting offset regions are handled.
        #[serde(default)]
        trim_policy: TrimPolicyKind,
    },
    /// Apply a draft (taper) angle to faces relative to a neutral plane.
    Draft {
        /// Faces to draft, identified by stable name.
        faces: Vec<EntityName>,
        /// Draft angle in radians (resolved via parameter table).
        angle: ParameterId,
        /// Neutral plane face, identified by stable name.
        neutral_plane: EntityName,
    },
    /// Sweep a face along a path curve.
    Sweep {
        /// The face to sweep, identified by stable name.
        face: EntityName,
        /// The edge whose curve defines the sweep path.
        path_edge: EntityName,
        /// Twist angle applied along the sweep (resolved via parameter table).
        twist: ParameterId,
        /// Whether to create a new body, boss, or cut.
        #[serde(default)]
        mode: OperationMode,
    },
    /// Loft through multiple cross-section faces.
    Loft {
        /// Ordered cross-section faces to loft between.
        sections: Vec<EntityName>,
        /// Whether the loft wraps around to connect last section to first.
        closed: bool,
        /// Whether to create a new body, boss, or cut.
        #[serde(default)]
        mode: OperationMode,
    },
    /// Create a pipe (tube) along a path with constant radius.
    Pipe {
        /// The edge whose curve defines the pipe path.
        path_edge: EntityName,
        /// Pipe radius (resolved via parameter table).
        radius: ParameterId,
        /// Whether to create a new body, boss, or cut.
        #[serde(default)]
        mode: OperationMode,
    },
    /// Revolve a face around an axis.
    Revolve {
        /// The face to revolve.
        face: EntityName,
        /// The axis of revolution.
        axis: AxisDesc,
        /// How far (and where) the revolution terminates.
        end_condition: RevolveEnd,
        /// Whether to create a new body, boss, or cut.
        #[serde(default)]
        mode: OperationMode,
    },
    /// A named measurement reference. Evaluates a geometric query and
    /// exposes the result as read-only parameter(s).
    Reference {
        /// User-visible name for this measurement.
        name: String,
        /// The geometric query to evaluate.
        query: crate::measure_query::MeasureQuery,
    },
    /// Replicate one or more seed features in a rectangular (linear) pattern.
    LinearPattern {
        /// Feature IDs of the operations to replicate.
        seed_features: Vec<crate::node::FeatureId>,
        /// Primary direction axis (controls count and spacing).
        primary: crate::pattern::PatternDirection,
        /// Optional secondary direction axis for a 2-D grid pattern.
        secondary: Option<crate::pattern::PatternDirection>,
        /// When `true`, replicate geometry directly instead of replaying
        /// parametric history.
        #[serde(default)]
        geometry_pattern: bool,
        /// Per-instance overrides (suppression and parameter values).
        #[serde(default)]
        instance_overrides: std::collections::HashMap<usize, crate::pattern::InstanceOverride>,
    },
    /// Replicate one or more seed features around a rotation axis.
    CircularPattern {
        /// Feature IDs of the operations to replicate.
        seed_features: Vec<crate::node::FeatureId>,
        /// The axis of rotation.
        axis: AxisDesc,
        /// Number of instances (including the seed).
        count: ParameterId,
        /// Total angular span of the pattern in radians.
        total_angle: ParameterId,
        /// When `true`, replicate geometry directly instead of replaying
        /// parametric history.
        #[serde(default)]
        geometry_pattern: bool,
        /// Per-instance overrides (suppression and parameter values).
        #[serde(default)]
        instance_overrides: std::collections::HashMap<usize, crate::pattern::InstanceOverride>,
    },
    /// Replicate one or more seed features along a curve.
    CurvePattern {
        /// Feature IDs of the operations to replicate.
        seed_features: Vec<crate::node::FeatureId>,
        /// The edge whose underlying curve defines the path.
        curve: EntityName,
        /// Number of instances (including the seed).
        count: ParameterId,
        /// When `true`, instances are spaced equally along the curve arc
        /// length.
        #[serde(default = "default_true")]
        equal_spacing: bool,
        /// When `true`, each instance is oriented tangent to the curve.
        #[serde(default)]
        align_to_curve: bool,
        /// When `true`, replicate geometry directly instead of replaying
        /// parametric history.
        #[serde(default)]
        geometry_pattern: bool,
        /// Per-instance overrides (suppression and parameter values).
        #[serde(default)]
        instance_overrides: std::collections::HashMap<usize, crate::pattern::InstanceOverride>,
    },
    /// Mirror one or more seed features about a plane.
    Mirror {
        /// Feature IDs of the operations to mirror.
        seed_features: Vec<crate::node::FeatureId>,
        /// The mirror plane.
        plane: crate::pattern::MirrorPlaneDesc,
        /// When `true`, the mirrored geometry is merged with the host body.
        #[serde(default = "default_true")]
        merge_result: bool,
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

/// How gap edges of an offset skin are closed (PK-328).
///
/// Mirrors [`popkern_shell::types::JoinStyle`] with serde; converted
/// at the replay boundary.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub enum JoinStyleKind {
    /// Rolling-ball fills (Minkowski-faithful; the default).
    #[default]
    Arc,
    /// Extend and intersect the adjacent offsets; sharp result.
    Intersection,
}

/// How globally self-intersecting offset regions are handled (PK-328).
///
/// Mirrors [`popkern_shell::types::TrimPolicy`] with serde; converted
/// at the replay boundary.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub enum TrimPolicyKind {
    /// Trim the colliding regions and re-close (the default).
    #[default]
    Trim,
    /// Decline with diagnostics instead of trimming.
    Decline,
}

/// Controls whether an operation creates new geometry or merges with the
/// active body.
///
/// `NewBody` is the default and creates an independent solid.  `Boss` and
/// `Cut` automatically union or subtract from the nearest host body,
/// mirroring the behaviour of direct-modelling "Boss" and "Cut" operations.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub enum OperationMode {
    /// Create a new, independent body (default).
    #[default]
    NewBody,
    /// Unite the result with the host body (additive).
    Boss,
    /// Subtract the result from the host body (subtractive).
    Cut,
}

/// End condition for an extrude operation.
///
/// Determines where (and how far) the extrusion terminates.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ExtrudeEnd {
    /// Extrude by an explicit distance.
    Blind {
        /// Extrusion distance (resolved via parameter table).
        distance: ParameterId,
    },
    /// Extrude until the profile exits all intersecting geometry.
    ThroughAll,
    /// Extrude until the profile touches the given face.
    UpToSurface {
        /// Target face, identified by stable name.
        face: EntityName,
    },
    /// Extrude until the profile touches the given face, then step back by an
    /// offset distance.
    OffsetFromSurface {
        /// Target face, identified by stable name.
        face: EntityName,
        /// Signed offset from the surface (resolved via parameter table).
        offset: ParameterId,
    },
    /// Extrude until the profile reaches the given vertex.
    UpToVertex {
        /// Target vertex, identified by stable name.
        vertex: EntityName,
    },
    /// Extrude symmetrically about the sketch plane by ±`half_distance`.
    MidPlane {
        /// Half the total extrusion distance (resolved via parameter table).
        half_distance: ParameterId,
    },
}

/// Additional direction for a two-directional extrude.
///
/// When present, the extrusion grows in both the primary direction (stored on
/// the parent feature) and this secondary direction simultaneously.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExtrudeSecondDir {
    /// The secondary extrusion direction vector.
    pub direction: Vec3Desc,
    /// End condition for the secondary direction.
    pub end_condition: ExtrudeEnd,
    /// Optional draft angle applied along the secondary direction (in radians,
    /// resolved via parameter table).
    pub draft_angle: Option<ParameterId>,
}

/// Identifies a rotation axis for revolve and circular-pattern operations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum AxisDesc {
    /// An existing model edge whose underlying curve defines the axis.
    Edge(EntityName),
    /// An explicitly specified origin point and direction vector.
    Literal {
        /// Axis origin in world space.
        origin: [f64; 3],
        /// Axis unit direction in world space.
        direction: [f64; 3],
    },
}

/// End condition for a revolve operation.
///
/// Determines the angular extent of the revolution.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RevolveEnd {
    /// Revolve by an explicit angle.
    Blind {
        /// Revolution angle in radians (resolved via parameter table).
        angle: ParameterId,
    },
    /// Revolve a full 360°.
    FullCircle,
    /// Revolve symmetrically about the sketch plane by ±`half_angle`.
    MidPlane {
        /// Half the total revolution angle in radians (resolved via parameter
        /// table).
        half_angle: ParameterId,
    },
    /// Revolve until the swept surface reaches the given face.
    UpToSurface {
        /// Target face, identified by stable name.
        face: EntityName,
    },
}

/// Returns `true`. Used as the `#[serde(default = "default_true")]` target for
/// boolean fields that should default to `true` rather than `false`.
fn default_true() -> bool {
    true
}
