//! The [`Feature`] enum and its supporting parameter types.
//!
//! Every variant stores the operation's parameters in a serialisable form.
//! Entity references use [`EntityName`] (stable names) rather than raw handles
//! so that references survive rollback and tree reordering.

use crate::feature_blend::{ChamferKindDesc, FilletRecipeDesc};
use crate::feature_draft::NeutralPlaneDesc;
use crate::feature_sketch::SketchPayload;
use crate::feature_snapshot::SnapshotOrigin;
use crate::feature_transform::RigidMotionDesc;
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
    /// A 2D sketch placed on a plane (or face): a constraint sketch
    /// (JSON, re-solved on replay) or an analytic profile (compiled
    /// exactly on replay) — see [`SketchPayload`].
    Sketch {
        /// Description of the plane on which the sketch lives.
        plane: SketchPlaneDesc,
        /// The sketch itself, in either dialect.
        sketch_data: SketchPayload,
    },
    /// A body recorded as a geometry snapshot rather than a recipe
    /// (PK-367): an import, a polyhedron, a smooth-mesh solid, an
    /// engraved result. Replay materialises the snapshot as a new body
    /// and mints every entity a stable name under this feature's
    /// operation — see [`crate::feature_snapshot`].
    RawBody {
        /// The body itself.
        snapshot: Box<popkern_core::topology::BRep>,
        /// What produced it.
        origin: SnapshotOrigin,
    },
    /// Replace the body (the design space) with the stiffest design at
    /// a volume fraction under the recipe's supports and loads (PK-956):
    /// a planar-faced body, its kept faces modified from the design
    /// space's.
    Generative {
        /// The recipe.
        spec: Box<crate::feature_generative::GenerativeDesc>,
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
    /// Round one or more edges. `radius` alone is a constant-radius
    /// `G1` fillet; `recipe` carries the full blend recipe (radius
    /// law, continuity, caps, profile, interference policy).
    Fillet {
        /// The edges to fillet, identified by stable name.
        edges: Vec<EntityName>,
        /// Fillet radius, or the start radius of a variable law
        /// (resolved via parameter table).
        radius: ParameterId,
        /// The blend recipe (PK-361). Defaults reproduce the
        /// constant-radius `G1` fillet trees saved before the recipe
        /// existed.
        #[serde(default)]
        recipe: FilletRecipeDesc,
    },
    /// Chamfer one or more edges. `distance` alone is an equal-distance
    /// bevel; `kind` selects two-distance or distance-angle geometry.
    Chamfer {
        /// The edges to chamfer, identified by stable name.
        edges: Vec<EntityName>,
        /// Chamfer distance, or the first distance of a two-distance /
        /// distance-angle kind (resolved via parameter table).
        distance: ParameterId,
        /// The chamfer geometry kind (PK-361). Defaults to the equal
        /// bevel trees saved before the kind existed.
        #[serde(default)]
        kind: ChamferKindDesc,
    },
    /// Move the whole current body by a rigid (similarity) motion
    /// (PK-362). Topology and stable names are preserved; only
    /// geometry moves.
    Transform {
        /// The motion to apply.
        motion: RigidMotionDesc,
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
        /// The neutral plane: a planar model face by stable name, or a
        /// literal `{origin, normal}` plane (PK-363). Its normal is the
        /// pull direction.
        neutral_plane: NeutralPlaneDesc,
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
    /// Sweep a face along a helical path about an axis (springs, coils,
    /// auger flights, thread forms).
    HelicalSweep {
        /// The face to sweep, identified by stable name.
        face: EntityName,
        /// The helix axis.
        axis: AxisDesc,
        /// Helix radius at the start (resolved via parameter table).
        start_radius: ParameterId,
        /// Helix radius at the end — equal to `start_radius` for a
        /// cylindrical helix, different for a conical one.
        end_radius: ParameterId,
        /// Axial advance per turn (resolved via parameter table).
        pitch: ParameterId,
        /// Number of turns, fractional allowed (resolved via parameter
        /// table).
        turns: ParameterId,
        /// Phase: angle about the axis at which the helix starts
        /// (resolved via parameter table) — what phases multiple
        /// helical bodies apart (propeller blades, multi-start threads).
        start_angle: ParameterId,
        /// Winding sense. `false` = right-hand.
        #[serde(default)]
        left_handed: bool,
        /// Whether to create a new body, boss, or cut.
        #[serde(default)]
        mode: OperationMode,
    },
    /// A generated fastener body with a real modeled thread
    /// (PK-460): bolt/screw/rod/nut, built in the insertion frame
    /// (axis +Z, origin at the bearing plane). Always a new body.
    Fastener {
        /// Which fastener to build.
        kind: crate::feature_fastener::FastenerKind,
        /// The thread.
        thread: crate::feature_fastener::ThreadDesc,
        /// Under-head / overall / rod length; a nut's thickness (≤ 0 = nominal).
        length: ParameterId,
        /// Thread only the last `thread_length` of the shank
        /// (`None` = full thread).
        thread_length: Option<ParameterId>,
        /// Internal fit clearance (nut only; `None` = exact basic
        /// profile).
        clearance: Option<ParameterId>,
        /// 45° lead-in chamfer on the threaded tip (screws/rods).
        #[serde(default = "crate::feature_fastener::default_true")]
        tip_chamfer: bool,
        /// 30° hex corner chamfers on both nut faces (nuts only; `true` before PK-686).
        #[serde(default = "crate::feature_fastener::default_true")]
        corner_chamfer: bool,
        /// 45° thread-mouth lead-in chamfers on both nut faces (nuts only; ditto).
        #[serde(default = "crate::feature_fastener::default_true")]
        mouth_chamfer: bool,
        /// Cosmetic thread (screws/rods): a root-diameter cylinder, no helix (PK-847).
        #[serde(default)]
        cosmetic_thread: bool,
    },
    /// A tapped (threaded) hole cut into the body: bore at the
    /// internal minor diameter plus a real modeled internal thread
    /// (PK-460). The axis direction points INTO the material.
    TappedHole {
        /// The hole axis (origin on the entry surface, direction
        /// into the material).
        axis: AxisDesc,
        /// The thread.
        thread: crate::feature_fastener::ThreadDesc,
        /// Blind depth along the axis; `None` = through.
        depth: Option<ParameterId>,
        /// Internal fit clearance (0 = exact basic profile).
        clearance: ParameterId,
        /// 45° lead-in chamfer at the entry.
        #[serde(default)]
        entry_chamfer: bool,
    },
    /// A plain clearance hole for a screw, with optional countersink
    /// or counterbore (PK-460).
    ClearanceHole {
        /// The hole axis (origin on the entry surface, direction
        /// into the material).
        axis: AxisDesc,
        /// The screw thread the hole clears.
        thread: crate::feature_fastener::ThreadDesc,
        /// Hole diameter series or explicit diameter.
        fit: crate::feature_fastener::ClearanceFitDesc,
        /// Blind depth along the axis; `None` = through.
        depth: Option<ParameterId>,
        /// Sink a 90° countersink sized from the countersunk-head
        /// standard plus this extra radial room.
        countersink: Option<ParameterId>,
        /// Sink a counterbore sized from the socket-head standard
        /// plus this extra radial/axial room.
        counterbore: Option<ParameterId>,
        /// 45° entry chamfer (skipped under a countersink or
        /// counterbore).
        #[serde(default)]
        entry_chamfer: bool,
    },
    /// One body of a hinge (PK-693); components share `spec`.
    Hinge {
        /// Butt or print-in-place.
        kind: crate::feature_hinge::HingeKindDesc,
        /// Leaf A, leaf B or the pin.
        component: crate::feature_hinge::HingeComponentDesc,
        /// The shared dimensions.
        spec: crate::feature_hinge::HingeDesc,
        /// Whether this feature owns the shared parameters.
        owner: bool,
    },
    /// A hinge grown on two bodies along their seam (PK-702): `LeafA` /
    /// `LeafB` knuckle the receiving body, `Pin` is a new body.
    SeamHinge {
        /// Butt or print-in-place.
        kind: crate::feature_hinge::HingeKindDesc,
        /// Which body's knuckles, or the pin.
        component: crate::feature_hinge::HingeComponentDesc,
        /// The seam, resolved at record time.
        seam: crate::feature_hinge::SeamDesc,
        /// The shared dimensions.
        spec: crate::feature_hinge::SeamHingeDesc,
        /// Whether this feature owns the shared parameters.
        owner: bool,
    },
    /// A living hinge body (PK-698); always a new body.
    LivingHinge {
        /// Straight web or slotted serpentine.
        style: crate::feature_hinge::LivingHingeStyleDesc,
        /// The dimensions.
        spec: crate::feature_hinge::LivingHingeDesc,
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

pub use crate::feature_shell_kinds::{JoinStyleKind, ShellDirectionKind, TrimPolicyKind};

pub use crate::feature_kinds::{BoolOpKind, OperationMode};

pub use crate::feature_extrude::{AxisDesc, ExtrudeEnd, ExtrudeSecondDir, RevolveEnd};

/// Returns `true`. Used as the `#[serde(default = "default_true")]` target for
/// boolean fields that should default to `true` rather than `false`.
fn default_true() -> bool {
    true
}
