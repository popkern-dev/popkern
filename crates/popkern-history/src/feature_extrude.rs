//! Extrude and revolve termination types: end conditions, second
//! directions, and rotation axes ([`AxisDesc`] also serves transforms,
//! patterns, and helical sweeps).

use serde::{Deserialize, Serialize};

use crate::feature::Vec3Desc;
use crate::parameter::ParameterId;
use popkern_core::naming::EntityName;

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
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
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
