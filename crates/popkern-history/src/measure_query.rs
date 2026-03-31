//! Measure query types for the expression evaluator.
//!
//! [`MeasureQuery`] variants encode entity type and property. They live in
//! popkern-history because they reference [`EntityName`] — the compute layer
//! (popkern-measure) uses typed handles instead.

use std::collections::HashMap;

use popkern_core::naming::EntityName;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use serde::{Deserialize, Serialize};

use crate::node::BodyId;
use crate::unit::DimensionKind;

/// Component accessor for vector-valued measurements.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Component {
    /// X component.
    X,
    /// Y component.
    Y,
    /// Z component.
    Z,
}

/// Bounding box extent accessor.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BoundExtent {
    /// Minimum X.
    MinX,
    /// Minimum Y.
    MinY,
    /// Minimum Z.
    MinZ,
    /// Maximum X.
    MaxX,
    /// Maximum Y.
    MaxY,
    /// Maximum Z.
    MaxZ,
    /// Size along X.
    SizeX,
    /// Size along Y.
    SizeY,
    /// Size along Z.
    SizeZ,
}

/// A geometric query that produces a scalar f64 result.
///
/// Each variant encodes the entity type and property. Invalid combinations
/// (e.g., `FaceRadius` on a vertex) are unrepresentable.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MeasureQuery {
    // -- Edge ---------------------------------------------------
    /// Arc length of an edge.
    EdgeLength(EntityName),
    /// Component of the tangent at the edge midpoint.
    EdgeTangent(EntityName, Component),

    // -- Face ---------------------------------------------------
    /// Surface area of a face.
    FaceArea(EntityName),
    /// Perimeter (outer loop length) of a face.
    FacePerimeter(EntityName),
    /// Component of the outward normal at the face centroid.
    FaceNormal(EntityName, Component),
    /// Component of the face centroid.
    FaceCentroid(EntityName, Component),
    /// Radius of a cylindrical, spherical, or toroidal face.
    FaceRadius(EntityName),
    /// Minimum principal curvature at the face centroid.
    FaceMinCurvature(EntityName),
    /// Maximum principal curvature at the face centroid.
    FaceMaxCurvature(EntityName),

    // -- Vertex -------------------------------------------------
    /// Component of a vertex position.
    VertexPosition(EntityName, Component),

    // -- Between entities ---------------------------------------
    /// Minimum distance between two entities.
    MinDistance(EntityName, EntityName),
    /// Dihedral angle between two planar faces.
    DihedralAngle(EntityName, EntityName),
    /// Angle between two edges (tangent directions).
    EdgeAngle(EntityName, EntityName),

    // -- Body (whole solid) -------------------------------------
    /// Enclosed volume.
    BodyVolume,
    /// Total surface area.
    BodySurfaceArea,
    /// Component of body centroid.
    BodyCentroid(Component),
    /// Bounding box extent.
    BodyBound(BoundExtent),
    /// Moment of inertia about a principal axis.
    BodyMomentOfInertia(popkern_measure::Axis),
}

impl MeasureQuery {
    /// Return the [`DimensionKind`] of this query's scalar result.
    pub fn dimension_kind(&self) -> DimensionKind {
        match self {
            Self::EdgeLength(_)
            | Self::FacePerimeter(_)
            | Self::FaceRadius(_)
            | Self::MinDistance(_, _)
            | Self::VertexPosition(_, _)
            | Self::FaceCentroid(_, _)
            | Self::FaceNormal(_, _)
            | Self::EdgeTangent(_, _)
            | Self::BodyCentroid(_)
            | Self::BodyBound(_) => DimensionKind::Length,

            Self::FaceArea(_) | Self::BodySurfaceArea => DimensionKind::Area,

            Self::BodyVolume => DimensionKind::Volume,

            Self::DihedralAngle(_, _) | Self::EdgeAngle(_, _) => DimensionKind::Angle,

            Self::FaceMinCurvature(_) | Self::FaceMaxCurvature(_) => DimensionKind::InverseLength,

            Self::BodyMomentOfInertia(_) => DimensionKind::MomentOfArea,
        }
    }

    /// Extract all [`EntityName`] references from this query.
    pub fn entity_names(&self) -> Vec<&EntityName> {
        match self {
            Self::EdgeLength(n)
            | Self::EdgeTangent(n, _)
            | Self::FaceArea(n)
            | Self::FacePerimeter(n)
            | Self::FaceNormal(n, _)
            | Self::FaceCentroid(n, _)
            | Self::FaceRadius(n)
            | Self::FaceMinCurvature(n)
            | Self::FaceMaxCurvature(n)
            | Self::VertexPosition(n, _) => vec![n],

            Self::MinDistance(a, b)
            | Self::DihedralAngle(a, b)
            | Self::EdgeAngle(a, b) => vec![a, b],

            Self::BodyVolume
            | Self::BodySurfaceArea
            | Self::BodyCentroid(_)
            | Self::BodyBound(_)
            | Self::BodyMomentOfInertia(_) => vec![],
        }
    }
}

/// Typed context for expression evaluation with BRep access.
///
/// Passed to `evaluate()` when MEASURE expressions may be present.
pub struct MeasureContext<'a> {
    /// All body BReps in the model.
    pub bodies: &'a HashMap<BodyId, BRep>,
    /// The currently active body.
    pub active_body: BodyId,
    /// Geometric tolerance.
    pub tolerance: &'a Tolerance,
}

impl MeasureContext<'_> {
    /// Create an empty context for non-model evaluation.
    ///
    /// Any MEASURE expression evaluated against this context will produce
    /// a descriptive error.
    pub fn empty() -> MeasureContext<'static> {
        static EMPTY_BODIES: std::sync::LazyLock<HashMap<BodyId, BRep>> =
            std::sync::LazyLock::new(HashMap::new);
        static DEFAULT_TOL: std::sync::LazyLock<Tolerance> =
            std::sync::LazyLock::new(|| Tolerance::default());
        MeasureContext {
            bodies: &EMPTY_BODIES,
            active_body: BodyId(0),
            tolerance: &DEFAULT_TOL,
        }
    }
}
