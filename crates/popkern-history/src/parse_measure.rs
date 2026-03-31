//! MEASURE expression parsing helpers.
//!
//! Contains the match tables that map dotted paths like `"body.volume"` or
//! `"face_0.area"` to [`MeasureQuery`] variants. Called from
//! [`Parser::parse_measure`](super::Parser::parse_measure).

use popkern_core::naming::EntityName;

use super::{EntityResolver, ParseError};
use crate::expression::Expression;
use crate::measure_query::{BoundExtent, Component, MeasureQuery};

// ---------------------------------------------------------------------------
// Body-level queries (no entity resolution needed)
// ---------------------------------------------------------------------------

/// Build a `MeasureQuery` from a single dotted path.
///
/// Body-level paths (e.g. `"body.volume"`) are handled directly. Entity-level
/// paths (e.g. `"face_0.area"`) require a resolver.
pub(super) fn build_single_measure(
    path: &str,
    resolver: Option<&dyn EntityResolver>,
    col: usize,
) -> Result<Expression, ParseError> {
    // Try body-level first (no resolver needed).
    if let Some(expr) = try_body_measure(path) {
        return Ok(expr);
    }

    // Entity-level: split into entity name + property suffix.
    let (entity_str, property) = split_entity_property(path)
        .ok_or_else(|| ParseError::new(col, format!("invalid MEASURE path: '{path}'")))?;

    let entity = resolve_entity_or_err(resolver, entity_str, col)?;
    build_entity_measure(entity, property, col)
}

/// Try to parse as a body-level MEASURE query.
fn try_body_measure(path: &str) -> Option<Expression> {
    let query = match path {
        "body.volume" => MeasureQuery::BodyVolume,
        "body.surface_area" => MeasureQuery::BodySurfaceArea,

        "body.centroid.x" => MeasureQuery::BodyCentroid(Component::X),
        "body.centroid.y" => MeasureQuery::BodyCentroid(Component::Y),
        "body.centroid.z" => MeasureQuery::BodyCentroid(Component::Z),

        "body.bound.min_x" => MeasureQuery::BodyBound(BoundExtent::MinX),
        "body.bound.min_y" => MeasureQuery::BodyBound(BoundExtent::MinY),
        "body.bound.min_z" => MeasureQuery::BodyBound(BoundExtent::MinZ),
        "body.bound.max_x" => MeasureQuery::BodyBound(BoundExtent::MaxX),
        "body.bound.max_y" => MeasureQuery::BodyBound(BoundExtent::MaxY),
        "body.bound.max_z" => MeasureQuery::BodyBound(BoundExtent::MaxZ),
        "body.bound.size_x" => MeasureQuery::BodyBound(BoundExtent::SizeX),
        "body.bound.size_y" => MeasureQuery::BodyBound(BoundExtent::SizeY),
        "body.bound.size_z" => MeasureQuery::BodyBound(BoundExtent::SizeZ),

        "body.moment_of_inertia.x" => MeasureQuery::BodyMomentOfInertia(popkern_measure::Axis::X),
        "body.moment_of_inertia.y" => MeasureQuery::BodyMomentOfInertia(popkern_measure::Axis::Y),
        "body.moment_of_inertia.z" => MeasureQuery::BodyMomentOfInertia(popkern_measure::Axis::Z),
        _ => return None,
    };
    Some(Expression::Measure(query))
}

// ---------------------------------------------------------------------------
// Entity-level queries
// ---------------------------------------------------------------------------

/// Known property suffixes for single-entity MEASURE queries, longest first.
const ENTITY_SUFFIXES: &[&str] = &[
    // vector-component (3-segment)
    ".normal.x",
    ".normal.y",
    ".normal.z",
    ".centroid.x",
    ".centroid.y",
    ".centroid.z",
    ".position.x",
    ".position.y",
    ".position.z",
    ".tangent.x",
    ".tangent.y",
    ".tangent.z",
    // scalar (2-segment)
    ".min_curvature",
    ".max_curvature",
    ".perimeter",
    ".length",
    ".area",
    ".radius",
];

/// Split a dotted path into `(entity_name, property_suffix)`.
///
/// Tries the longest known suffix first. Returns `None` if no known property
/// suffix matches.
fn split_entity_property(path: &str) -> Option<(&str, &str)> {
    for suffix in ENTITY_SUFFIXES {
        if let Some(entity) = path.strip_suffix(suffix)
            && !entity.is_empty()
        {
            // Return suffix without the leading '.'
            return Some((entity, &suffix[1..]));
        }
    }
    None
}

/// Build a `MeasureQuery` from a resolved entity and a property string.
fn build_entity_measure(
    entity: EntityName,
    property: &str,
    col: usize,
) -> Result<Expression, ParseError> {
    let query = match property {
        // Edge
        "length" => MeasureQuery::EdgeLength(entity),
        "tangent.x" => MeasureQuery::EdgeTangent(entity, Component::X),
        "tangent.y" => MeasureQuery::EdgeTangent(entity, Component::Y),
        "tangent.z" => MeasureQuery::EdgeTangent(entity, Component::Z),

        // Face scalar
        "area" => MeasureQuery::FaceArea(entity),
        "perimeter" => MeasureQuery::FacePerimeter(entity),
        "radius" => MeasureQuery::FaceRadius(entity),
        "min_curvature" => MeasureQuery::FaceMinCurvature(entity),
        "max_curvature" => MeasureQuery::FaceMaxCurvature(entity),

        // Face vector-component
        "normal.x" => MeasureQuery::FaceNormal(entity, Component::X),
        "normal.y" => MeasureQuery::FaceNormal(entity, Component::Y),
        "normal.z" => MeasureQuery::FaceNormal(entity, Component::Z),
        "centroid.x" => MeasureQuery::FaceCentroid(entity, Component::X),
        "centroid.y" => MeasureQuery::FaceCentroid(entity, Component::Y),
        "centroid.z" => MeasureQuery::FaceCentroid(entity, Component::Z),

        // Vertex
        "position.x" => MeasureQuery::VertexPosition(entity, Component::X),
        "position.y" => MeasureQuery::VertexPosition(entity, Component::Y),
        "position.z" => MeasureQuery::VertexPosition(entity, Component::Z),

        _ => {
            return Err(ParseError::new(
                col,
                format!("unknown entity property '{property}'"),
            ));
        }
    };
    Ok(Expression::Measure(query))
}

// ---------------------------------------------------------------------------
// Pair queries (two entities)
// ---------------------------------------------------------------------------

/// Known property suffixes for pair MEASURE queries, longest first.
const PAIR_SUFFIXES: &[&str] = &[".dihedral_angle", ".min_distance", ".edge_angle"];

/// Build a `MeasureQuery` from a pair of entities.
///
/// `second` is the second dotted identifier including the property suffix
/// (e.g. `"face_1.min_distance"`).
pub(super) fn build_pair_measure(
    first: &str,
    second: &str,
    resolver: Option<&dyn EntityResolver>,
    col: usize,
) -> Result<Expression, ParseError> {
    let (entity_b_str, property) = split_pair_property(second)
        .ok_or_else(|| ParseError::new(col, format!("invalid pair MEASURE syntax: '{second}'")))?;

    let entity_a = resolve_entity_or_err(resolver, first, col)?;
    let entity_b = resolve_entity_or_err(resolver, entity_b_str, col)?;

    let query = match property {
        "min_distance" => MeasureQuery::MinDistance(entity_a, entity_b),
        "dihedral_angle" => MeasureQuery::DihedralAngle(entity_a, entity_b),
        "edge_angle" => MeasureQuery::EdgeAngle(entity_a, entity_b),
        _ => {
            return Err(ParseError::new(
                col,
                format!("unknown pair property '{property}'"),
            ));
        }
    };
    Ok(Expression::Measure(query))
}

/// Split a pair second-argument path into `(entity_name, property)`.
fn split_pair_property(path: &str) -> Option<(&str, &str)> {
    for suffix in PAIR_SUFFIXES {
        if let Some(entity) = path.strip_suffix(suffix)
            && !entity.is_empty()
        {
            return Some((entity, &suffix[1..]));
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Entity resolution helper
// ---------------------------------------------------------------------------

/// Resolve an entity name string or produce a parse error.
fn resolve_entity_or_err(
    resolver: Option<&dyn EntityResolver>,
    name: &str,
    col: usize,
) -> Result<EntityName, ParseError> {
    match resolver {
        Some(r) => r
            .resolve_entity(name)
            .ok_or_else(|| ParseError::new(col, format!("unknown entity '{name}'"))),
        None => Err(ParseError::new(
            col,
            format!(
                "entity resolution not available for '{name}'; \
                 MEASURE with entity references requires a resolver"
            ),
        )),
    }
}
