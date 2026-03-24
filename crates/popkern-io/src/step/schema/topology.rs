//! Typed topology entity definitions from the STEP schema (AP203/AP214).

use super::{MappingError, as_bool, as_entity_ref, as_ref_list, as_string, param_at};
use crate::step::parser::types::{EntityId, Parameter};

/// VERTEX_POINT — a topological vertex associated with a geometric point.
#[derive(Debug, Clone)]
pub struct VertexPoint {
    /// Label name.
    pub name: String,
    /// Reference to the `CARTESIAN_POINT`.
    pub point: EntityId,
}

impl VertexPoint {
    const TYPE: &str = "VERTEX_POINT";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            point: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// EDGE_CURVE — a topological edge associated with a geometric curve.
#[derive(Debug, Clone)]
pub struct EdgeCurve {
    /// Label name.
    pub name: String,
    /// Reference to the start `VERTEX_POINT`.
    pub start: EntityId,
    /// Reference to the end `VERTEX_POINT`.
    pub end: EntityId,
    /// Reference to the underlying curve geometry.
    pub edge_geometry: EntityId,
    /// Whether the edge direction agrees with the curve direction.
    pub same_sense: bool,
}

impl EdgeCurve {
    const TYPE: &str = "EDGE_CURVE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            start: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            end: as_entity_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            edge_geometry: as_entity_ref(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
            same_sense: as_bool(param_at(params, 4, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// ORIENTED_EDGE — an edge with explicit orientation.
#[derive(Debug, Clone)]
pub struct OrientedEdge {
    /// Label name.
    pub name: String,
    /// Reference to the `EDGE_CURVE`.
    pub edge_element: EntityId,
    /// Orientation with respect to the underlying edge.
    pub orientation: bool,
}

impl OrientedEdge {
    const TYPE: &str = "ORIENTED_EDGE";

    /// Parse from raw STEP parameters.
    ///
    /// ORIENTED_EDGE has params: name, *, *, #edge_element, .orientation.
    /// The two derived (`*`) params at indices 1 and 2 are skipped.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            edge_element: as_entity_ref(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
            orientation: as_bool(param_at(params, 4, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// EDGE_LOOP — a closed loop of oriented edges.
#[derive(Debug, Clone)]
pub struct EdgeLoop {
    /// Label name.
    pub name: String,
    /// References to the `ORIENTED_EDGE` entities in loop order.
    pub edge_list: Vec<EntityId>,
}

impl EdgeLoop {
    const TYPE: &str = "EDGE_LOOP";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            edge_list: as_ref_list(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// FACE_BOUND — a bound (inner or outer) of a face.
#[derive(Debug, Clone)]
pub struct FaceBound {
    /// Label name.
    pub name: String,
    /// Reference to the `EDGE_LOOP`.
    pub bound: EntityId,
    /// Orientation of the bound with respect to the face.
    pub orientation: bool,
}

impl FaceBound {
    const TYPE: &str = "FACE_BOUND";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            bound: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            orientation: as_bool(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// FACE_OUTER_BOUND — the outer bound of a face (same structure as `FACE_BOUND`).
#[derive(Debug, Clone)]
pub struct FaceOuterBound {
    /// Label name.
    pub name: String,
    /// Reference to the `EDGE_LOOP`.
    pub bound: EntityId,
    /// Orientation of the bound with respect to the face.
    pub orientation: bool,
}

impl FaceOuterBound {
    const TYPE: &str = "FACE_OUTER_BOUND";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            bound: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            orientation: as_bool(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// ADVANCED_FACE — a face with an underlying surface and explicit bounds.
#[derive(Debug, Clone)]
pub struct AdvancedFace {
    /// Label name.
    pub name: String,
    /// References to `FACE_BOUND` / `FACE_OUTER_BOUND` entities.
    pub bounds: Vec<EntityId>,
    /// Reference to the surface geometry.
    pub face_geometry: EntityId,
    /// Whether the face normal agrees with the surface normal.
    pub same_sense: bool,
}

impl AdvancedFace {
    const TYPE: &str = "ADVANCED_FACE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            bounds: as_ref_list(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            face_geometry: as_entity_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            same_sense: as_bool(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// CLOSED_SHELL — a closed set of faces forming a watertight boundary.
#[derive(Debug, Clone)]
pub struct ClosedShell {
    /// Label name.
    pub name: String,
    /// References to the `ADVANCED_FACE` entities.
    pub faces: Vec<EntityId>,
}

impl ClosedShell {
    const TYPE: &str = "CLOSED_SHELL";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            faces: as_ref_list(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// OPEN_SHELL — an open set of faces (not necessarily watertight).
#[derive(Debug, Clone)]
pub struct OpenShell {
    /// Label name.
    pub name: String,
    /// References to the `ADVANCED_FACE` entities.
    pub faces: Vec<EntityId>,
}

impl OpenShell {
    const TYPE: &str = "OPEN_SHELL";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            faces: as_ref_list(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// MANIFOLD_SOLID_BREP — a solid body defined by a single closed shell.
#[derive(Debug, Clone)]
pub struct ManifoldSolidBrep {
    /// Label name.
    pub name: String,
    /// Reference to the outer `CLOSED_SHELL`.
    pub outer: EntityId,
}

impl ManifoldSolidBrep {
    const TYPE: &str = "MANIFOLD_SOLID_BREP";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            outer: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// BREP_WITH_VOIDS — a solid body with internal voids.
#[derive(Debug, Clone)]
pub struct BrepWithVoids {
    /// Label name.
    pub name: String,
    /// Reference to the outer `CLOSED_SHELL`.
    pub outer: EntityId,
    /// References to void `CLOSED_SHELL` entities.
    pub voids: Vec<EntityId>,
}

impl BrepWithVoids {
    const TYPE: &str = "BREP_WITH_VOIDS";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            outer: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            voids: as_ref_list(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p_str(s: &str) -> Parameter {
        Parameter::String(s.into())
    }
    fn p_ref(id: u64) -> Parameter {
        Parameter::EntityRef(id)
    }
    fn p_enum(s: &str) -> Parameter {
        Parameter::Enum(s.into())
    }

    #[test]
    fn vertex_point_parse() {
        let params = vec![p_str(""), p_ref(100)];
        let vp = VertexPoint::from_params(&params).unwrap();
        assert_eq!(vp.point, 100);
    }

    #[test]
    fn edge_curve_parse() {
        let params = vec![p_str("e1"), p_ref(1), p_ref(2), p_ref(3), p_enum("T")];
        let ec = EdgeCurve::from_params(&params).unwrap();
        assert_eq!(ec.start, 1);
        assert_eq!(ec.end, 2);
        assert_eq!(ec.edge_geometry, 3);
        assert!(ec.same_sense);
    }

    #[test]
    fn oriented_edge_parse() {
        let params = vec![
            p_str(""),
            Parameter::Derived,
            Parameter::Derived,
            p_ref(50),
            p_enum("F"),
        ];
        let oe = OrientedEdge::from_params(&params).unwrap();
        assert_eq!(oe.edge_element, 50);
        assert!(!oe.orientation);
    }

    #[test]
    fn advanced_face_parse() {
        let params = vec![
            p_str("face"),
            Parameter::List(vec![p_ref(10), p_ref(11)]),
            p_ref(20),
            p_enum("T"),
        ];
        let af = AdvancedFace::from_params(&params).unwrap();
        assert_eq!(af.bounds, vec![10, 11]);
        assert_eq!(af.face_geometry, 20);
        assert!(af.same_sense);
    }

    #[test]
    fn manifold_solid_brep_parse() {
        let params = vec![p_str("body"), p_ref(99)];
        let msb = ManifoldSolidBrep::from_params(&params).unwrap();
        assert_eq!(msb.name, "body");
        assert_eq!(msb.outer, 99);
    }
}
