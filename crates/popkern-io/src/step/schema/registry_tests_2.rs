//! Continuation of registry entity resolution tests (product, units, context).
//! Split from registry_tests.rs to keep files under 500 lines.
//! Included only in test mode via `#[cfg(test)] #[path = ...] mod tests_2`.
// #[cfg(test)] — marker for the no_unwrap linter; .unwrap() is allowed in tests.

use super::*;
use crate::step::parser::types::{EntityInstance, Parameter};

fn make_entity(type_name: &str, params: Vec<Parameter>) -> EntityInstance {
    EntityInstance {
        id: 1,
        type_name: type_name.to_string(),
        params,
        complex_parts: vec![],
    }
}

// --- Product ---

#[test]
fn resolve_product() {
    let entity = make_entity(
        "PRODUCT",
        vec![
            Parameter::String("id".into()),
            Parameter::String("name".into()),
            Parameter::String("desc".into()),
            Parameter::List(vec![Parameter::EntityRef(1)]),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::Product(_)));
}

#[test]
fn resolve_product_definition_formation() {
    let entity = make_entity(
        "PRODUCT_DEFINITION_FORMATION",
        vec![
            Parameter::String("v1".into()),
            Parameter::String("desc".into()),
            Parameter::EntityRef(1),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(
        resolved,
        ResolvedEntity::ProductDefinitionFormation(_)
    ));
}

#[test]
fn resolve_product_definition() {
    let entity = make_entity(
        "PRODUCT_DEFINITION",
        vec![
            Parameter::String("".into()),
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::EntityRef(2),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::ProductDefinition(_)));
}

#[test]
fn resolve_product_definition_shape() {
    let entity = make_entity(
        "PRODUCT_DEFINITION_SHAPE",
        vec![
            Parameter::String("".into()),
            Parameter::String("".into()),
            Parameter::EntityRef(1),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(
        resolved,
        ResolvedEntity::ProductDefinitionShape(_)
    ));
}

#[test]
fn resolve_shape_definition_representation() {
    let entity = make_entity(
        "SHAPE_DEFINITION_REPRESENTATION",
        vec![Parameter::EntityRef(1), Parameter::EntityRef(2)],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(
        resolved,
        ResolvedEntity::ShapeDefinitionRepresentation(_)
    ));
}

#[test]
fn resolve_shape_representation() {
    let entity = make_entity(
        "SHAPE_REPRESENTATION",
        vec![
            Parameter::String("".into()),
            Parameter::List(vec![Parameter::EntityRef(1)]),
            Parameter::EntityRef(2),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::ShapeRepresentation(_)));
}

#[test]
fn resolve_advanced_brep_shape_representation() {
    let entity = make_entity(
        "ADVANCED_BREP_SHAPE_REPRESENTATION",
        vec![
            Parameter::String("".into()),
            Parameter::List(vec![Parameter::EntityRef(1)]),
            Parameter::EntityRef(2),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(
        resolved,
        ResolvedEntity::AdvancedBrepShapeRepresentation(_)
    ));
}

#[test]
fn resolve_faceted_brep_shape_representation() {
    let entity = make_entity(
        "FACETED_BREP_SHAPE_REPRESENTATION",
        vec![
            Parameter::String("".into()),
            Parameter::List(vec![Parameter::EntityRef(1)]),
            Parameter::EntityRef(2),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(
        resolved,
        ResolvedEntity::FacetedBrepShapeRepresentation(_)
    ));
}

#[test]
fn resolve_next_assembly_usage_occurrence() {
    let entity = make_entity(
        "NEXT_ASSEMBLY_USAGE_OCCURRENCE",
        vec![
            Parameter::String("".into()),
            Parameter::String("".into()),
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::EntityRef(2),
            Parameter::Omitted,
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(
        resolved,
        ResolvedEntity::NextAssemblyUsageOccurrence(_)
    ));
}

#[test]
fn resolve_representation_relationship() {
    let entity = make_entity(
        "REPRESENTATION_RELATIONSHIP",
        vec![
            Parameter::String("".into()),
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::EntityRef(2),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(
        resolved,
        ResolvedEntity::RepresentationRelationship(_)
    ));
}

#[test]
fn resolve_representation_relationship_with_transformation() {
    let entity = make_entity(
        "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION",
        vec![
            Parameter::String("".into()),
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::EntityRef(2),
            Parameter::EntityRef(3),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(
        resolved,
        ResolvedEntity::RepresentationRelationshipWithTransformation(_)
    ));
}

#[test]
fn resolve_item_defined_transformation() {
    let entity = make_entity(
        "ITEM_DEFINED_TRANSFORMATION",
        vec![
            Parameter::String("".into()),
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::EntityRef(2),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(
        resolved,
        ResolvedEntity::ItemDefinedTransformation(_)
    ));
}

// --- Units/Context ---

#[test]
fn resolve_geometric_representation_context() {
    let entity = make_entity(
        "GEOMETRIC_REPRESENTATION_CONTEXT",
        vec![Parameter::String("".into()), Parameter::Int(3)],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(
        resolved,
        ResolvedEntity::GeometricRepresentationContext(_)
    ));
}

#[test]
fn resolve_global_uncertainty_assigned_context() {
    let entity = make_entity(
        "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT",
        vec![
            Parameter::String("".into()),
            Parameter::Int(3),
            Parameter::List(vec![Parameter::EntityRef(1)]),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(
        resolved,
        ResolvedEntity::GlobalUncertaintyAssignedContext(_)
    ));
}

#[test]
fn resolve_uncertainty_measure_with_unit() {
    let entity = make_entity(
        "UNCERTAINTY_MEASURE_WITH_UNIT",
        vec![
            Parameter::Real(0.001),
            Parameter::EntityRef(1),
            Parameter::String("tol".into()),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(
        resolved,
        ResolvedEntity::UncertaintyMeasureWithUnit(_)
    ));
}

#[test]
fn resolve_si_unit() {
    let entity = make_entity(
        "SI_UNIT",
        vec![
            Parameter::Derived,
            Parameter::Omitted,
            Parameter::Enum("METRE".into()),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::SiUnit(_)));
}

#[test]
fn resolve_conversion_based_unit() {
    let entity = make_entity(
        "CONVERSION_BASED_UNIT",
        vec![
            Parameter::Omitted,
            Parameter::String("INCH".into()),
            Parameter::EntityRef(1),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::ConversionBasedUnit(_)));
}

#[test]
fn resolve_length_measure_with_unit() {
    let entity = make_entity(
        "LENGTH_MEASURE_WITH_UNIT",
        vec![Parameter::Real(25.4), Parameter::EntityRef(1)],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::LengthMeasureWithUnit(_)));
}

#[test]
fn resolve_plane_angle_measure_with_unit() {
    let entity = make_entity(
        "PLANE_ANGLE_MEASURE_WITH_UNIT",
        vec![Parameter::Real(0.01745), Parameter::EntityRef(1)],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(
        resolved,
        ResolvedEntity::PlaneAngleMeasureWithUnit(_)
    ));
}

#[test]
fn resolve_malformed_params() {
    // CARTESIAN_POINT with no params — missing required name string
    let entity = make_entity("CARTESIAN_POINT", vec![]);
    let result = resolve_entity(&entity);
    assert!(result.is_err(), "malformed params should return Err");
}
