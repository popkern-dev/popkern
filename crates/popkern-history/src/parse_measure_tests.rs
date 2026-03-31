use crate::expression::Expression;
use crate::measure_query::MeasureQuery;
use crate::parameter::ParameterTable;

use super::parse_expression;

#[test]
fn test_measure_body_volume() {
    let t = ParameterTable::new();
    let expr = parse_expression("MEASURE(body.volume)", &t).unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::BodyVolume)
    ));
}

#[test]
fn test_measure_body_surface_area() {
    let t = ParameterTable::new();
    let expr = parse_expression("MEASURE(body.surface_area)", &t).unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::BodySurfaceArea)
    ));
}

#[test]
fn test_measure_body_centroid_x() {
    let t = ParameterTable::new();
    let expr = parse_expression("MEASURE(body.centroid.x)", &t).unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::BodyCentroid(
            crate::measure_query::Component::X
        ))
    ));
}

#[test]
fn test_measure_body_centroid_y() {
    let t = ParameterTable::new();
    let expr = parse_expression("MEASURE(body.centroid.y)", &t).unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::BodyCentroid(
            crate::measure_query::Component::Y
        ))
    ));
}

#[test]
fn test_measure_body_centroid_z() {
    let t = ParameterTable::new();
    let expr = parse_expression("MEASURE(body.centroid.z)", &t).unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::BodyCentroid(
            crate::measure_query::Component::Z
        ))
    ));
}

#[test]
fn test_measure_body_bound_min_x() {
    let t = ParameterTable::new();
    let expr = parse_expression("MEASURE(body.bound.min_x)", &t).unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::BodyBound(
            crate::measure_query::BoundExtent::MinX
        ))
    ));
}

#[test]
fn test_measure_body_bound_size_x() {
    let t = ParameterTable::new();
    let expr = parse_expression("MEASURE(body.bound.size_x)", &t).unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::BodyBound(
            crate::measure_query::BoundExtent::SizeX
        ))
    ));
}

#[test]
fn test_measure_body_bound_max_z() {
    let t = ParameterTable::new();
    let expr = parse_expression("MEASURE(body.bound.max_z)", &t).unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::BodyBound(
            crate::measure_query::BoundExtent::MaxZ
        ))
    ));
}

#[test]
fn test_measure_body_moment_of_inertia_z() {
    let t = ParameterTable::new();
    let expr = parse_expression("MEASURE(body.moment_of_inertia.z)", &t).unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::BodyMomentOfInertia(popkern_measure::Axis::Z))
    ));
}

#[test]
fn test_measure_body_moment_of_inertia_x() {
    let t = ParameterTable::new();
    let expr = parse_expression("MEASURE(body.moment_of_inertia.x)", &t).unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::BodyMomentOfInertia(popkern_measure::Axis::X))
    ));
}

#[test]
fn test_measure_in_expression() {
    let t = ParameterTable::new();
    let expr = parse_expression("MEASURE(body.volume) * 2", &t).unwrap();
    assert!(matches!(expr, Expression::BinOp(..)));
}

#[test]
fn test_measure_in_addition() {
    let t = ParameterTable::new();
    let expr = parse_expression(
        "MEASURE(body.bound.size_x) + MEASURE(body.bound.size_y)",
        &t,
    )
    .unwrap();
    assert!(matches!(expr, Expression::BinOp(..)));
}

#[test]
fn test_measure_entity_without_resolver() {
    let t = ParameterTable::new();
    let result = parse_expression("MEASURE(face_0.area)", &t);
    assert!(result.is_err());
    assert!(result.unwrap_err().message.contains("entity resolution"));
}

#[test]
fn test_measure_entity_with_resolver() {
    use super::parse_expression_with_resolver;
    use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};

    struct TestResolver;
    impl super::EntityResolver for TestResolver {
        fn resolve_entity(&self, name: &str) -> Option<EntityName> {
            if name == "face_0" {
                Some(EntityName {
                    operation: OperationId(1),
                    role: EntityRole {
                        category: RoleCategory::Generated,
                        index: 0,
                    },
                    ancestors: vec![],
                })
            } else {
                None
            }
        }
    }

    let t = ParameterTable::new();
    let resolver = TestResolver;
    let expr = parse_expression_with_resolver("MEASURE(face_0.area)", &t, Some(&resolver)).unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::FaceArea(_))
    ));
}

#[test]
fn test_measure_entity_normal_with_resolver() {
    use super::parse_expression_with_resolver;
    use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};

    struct TestResolver;
    impl super::EntityResolver for TestResolver {
        fn resolve_entity(&self, name: &str) -> Option<EntityName> {
            if name == "face_0" {
                Some(EntityName {
                    operation: OperationId(1),
                    role: EntityRole {
                        category: RoleCategory::Generated,
                        index: 0,
                    },
                    ancestors: vec![],
                })
            } else {
                None
            }
        }
    }

    let t = ParameterTable::new();
    let resolver = TestResolver;
    let expr =
        parse_expression_with_resolver("MEASURE(face_0.normal.x)", &t, Some(&resolver)).unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::FaceNormal(
            _,
            crate::measure_query::Component::X
        ))
    ));
}

#[test]
fn test_measure_entity_unknown() {
    use super::parse_expression_with_resolver;
    use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};

    struct TestResolver;
    impl super::EntityResolver for TestResolver {
        fn resolve_entity(&self, _name: &str) -> Option<EntityName> {
            None
        }
    }

    let t = ParameterTable::new();
    let resolver = TestResolver;
    let result = parse_expression_with_resolver("MEASURE(face_99.area)", &t, Some(&resolver));
    assert!(result.is_err());
    assert!(result.unwrap_err().message.contains("unknown entity"));
}

#[test]
fn test_measure_pair_min_distance() {
    use super::parse_expression_with_resolver;
    use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};

    struct TestResolver;
    impl super::EntityResolver for TestResolver {
        fn resolve_entity(&self, name: &str) -> Option<EntityName> {
            match name {
                "face_0" => Some(EntityName {
                    operation: OperationId(1),
                    role: EntityRole {
                        category: RoleCategory::Generated,
                        index: 0,
                    },
                    ancestors: vec![],
                }),
                "face_1" => Some(EntityName {
                    operation: OperationId(1),
                    role: EntityRole {
                        category: RoleCategory::Generated,
                        index: 1,
                    },
                    ancestors: vec![],
                }),
                _ => None,
            }
        }
    }

    let t = ParameterTable::new();
    let resolver = TestResolver;
    let expr =
        parse_expression_with_resolver("MEASURE(face_0, face_1.min_distance)", &t, Some(&resolver))
            .unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::MinDistance(_, _))
    ));
}

#[test]
fn test_measure_pair_dihedral_angle() {
    use super::parse_expression_with_resolver;
    use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};

    struct TestResolver;
    impl super::EntityResolver for TestResolver {
        fn resolve_entity(&self, name: &str) -> Option<EntityName> {
            match name {
                "face_0" => Some(EntityName {
                    operation: OperationId(1),
                    role: EntityRole {
                        category: RoleCategory::Generated,
                        index: 0,
                    },
                    ancestors: vec![],
                }),
                "face_1" => Some(EntityName {
                    operation: OperationId(1),
                    role: EntityRole {
                        category: RoleCategory::Generated,
                        index: 1,
                    },
                    ancestors: vec![],
                }),
                _ => None,
            }
        }
    }

    let t = ParameterTable::new();
    let resolver = TestResolver;
    let expr = parse_expression_with_resolver(
        "MEASURE(face_0, face_1.dihedral_angle)",
        &t,
        Some(&resolver),
    )
    .unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::DihedralAngle(_, _))
    ));
}

#[test]
fn test_measure_invalid_path() {
    let t = ParameterTable::new();
    let result = parse_expression("MEASURE(body.invalid)", &t);
    assert!(result.is_err());
    let msg = result.unwrap_err().message;
    assert!(msg.contains("invalid MEASURE path"));
}

#[test]
fn test_measure_missing_rparen() {
    let t = ParameterTable::new();
    let result = parse_expression("MEASURE(body.volume", &t);
    assert!(result.is_err());
    assert!(result.unwrap_err().message.contains("expected ')'"));
}

#[test]
fn test_multi_dot_lexer() {
    // Ensure the lexer produces a single Ident token for multi-dot paths
    let t = ParameterTable::new();
    // "body.centroid.x" is a valid body measure path, so this should parse
    let expr = parse_expression("MEASURE(body.centroid.x)", &t).unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::BodyCentroid(
            crate::measure_query::Component::X
        ))
    ));
}

#[test]
fn test_measure_vertex_position() {
    use super::parse_expression_with_resolver;
    use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};

    struct TestResolver;
    impl super::EntityResolver for TestResolver {
        fn resolve_entity(&self, name: &str) -> Option<EntityName> {
            if name == "vertex_0" {
                Some(EntityName {
                    operation: OperationId(1),
                    role: EntityRole {
                        category: RoleCategory::Generated,
                        index: 0,
                    },
                    ancestors: vec![],
                })
            } else {
                None
            }
        }
    }

    let t = ParameterTable::new();
    let resolver = TestResolver;
    let expr = parse_expression_with_resolver("MEASURE(vertex_0.position.z)", &t, Some(&resolver))
        .unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::VertexPosition(
            _,
            crate::measure_query::Component::Z
        ))
    ));
}

#[test]
fn test_measure_edge_length() {
    use super::parse_expression_with_resolver;
    use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};

    struct TestResolver;
    impl super::EntityResolver for TestResolver {
        fn resolve_entity(&self, name: &str) -> Option<EntityName> {
            if name == "edge_0" {
                Some(EntityName {
                    operation: OperationId(1),
                    role: EntityRole {
                        category: RoleCategory::Generated,
                        index: 0,
                    },
                    ancestors: vec![],
                })
            } else {
                None
            }
        }
    }

    let t = ParameterTable::new();
    let resolver = TestResolver;
    let expr =
        parse_expression_with_resolver("MEASURE(edge_0.length)", &t, Some(&resolver)).unwrap();
    assert!(matches!(
        expr,
        Expression::Measure(MeasureQuery::EdgeLength(_))
    ));
}
