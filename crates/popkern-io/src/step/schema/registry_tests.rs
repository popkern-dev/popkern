// #[cfg(test)] — included only in test mode; .unwrap() allowed in tests.
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

#[test]
fn resolve_cartesian_point() {
    let entity = make_entity(
        "CARTESIAN_POINT",
        vec![
            Parameter::String("origin".into()),
            Parameter::List(vec![
                Parameter::Real(1.0),
                Parameter::Real(2.0),
                Parameter::Real(3.0),
            ]),
        ],
    );
    let resolved = resolve_entity(&entity).expect("should resolve");
    assert!(matches!(resolved, ResolvedEntity::CartesianPoint(p) if p.coords == [1.0, 2.0, 3.0]));
}

#[test]
fn resolve_manifold_solid_brep() {
    let entity = make_entity(
        "MANIFOLD_SOLID_BREP",
        vec![Parameter::String("body".into()), Parameter::EntityRef(42)],
    );
    let resolved = resolve_entity(&entity).expect("should resolve");
    assert!(matches!(resolved, ResolvedEntity::ManifoldSolidBrep(b) if b.outer == 42));
}

#[test]
fn resolve_unknown_type() {
    let entity = make_entity("FOOBAR", vec![]);
    let resolved = resolve_entity(&entity).expect("should not error for unknown types");
    assert!(matches!(resolved, ResolvedEntity::Unknown(_)));
}

// --- Geometry: points, directions, placements ---

#[test]
fn resolve_direction() {
    let entity = make_entity(
        "DIRECTION",
        vec![
            Parameter::String("z".into()),
            Parameter::List(vec![
                Parameter::Real(0.0),
                Parameter::Real(0.0),
                Parameter::Real(1.0),
            ]),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::Direction(_)));
}

#[test]
fn resolve_vector() {
    let entity = make_entity(
        "VECTOR",
        vec![
            Parameter::String("v".into()),
            Parameter::EntityRef(2),
            Parameter::Real(1.0),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::Vector(_)));
}

#[test]
fn resolve_axis1_placement() {
    let entity = make_entity(
        "AXIS1_PLACEMENT",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::EntityRef(2),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::Axis1Placement(_)));
}

#[test]
fn resolve_axis2_placement_3d() {
    let entity = make_entity(
        "AXIS2_PLACEMENT_3D",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::Omitted,
            Parameter::Omitted,
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::Axis2Placement3d(_)));
}

// --- Geometry: curves ---

#[test]
fn resolve_line() {
    let entity = make_entity(
        "LINE",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::EntityRef(2),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::Line(_)));
}

#[test]
fn resolve_circle() {
    let entity = make_entity(
        "CIRCLE",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::Real(5.0),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::Circle(_)));
}

#[test]
fn resolve_ellipse() {
    let entity = make_entity(
        "ELLIPSE",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::Real(10.0),
            Parameter::Real(5.0),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::Ellipse(_)));
}

#[test]
fn resolve_bspline_curve() {
    let entity = make_entity(
        "B_SPLINE_CURVE_WITH_KNOTS",
        vec![
            Parameter::String("".into()),
            Parameter::Int(3),
            Parameter::List(vec![Parameter::EntityRef(1), Parameter::EntityRef(2)]),
            Parameter::Enum("UNSPECIFIED".into()),
            Parameter::Enum("F".into()),
            Parameter::Enum("F".into()),
            Parameter::List(vec![Parameter::Int(2), Parameter::Int(2)]),
            Parameter::List(vec![Parameter::Real(0.0), Parameter::Real(1.0)]),
            Parameter::Enum("UNSPECIFIED".into()),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::BSplineCurveWithKnots(_)));
}

#[test]
fn resolve_trimmed_curve() {
    let entity = make_entity(
        "TRIMMED_CURVE",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::List(vec![Parameter::Real(0.0)]),
            Parameter::List(vec![Parameter::Real(1.0)]),
            Parameter::Enum("T".into()),
            Parameter::Enum("PARAMETER".into()),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::TrimmedCurve(_)));
}

#[test]
fn resolve_composite_curve() {
    let entity = make_entity(
        "COMPOSITE_CURVE",
        vec![
            Parameter::String("".into()),
            Parameter::List(vec![Parameter::EntityRef(1)]),
            Parameter::Enum("F".into()),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::CompositeCurve(_)));
}

#[test]
fn resolve_surface_curve() {
    let entity = make_entity(
        "SURFACE_CURVE",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::List(vec![Parameter::EntityRef(2)]),
            Parameter::Enum("PCURVE_S1".into()),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::SurfaceCurve(_)));
}

#[test]
fn resolve_seam_curve() {
    let entity = make_entity(
        "SEAM_CURVE",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::List(vec![Parameter::EntityRef(2)]),
            Parameter::Enum("PCURVE_S1".into()),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::SeamCurve(_)));
}

#[test]
fn resolve_pcurve() {
    let entity = make_entity(
        "PCURVE",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::EntityRef(2),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::Pcurve(_)));
}

// --- Geometry: surfaces ---

#[test]
fn resolve_plane() {
    let entity = make_entity(
        "PLANE",
        vec![Parameter::String("".into()), Parameter::EntityRef(1)],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::Plane(_)));
}

#[test]
fn resolve_cylindrical_surface() {
    let entity = make_entity(
        "CYLINDRICAL_SURFACE",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::Real(5.0),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::CylindricalSurface(_)));
}

#[test]
fn resolve_conical_surface() {
    let entity = make_entity(
        "CONICAL_SURFACE",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::Real(3.0),
            Parameter::Real(0.5),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::ConicalSurface(_)));
}

#[test]
fn resolve_spherical_surface() {
    let entity = make_entity(
        "SPHERICAL_SURFACE",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::Real(10.0),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::SphericalSurface(_)));
}

#[test]
fn resolve_toroidal_surface() {
    let entity = make_entity(
        "TOROIDAL_SURFACE",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::Real(10.0),
            Parameter::Real(2.0),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::ToroidalSurface(_)));
}

#[test]
fn resolve_bspline_surface() {
    let entity = make_entity(
        "B_SPLINE_SURFACE_WITH_KNOTS",
        vec![
            Parameter::String("".into()),
            Parameter::Int(2),
            Parameter::Int(2),
            Parameter::List(vec![Parameter::List(vec![
                Parameter::EntityRef(1),
                Parameter::EntityRef(2),
            ])]),
            Parameter::Enum("UNSPECIFIED".into()),
            Parameter::Enum("F".into()),
            Parameter::Enum("F".into()),
            Parameter::Enum("F".into()),
            Parameter::List(vec![Parameter::Int(3), Parameter::Int(3)]),
            Parameter::List(vec![Parameter::Int(3), Parameter::Int(3)]),
            Parameter::List(vec![Parameter::Real(0.0), Parameter::Real(1.0)]),
            Parameter::List(vec![Parameter::Real(0.0), Parameter::Real(1.0)]),
            Parameter::Enum("UNSPECIFIED".into()),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(
        resolved,
        ResolvedEntity::BSplineSurfaceWithKnots(_)
    ));
}

// --- Topology ---

#[test]
fn resolve_vertex_point() {
    let entity = make_entity(
        "VERTEX_POINT",
        vec![Parameter::String("".into()), Parameter::EntityRef(1)],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::VertexPoint(_)));
}

#[test]
fn resolve_edge_curve() {
    let entity = make_entity(
        "EDGE_CURVE",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::EntityRef(2),
            Parameter::EntityRef(3),
            Parameter::Enum("T".into()),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::EdgeCurve(_)));
}

#[test]
fn resolve_oriented_edge() {
    let entity = make_entity(
        "ORIENTED_EDGE",
        vec![
            Parameter::String("".into()),
            Parameter::Omitted,
            Parameter::Omitted,
            Parameter::EntityRef(1),
            Parameter::Enum("T".into()),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::OrientedEdge(_)));
}

#[test]
fn resolve_edge_loop() {
    let entity = make_entity(
        "EDGE_LOOP",
        vec![
            Parameter::String("".into()),
            Parameter::List(vec![Parameter::EntityRef(1), Parameter::EntityRef(2)]),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::EdgeLoop(_)));
}

#[test]
fn resolve_face_bound() {
    let entity = make_entity(
        "FACE_BOUND",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::Enum("T".into()),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::FaceBound(_)));
}

#[test]
fn resolve_face_outer_bound() {
    let entity = make_entity(
        "FACE_OUTER_BOUND",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::Enum("T".into()),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::FaceOuterBound(_)));
}

#[test]
fn resolve_advanced_face() {
    let entity = make_entity(
        "ADVANCED_FACE",
        vec![
            Parameter::String("".into()),
            Parameter::List(vec![Parameter::EntityRef(1)]),
            Parameter::EntityRef(2),
            Parameter::Enum("T".into()),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::AdvancedFace(_)));
}

#[test]
fn resolve_closed_shell() {
    let entity = make_entity(
        "CLOSED_SHELL",
        vec![
            Parameter::String("".into()),
            Parameter::List(vec![Parameter::EntityRef(1)]),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::ClosedShell(_)));
}

#[test]
fn resolve_open_shell() {
    let entity = make_entity(
        "OPEN_SHELL",
        vec![
            Parameter::String("".into()),
            Parameter::List(vec![Parameter::EntityRef(1)]),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::OpenShell(_)));
}

#[test]
fn resolve_brep_with_voids() {
    let entity = make_entity(
        "BREP_WITH_VOIDS",
        vec![
            Parameter::String("".into()),
            Parameter::EntityRef(1),
            Parameter::List(vec![Parameter::EntityRef(2)]),
        ],
    );
    let resolved = resolve_entity(&entity).unwrap();
    assert!(matches!(resolved, ResolvedEntity::BrepWithVoids(_)));
}
