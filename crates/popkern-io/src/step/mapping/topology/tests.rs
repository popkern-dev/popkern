use super::*;
use crate::step::parser::types::{EntityInstance, Parameter};
use popkern_core::geometry::Vec3;

fn make_entity(id: EntityId, type_name: &str, params: Vec<Parameter>) -> EntityInstance {
    EntityInstance {
        id,
        type_name: type_name.into(),
        params,
        complex_parts: vec![],
    }
}

fn p_str(s: &str) -> Parameter {
    Parameter::String(s.into())
}

fn p_ref(id: u64) -> Parameter {
    Parameter::EntityRef(id)
}

fn p_enum(s: &str) -> Parameter {
    Parameter::Enum(s.into())
}

fn p_real(v: f64) -> Parameter {
    Parameter::Real(v)
}

fn p_list(items: Vec<Parameter>) -> Parameter {
    Parameter::List(items)
}

#[test]
fn recovers_edge_curve_trim_from_topological_vertices() {
    let curve = CurveOnEdge {
        geometry: CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        },
        t_start: 0.0,
        t_end: std::f64::consts::TAU,
    };
    let start = Point3::new(1.0, 0.0, 0.0);
    let end = Point3::new(0.0, 1.0, 0.0);
    let mut log = DiagnosticLog::new();

    let trimmed =
        recover_edge_curve_trim(99, curve, start, end, true, Tolerance::default(), &mut log);

    assert!(
        trimmed
            .geometry
            .point_at(trimmed.t_start)
            .distance_to(&start)
            <= Tolerance::default().linear()
    );
    assert!(
        trimmed.geometry.point_at(trimmed.t_end).distance_to(&end) <= Tolerance::default().linear()
    );
    assert!(
        (trimmed.t_end - trimmed.t_start).abs() < std::f64::consts::TAU,
        "trim recovery must not leave a full circle parameter range"
    );
    assert!(
        log.errors().next().is_none(),
        "trim recovery should not emit errors: {:?}",
        log.entries()
    );
}

#[test]
fn recovers_periodic_trim_across_zero_angle_with_same_sense() {
    let curve = CurveOnEdge {
        geometry: CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        },
        t_start: 0.0,
        t_end: std::f64::consts::TAU,
    };
    let start_angle = 7.0 * std::f64::consts::TAU / 8.0;
    let end_angle = std::f64::consts::TAU / 8.0;
    let start = curve.geometry.point_at(start_angle);
    let end = curve.geometry.point_at(end_angle);
    let mut log = DiagnosticLog::new();

    let trimmed =
        recover_edge_curve_trim(99, curve, start, end, true, Tolerance::default(), &mut log);

    assert!(
        trimmed.t_end > trimmed.t_start,
        "same-sense periodic trims crossing zero must unwrap forward"
    );
    assert!(
        (trimmed.t_end - trimmed.t_start) < std::f64::consts::PI,
        "same-sense periodic trim should take the short forward span"
    );
    assert!(
        trimmed
            .geometry
            .point_at(trimmed.t_start)
            .distance_to(&start)
            <= Tolerance::default().linear()
    );
    assert!(
        trimmed.geometry.point_at(trimmed.t_end).distance_to(&end) <= Tolerance::default().linear()
    );
}

/// Build a minimal triangle solid from STEP entities and verify topology.
#[test]
fn build_triangle_solid() {
    let entities = vec![
        make_entity(
            1,
            "CARTESIAN_POINT",
            vec![
                p_str(""),
                p_list(vec![p_real(0.0), p_real(0.0), p_real(0.0)]),
            ],
        ),
        make_entity(
            2,
            "CARTESIAN_POINT",
            vec![
                p_str(""),
                p_list(vec![p_real(1.0), p_real(0.0), p_real(0.0)]),
            ],
        ),
        make_entity(
            3,
            "CARTESIAN_POINT",
            vec![
                p_str(""),
                p_list(vec![p_real(0.0), p_real(1.0), p_real(0.0)]),
            ],
        ),
        make_entity(
            4,
            "DIRECTION",
            vec![
                p_str(""),
                p_list(vec![p_real(0.0), p_real(0.0), p_real(1.0)]),
            ],
        ),
        make_entity(
            5,
            "DIRECTION",
            vec![
                p_str(""),
                p_list(vec![p_real(1.0), p_real(0.0), p_real(0.0)]),
            ],
        ),
        make_entity(
            6,
            "AXIS2_PLACEMENT_3D",
            vec![p_str(""), p_ref(1), p_ref(4), p_ref(5)],
        ),
        make_entity(7, "PLANE", vec![p_str(""), p_ref(6)]),
        make_entity(
            8,
            "DIRECTION",
            vec![
                p_str(""),
                p_list(vec![p_real(1.0), p_real(0.0), p_real(0.0)]),
            ],
        ),
        make_entity(9, "VECTOR", vec![p_str(""), p_ref(8), p_real(1.0)]),
        make_entity(
            10,
            "DIRECTION",
            vec![
                p_str(""),
                p_list(vec![p_real(-1.0), p_real(1.0), p_real(0.0)]),
            ],
        ),
        make_entity(11, "VECTOR", vec![p_str(""), p_ref(10), p_real(1.0)]),
        make_entity(
            12,
            "DIRECTION",
            vec![
                p_str(""),
                p_list(vec![p_real(0.0), p_real(-1.0), p_real(0.0)]),
            ],
        ),
        make_entity(13, "VECTOR", vec![p_str(""), p_ref(12), p_real(1.0)]),
        make_entity(14, "LINE", vec![p_str(""), p_ref(1), p_ref(9)]),
        make_entity(15, "LINE", vec![p_str(""), p_ref(2), p_ref(11)]),
        make_entity(16, "LINE", vec![p_str(""), p_ref(3), p_ref(13)]),
        make_entity(20, "VERTEX_POINT", vec![p_str(""), p_ref(1)]),
        make_entity(21, "VERTEX_POINT", vec![p_str(""), p_ref(2)]),
        make_entity(22, "VERTEX_POINT", vec![p_str(""), p_ref(3)]),
        make_entity(
            30,
            "EDGE_CURVE",
            vec![p_str(""), p_ref(20), p_ref(21), p_ref(14), p_enum("T")],
        ),
        make_entity(
            31,
            "EDGE_CURVE",
            vec![p_str(""), p_ref(21), p_ref(22), p_ref(15), p_enum("T")],
        ),
        make_entity(
            32,
            "EDGE_CURVE",
            vec![p_str(""), p_ref(22), p_ref(20), p_ref(16), p_enum("T")],
        ),
        make_entity(
            40,
            "ORIENTED_EDGE",
            vec![
                p_str(""),
                Parameter::Derived,
                Parameter::Derived,
                p_ref(30),
                p_enum("T"),
            ],
        ),
        make_entity(
            41,
            "ORIENTED_EDGE",
            vec![
                p_str(""),
                Parameter::Derived,
                Parameter::Derived,
                p_ref(31),
                p_enum("T"),
            ],
        ),
        make_entity(
            42,
            "ORIENTED_EDGE",
            vec![
                p_str(""),
                Parameter::Derived,
                Parameter::Derived,
                p_ref(32),
                p_enum("T"),
            ],
        ),
        make_entity(
            50,
            "EDGE_LOOP",
            vec![p_str(""), p_list(vec![p_ref(40), p_ref(41), p_ref(42)])],
        ),
        make_entity(
            60,
            "FACE_OUTER_BOUND",
            vec![p_str(""), p_ref(50), p_enum("T")],
        ),
        make_entity(
            70,
            "ADVANCED_FACE",
            vec![p_str(""), p_list(vec![p_ref(60)]), p_ref(7), p_enum("T")],
        ),
        make_entity(80, "CLOSED_SHELL", vec![p_str(""), p_list(vec![p_ref(70)])]),
        make_entity(90, "MANIFOLD_SOLID_BREP", vec![p_str("body"), p_ref(80)]),
    ];

    let mut log = DiagnosticLog::new();
    let graph = EntityGraph::build(entities, &mut log);
    assert!(!log.has_errors(), "graph build errors: {:?}", log.entries());

    let mut geo_cache = GeometryCache::new(&graph);
    let tolerance = Tolerance::default();

    let solids = map_brep_solids(&graph, &mut geo_cache, &mut log, tolerance);
    assert_eq!(solids.len(), 1, "expected 1 solid, got {}", solids.len());

    let (entity_id, brep) = &solids[0];
    assert_eq!(*entity_id, 90);
    assert_eq!(brep.vertices.len(), 3);
    assert_eq!(brep.edges.len(), 3);
    assert_eq!(brep.faces.len(), 1);
    assert_eq!(brep.solids.len(), 1);
    assert_eq!(brep.shells.len(), 1);
}
