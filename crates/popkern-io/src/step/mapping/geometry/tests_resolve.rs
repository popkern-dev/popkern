//! Extended resolution tests for surfaces, curves, and placements.

use super::*;
use crate::step::diagnostics::DiagnosticLog;
use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::{EntityId, EntityInstance, Parameter};
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::surface::SurfaceGeometry;

fn ent(id: EntityId, ty: &str, p: Vec<Parameter>) -> EntityInstance {
    EntityInstance {
        id,
        type_name: ty.into(),
        params: p,
        complex_parts: vec![],
    }
}
fn s(v: &str) -> Parameter {
    Parameter::String(v.into())
}
fn r(v: f64) -> Parameter {
    Parameter::Real(v)
}
fn er(id: u64) -> Parameter {
    Parameter::EntityRef(id)
}
fn l(items: Vec<Parameter>) -> Parameter {
    Parameter::List(items)
}

fn graph(entities: Vec<EntityInstance>) -> EntityGraph {
    EntityGraph::build(entities, &mut DiagnosticLog::new())
}

fn placement_entities(px: f64, py: f64, pz: f64) -> Vec<EntityInstance> {
    vec![
        ent(
            1,
            "CARTESIAN_POINT",
            vec![s(""), l(vec![r(px), r(py), r(pz)])],
        ),
        ent(2, "DIRECTION", vec![s(""), l(vec![r(0.0), r(0.0), r(1.0)])]),
        ent(3, "DIRECTION", vec![s(""), l(vec![r(1.0), r(0.0), r(0.0)])]),
        ent(4, "AXIS2_PLACEMENT_3D", vec![s(""), er(1), er(2), er(3)]),
    ]
}

fn resolve_surf(ty: &str, params: Vec<Parameter>) -> SurfaceGeometry {
    let mut ents = placement_entities(1.0, 2.0, 3.0);
    ents.push(ent(5, ty, params));
    let g = graph(ents);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    c.resolve_surface(5, true, &mut log).unwrap().geometry
}

#[test]
fn resolve_cylindrical_surface() {
    match resolve_surf("CYLINDRICAL_SURFACE", vec![s(""), er(4), r(7.5)]) {
        SurfaceGeometry::Cylinder { origin, radius, .. } => {
            assert!((origin.x - 1.0).abs() < 1e-10);
            assert!((radius - 7.5).abs() < 1e-10);
        }
        other => panic!("expected Cylinder, got {other:?}"),
    }
}

#[test]
fn resolve_conical_surface() {
    match resolve_surf("CONICAL_SURFACE", vec![s(""), er(4), r(2.0), r(0.5)]) {
        SurfaceGeometry::Cone { half_angle, .. } => {
            assert!((half_angle - 0.5).abs() < 1e-10);
        }
        other => panic!("expected Cone, got {other:?}"),
    }
}

#[test]
fn resolve_spherical_surface() {
    match resolve_surf("SPHERICAL_SURFACE", vec![s(""), er(4), r(10.0)]) {
        SurfaceGeometry::Sphere { center, radius } => {
            assert!((center.x - 1.0).abs() < 1e-10);
            assert!((radius - 10.0).abs() < 1e-10);
        }
        other => panic!("expected Sphere, got {other:?}"),
    }
}

#[test]
fn resolve_toroidal_surface() {
    match resolve_surf("TOROIDAL_SURFACE", vec![s(""), er(4), r(8.0), r(2.0)]) {
        SurfaceGeometry::Torus {
            major_radius,
            minor_radius,
            ..
        } => {
            assert!((major_radius - 8.0).abs() < 1e-10);
            assert!((minor_radius - 2.0).abs() < 1e-10);
        }
        other => panic!("expected Torus, got {other:?}"),
    }
}

#[test]
fn resolve_surface_unsupported_type() {
    let g = graph(vec![ent(1, "OFFSET_SURFACE", vec![s(""), er(99)])]);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    assert!(c.resolve_surface(1, true, &mut log).is_none());
}

#[test]
fn resolve_line_curve() {
    let ents = vec![
        ent(
            1,
            "CARTESIAN_POINT",
            vec![s(""), l(vec![r(0.0), r(0.0), r(0.0)])],
        ),
        ent(2, "DIRECTION", vec![s(""), l(vec![r(1.0), r(0.0), r(0.0)])]),
        ent(3, "VECTOR", vec![s(""), er(2), r(1.0)]),
        ent(4, "LINE", vec![s(""), er(1), er(3)]),
    ];
    let g = graph(ents);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    assert!(matches!(
        c.resolve_curve(4, &mut log).unwrap().geometry,
        CurveGeometry::Line { .. }
    ));
}

#[test]
fn resolve_circle_curve() {
    let mut ents = placement_entities(0.0, 0.0, 0.0);
    ents.push(ent(5, "CIRCLE", vec![s(""), er(4), r(5.0)]));
    let g = graph(ents);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    match c.resolve_curve(5, &mut log).unwrap().geometry {
        CurveGeometry::Circle { radius, .. } => assert!((radius - 5.0).abs() < 1e-10),
        other => panic!("expected Circle, got {other:?}"),
    }
}

#[test]
fn resolve_curve_unsupported_and_missing() {
    let mut log = DiagnosticLog::new();
    let g1 = graph(vec![ent(1, "ELLIPSE", vec![s(""), er(99), r(5.0), r(3.0)])]);
    assert!(GeometryCache::new(&g1).resolve_curve(1, &mut log).is_none());
    let g2 = graph(vec![]);
    assert!(
        GeometryCache::new(&g2)
            .resolve_curve(999, &mut log)
            .is_none()
    );
}

#[test]
fn resolve_placement_missing_axis_defaults() {
    let ents = vec![
        ent(
            1,
            "CARTESIAN_POINT",
            vec![s(""), l(vec![r(0.0), r(0.0), r(0.0)])],
        ),
        ent(
            2,
            "AXIS2_PLACEMENT_3D",
            vec![s(""), er(1), Parameter::Omitted, Parameter::Omitted],
        ),
    ];
    let g = graph(ents);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    let (_o, z, x) = c.resolve_placement(2, &mut log).unwrap();
    assert!((z.z - 1.0).abs() < 1e-10);
    assert!((x.x - 1.0).abs() < 1e-10);
}

#[test]
fn resolve_missing_entities_return_none() {
    let g = graph(vec![]);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    assert!(c.resolve_point(999, &mut log).is_none());
    assert!(c.resolve_direction(999, &mut log).is_none());
    assert!(c.resolve_placement(999, &mut log).is_none());
}
