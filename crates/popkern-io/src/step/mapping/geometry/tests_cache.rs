use super::*;
use crate::step::parser::types::{EntityInstance, Parameter};
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

/// Build the standard origin-Z-X placement (point #1, dir #2 = Z, dir #3 = X, a2p #4).
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

#[test]
fn resolve_cartesian_point() {
    let g = graph(vec![ent(
        1,
        "CARTESIAN_POINT",
        vec![s("p"), l(vec![r(1.0), r(2.0), r(3.0)])],
    )]);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    let pt = c.resolve_point(1, &mut log).map(|p| (p.x, p.y, p.z));
    assert_eq!(pt, Some((1.0, 2.0, 3.0)));
    assert!(c.resolve_point(1, &mut log).is_some()); // cached
}

#[test]
fn resolve_direction_normalized() {
    let g = graph(vec![ent(
        1,
        "DIRECTION",
        vec![s("d"), l(vec![r(0.0), r(0.0), r(2.0)])],
    )]);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    let d = c.resolve_direction(1, &mut log);
    assert!(d.is_some());
    let v = d.unwrap();
    assert!(v.x.abs() < 1e-10 && v.y.abs() < 1e-10 && (v.z - 1.0).abs() < 1e-10);
}

#[test]
fn resolve_axis2_placement() {
    let ents = placement_entities(10.0, 20.0, 30.0);
    // Override point to (10, 20, 30)
    let g = graph(ents);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    let (origin, z, x) = c.resolve_placement(4, &mut log).unwrap();
    assert!((origin.x - 10.0).abs() < 1e-10);
    assert!((z.z - 1.0).abs() < 1e-10);
    assert!((x.x - 1.0).abs() < 1e-10);
}

#[test]
fn resolve_plane_surface() {
    let mut ents = placement_entities(0.0, 0.0, 0.0);
    ents.push(ent(5, "PLANE", vec![s(""), er(4)]));
    let g = graph(ents);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    let sof = c.resolve_surface(5, true, &mut log).unwrap();
    assert!(sof.orientation);
    assert!(matches!(sof.geometry, SurfaceGeometry::Plane { .. }));
}

#[test]
fn expand_knots_works() {
    assert_eq!(
        expand_knots(&[0.0, 1.0], &[3, 3]),
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    );
}
