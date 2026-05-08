//! Extended resolution tests for surfaces, curves, and placements.
//! Included only in test mode; .unwrap() is allowed in tests.
// #[cfg(test)] — marker for the no_unwrap linter.

use super::*;
use crate::step::diagnostics::DiagnosticLog;
use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::{EntityId, EntityInstance, Parameter};
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::surface::SurfaceGeometry;

#[path = "tests_resolve_curves.rs"]
mod tests_resolve_curves;

pub(super) fn ent(id: EntityId, ty: &str, p: Vec<Parameter>) -> EntityInstance {
    EntityInstance {
        id,
        type_name: ty.into(),
        params: p,
        complex_parts: vec![],
    }
}
fn ent_complex(
    id: EntityId,
    ty: &str,
    p: Vec<Parameter>,
    complex_parts: Vec<(&str, Vec<Parameter>)>,
) -> EntityInstance {
    EntityInstance {
        id,
        type_name: ty.into(),
        params: p,
        complex_parts: complex_parts
            .into_iter()
            .map(|(name, params)| (name.into(), params))
            .collect(),
    }
}
pub(super) fn s(v: &str) -> Parameter {
    Parameter::String(v.into())
}
pub(super) fn r(v: f64) -> Parameter {
    Parameter::Real(v)
}
pub(super) fn er(id: u64) -> Parameter {
    Parameter::EntityRef(id)
}
pub(super) fn l(items: Vec<Parameter>) -> Parameter {
    Parameter::List(items)
}

pub(super) fn graph(entities: Vec<EntityInstance>) -> EntityGraph {
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

// ---------------------------------------------------------------------------
// B-spline surface
// ---------------------------------------------------------------------------

fn bspline_surface_params(cp_ids: Vec<Vec<u64>>) -> Vec<Parameter> {
    // B_SPLINE_SURFACE_WITH_KNOTS:
    // name, u_degree, v_degree, [[cp...]], form, u_closed, v_closed,
    // self_intersect, u_mults, v_mults, u_knots, v_knots, knot_spec
    let cp_list = Parameter::List(
        cp_ids
            .iter()
            .map(|row| Parameter::List(row.iter().map(|&id| er(id)).collect()))
            .collect(),
    );
    vec![
        s(""),
        Parameter::Int(2),
        Parameter::Int(2),
        cp_list,
        Parameter::Enum("UNSPECIFIED".into()),
        Parameter::Enum("F".into()),
        Parameter::Enum("F".into()),
        Parameter::Enum("F".into()),
        l(vec![Parameter::Int(3), Parameter::Int(3)]),
        l(vec![Parameter::Int(3), Parameter::Int(3)]),
        l(vec![r(0.0), r(1.0)]),
        l(vec![r(0.0), r(1.0)]),
        Parameter::Enum("UNSPECIFIED".into()),
    ]
}

#[test]
fn resolve_bspline_surface() {
    // 2x2 grid of control points (IDs 10..=13), placement at IDs 1-4, surface at 5
    let mut ents = placement_entities(0.0, 0.0, 0.0);
    // 4 control points
    for id in 10u64..=13 {
        ents.push(ent(
            id,
            "CARTESIAN_POINT",
            vec![s(""), l(vec![r(id as f64), r(0.0), r(0.0)])],
        ));
    }
    ents.push(ent(
        20,
        "B_SPLINE_SURFACE_WITH_KNOTS",
        bspline_surface_params(vec![vec![10, 11], vec![12, 13]]),
    ));
    let g = graph(ents);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    let sof = c.resolve_surface(20, false, &mut log).unwrap();
    assert!(!sof.orientation);
    match sof.geometry {
        SurfaceGeometry::Nurbs(ns) => {
            assert_eq!(ns.degree_u, 2);
            assert_eq!(ns.degree_v, 2);
            assert_eq!(ns.num_u, 2);
            assert_eq!(ns.num_v, 2);
            assert_eq!(ns.control_points.len(), 4);
            // default weights all 1.0
            assert!(ns.weights.iter().all(|&w| (w - 1.0).abs() < 1e-10));
        }
        other => panic!("expected Nurbs surface, got {other:?}"),
    }
}

#[test]
fn resolve_rational_bspline_surface_weights_from_complex_part() {
    let mut ents = placement_entities(0.0, 0.0, 0.0);
    for id in 10u64..=13 {
        ents.push(ent(
            id,
            "CARTESIAN_POINT",
            vec![s(""), l(vec![r(id as f64), r(0.0), r(0.0)])],
        ));
    }
    ents.push(ent_complex(
        20,
        "B_SPLINE_SURFACE_WITH_KNOTS",
        bspline_surface_params(vec![vec![10, 11], vec![12, 13]]),
        vec![(
            "RATIONAL_B_SPLINE_SURFACE",
            vec![l(vec![l(vec![r(1.0), r(0.5)]), l(vec![r(0.75), r(1.25)])])],
        )],
    ));
    let g = graph(ents);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    let sof = c.resolve_surface(20, true, &mut log).unwrap();
    match sof.geometry {
        SurfaceGeometry::Nurbs(ns) => {
            assert_eq!(ns.weights, vec![1.0, 0.5, 0.75, 1.25]);
        }
        other => panic!("expected Nurbs surface, got {other:?}"),
    }
}

#[test]
fn resolve_bspline_surface_bad_params_returns_none() {
    // Give the entity wrong params (missing fields) so from_params fails
    let g = graph(vec![ent(
        1,
        "B_SPLINE_SURFACE_WITH_KNOTS",
        vec![s("only_one_param")],
    )]);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    assert!(c.resolve_surface(1, true, &mut log).is_none());
}
