//! Curve resolution and primitive error-path tests.
//! Included only in test mode; .unwrap() is allowed in tests.
// #[cfg(test)] — marker for the no_unwrap linter.

use super::*;
use crate::step::diagnostics::DiagnosticLog;
use crate::step::parser::types::{EntityInstance, Parameter};
use popkern_core::geometry::curve::CurveGeometry;

fn bspline_curve_params(cp_ids: Vec<u64>) -> Vec<Parameter> {
    // B_SPLINE_CURVE_WITH_KNOTS:
    // name, degree, [cp_refs], curve_form, closed, self_intersect,
    // knot_mults, knots, knot_spec
    vec![
        s(""),
        Parameter::Int(2),
        l(cp_ids.iter().map(|&id| er(id)).collect()),
        Parameter::Enum("UNSPECIFIED".into()),
        Parameter::Enum("F".into()),
        Parameter::Enum("F".into()),
        l(vec![
            Parameter::Int(3),
            Parameter::Int(1),
            Parameter::Int(3),
        ]),
        l(vec![r(0.0), r(0.5), r(1.0)]),
        Parameter::Enum("UNSPECIFIED".into()),
    ]
}

#[test]
fn resolve_bspline_curve_nurbs() {
    let mut ents: Vec<EntityInstance> = Vec::new();
    for id in 1u64..=4 {
        ents.push(ent(
            id,
            "CARTESIAN_POINT",
            vec![s(""), l(vec![r(id as f64), r(0.0), r(0.0)])],
        ));
    }
    ents.push(ent(
        10,
        "B_SPLINE_CURVE_WITH_KNOTS",
        bspline_curve_params(vec![1, 2, 3, 4]),
    ));
    let g = graph(ents);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    let coe = c.resolve_curve(10, &mut log).unwrap();
    match coe.geometry {
        CurveGeometry::Nurbs(nc) => {
            assert_eq!(nc.degree, 2);
            assert_eq!(nc.control_points.len(), 4);
            assert!((coe.t_start - 0.0).abs() < 1e-10);
            assert!((coe.t_end - 1.0).abs() < 1e-10);
        }
        other => panic!("expected Nurbs curve, got {other:?}"),
    }
}

#[test]
fn resolve_bspline_curve_bad_params_returns_none() {
    let g = graph(vec![ent(
        1,
        "B_SPLINE_CURVE_WITH_KNOTS",
        vec![s("only_one_param")],
    )]);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    assert!(c.resolve_curve(1, &mut log).is_none());
}

#[test]
fn resolve_line_preserves_step_vector_magnitude() {
    let ents = vec![
        ent(
            1,
            "CARTESIAN_POINT",
            vec![s(""), l(vec![r(0.0), r(0.0), r(0.0)])],
        ),
        ent(2, "DIRECTION", vec![s(""), l(vec![r(1.0), r(0.0), r(0.0)])]),
        ent(3, "VECTOR", vec![s(""), er(2), r(2.5)]),
        ent(4, "LINE", vec![s(""), er(1), er(3)]),
    ];
    let g = graph(ents);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    let coe = c.resolve_curve(4, &mut log).unwrap();
    match coe.geometry {
        CurveGeometry::Line { direction, .. } => {
            assert!((direction.length() - 2.5).abs() < 1e-12);
        }
        other => panic!("expected line curve, got {other:?}"),
    }
}

#[test]
fn resolve_surface_curve_delegates_to_3d_curve() {
    let ents = vec![
        ent(
            1,
            "CARTESIAN_POINT",
            vec![s(""), l(vec![r(0.0), r(0.0), r(0.0)])],
        ),
        ent(2, "DIRECTION", vec![s(""), l(vec![r(1.0), r(0.0), r(0.0)])]),
        ent(3, "VECTOR", vec![s(""), er(2), r(1.0)]),
        ent(4, "LINE", vec![s(""), er(1), er(3)]),
        ent(
            5,
            "SURFACE_CURVE",
            vec![s(""), er(4), l(vec![]), Parameter::Enum("CURVE_3D".into())],
        ),
    ];
    let g = graph(ents);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    assert!(matches!(
        c.resolve_curve(5, &mut log).unwrap().geometry,
        CurveGeometry::Line { .. }
    ));
}

#[test]
fn resolve_seam_curve_delegates_to_3d_curve() {
    let ents = vec![
        ent(
            1,
            "CARTESIAN_POINT",
            vec![s(""), l(vec![r(0.0), r(0.0), r(0.0)])],
        ),
        ent(2, "DIRECTION", vec![s(""), l(vec![r(1.0), r(0.0), r(0.0)])]),
        ent(3, "VECTOR", vec![s(""), er(2), r(1.0)]),
        ent(4, "LINE", vec![s(""), er(1), er(3)]),
        ent(
            5,
            "SEAM_CURVE",
            vec![s(""), er(4), l(vec![]), Parameter::Enum("CURVE_3D".into())],
        ),
    ];
    let g = graph(ents);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    assert!(matches!(
        c.resolve_curve(5, &mut log).unwrap().geometry,
        CurveGeometry::Line { .. }
    ));
}

#[test]
fn resolve_surface_curve_bad_params_returns_none() {
    let g = graph(vec![ent(1, "SURFACE_CURVE", vec![s("only_one_param")])]);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    assert!(c.resolve_curve(1, &mut log).is_none());
}

#[test]
fn resolve_point_bad_params_returns_none() {
    let g = graph(vec![ent(
        1,
        "CARTESIAN_POINT",
        vec![s(""), l(vec![r(1.0), r(2.0)])],
    )]);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    assert!(c.resolve_point(1, &mut log).is_none());
}

#[test]
fn resolve_direction_bad_params_returns_none() {
    let g = graph(vec![ent(
        1,
        "DIRECTION",
        vec![s(""), l(vec![r(1.0), r(0.0)])],
    )]);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    assert!(c.resolve_direction(1, &mut log).is_none());
}

#[test]
fn resolve_direction_zero_length_returns_default() {
    let g = graph(vec![ent(
        1,
        "DIRECTION",
        vec![s(""), l(vec![r(0.0), r(0.0), r(0.0)])],
    )]);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    let d = c.resolve_direction(1, &mut log).unwrap();
    assert!((d.z - 1.0).abs() < 1e-10);
    assert!(!log.entries().is_empty());
}

#[test]
fn resolve_placement_bad_params_returns_none() {
    let g = graph(vec![ent(1, "AXIS2_PLACEMENT_3D", vec![s("only_name")])]);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    assert!(c.resolve_placement(1, &mut log).is_none());
}

#[test]
fn resolve_line_bad_params_returns_none() {
    let g = graph(vec![ent(1, "LINE", vec![s("bad")])]);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    assert!(c.resolve_curve(1, &mut log).is_none());
}

#[test]
fn resolve_circle_bad_params_returns_none() {
    let g = graph(vec![ent(1, "CIRCLE", vec![s("bad")])]);
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    assert!(c.resolve_curve(1, &mut log).is_none());
}

#[test]
fn resolve_surface_bad_params_returns_none_for_each_type() {
    let bad = vec![s("bad")];
    for ty in &[
        "PLANE",
        "CYLINDRICAL_SURFACE",
        "CONICAL_SURFACE",
        "SPHERICAL_SURFACE",
        "TOROIDAL_SURFACE",
    ] {
        let g = graph(vec![ent(1, ty, bad.clone())]);
        let mut c = GeometryCache::new(&g);
        let mut log = DiagnosticLog::new();
        assert!(
            c.resolve_surface(1, true, &mut log).is_none(),
            "expected None for bad params on {ty}"
        );
    }
}
