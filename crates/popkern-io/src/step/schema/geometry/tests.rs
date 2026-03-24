use super::*;
use crate::step::parser::types::Parameter;

fn p_str(s: &str) -> Parameter {
    Parameter::String(s.into())
}
fn p_real(v: f64) -> Parameter {
    Parameter::Real(v)
}
fn p_ref(id: u64) -> Parameter {
    Parameter::EntityRef(id)
}
fn p_enum(s: &str) -> Parameter {
    Parameter::Enum(s.into())
}
fn p_int(v: i64) -> Parameter {
    Parameter::Int(v)
}

#[test]
fn cartesian_point_parse() {
    let params = vec![
        p_str("origin"),
        Parameter::List(vec![p_real(1.0), p_real(2.0), p_real(3.0)]),
    ];
    let pt = CartesianPoint::from_params(&params).unwrap();
    assert_eq!(pt.name, "origin");
    assert!((pt.coords[0] - 1.0).abs() < f64::EPSILON);
    assert!((pt.coords[2] - 3.0).abs() < f64::EPSILON);
}

#[test]
fn direction_parse() {
    let params = vec![
        p_str("z"),
        Parameter::List(vec![p_real(0.0), p_real(0.0), p_real(1.0)]),
    ];
    let d = Direction::from_params(&params).unwrap();
    assert_eq!(d.name, "z");
    assert!((d.ratios[2] - 1.0).abs() < f64::EPSILON);
}

#[test]
fn axis2_placement_3d_parse() {
    let params = vec![p_str(""), p_ref(10), p_ref(11), Parameter::Omitted];
    let a = Axis2Placement3d::from_params(&params).unwrap();
    assert_eq!(a.location, 10);
    assert_eq!(a.axis, Some(11));
    assert_eq!(a.ref_direction, None);
}

#[test]
fn circle_parse() {
    let params = vec![p_str("c1"), p_ref(5), p_real(12.5)];
    let c = StepCircle::from_params(&params).unwrap();
    assert_eq!(c.position, 5);
    assert!((c.radius - 12.5).abs() < f64::EPSILON);
}

#[test]
fn vector_parse() {
    let params = vec![p_str("v1"), p_ref(8), p_real(2.5)];
    let v = StepVector::from_params(&params).unwrap();
    assert_eq!(v.name, "v1");
    assert_eq!(v.orientation, 8);
    assert!((v.magnitude - 2.5).abs() < f64::EPSILON);
}

#[test]
fn axis1_placement_parse() {
    let params = vec![p_str("a1"), p_ref(10), p_ref(11)];
    let a = Axis1Placement::from_params(&params).unwrap();
    assert_eq!(a.name, "a1");
    assert_eq!(a.location, 10);
    assert_eq!(a.axis, Some(11));
}

#[test]
fn axis1_placement_no_axis() {
    let params = vec![p_str("a1"), p_ref(10), Parameter::Omitted];
    let a = Axis1Placement::from_params(&params).unwrap();
    assert_eq!(a.axis, None);
}

#[test]
fn ellipse_parse() {
    let params = vec![p_str("e1"), p_ref(5), p_real(10.0), p_real(6.0)];
    let e = StepEllipse::from_params(&params).unwrap();
    assert_eq!(e.name, "e1");
    assert_eq!(e.position, 5);
    assert!((e.semi_axis_1 - 10.0).abs() < f64::EPSILON);
    assert!((e.semi_axis_2 - 6.0).abs() < f64::EPSILON);
}

#[test]
fn trimmed_curve_parse() {
    let params = vec![
        p_str("tc"),
        p_ref(100),
        Parameter::List(vec![p_ref(50)]),
        Parameter::List(vec![p_real(1.0)]),
        p_enum("T"),
        p_enum("PARAMETER"),
    ];
    let tc = TrimmedCurve::from_params(&params).unwrap();
    assert_eq!(tc.name, "tc");
    assert_eq!(tc.basis_curve, 100);
    assert_eq!(tc.trim_1.len(), 1);
    assert_eq!(tc.trim_2.len(), 1);
    assert!(tc.sense);
    assert_eq!(tc.master_representation, "PARAMETER");
}

#[test]
fn composite_curve_parse() {
    let params = vec![
        p_str("cc"),
        Parameter::List(vec![p_ref(1), p_ref(2), p_ref(3)]),
        p_enum("F"),
    ];
    let cc = CompositeCurve::from_params(&params).unwrap();
    assert_eq!(cc.name, "cc");
    assert_eq!(cc.segments, vec![1, 2, 3]);
    assert!(!cc.self_intersect);
}

#[test]
fn surface_curve_parse() {
    let params = vec![
        p_str("sc"),
        p_ref(10),
        Parameter::List(vec![p_ref(20), p_ref(21)]),
        p_enum("PCURVE_S1"),
    ];
    let sc = SurfaceCurve::from_params(&params).unwrap();
    assert_eq!(sc.name, "sc");
    assert_eq!(sc.curve_3d, 10);
    assert_eq!(sc.associated_geometry, vec![20, 21]);
    assert_eq!(sc.master_representation, "PCURVE_S1");
}

#[test]
fn seam_curve_parse() {
    let params = vec![
        p_str("seam"),
        p_ref(15),
        Parameter::List(vec![p_ref(25)]),
        p_enum("PCURVE_S2"),
    ];
    let sc = SeamCurve::from_params(&params).unwrap();
    assert_eq!(sc.name, "seam");
    assert_eq!(sc.curve_3d, 15);
    assert_eq!(sc.associated_geometry, vec![25]);
    assert_eq!(sc.master_representation, "PCURVE_S2");
}

#[test]
fn pcurve_parse() {
    let params = vec![p_str("pc"), p_ref(40), p_ref(41)];
    let pc = StepPcurve::from_params(&params).unwrap();
    assert_eq!(pc.name, "pc");
    assert_eq!(pc.basis_surface, 40);
    assert_eq!(pc.reference_curve, 41);
}

#[test]
fn bspline_curve_parse() {
    let params = vec![
        p_str("curve"),
        p_int(3),
        Parameter::List(vec![p_ref(1), p_ref(2), p_ref(3), p_ref(4)]),
        p_enum("UNSPECIFIED"),
        p_enum("F"),
        p_enum("F"),
        Parameter::List(vec![p_int(4), p_int(4)]),
        Parameter::List(vec![p_real(0.0), p_real(1.0)]),
        p_enum("UNSPECIFIED"),
    ];
    let bs = BSplineCurveWithKnots::from_params(&params).unwrap();
    assert_eq!(bs.degree, 3);
    assert_eq!(bs.control_points.len(), 4);
    assert!(!bs.closed);
    assert_eq!(bs.knot_multiplicities, vec![4, 4]);
    assert!(bs.weights.is_none());
}
