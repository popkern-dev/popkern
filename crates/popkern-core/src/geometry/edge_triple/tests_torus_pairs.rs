use super::tests_support::*;
use crate::geometry::edge_triple::{
    AnalyticPair, EdgeAnchors, EdgeTripleOptions, EdgeTripleProvenance, ParameterMap,
    compute_edge_geometry,
};
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

#[test]
fn cone_torus_axial_analytic_circle() {
    // Stage 2b.28: cone apex on torus axis, axes coincident. Octic
    // intersection factors into ≤ 2 horizontal circles via the
    // axisymmetric collapse argument.
    //
    // Torus R=5, r=1. Cone apex (0,0,0), axis (0,0,1), half-angle
    // arctan(1/2) ≈ 0.4636 (so tan θ = 0.5). With h_axis = 0
    // (apex at torus centre):
    //   cos(v + θ) = (sin θ · 0 − R · cos θ) / r
    //              = -5 · cos(arctan(1/2)) = -5 · (2/√5) = -2√5 ≈ -4.47
    // |rhs| > 1 ⇒ no intersection. So shift apex.
    //
    // Take h_axis = 8 (apex below torus): cone apex (0, 0, -8). Then
    //   sin θ · 8 − R · cos θ = (8/√5) − 5·(2/√5) = (8 − 10)/√5
    //                         = -2/√5
    //   rhs = -2/√5 / 1 ≈ -0.894 — within [-1, 1] ✓
    //
    //   v_± = -θ ± acos(-2/√5)
    //   acos(-2/√5) ≈ 2.678.
    //   θ ≈ 0.4636.
    //   v_+ ≈ -0.4636 + 2.678 ≈ 2.215.
    //   v_- ≈ -0.4636 - 2.678 ≈ -3.142 (= -π).
    //
    // At v = -π (= π): cos = -1, sin = 0. Circle: radius R + r·(-1) = 4
    // at axial height 0 from torus centre. Cone v = h_axis + 0 = 8.
    // Cone radius at v=8 = 8·tan θ = 4 ✓.
    let torus = nurbs_torus_dummy(); // R=5, r=1
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, -8.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: (1.0_f64 / 2.0).atan(),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    // Inner-tube (v=π) circle: radius 4 at z=0.
    let start = Point3::new(4.0, 0.0, 0.0);
    let end = Point3::new(0.0, 4.0, 0.0);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&cone, &torus, anchors, None, &tol, &opts)
        .expect("cone×torus axial dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::ConeTorus,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Arc { .. }
    ));
    if let crate::geometry::curve::CurveGeometry::Arc { radius, .. } = triple.curve_3d {
        assert!(
            (radius - 4.0).abs() < 1e-9,
            "cone×torus circle radius should be 4, got {radius}"
        );
    }
    assert_guarantees(&triple, &cone, &torus, start, end, &tol);
}

#[test]
fn sphere_torus_axial_analytic_circle() {
    // Stage 2b.27: sphere centred on torus axis. The octic
    // intersection factors into ≤ 2 horizontal circles. With torus
    // R=5, r=1, sphere centre (0,0,0) (= torus centre), radius
    // ρ = √26: the constraint
    //   R cos v − h sin v = (ρ² − R² − r² − h²) / (2r)
    // with h=0 reduces to cos v = (26 − 25 − 1)/2 = 0, so v = ±π/2.
    // At v = π/2: axial height = sin(π/2) = 1, circle radius
    // = R + r·cos(π/2) = 5. Sphere check: distance from origin to
    // (5, 0, 1) = √26 ✓.
    let torus = nurbs_torus_dummy(); // R=5, r=1
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 26.0_f64.sqrt(),
        axis: crate::geometry::Vec3::Z,
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    // Anchors on the +z circle (centre (0, 0, 1), radius 5, in
    // xy-plane translated up by 1).
    let start = Point3::new(5.0, 0.0, 1.0);
    let end = Point3::new(0.0, 5.0, 1.0);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&sphere, &torus, anchors, None, &tol, &opts)
        .expect("sphere×torus axial dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::SphereTorus,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Arc { .. }
    ));
    if let crate::geometry::curve::CurveGeometry::Arc { radius, .. } = triple.curve_3d {
        assert!(
            (radius - 5.0).abs() < 1e-9,
            "sphere×torus circle radius should be R = 5, got {radius}"
        );
    }
    assert_guarantees(&triple, &sphere, &torus, start, end, &tol);
}

#[test]
fn plane_torus_oblique_quartic_curve() {
    // Stage 2b.26: torus R=5, r=1; plane through torus centre tilted
    // by α with sin α = 0.1 (TIGHTER than the Villarceau angle of
    // arcsin(r/R) = arcsin(0.2) ≈ 0.2014). The intersection is a
    // 4th-order curve covering all u ∈ [0, TAU] (no u-range
    // restrictions since R·sin α / r = 0.5 < 1).
    //
    // Closed form: cos(v − φ(u)) = −R·sin α·cos u / M(u) where
    // M(u) = √(sin²α·cos²u + cos²α). At u = 0: cos(v − φ) = −0.5,
    // so v = φ(0) + acos(−0.5).
    let torus = nurbs_torus_dummy(); // R=5, r=1
    let sin_a = 0.1_f64;
    let cos_a = (1.0 - sin_a * sin_a).sqrt();
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(sin_a, 0.0, cos_a),
        u_axis: Vec3::new(cos_a, 0.0, -sin_a),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let big_r = 5.0;
    let r = 1.0;
    let v_at = |u: f64, branch_sign: f64| -> f64 {
        let cu = u.cos();
        let x = sin_a * cu;
        let m = (x * x + cos_a * cos_a).sqrt();
        let phi = cos_a.atan2(x);
        let arg = (-big_r * x) / (r * m);
        phi + branch_sign * arg.acos()
    };
    let pt_at = |u: f64| -> Point3 {
        let v = v_at(u, 1.0);
        let cv = v.cos();
        let sv = v.sin();
        Point3::new(
            (big_r + r * cv) * u.cos(),
            (big_r + r * cv) * u.sin(),
            r * sv,
        )
    };
    let start = pt_at(0.0);
    let end = pt_at(std::f64::consts::FRAC_PI_3);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&plane, &torus, anchors, None, &tol, &opts)
        .expect("plane×torus oblique non-Villarceau dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::PlaneTorus,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Nurbs(_)
    ));
    assert_guarantees(&triple, &plane, &torus, start, end, &tol);
}

#[test]
fn plane_torus_villarceau_analytic_circle() {
    // Stage 2b.22: Villarceau plane bitangent to a torus, factoring
    // the quartic into 2 circles of major radius. Torus R=5, r=1.
    // Villarceau angle α = arcsin(1/5) ≈ 0.2014 rad.
    // Plane normal: (sin α, 0, cos α) = (1/5, 0, √24/5).
    let torus = nurbs_torus_dummy(); // R=5, r=1
    let s_a = 0.2_f64; // sin α = 1/5
    let c_a = (1.0_f64 - s_a * s_a).sqrt(); // cos α = √24/5
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(s_a, 0.0, c_a),
        u_axis: Vec3::new(c_a, 0.0, -s_a),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    // Anchors on the +y Villarceau circle (centre (0, 1, 0), R=5).
    // At θ=0 in arc-basis frame: (R, 0, 0) projected — but the arc
    // basis depends on plane normal. Just pick two known points on
    // the +y circle: outer-equator point (0, 6, 0) (= (0, 1, 0) +
    // 5·(0, 1, 0))? Wait, circle centre is (0, 1, 0) and radius 5,
    // lies in the Villarceau plane. So circle points satisfy
    // x·sin α + z·cos α = 0. At (0, 6, 0): 0+0=0 ✓ and distance
    // from centre = (0, 5, 0) magnitude 5 ✓.
    let start = Point3::new(0.0, 6.0, 0.0);
    // Other known point: where x·s_a + z·c_a = 0, distance 5 from
    // (0, 1, 0). Try x = 5·c_a, z = -5·s_a, y = 1: distance from
    // centre = (5·c_a, 0, -5·s_a), magnitude = 5·sqrt(c_a²+s_a²) = 5 ✓.
    let end = Point3::new(5.0 * c_a, 1.0, -5.0 * s_a);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&plane, &torus, anchors, None, &tol, &opts)
        .expect("plane×torus Villarceau analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::PlaneTorus,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Arc { .. }
    ));
    // Major radius circle.
    if let crate::geometry::curve::CurveGeometry::Arc { radius, .. } = triple.curve_3d {
        assert!(
            (radius - 5.0).abs() < 1e-9,
            "Villarceau circle should have radius = R = 5, got {radius}"
        );
    }
    assert_guarantees(&triple, &plane, &torus, start, end, &tol);
}

#[test]
fn cyl_torus_axial_analytic_circle() {
    // Stage 2b.29: cyl axis ∥ torus axis, cyl origin on torus axis.
    // Octic intersection factors into ≤ 2 horizontal circles.
    //
    // Torus R=5, r=1. Cylinder origin (0,0,0) (= torus centre),
    // axis (0,0,1), radius 5 (= R, mid-tube). Then radial_diff = 0,
    // sqrt_radical = 1, h_offset = 0, so v_cyl_± = ±1: two circles
    // at z = ±1, each radius 5. v_torus = ±π/2.
    let torus = nurbs_torus_dummy(); // R=5, r=1
    let cyl = z_axis_cyl(5.0);
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    // Anchors on the +z circle (centre (0, 0, 1), radius 5).
    let start = Point3::new(5.0, 0.0, 1.0);
    let end = Point3::new(0.0, 5.0, 1.0);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&cyl, &torus, anchors, None, &tol, &opts)
        .expect("cyl×torus axial dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::CylinderTorus,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Arc { .. }
    ));
    if let crate::geometry::curve::CurveGeometry::Arc { radius, .. } = triple.curve_3d {
        assert!(
            (radius - 5.0).abs() < 1e-9,
            "cyl×torus circle radius should be 5, got {radius}"
        );
    }
    assert_guarantees(&triple, &cyl, &torus, start, end, &tol);
}

#[test]
fn torus_torus_coaxial_analytic_circle() {
    // Stage 2b.30: two coaxial tori, axes coincident, centres offset
    // axially. The octic intersection collapses to ≤ 2 horizontal
    // circles via the axisymmetric collapse argument.
    //
    // Two unit tubes at R=5, vertical separation h=√3:
    //   C1 = (5, 0), r1 = 1; C2 = (5, √3), r2 = 1
    //   a = √3/2, h_off = 1/2 (in (X, Y) two-circle intersection)
    //   solutions: (X, Y) = (5 ± 1/2, √3/2)
    // → two horizontal circles at z = √3/2, radii 4.5 and 5.5.
    let torus_a = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 5.0,
        minor_radius: 1.0,
    };
    let torus_b = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 3.0_f64.sqrt()),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 5.0,
        minor_radius: 1.0,
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let z = 3.0_f64.sqrt() / 2.0;
    // Anchor on the outer (X = 5.5) circle.
    let start = Point3::new(5.5, 0.0, z);
    let end = Point3::new(0.0, 5.5, z);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&torus_a, &torus_b, anchors, None, &tol, &opts)
        .expect("torus×torus coaxial dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::TorusTorus,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Arc { .. }
    ));
    if let crate::geometry::curve::CurveGeometry::Arc { radius, .. } = triple.curve_3d {
        assert!(
            (radius - 5.5).abs() < 1e-9,
            "torus×torus outer circle radius should be 5.5, got {radius}"
        );
    }
    assert_guarantees(&triple, &torus_a, &torus_b, start, end, &tol);
}
