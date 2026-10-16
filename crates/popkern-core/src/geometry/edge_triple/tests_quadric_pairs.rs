use super::tests_support::*;
use crate::geometry::edge_triple::{
    AnalyticPair, EdgeAnchors, EdgeTripleOptions, EdgeTripleProvenance, ParameterMap,
    compute_edge_geometry,
};
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

#[test]
fn sphere_cone_axial_analytic_circle() {
    // Stage 2b.18: unit sphere centred at apex of cone (axis (0,0,1),
    // half-angle 30°). Solving v²·sec²θ − 2dv + (d² − 1) = 0 with
    // d = 0 gives v = cos(30°) = √3/2; intersection circle radius
    // = v·tan(30°) = 1/2.
    let sphere = unit_sphere_at(Point3::new(0.0, 0.0, 0.0));
    let cone = z_axis_cone(std::f64::consts::FRAC_PI_6);
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let v = 3.0_f64.sqrt() / 2.0;
    let r_circ = 0.5;
    let start = Point3::new(r_circ, 0.0, v);
    let end = Point3::new(0.0, r_circ, v);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&sphere, &cone, anchors, None, &tol, &opts)
        .expect("sphere×cone axial analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::SphereCone,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Arc { .. }
    ));
    assert_guarantees(&triple, &sphere, &cone, start, end, &tol);
}

#[test]
fn cone_cone_coaxial_analytic_circle() {
    // Stage 2b.19: cone A apex at origin, half-angle 30°; cone B
    // apex at (0, 0, 1), half-angle 60°. Solving v_a·tan30 =
    // (v_a−1)·tan60 gives v_a · (1/√3 − √3) = -√3, so
    // v_a = √3/(√3 − 1/√3) = √3/(2/√3) = 3/2. circle_radius =
    // (3/2)·tan(30°) = √3/2. Centre at (0, 0, 1.5).
    let cone_a = z_axis_cone(std::f64::consts::FRAC_PI_6); // apex (0,0,0), 30°
    let cone_b = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 1.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: std::f64::consts::FRAC_PI_3, // 60°
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let r_circ = 3.0_f64.sqrt() / 2.0;
    let start = Point3::new(r_circ, 0.0, 1.5);
    let end = Point3::new(0.0, r_circ, 1.5);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&cone_a, &cone_b, anchors, None, &tol, &opts)
        .expect("cone×cone coaxial analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::ConeCone,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Arc { .. }
    ));
    assert_guarantees(&triple, &cone_a, &cone_b, start, end, &tol);
}

#[test]
fn cyl_cone_offset_apex_quartic_curve() {
    // Stage 2b.25: cylinder along z-axis, cone with parallel axis
    // but apex offset off the cyl axis. Coaxial dispatcher returns
    // None; general dispatcher samples the closed-form quadratic in
    // v at each cone-u angle.
    //
    // Cyl: origin (0,0,0), axis (0,0,1), radius 1.
    // Cone: apex (0.5, 0, 0), axis (0,0,1), half-angle 30°.
    //
    // Cylinder constraint at cone point (0.5 + v/√3·cos u, v/√3·sin u, v):
    //   (0.5 + v/√3·cos u)² + (v/√3·sin u)² = 1
    //   ⇒ v² + v·√3·cos u − 2.25 = 0
    //   ⇒ v = (−√3·cos u + √(3·cos²u + 9)) / 2  (positive root)
    let cyl = z_axis_cyl(1.0);
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.5, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: std::f64::consts::FRAC_PI_6,
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let v_at = |u: f64| -> f64 {
        let cu = u.cos();
        (-3.0_f64.sqrt() * cu + (3.0 * cu * cu + 9.0).sqrt()) / 2.0
    };
    let pt_at = |u: f64| -> Point3 {
        let v = v_at(u);
        let r = v / 3.0_f64.sqrt();
        Point3::new(0.5 + r * u.cos(), r * u.sin(), v)
    };
    let start = pt_at(0.0);
    let end = pt_at(std::f64::consts::FRAC_PI_2);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&cyl, &cone, anchors, None, &tol, &opts)
        .expect("cyl×cone non-coaxial dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::CylinderCone,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Nurbs(_)
    ));
    assert_guarantees(&triple, &cyl, &cone, start, end, &tol);
}

#[test]
fn cyl_cone_coaxial_analytic_circle() {
    // Stage 2b.17: cylinder of radius 1 along z-axis, cone with apex
    // at origin, axis (0,0,1), half-angle 30°. The intersection is a
    // single circle at height v_cone = 1 / tan(30°) = √3, radius 1.
    let cyl = z_axis_cyl(1.0);
    let cone = z_axis_cone(std::f64::consts::FRAC_PI_6);
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let h = 3.0_f64.sqrt();
    let start = Point3::new(1.0, 0.0, h);
    let end = Point3::new(0.0, 1.0, h);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&cyl, &cone, anchors, None, &tol, &opts)
        .expect("cyl×cone coaxial analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::CylinderCone,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Arc { .. }
    ));
    assert_guarantees(&triple, &cyl, &cone, start, end, &tol);
}

#[test]
fn sphere_cone_non_axial_quartic_curve() {
    // Stage 2b.24: sphere centred OFF the cone axis. The intersection
    // is a genuine 4th-degree curve. The dispatcher samples via the
    // closed-form quadratic in v at each cone-u angle and emits a
    // NURBS curve_3d.
    //
    // Cone: apex (0,0,0), axis (0,0,1), half-angle 30°.
    // Sphere: centre (1, 0, 0), radius √2. Solving:
    //   |(v/√3·cos u − 1, v/√3·sin u, v) − (0,0,0)|² = 2
    //   ⇒ 4v²/3 − (2v/√3)·cos u − 1 = 0
    //   ⇒ v = (√3/4)·(cos u + √(cos²u + 4))     (positive branch)
    let cone = z_axis_cone(std::f64::consts::FRAC_PI_6);
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(1.0, 0.0, 0.0),
        radius: 2.0_f64.sqrt(),
        axis: crate::geometry::Vec3::Z,
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let v_at = |u: f64| -> f64 {
        let cu = u.cos();
        (3.0_f64.sqrt() / 4.0) * (cu + (cu * cu + 4.0).sqrt())
    };
    let pt_at = |u: f64| -> Point3 {
        let v = v_at(u);
        let r = v / 3.0_f64.sqrt();
        Point3::new(r * u.cos(), r * u.sin(), v)
    };
    let start = pt_at(0.0);
    let end = pt_at(std::f64::consts::FRAC_PI_3);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&sphere, &cone, anchors, None, &tol, &opts)
        .expect("sphere×cone non-axial dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::SphereCone,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    // Genuinely quartic ⇒ NURBS curve_3d.
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Nurbs(_)
    ));
    assert_guarantees(&triple, &sphere, &cone, start, end, &tol);
}

#[test]
fn cyl_cyl_perpendicular_viviani_curve() {
    // Stage 2b.23: two unit cylinders with perpendicular intersecting
    // axes — the classical Viviani's curve. The intersection is
    // genuinely 4th-order and has no analytic Line/Arc/Ellipse/Para/
    // Hyper representation; the dispatcher emits a NURBS curve_3d
    // sampled along the implicit-pair solution.
    let cyl_a = z_axis_cyl(1.0);
    let cyl_b = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(1.0, 0.0, 0.0),
        radius: 1.0,
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    // Anchors on the upper-z branch: at u_cyl_a = 0 and u_cyl_a = π/4.
    // x = cos u, y = sin u, z = +sqrt(1 − sin²u) = +cos u (Viviani:
    // z = ±x for equal radii). At u=0: (1, 0, 1). At u=π/4:
    // (√2/2, √2/2, +sqrt(1 − 1/2)) = (√2/2, √2/2, √2/2).
    let r2 = 2.0_f64.sqrt() / 2.0;
    let start = Point3::new(1.0, 0.0, 1.0);
    let end = Point3::new(r2, r2, r2);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&cyl_a, &cyl_b, anchors, None, &tol, &opts)
        .expect("cyl×cyl perpendicular dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::CylinderCylinder,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    // 4th-order quartic ⇒ NURBS curve_3d (no analytic primitive).
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Nurbs(_)
    ));
    assert_guarantees(&triple, &cyl_a, &cyl_b, start, end, &tol);
}

#[test]
fn cyl_cyl_parallel_analytic_line() {
    // Two unit cylinders along z-axis, centres at (0,0,0) and (1,0,0)
    // (axis offset 1.0 < r1+r2 = 2.0). Two intersection lines at
    // x = 0.5, y = ±√(1 - 0.25) = ±√0.75.
    let cyl_a = z_axis_cyl(1.0);
    let cyl_b = SurfaceGeometry::Cylinder {
        origin: Point3::new(1.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let h = (1.0_f64 - 0.25).sqrt();
    let start = Point3::new(0.5, h, 0.0);
    let end = Point3::new(0.5, h, 2.0);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&cyl_a, &cyl_b, anchors, None, &tol, &opts)
        .expect("cyl×cyl parallel analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::CylinderCylinder,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Line { .. }
    ));
    assert_guarantees(&triple, &cyl_a, &cyl_b, start, end, &tol);
}

#[test]
fn cyl_sphere_axial_analytic_circle() {
    // Z-axis cylinder of radius 1; sphere at origin, radius √2.
    // Intersection: two circles at z = ±1, each of radius 1.
    let cyl = z_axis_cyl(1.0);
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 2.0_f64.sqrt(),
        axis: crate::geometry::Vec3::Z,
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    // Anchors on the z=+1 circle.
    let start = Point3::new(1.0, 0.0, 1.0);
    let end = Point3::new(-1.0, 0.0, 1.0);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&cyl, &sphere, anchors, None, &tol, &opts)
        .expect("cyl×sphere axial analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::CylinderSphere,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert_guarantees(&triple, &cyl, &sphere, start, end, &tol);
}

#[test]
fn cone_cone_parallel_offset_levin_analytic_hyperbola() {
    // Stage 2b.34b: cone × cone with parallel-but-offset axes.
    //
    // cone_a: apex (0, 0, 0), axis (0, 0, 1), half-angle 45° → x² + y² = z²
    // cone_b: apex (1, 0, 0), axis (0, 0, 1), half-angle 45° → (x-1)² + y² = z²
    // Subtracting: 2x − 1 = 0 → x = 1/2 (single plane).
    // Plane x = 1/2 cuts each cone in z² − y² = 1/4 — a hyperbola.
    let cone_a = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: std::f64::consts::FRAC_PI_4,
    };
    let cone_b = SurfaceGeometry::Cone {
        apex: Point3::new(1.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: std::f64::consts::FRAC_PI_4,
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    // Anchors on the upper-nappe branch of the hyperbola in x = 1/2.
    let start = Point3::new(0.5, 0.0, 0.5);
    let end = Point3::new(0.5, 1.0, 1.25_f64.sqrt());
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&cone_a, &cone_b, anchors, None, &tol, &opts)
        .expect("cone × cone parallel-offset via Levin");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::ConeCone,
            ..
        }
    ));
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Hyperbola { .. }
    ));
    assert_guarantees(&triple, &cone_a, &cone_b, start, end, &tol);
}

#[test]
fn cyl_cyl_skew_levin_analytic_ellipse() {
    // Stage 2b.34a: Levin's algorithm wired in, decompose_2x2 fixed
    // to return parameterisation-preserving phase.
    //
    // Two unit cylinders through the origin with axes at 45° (skew —
    // not parallel, not perpendicular):
    //   cyl A: axis (0, 0, 1), radius 1 → x² + y² = 1
    //   cyl B: axis (1, 0, 1)/√2, radius 1 → (x−z)²/2 + y² = 1
    //
    // Subtracting: x² + 2xz − z² = 0 → x = (√2 − 1)·z OR x = −(√2 + 1)·z.
    // Two real plane fibres in the pencil; Levin picks the one
    // through the anchors (z=±1 on plane 1).
    let cyl_z = z_axis_cyl(1.0);
    let inv_root2 = 1.0 / 2.0_f64.sqrt();
    let cyl_45 = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(inv_root2, 0.0, inv_root2),
        radius: 1.0,
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let s2m1 = 2.0_f64.sqrt() - 1.0;
    let y_at = (2.0 * 2.0_f64.sqrt() - 2.0).sqrt();
    let start = Point3::new(s2m1, y_at, 1.0);
    let end = Point3::new(-s2m1, y_at, -1.0);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&cyl_z, &cyl_45, anchors, None, &tol, &opts)
        .expect("cyl × cyl skew via Levin");
    let EdgeTripleProvenance::AnalyticIntersection {
        kind, cardinality, ..
    } = &triple.provenance
    else {
        panic!("expected AnalyticIntersection, got {:?}", triple.provenance);
    };
    assert_eq!(*kind, AnalyticPair::CylinderCylinder);
    // Stage 2b.37: Levin populates the cardinality report — 2 plane
    // fibres in the pencil, both yielding ellipses.
    let card = cardinality
        .as_ref()
        .expect("Levin should populate cardinality");
    assert_eq!(card.total, 2);
    use crate::geometry::edge_triple::types::BranchKind;
    for k in &card.kinds {
        assert_eq!(*k, BranchKind::Ellipse);
    }
    // Stage 2b.42: skew cylinders cross transversely, so the Levin
    // path must NOT tag this as tangent.
    assert!(!card.tangent, "skew cyl×cyl is a transverse intersection");
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Ellipse { .. }
    ));
    assert_guarantees(&triple, &cyl_z, &cyl_45, start, end, &tol);
}

#[test]
fn sphere_sphere_polar_analytic_circle() {
    // Two unit spheres at z=0 and z=1, centres along the polar axis.
    // Intersection: a circle at z=0.5 with radius √(1 - 0.25) ≈ 0.866.
    let a = unit_sphere_at(Point3::new(0.0, 0.0, 0.0));
    let b = unit_sphere_at(Point3::new(0.0, 0.0, 1.0));
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let r = (1.0_f64 - 0.25).sqrt();
    let start = Point3::new(r, 0.0, 0.5);
    let end = Point3::new(-r, 0.0, 0.5);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&a, &b, anchors, None, &tol, &opts)
        .expect("sphere×sphere analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::SphereSphere,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert_guarantees(&triple, &a, &b, start, end, &tol);
}
