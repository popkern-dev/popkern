use crate::plane::SketchPlane;
use crate::profile::{AnalyticProfile, PathSegment, ProfileError, ProfileShape};
use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use std::f64::consts::PI;

fn circle(r: f64) -> AnalyticProfile<f64> {
    AnalyticProfile::new(ProfileShape::Circle { radius: r })
}

/// A rectangle whose front side bulges as an arc through `via`
/// (the mixed-path class of the part API), CCW or CW by `via`'s
/// side.
fn bulged(via: [f64; 2]) -> AnalyticProfile<f64> {
    AnalyticProfile::new(ProfileShape::Path {
        start: [-30.0, 0.0],
        segments: vec![
            PathSegment::Line { to: [30.0, 0.0] },
            PathSegment::Arc {
                via,
                to: [-30.0, 0.0],
            },
        ],
    })
}

fn max_radial_error(curve: &popkern_core::geometry::NurbsCurve, c: Point3, r: f64) -> f64 {
    let (t0, t1) = curve.parameter_range();
    (0..=500)
        .map(|i| {
            let p = curve.point_at(t0 + (t1 - t0) * i as f64 / 500.0);
            (p.distance_to(&c) - r).abs()
        })
        .fold(0.0, f64::max)
}

#[test]
fn circle_wire_is_the_standard_nine_point_form_lifted_to_the_plane() {
    let plane = SketchPlane::new(
        Point3::new(0.0, 0.0, 2.5),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
    )
    .unwrap();
    let w = circle(0.25).local_wire(&plane).unwrap();
    assert_eq!((w.degree, w.control_points.len()), (2, 9));
    assert!(max_radial_error(&w, Point3::new(0.0, 0.0, 2.5), 0.25) < 1e-12);
    // Placement moves the wire; local keeps it at the origin.
    let placed = AnalyticProfile {
        center: [1.0, 0.0],
        ..circle(0.25)
    };
    let pw = placed.placed_wire(&plane).unwrap();
    assert!(max_radial_error(&pw, Point3::new(1.0, 0.0, 2.5), 0.25) < 1e-12);
    let lw = placed.local_wire(&plane).unwrap();
    assert!(max_radial_error(&lw, Point3::new(0.0, 0.0, 2.5), 0.25) < 1e-12);
}

#[test]
fn rect_polygon_and_path_wires_are_closed_and_exact() {
    let plane = SketchPlane::xy();
    let rect = AnalyticProfile::new(ProfileShape::Rect {
        width: 2.0,
        height: 1.0,
    })
    .local_wire(&plane)
    .unwrap();
    assert_eq!((rect.degree, rect.control_points.len()), (1, 5));
    assert!(
        rect.point_at(0.0)
            .distance_to(&Point3::new(-1.0, -0.5, 0.0))
            < 1e-12
    );

    let tri = AnalyticProfile::new(ProfileShape::Polygon {
        points: vec![[0.0, 0.0], [4.0, 0.0], [0.0, 3.0]],
    })
    .local_wire(&plane)
    .unwrap();
    assert!(tri.point_at(0.0).distance_to(&tri.point_at(1.0)) < 1e-12);

    for via in [[0.0, 18.0], [0.0, -18.0]] {
        let w = bulged(via).local_wire(&plane).unwrap();
        assert_eq!(w.degree, 2);
        assert!(
            w.point_at(0.0).distance_to(&w.point_at(1.0)) < 1e-12,
            "closed"
        );
        // The bulge passes through `via` and the arc is circular about
        // its circumcentre (0, (18² − 30²)/36).
        let cy = (18.0_f64.powi(2) - 30.0_f64.powi(2)) / 36.0 * via[1].signum();
        let c = Point3::new(0.0, cy, 0.0);
        let r = (30.0_f64.powi(2) + cy * cy).sqrt();
        let (t0, t1) = w.parameter_range();
        let arc_err = (0..=200)
            .map(|i| {
                // The arc occupies the second half of the parameter
                // range (line span + arc sub-spans, uniform knots).
                let t = t0 + (t1 - t0) * (0.5 + 0.5 * i as f64 / 200.0);
                (w.point_at(t).distance_to(&c) - r).abs()
            })
            .fold(0.0, f64::max);
        assert!(
            arc_err < 1e-9,
            "arc through {via:?}: radial error {arc_err:.3e}"
        );
        let hits_via = (0..=2000).any(|i| {
            w.point_at(t0 + (t1 - t0) * i as f64 / 2000.0)
                .distance_to(&Point3::new(via[0], via[1], 0.0))
                < 0.05
        });
        assert!(hits_via, "bulge must pass through via {via:?}");
    }
}

#[test]
fn uncompilable_descriptions_decline_with_the_named_defect() {
    let plane = SketchPlane::xy();
    let err = AnalyticProfile::new(ProfileShape::Polygon {
        points: vec![[0.0, 0.0], [1.0, 0.0]],
    })
    .local_wire(&plane)
    .unwrap_err();
    assert_eq!(err, ProfileError::PolygonTooShort { got: 2 });
    let err = AnalyticProfile::new(ProfileShape::Path {
        start: [0.0, 0.0],
        segments: vec![PathSegment::Line { to: [1.0, 0.0] }],
    })
    .local_wire(&plane)
    .unwrap_err();
    assert_eq!(err, ProfileError::PathTooShort);
    let err = AnalyticProfile::new(ProfileShape::Path {
        start: [0.0, 0.0],
        segments: vec![
            PathSegment::Line { to: [1.0, 0.0] },
            PathSegment::Arc {
                via: [2.0, 0.0],
                to: [3.0, 0.0],
            },
        ],
    })
    .local_wire(&plane)
    .unwrap_err();
    assert_eq!(err, ProfileError::CollinearArc);
    assert!(matches!(
        circle(0.0).local_wire(&plane).unwrap_err(),
        ProfileError::NonPositiveDimension {
            what: "circle radius",
            ..
        }
    ));
    assert!(matches!(
        AnalyticProfile::new(ProfileShape::Rect {
            width: 1.0,
            height: f64::NAN
        })
        .local_wire(&plane)
        .unwrap_err(),
        ProfileError::NonPositiveDimension {
            what: "rect height",
            ..
        }
    ));
}

/// The face path: exact `Arc` edges, a lamina that validates, and —
/// through the extrude operator's promotion (PK-377) — true
/// `Cylinder` walls. The circle case is the one the constraint bridge
/// used to split into two IDENTICAL half-arcs (validate red).
#[test]
fn circle_face_extrudes_to_cylinder_walls() {
    let placed = AnalyticProfile {
        center: [1.0, 0.0],
        ..circle(0.25)
    };
    let result = placed
        .to_brep(SketchPlane::xy(), Tolerance::default(), None)
        .unwrap();
    let mut brep = result.brep;
    brep.validate().expect("circle lamina validates");
    let arcs: Vec<_> = brep
        .edges
        .values()
        .filter_map(|e| match e.curve.as_ref().map(|c| &c.geometry) {
            Some(CurveGeometry::Arc {
                center,
                radius,
                start_angle,
                end_angle,
                ..
            }) => Some((*center, *radius, end_angle - start_angle)),
            _ => None,
        })
        .collect();
    assert_eq!(arcs.len(), 2, "two half-arc edges");
    for (c, r, sweep) in &arcs {
        assert!(c.distance_to(&Point3::new(1.0, 0.0, 0.0)) < 1e-12);
        assert!((r - 0.25).abs() < 1e-12);
        assert!(
            (sweep.abs() - PI).abs() < 1e-9,
            "each half sweeps π, got {sweep}"
        );
    }
    // The two halves cover DIFFERENT halves of the circle.
    let mids: Vec<Point3> = brep
        .edges
        .values()
        .filter_map(|e| e.curve.as_ref())
        .map(|c| c.geometry.point_at((c.t_start + c.t_end) / 2.0))
        .collect();
    assert!(
        mids[0].distance_to(&mids[1]) > 0.49,
        "halves must be distinct: {mids:?}"
    );

    let face = result.regions[0].face;
    let ext = brep
        .extrude(face, Vec3::new(0.0, 0.0, 1.0), 0.5, None)
        .unwrap();
    brep.validate().expect("extruded circle validates");
    assert!(brep.validate_all().is_empty(), "{:?}", brep.validate_all());
    for &f in &ext.side_faces {
        let SurfaceGeometry::Cylinder { origin, radius, .. } =
            &brep.faces[f].surface.as_ref().unwrap().geometry
        else {
            panic!("circle walls must be true Cylinders");
        };
        assert!((origin.x - 1.0).abs() < 1e-12 && origin.y.abs() < 1e-12);
        assert!((radius - 0.25).abs() < 1e-12);
    }
}

/// Arcs that cross the ±π branch of `atan2`, in both senses, come
/// through the bridge with the right sweep (the branch-cut defect
/// sent them the long way round).
#[test]
fn branch_crossing_and_clockwise_arcs_serve_exact_faces() {
    for via in [[-32.0, 0.0], [-28.0, 0.0]] {
        // A "C"/"D" outline: left side an arc from (-30, 20) through
        // `via` to (-30, -20) — centred near x = −30, crossing angle π.
        let p = AnalyticProfile::new(ProfileShape::Path {
            start: [-30.0, 20.0],
            segments: vec![
                PathSegment::Line { to: [0.0, 20.0] },
                PathSegment::Line { to: [0.0, -20.0] },
                PathSegment::Line { to: [-30.0, -20.0] },
                PathSegment::Arc {
                    via,
                    to: [-30.0, 20.0],
                },
            ],
        });
        let result = p
            .to_brep(SketchPlane::xy(), Tolerance::default(), None)
            .unwrap();
        let mut brep = result.brep;
        brep.validate().expect("outline validates");
        let (c, _) = brep
            .edges
            .values()
            .filter_map(|e| e.curve.as_ref())
            .find_map(|c| match c.geometry {
                CurveGeometry::Arc { .. } => Some((
                    c.geometry.point_at((c.t_start + c.t_end) / 2.0),
                    c.geometry.point_at(c.t_start),
                )),
                _ => None,
            })
            .expect("one arc edge");
        assert!(
            c.distance_to(&Point3::new(via[0], via[1], 0.0)) < 1e-9,
            "arc midpoint must be via {via:?}, got {c:?}"
        );
        let face = result.regions[0].face;
        let ext = brep
            .extrude(face, Vec3::new(0.0, 0.0, 1.0), 5.0, None)
            .unwrap();
        assert!(brep.validate_all().is_empty(), "{:?}", brep.validate_all());
        let cylinders = ext
            .side_faces
            .iter()
            .filter(|&&f| {
                matches!(
                    brep.faces[f].surface.as_ref().unwrap().geometry,
                    SurfaceGeometry::Cylinder { .. }
                )
            })
            .count();
        assert_eq!(cylinders, 1, "via {via:?}");
    }
}

#[test]
fn sketch_container_carries_zero_constraints() {
    let sk = bulged([0.0, 18.0]).to_sketch_elements().unwrap();
    assert_eq!(sk.constraint_count(), 0);
    assert_eq!(sk.non_construction_curves().len(), 2);
}
