use super::*;
use crate::geometry::curve::CurveEval;
use crate::geometry::edge_triple::{BoundarySide, SurfaceSide};
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::{Point3, Vec3};

fn z_axis_cyl(radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius,
    }
}

#[test]
fn identical_cylinders_short_circuit_to_coincident() {
    let cyl = z_axis_cyl(1.0);
    let anchors = EdgeAnchors::Open {
        start: Point3::new(1.0, 0.0, 0.0),
        end: Point3::new(1.0, 0.0, 1.0),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let result = compute_edge_geometry(&cyl, &cyl, anchors, None, &tol, &opts);
    match result {
        Err(EdgeTripleError::CoincidentSurfaces { reason }) => {
            assert!(!reason.is_empty(), "reason must be populated");
        }
        other => panic!("expected CoincidentSurfaces, got {other:?}"),
    }
}

#[test]
fn nurbs_surface_boundary_hint_constructs_owned_boundary_edge() {
    let patch = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    });
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let triple = compute_edge_geometry(
        &patch,
        &plane,
        EdgeAnchors::Open {
            start: Point3::new(0.0, 0.0, 0.0),
            end: Point3::new(0.0, 1.0, 0.0),
        },
        Some(EdgeHint::SurfaceBoundary {
            side: SurfaceSide::A,
            boundary: BoundarySide::UMin,
        }),
        &Tolerance::default(),
        &EdgeTripleOptions::default(),
    )
    .expect("NURBS boundary hint constructs edge triple");

    assert!(matches!(
        triple.provenance,
        crate::geometry::edge_triple::types::EdgeTripleProvenance::SurfaceBoundary { .. }
    ));
    let quality = triple
        .verify(&patch, &plane, &Tolerance::default(), 8)
        .expect("NURBS boundary triple verifies");
    assert!(quality.max_positional_residual <= 1e-6);
}

#[test]
fn nurbs_boundary_hint_projects_onto_nurbs_target_surface() {
    let patch = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    });
    let target = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(-1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    });

    let triple = compute_edge_geometry(
        &patch,
        &target,
        EdgeAnchors::Open {
            start: Point3::new(0.0, 0.0, 0.0),
            end: Point3::new(0.0, 1.0, 0.0),
        },
        Some(EdgeHint::SurfaceBoundary {
            side: SurfaceSide::A,
            boundary: BoundarySide::UMin,
        }),
        &Tolerance::default(),
        &EdgeTripleOptions::default(),
    )
    .expect("NURBS boundary should project onto NURBS target");

    assert!(matches!(
        triple.provenance,
        crate::geometry::edge_triple::types::EdgeTripleProvenance::SurfaceBoundary { .. }
    ));
    let quality = triple
        .verify(&patch, &target, &Tolerance::default(), 8)
        .expect("NURBS-to-NURBS boundary triple verifies");
    assert!(quality.max_positional_residual <= 1e-6);
}

#[test]
fn coincident_unit_spheres_short_circuit() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let anchors = EdgeAnchors::Open {
        start: Point3::new(1.0, 0.0, 0.0),
        end: Point3::new(0.0, 1.0, 0.0),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let result = compute_edge_geometry(&sphere, &sphere, anchors, None, &tol, &opts);
    assert!(
        matches!(result, Err(EdgeTripleError::CoincidentSurfaces { .. })),
        "expected CoincidentSurfaces, got {result:?}"
    );
}

#[test]
fn cyl_inscribed_sphere_returns_unit_equator_circle() {
    let cyl = z_axis_cyl(1.0);
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let anchors = EdgeAnchors::Open {
        start: Point3::new(1.0, 0.0, 0.0),
        end: Point3::new(0.0, 1.0, 0.0),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let triple = compute_edge_geometry(&cyl, &sphere, anchors, None, &tol, &opts)
        .expect("tangent path produces a triple");
    let p_start = triple.curve_3d.point_at(triple.t_start);
    let p_end = triple.curve_3d.point_at(triple.t_end);
    assert!(
        p_start.distance_to(&Point3::new(1.0, 0.0, 0.0)) <= tol.linear(),
        "start anchor mismatch: {:?}",
        p_start
    );
    assert!(
        p_end.distance_to(&Point3::new(0.0, 1.0, 0.0)) <= tol.linear(),
        "end anchor mismatch: {:?}",
        p_end
    );
    let q = triple
        .verify(&cyl, &sphere, &tol, 32)
        .expect("verify on tangent intersection");
    assert!(
        q.max_positional_residual <= tol.linear(),
        "G2 residual {}",
        q.max_positional_residual
    );
    match &triple.provenance {
        crate::geometry::edge_triple::types::EdgeTripleProvenance::AnalyticIntersection {
            cardinality: Some(c),
            ..
        } => {
            assert!(c.tangent, "expected tangent: true on cyl-inscribed-sphere");
            assert_eq!(c.total, 1, "tangent intersection has a single branch");
        }
        other => panic!("expected AnalyticIntersection w/ cardinality, got {other:?}"),
    }
}

#[test]
fn coaxial_cylinders_with_different_radii_do_not_short_circuit() {
    let cyl_a = z_axis_cyl(1.0);
    let cyl_b = z_axis_cyl(2.0);
    let anchors = EdgeAnchors::Open {
        start: Point3::new(1.0, 0.0, 0.0),
        end: Point3::new(1.0, 0.0, 1.0),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let result = compute_edge_geometry(&cyl_a, &cyl_b, anchors, None, &tol, &opts);
    assert!(
        !matches!(result, Err(EdgeTripleError::CoincidentSurfaces { .. })),
        "must not classify distinct radii as coincident: {result:?}"
    );
}

#[test]
fn all_branches_returns_single_element_vec_on_tangent() {
    let cyl = z_axis_cyl(1.0);
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let anchors = EdgeAnchors::Open {
        start: Point3::new(1.0, 0.0, 0.0),
        end: Point3::new(0.0, 1.0, 0.0),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let branches = compute_edge_geometry_all_branches(&cyl, &sphere, anchors, &tol, &opts)
        .expect("all-branches returns Ok on tangent intersection");
    assert_eq!(branches.len(), 1, "single-branch path returns one triple");
}

#[test]
fn all_branches_propagates_coincident_surfaces_error() {
    let cyl = z_axis_cyl(1.0);
    let anchors = EdgeAnchors::Open {
        start: Point3::new(1.0, 0.0, 0.0),
        end: Point3::new(1.0, 0.0, 1.0),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let result = compute_edge_geometry_all_branches(&cyl, &cyl, anchors, &tol, &opts);
    assert!(
        matches!(result, Err(EdgeTripleError::CoincidentSurfaces { .. })),
        "expected CoincidentSurfaces, got {result:?}"
    );
}

#[test]
fn all_branches_returns_two_branches_for_steinmetz_cyl_cyl() {
    let cyl_z = z_axis_cyl(1.0);
    let cyl_x = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(1.0, 0.0, 0.0),
        radius: 1.0,
    };
    let anchors = EdgeAnchors::Open {
        start: Point3::new(1.0, 0.0, -1.0),
        end: Point3::new(0.0, 1.0, 0.0),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let branches = compute_edge_geometry_all_branches(&cyl_z, &cyl_x, anchors, &tol, &opts)
        .expect("Steinmetz produces branches");
    assert_eq!(
        branches.len(),
        2,
        "Steinmetz has 2 ellipse branches (one per plane fibre)"
    );
    for (i, t) in branches.iter().enumerate() {
        let q = t.verify(&cyl_z, &cyl_x, &tol, 32).unwrap_or_else(|e| {
            panic!("branch {} failed verify: {:?}", i, e);
        });
        assert!(
            q.max_positional_residual <= tol.linear(),
            "branch {} G2 residual {} > tol {}",
            i,
            q.max_positional_residual,
            tol.linear()
        );
    }
    for (i, t) in branches.iter().enumerate() {
        let crate::geometry::edge_triple::types::EdgeTripleProvenance::AnalyticIntersection {
            cardinality: Some(c),
            ..
        } = &t.provenance
        else {
            panic!("branch {} missing AnalyticIntersection cardinality", i);
        };
        assert_eq!(c.total, 2, "branch {} cardinality.total", i);
        assert!(!c.tangent, "Steinmetz is transverse, not tangent");
    }
}

#[test]
fn all_branches_propagates_dispatch_outcome() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 5.0,
        minor_radius: 1.0,
    };
    let anchors = EdgeAnchors::Open {
        start: Point3::new(6.0, 0.0, 0.0),
        end: Point3::new(0.0, 6.0, 0.0),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let _ = compute_edge_geometry_all_branches(&torus, &torus, anchors, &tol, &opts);
}
