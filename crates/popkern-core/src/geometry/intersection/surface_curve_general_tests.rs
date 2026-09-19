use super::*;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

// -----------------------------------------------------------------------
// Test 1: NURBS line piercing a plane (verify against analytic line-plane)
//
// Plane z=0, line along z-axis from (0,0,-1) to (0,0,1).
// Should find one intersection at t=0.5 (midpoint), point (0,0,0).
// -----------------------------------------------------------------------
#[test]
fn nurbs_line_piercing_plane() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };

    // Build a NURBS curve that represents a line from (0,0,-1) to (0,0,1).
    let nurbs_line = CurveGeometry::Nurbs(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, -1.0), Point3::new(0.0, 0.0, 1.0)],
        weights: vec![1.0, 1.0],
    });

    let tol = Tolerance::default();
    let results = intersect_surface_curve_general(
        &plane,
        (-10.0, 10.0),
        (-10.0, 10.0),
        &nurbs_line,
        (0.0, 1.0),
        &tol,
        false,
    );

    let points: Vec<_> = results
        .iter()
        .filter(|r| matches!(r, SurfaceCurveIntersection::Point { .. }))
        .collect();

    assert_eq!(
        points.len(),
        1,
        "Expected exactly 1 point intersection, got {}",
        points.len()
    );
    if let SurfaceCurveIntersection::Point { point, t, .. } = &points[0] {
        assert!(
            (point.z).abs() < 1e-5,
            "Point should be on z=0 plane, got z={}",
            point.z
        );
        assert!((*t - 0.5).abs() < 1e-5, "t should be ~0.5, got {}", t);
    }
}

// -----------------------------------------------------------------------
// Test 2: Circle (as NURBS) piercing a plane — expect 2 intersections.
//
// Plane z=0.  Circle in xz-plane centered at origin, radius 1.
// The circle crosses z=0 at (1,0,0) and (-1,0,0).
// -----------------------------------------------------------------------
#[test]
fn circle_piercing_plane_two_points() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };

    // Circle in the xz-plane: center=(0,0,0), normal=(0,1,0), radius=1.
    let circle = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 1.0, 0.0),
        radius: 1.0,
    };

    let tol = Tolerance::default();
    let results = intersect_surface_curve_general(
        &plane,
        (-10.0, 10.0),
        (-10.0, 10.0),
        &circle,
        (0.0, std::f64::consts::TAU),
        &tol,
        false,
    );

    let points: Vec<_> = results
        .iter()
        .filter(|r| matches!(r, SurfaceCurveIntersection::Point { .. }))
        .collect();

    assert!(
        points.len() >= 2,
        "Expected at least 2 point intersections for circle piercing plane, got {}",
        points.len()
    );

    // All intersection points should be on the z=0 plane.
    for p in &points {
        if let SurfaceCurveIntersection::Point { point, .. } = p {
            assert!(
                point.z.abs() < 1e-4,
                "Intersection point should be on z=0 plane, got z={}",
                point.z
            );
        }
    }
}

// -----------------------------------------------------------------------
// Test 3: Non-intersecting curve and surface (empty).
//
// Plane z=0.  Line from (0,0,1) to (0,0,2) — entirely above the plane.
// -----------------------------------------------------------------------
#[test]
fn non_intersecting_empty() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };

    let nurbs_line = CurveGeometry::Nurbs(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 1.0), Point3::new(0.0, 0.0, 2.0)],
        weights: vec![1.0, 1.0],
    });

    let tol = Tolerance::default();
    let results = intersect_surface_curve_general(
        &plane,
        (-10.0, 10.0),
        (-10.0, 10.0),
        &nurbs_line,
        (0.0, 1.0),
        &tol,
        false,
    );

    let points: Vec<_> = results
        .iter()
        .filter(|r| matches!(r, SurfaceCurveIntersection::Point { .. }))
        .collect();

    assert!(
        points.is_empty(),
        "Expected no intersections, got {}",
        points.len()
    );
}

// -----------------------------------------------------------------------
// Test 4: Invariant — curve.point_at(t) ~ surface.point_at(u, v)
//
// For every Point intersection returned, verify the invariant.
// -----------------------------------------------------------------------
#[test]
fn invariant_curve_surface_match() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };

    let nurbs_line = CurveGeometry::Nurbs(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, -1.0), Point3::new(0.0, 0.0, 1.0)],
        weights: vec![1.0, 1.0],
    });

    let tol = Tolerance::default();
    let results = intersect_surface_curve_general(
        &plane,
        (-10.0, 10.0),
        (-10.0, 10.0),
        &nurbs_line,
        (0.0, 1.0),
        &tol,
        false,
    );

    for r in &results {
        if let SurfaceCurveIntersection::Point { point, t, u, v } = r {
            let curve_pt = nurbs_line.point_at(*t);
            let surface_pt = plane.point_at(*u, *v);
            assert!(
                point.distance_to(&curve_pt) < 1e-4,
                "curve.point_at(t) should match intersection point"
            );
            assert!(
                point.distance_to(&surface_pt) < 1e-4,
                "surface.point_at(u,v) should match intersection point"
            );
        }
    }
}

// -----------------------------------------------------------------------
// padded_cell AABB containment: the box must contain its patch even
// when the radial extremum lies BETWEEN the centre and a corner
// (PK-460 tl=8: a 74°-azimuth helical-band cell phased by the classify
// window pad boxed x ≤ 2.95 while the patch bulged to 3.0625; the
// exit crossing was AABB-culled before Newton ran and an inside-tool
// sample classified Outside). A cylinder cell reproduces the blind
// spot exactly: azimuth span [-15°, +58°] puts the x-extremum at
// azimuth 0 — off-centre and outside the corner box.
// -----------------------------------------------------------------------
#[test]
fn padded_cell_contains_off_center_bulge() {
    let cylinder = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 3.0625,
    };
    let deg = std::f64::consts::PI / 180.0;
    let u_range = (-15.0 * deg, 58.0 * deg);
    let v_range = (0.0, 0.4);
    let cell = super::surface_curve_grid::padded_cell(&cylinder, u_range, v_range);
    for i in 0..=20 {
        for j in 0..=20 {
            let u = u_range.0 + (u_range.1 - u_range.0) * i as f64 / 20.0;
            let v = v_range.0 + (v_range.1 - v_range.0) * j as f64 / 20.0;
            let p = cylinder.point_at(u, v);
            assert!(
                p.x >= cell.min.x
                    && p.x <= cell.max.x
                    && p.y >= cell.min.y
                    && p.y <= cell.max.y
                    && p.z >= cell.min.z
                    && p.z <= cell.max.z,
                "patch point ({},{},{}) outside cell box {:?}..{:?} at uv ({u},{v})",
                p.x,
                p.y,
                p.z,
                (cell.min.x, cell.min.y, cell.min.z),
                (cell.max.x, cell.max.y, cell.max.z),
            );
        }
    }
}
