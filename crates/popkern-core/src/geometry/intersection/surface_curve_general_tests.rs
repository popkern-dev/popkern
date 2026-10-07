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

// --- PK-513: straight-segment culling and grid memoization -------------

#[test]
fn a_diagonal_segment_misses_a_box_its_aabb_contains() {
    use super::segment_meets_box;
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(10.0, 10.0, 10.0);
    // The segment's AABB is [0,10]^3 and contains this box, but the
    // diagonal passes far from it.
    let (lo, hi) = (Point3::new(8.0, 0.0, 0.0), Point3::new(9.0, 1.0, 1.0));
    assert!(!segment_meets_box(&a, &b, &lo, &hi, 1e-6));
    // A box on the diagonal is met, and the margin admits a near miss.
    let (lo, hi) = (Point3::new(4.0, 4.0, 4.0), Point3::new(5.0, 5.0, 5.0));
    assert!(segment_meets_box(&a, &b, &lo, &hi, 0.0));
    let (lo, hi) = (Point3::new(5.1, 5.0, 5.0), Point3::new(6.0, 5.05, 6.0));
    assert!(segment_meets_box(&a, &b, &lo, &hi, 0.06));
    // Axis-parallel segment: a zero direction component is a slab
    // containment test.
    let c = Point3::new(0.5, 0.5, -3.0);
    let d = Point3::new(0.5, 0.5, 3.0);
    let (lo, hi) = (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0));
    assert!(segment_meets_box(&c, &d, &lo, &hi, 0.0));
    let c2 = Point3::new(1.5, 0.5, -3.0);
    let d2 = Point3::new(1.5, 0.5, 3.0);
    assert!(!segment_meets_box(&c2, &d2, &lo, &hi, 0.0));
}

#[test]
fn the_memoized_grid_equals_a_fresh_build() {
    use super::surface_curve_grid::{build_surface_grid, surface_grid_cached};
    let surface = nurbs_plane_z(0.25);
    let fresh = build_surface_grid(&surface, (0.0, 1.0), (0.0, 1.0), 8, true);
    for _ in 0..2 {
        let cached = surface_grid_cached(&surface, (0.0, 1.0), (0.0, 1.0), 8, true);
        assert_eq!(cached.len(), fresh.len());
        for (c, f) in cached.iter().zip(&fresh) {
            assert_eq!(c.u_range, f.u_range);
            assert_eq!(c.v_range, f.v_range);
            assert_eq!(c.min, f.min);
            assert_eq!(c.max, f.max);
        }
    }
}

/// A bicubic NURBS patch over the unit square at height `z` (two spans
/// per axis, so the grid refines on interior knots).
fn nurbs_plane_z(z: f64) -> SurfaceGeometry {
    use crate::geometry::nurbs_surface::NurbsSurface;
    let n = 5;
    let mut control_points = Vec::new();
    for i in 0..n {
        for j in 0..n {
            control_points.push(Point3::new(
                i as f64 / (n - 1) as f64,
                j as f64 / (n - 1) as f64,
                z + 0.05 * ((i + j) % 2) as f64,
            ));
        }
    }
    let knots = vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0];
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 3,
        degree_v: 3,
        knots_u: knots.clone(),
        knots_v: knots,
        control_points,
        num_u: n,
        num_v: n,
        weights: vec![1.0; n * n],
    })
}
