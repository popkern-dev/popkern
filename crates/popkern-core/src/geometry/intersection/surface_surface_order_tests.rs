//! PK-252 solver-level regressions: operand-order equivalence and
//! single-emission of the intersection curve set.
//!
//! The marching stage consumes seeds by distance to the traced
//! polyline SEGMENTS. Before that fix, seeds landing mid-segment
//! (sample spacing up to 1000× linear tolerance vs a 10× consumption
//! radius) re-marched the identical cylinder×NURBS arc 18–20 times,
//! with the count varying by operand order.

use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::surface_surface::intersect_surfaces;
use crate::geometry::intersection::{SurfaceSurfaceIntersection, SurfaceSurfaceIntersectionResult};
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

fn quarter_cylinder() -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 0.5,
    }
}

/// Degree-1 bilinear quad at z = 0.5 crossing the cylinder wall.
fn crossing_quad() -> SurfaceGeometry {
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.1, 0.1, 0.5),
            Point3::new(0.1, 1.0, 0.5),
            Point3::new(1.0, 0.1, 0.5),
            Point3::new(1.0, 1.0, 0.5),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    })
}

fn curve_endpoints(result: &SurfaceSurfaceIntersectionResult) -> Vec<(Point3, Point3)> {
    result
        .intersections
        .iter()
        .filter_map(|ix| match ix {
            SurfaceSurfaceIntersection::Curve(c) => {
                let (t0, t1) = c.curve_3d.parameter_range();
                Some((c.curve_3d.point_at(t0), c.curve_3d.point_at(t1)))
            }
            SurfaceSurfaceIntersection::CoincidentRegion { .. } => None,
        })
        .collect()
}

#[test]
fn cylinder_nurbs_arc_emitted_once_in_both_operand_orders() {
    let tol = Tolerance::default();
    let cyl = quarter_cylinder();
    let quad = crossing_quad();
    let cyl_u = (-0.5, std::f64::consts::FRAC_PI_2 + 0.5);
    let cyl_v = (-0.5, 1.5);
    let quad_uv = (0.0, 1.0);

    let analytic_first =
        intersect_surfaces(&cyl, cyl_u, cyl_v, &quad, quad_uv, quad_uv, &tol).unwrap();
    let free_form_first =
        intersect_surfaces(&quad, quad_uv, quad_uv, &cyl, cyl_u, cyl_v, &tol).unwrap();

    let af = curve_endpoints(&analytic_first);
    let ff = curve_endpoints(&free_form_first);
    assert_eq!(
        af.len(),
        1,
        "analytic-first must emit the arc exactly once: {} curves",
        af.len()
    );
    assert_eq!(
        ff.len(),
        1,
        "free-form-first must emit the arc exactly once: {} curves",
        ff.len()
    );

    // The two orders trace the same arc (endpoint sets match within a
    // seed-placement band; the trace direction may differ).
    let (a0, a1) = af[0];
    let (b0, b1) = ff[0];
    let band = 5.0e-3;
    let direct = a0.distance_to(&b0) < band && a1.distance_to(&b1) < band;
    let reversed = a0.distance_to(&b1) < band && a1.distance_to(&b0) < band;
    assert!(
        direct || reversed,
        "orders must trace the same arc: af=({a0:?},{a1:?}) ff=({b0:?},{b1:?})"
    );

    // And it is the true intersection: on the cylinder and in the quad plane.
    for (p, q) in [af[0], ff[0]] {
        for pt in [p, q] {
            let radial = (pt.x * pt.x + pt.y * pt.y).sqrt();
            assert!(
                (radial - 0.5).abs() < 1.0e-3 && (pt.z - 0.5).abs() < 1.0e-6,
                "arc endpoint off the true intersection: {pt:?}"
            );
        }
    }
}
