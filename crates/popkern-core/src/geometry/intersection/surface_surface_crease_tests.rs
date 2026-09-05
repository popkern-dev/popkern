//! SSI marching across C0 creases of a faceted closed NURBS wall
//! (PK-442). The march wraps the geometrically-closed chart through
//! every crease, and the segment honesty gate keeps seam-corner wrap
//! jitter from emitting phantom records off the surface.

use crate::geometry::intersection::{SurfaceSurfaceIntersection, intersect_surfaces};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// A square prism as ONE degree-1 closed NURBS wall: u runs around
/// the square (a C0 crease at every corner, one doubling as the chart
/// seam), v runs along the height.
fn square_prism() -> SurfaceGeometry {
    use crate::geometry::nurbs_surface::NurbsSurface;
    let ring = [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0), (1.0, 0.0)];
    // Layout: index = i_u * num_v + i_v (u-major).
    let mut cps = Vec::new();
    for (x, y) in ring {
        for z in [0.0, 2.0] {
            cps.push(Point3::new(x, y, z));
        }
    }
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        num_u: 5,
        num_v: 2,
        control_points: cps,
        weights: vec![1.0; 10],
    })
}

/// A tilted plane through the prism: the intersection is one closed
/// quadrilateral-ish loop crossing all four creases (including the
/// seam corner). Every emitted intersection must be that closed loop
/// — traced whole through the creases — and every emission must lie
/// ON both surfaces (no extrapolated-facet phantoms, no open
/// half-traces).
#[test]
fn march_traces_a_closed_loop_through_c0_creases() {
    let prism = square_prism();
    let u_axis = Vec3::new(1.0, 0.0, -0.2).normalized().expect("axis");
    let v_axis = Vec3::new(0.0, 1.0, -0.1).normalized().expect("axis");
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 1.0),
        normal: u_axis.cross(&v_axis).normalized().expect("normal"),
        u_axis,
        v_axis,
    };
    let tol = Tolerance::default();
    let result = intersect_surfaces(
        &prism,
        (0.0, 1.0),
        (0.0, 1.0),
        &plane,
        (-3.0, 3.0),
        (-3.0, 3.0),
        &tol,
    )
    .expect("intersect");
    assert!(
        !result.intersections.is_empty(),
        "the crease-crossing loop must be traced"
    );
    // The loop is emitted as crease-bounded PIECES (the same contract
    // as seam splitting): each piece is smooth, fits the wall to
    // tolerance, and terminates exactly on a crease or seam iso line.
    // Before the crease split, one smooth fit rounded the corners by
    // 1.1e-2; before the segment honesty gate, seam-corner jitter
    // emitted phantom pieces on the EXTRAPOLATED facet.
    let crease_us = [0.0, 0.25, 0.5, 0.75, 1.0];
    let mut pieces = 0usize;
    for item in &result.intersections {
        let SurfaceSurfaceIntersection::Curve(c) = item else {
            panic!("unexpected coincident region on a transversal cut");
        };
        pieces += 1;
        use crate::geometry::curve::CurveEval;
        let (t0, t1) = c.curve_3d.parameter_range();
        for k in 0..=16 {
            let t = t0 + (t1 - t0) * k as f64 / 16.0;
            let p = c.curve_3d.point_at(t);
            let q = c.curve_uv1.point_at(t);
            let on1 = prism.point_at(q.x, q.y).distance_to(&p);
            assert!(
                on1 < 1.0e-6,
                "emitted record leaves the wall: station {k} deviates {on1:.3e}"
            );
        }
        for t in [t0, t1] {
            let q = c.curve_uv1.point_at(t);
            let nearest = crease_us
                .iter()
                .map(|&u| (q.x - u).abs())
                .fold(f64::INFINITY, f64::min);
            assert!(
                nearest < 1.0e-6,
                "piece endpoint u {:.6} is not on a crease/seam iso line",
                q.x
            );
        }
    }
    assert!(
        pieces >= 4,
        "the loop crosses four creases; expected at least four pieces, got {pieces}"
    );
}
