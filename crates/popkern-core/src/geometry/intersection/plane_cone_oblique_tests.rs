//! Tests for the oblique plane-cone section.

use super::*;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use std::f64::consts::TAU;

/// A near-parabolic ellipse turns most of a turn round its far end within
/// one quarter of cone angle: its section must still be one simple loop
/// on both surfaces, its tangent turning one whole turn (a box face cut
/// from a cone: the quarter there came out as the complementary arc, the
/// loop retraced its next quarter and turned past a whole turn, and the
/// box face's edge met the record twice at one point; the generated
/// corpus's box x cone crossings). Failing state: 12.49 rad of absolute
/// turning, two turns.
#[test]
fn a_near_parabolic_section_is_one_simple_loop() {
    use crate::geometry::curve::CurveEval;
    use crate::geometry::surface::SurfaceEval;
    use popkern_math::PortableMath;
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.46530809839046894, -2.5194065453285526, 0.586118127026159),
        normal: Vec3::new(-0.4850328944979293, 0.23509721596986988, -0.842301840374445),
        u_axis: Vec3::new(
            0.45336492525703553,
            0.8912399815474512,
            -0.012310151826100527,
        ),
        v_axis: Vec3::new(
            0.7477989942502219,
            -0.38784093947714626,
            -0.5388655396885585,
        ),
    };
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(1.9227099260235259, -4.2315687546604694, 0.02043632002898843),
        axis: Vec3::new(0.321645266816798, 0.684646403347006, -0.6540669879441707),
        half_angle: 0.4636476090008061,
    };
    let tol = Tolerance::default();
    let result = super::super::try_plane_cone(
        &plane,
        (-1.0, 2.9),
        (-1.3, 3.8),
        &cone,
        (0.0, TAU),
        (-1.1, 3.3),
        &tol,
    )
    .expect("an exact section");
    let SurfaceSurfaceIntersection::Curve(c) = &result[0] else {
        panic!("a curve");
    };
    let (t0, t1) = c.curve_uv1.parameter_range();
    let n = 4000;
    // The sampled polygon's absolute turning: one whole turn for a simple
    // convex loop; a retrace adds a half turn at each reversal.
    let mut pts = Vec::with_capacity(n + 1);
    for k in 0..=n {
        let t = t0 + (t1 - t0) * k as f64 / n as f64;
        let q = c.curve_3d.point_at(t);
        for s in [&plane, &cone] {
            let (u, v) = s.closest_point(&q).expect("projects");
            assert!(
                s.point_at(u, v).distance_to(&q) < 1e-9,
                "off a surface at {t}"
            );
        }
        let p = c.curve_uv1.point_at(t);
        pts.push((p.x, p.y));
    }
    let mut total = 0.0;
    for k in 0..n {
        let (a, b, c) = (pts[k], pts[k + 1], pts[(k + 2) % n]);
        let (d0, d1) = ((b.0 - a.0, b.1 - a.1), (c.0 - b.0, c.1 - b.1));
        total += (d0.0 * d1.1 - d0.1 * d1.0)
            .patan2(d0.0 * d1.0 + d0.1 * d1.1)
            .abs();
    }
    assert!((total - TAU).abs() < 1e-2, "tangent turns {total}");
}
