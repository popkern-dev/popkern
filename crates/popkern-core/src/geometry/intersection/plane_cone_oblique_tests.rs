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

/// The section's pcurve on the cone is the curve's own, at every
/// parameter where the section passes the cone's window: `cone(uv(t))`
/// within the linear tolerance of `curve(t)`, for an oblique ellipse, a
/// perpendicular circle and a near-parabolic ellipse at scale 1000. The
/// 3D control points projected to the chart with the curve's weights
/// strayed 1.3e-2 (a box face's ellipse on a frustum, corpus PF-X3 379):
/// a clip crossing read on the pcurve landed 7.6e-3 from the junction its
/// record had to reach, the record was dropped, and the union's loop
/// broke. The shared degree-1 chart fit, tried next, reached its sample
/// cap on the near-parabolic one (PF-X 411) and returned a pcurve 3.4
/// off; the Hermite fit, doubling uniformly to its cap, could not reach
/// the far end of it.
#[test]
fn the_cone_pcurve_is_co_parameterized_with_the_section() {
    use crate::geometry::curve::CurveEval;
    use crate::geometry::surface::SurfaceEval;
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(1.3070545268287006, -0.11351508012919986, 5.652775099718006),
        axis: Vec3::new(0.2845865245180198, -0.3432775210217751, -0.8950815904842975),
        half_angle: 0.30587887140485215,
    };
    let oblique = SurfaceGeometry::Plane {
        origin: Point3::new(1.6808273084955485, -0.5042598616498635, 2.2583413476934675),
        normal: Vec3::new(-0.0860508000801483, 0.08083691267671221, 0.9930058677341558),
        u_axis: Vec3::new(
            0.9784242837598383,
            -0.18105378305915348,
            0.09952612013444798,
        ),
        v_axis: Vec3::new(
            0.18780574898797478,
            0.9803808624227346,
            -0.05975935018442719,
        ),
    };
    let axis = Vec3::new(0.2845865245180198, -0.3432775210217751, -0.8950815904842975);
    let (e1, e2) = crate::geometry::analytic_surfaces::orthonormal_basis(axis);
    let square = SurfaceGeometry::Plane {
        origin: Point3::new(1.3070545268287006, -0.11351508012919986, 5.652775099718006)
            + axis * 3.0,
        normal: axis,
        u_axis: e1,
        v_axis: e2,
    };
    let frustum = SurfaceGeometry::Cone {
        apex: Point3::new(1361.850193720257, -108.59607800301, 129.04687950558895),
        axis: Vec3::new(
            -0.7185688548503939,
            0.6900844510627112,
            -0.08626848347148491,
        ),
        half_angle: 0.30805278102377637,
    };
    let near_parabolic = SurfaceGeometry::Plane {
        origin: Point3::new(-3142.871087794154, 3292.286877229266, -40.88422173924586),
        normal: Vec3::new(-0.8184363057636568, -0.3447470139290311, 0.4596863167367009),
        u_axis: Vec3::new(0.270826265541727, -0.9370232870111803, -0.22054589883170625),
        v_axis: Vec3::new(0.5067693235592112, -0.05600764220871479, 0.8602604237750102),
    };
    let tol = Tolerance::default();
    let cases = [
        (&oblique, &cone, (-10.0, 10.0), 20.0),
        (&square, &cone, (-10.0, 10.0), 20.0),
        (&near_parabolic, &frustum, (-1500.0, 9000.0), 20000.0),
    ];
    for (plane, cone, heights, reach) in cases {
        let result = super::super::try_plane_cone(
            plane,
            (-reach, reach),
            (-reach, reach),
            cone,
            (0.0, TAU),
            heights,
            &tol,
        )
        .expect("an exact section");
        let SurfaceSurfaceIntersection::Curve(c) = &result[0] else {
            panic!("a curve");
        };
        let (t0, t1) = c.curve_3d.parameter_range();
        let height = |q: Point3| {
            let SurfaceGeometry::Cone { apex, axis, .. } = cone else {
                unreachable!()
            };
            (q - *apex).dot(axis)
        };
        let worst = (0..=20000)
            .map(|k| t0 + (t1 - t0) * f64::from(k) / 20000.0)
            .filter(|&t| {
                let h = height(c.curve_3d.point_at(t));
                h >= heights.0 && h <= heights.1
            })
            .map(|t| {
                let uv = c.curve_uv2.point_at(t);
                cone.point_at(uv.x, uv.y)
                    .distance_to(&c.curve_3d.point_at(t))
            })
            .fold(0.0_f64, f64::max);
        assert!(
            worst <= tol.linear(),
            "the cone pcurve strays {worst:.3e} from the section"
        );
    }
}
