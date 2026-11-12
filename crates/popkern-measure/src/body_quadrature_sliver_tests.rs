use super::*;
use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::Point3;
use popkern_core::geometry::Transform3;
use popkern_core::geometry::pcurve::PcurveGeometry;

/// PK-938: `locate` refines by Newton on `|pc(t) − q|²`, a step of
/// `(pc − q)·pc′ / |pc′|²`. A NURBS trim's `tangent_at` was its unit
/// tangent, so on a pcurve moving 6 units of u per unit t each step
/// overshot six-fold and the refinement oscillated off the sample it
/// started from. It now converges to the foot.
#[test]
fn locate_converges_on_a_fast_nurbs_trim() {
    let pcurve = Pcurve {
        geometry: PcurveGeometry::Nurbs(NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(3.0, 2.0, 0.0),
                Point3::new(6.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        }),
        t_start: 0.0,
        t_end: 1.0,
    };
    let samples = (0..=LOCATE_SAMPLES)
        .map(|k| {
            let t = k as f64 / LOCATE_SAMPLES as f64;
            let p = pcurve.point_at(t);
            (t, (p.x, p.y))
        })
        .collect();
    let trim = Trim {
        pcurve: &pcurve,
        edge: None,
        lo: 0.0,
        hi: 1.0,
        closed: false,
        samples,
    };
    for foot in [0.1037, 0.304, 0.77] {
        let (p, d) = (pcurve.point_at(foot), pcurve.tangent_at(foot));
        let n = (-d.y, d.x);
        let len = (n.0 * n.0 + n.1 * n.1).sqrt();
        let q = (p.x + 0.01 * n.0 / len, p.y + 0.01 * n.1 / len);
        let (t, distance) = locate(&trim, q, (None, None));
        assert!((t - foot).abs() < 1e-10, "foot {foot}: located {t}");
        assert!((distance - 0.01).abs() < 1e-10, "foot {foot}: {distance}");
    }
}

/// Viviani's common of a ball of radius 10 (poles on z, as built, then
/// turned about z) and a cylinder of radius 5 on x = 5, the ball's poles
/// on the caps' rims, has the
/// volume (2/9)(3π − 4)R³: a boundary chord from a pole vertex stands
/// under its other end (the pole row is the pole), the trim reaching
/// the pole at its own longitude, and its sliver is found on the trim
/// in space and fanned along the pole row. Located in (u, v) the chord's
/// pole end lay a longitude step off the trim, no trim carried it, and
/// the triangle between was lost: 5.7e-7 of the volume (PK-1169).
#[test]
fn a_chord_from_a_pole_vertex_keeps_its_sliver() {
    let want = 2.0 / 9.0 * (3.0 * std::f64::consts::PI - 4.0) * 1000.0;
    let tol = popkern_core::tolerance::Tolerance::default();
    for turn in [0.0, 0.7] {
        let ball = BRep::sphere(10.0, tol)
            .and_then(|b| b.transform(&Transform3::rotation_z(turn)))
            .expect("ball");
        let rod = BRep::cylinder(5.0, 24.0, tol)
            .and_then(|b| b.transform(&Transform3::rotation_z(std::f64::consts::PI)))
            .and_then(|b| b.transform(&Transform3::translation(5.0, 0.0, -12.0)))
            .expect("rod");
        let body = popkern_bool::boolean(
            popkern_bool::BoolOperation::Intersection,
            ball,
            rod,
            tol,
            None,
        )
        .expect("Viviani's common");
        let v = crate::body_volume(&body, &tol).expect("a volume");
        assert!(
            (v - want).abs() < 1e-9 * want,
            "turn {turn}: {v} for {want}"
        );
    }
}
