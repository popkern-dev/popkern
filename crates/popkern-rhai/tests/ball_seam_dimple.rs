//! PK-1194: a ball dimpled across its seam keeps the seam's two sides a
//! period apart in its chart. Split by the dimple, both came back on
//! u = 0, the face's loop enclosed two slivers beside the seam, the
//! classify ray rejected every hit on the ball, and the next dimple kept
//! its outer cap (a golf ball read 4174.7 for 4124.0).

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;

#[test]
fn a_ball_dimpled_across_its_seam_reads_inside_and_its_volume() {
    let source = "let m = model();\nlet phi = (1.0 + sqrt(5.0)) / 2.0;\nlet ball = m.sphere(10.0);\n\
        for d in [[0.0, 1.0, phi], [phi, 0.0, 1.0], [0.0, 1.0, -phi]] {\n\
            let l = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);\n\
            let k = 10.6 / l;\n\
            ball = ball.subtract(m.tool(|m| m.sphere(1.2).translate(d[0] * k, d[1] * k, d[2] * k)));\n\
        }\nm";
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host
        .run(source, "dimples")
        .unwrap_or_else(|e| panic!("{e}"));
    let body = run.model.expect("model").to_parts_all().expect("parts")[0]
        .brep()
        .clone();
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
    let boxes: Vec<_> = body
        .faces
        .keys()
        .map(|f| {
            let far = 100.0;
            (
                f,
                Aabb {
                    min: Point3::new(-far, -far, -far),
                    max: Point3::new(far, far, far),
                },
            )
        })
        .collect();
    let at = popkern_bool::point_in_solid_final(
        &Point3::new(0.0, 0.0, 0.0),
        &body,
        &boxes,
        &body.tolerance,
    );
    assert!(
        matches!(at, Ok(popkern_bool::error::PointClassification::Inside)),
        "{at:?}"
    );
    // The ball less three lenses with spheres of radius 1.2 at 10.6.
    let (r, s, d) = (10.0_f64, 1.2_f64, 10.6_f64);
    let lens = std::f64::consts::PI
        * (r + s - d).powi(2)
        * (d * d + 2.0 * d * s - 3.0 * s * s + 2.0 * d * r + 6.0 * s * r - 3.0 * r * r)
        / (12.0 * d);
    let want = 4.0 / 3.0 * std::f64::consts::PI * r.powi(3) - 3.0 * lens;
    let v = popkern_measure::body_volume(&body, &Tolerance::default()).unwrap();
    assert!((v - want).abs() < 1e-6 * want, "{v} want {want}");
}
