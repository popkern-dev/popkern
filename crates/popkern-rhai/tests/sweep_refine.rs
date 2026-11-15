//! A pipe along a cubic holds its circle between frames: the frames
//! double until the skinned wall meets the profile placed at every
//! midpoint within the linear tolerance. A tube of radius 20 along a
//! planar cubic, skinned through seventeen fixed frames, bulged between
//! them (its volume read πr²L 1.2e-6 short); every point of its wall now
//! lies 20 from the path within 1e-5.

use popkern_core::geometry::curve::closest_point_in_window;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{CurveEval, Point3};

#[test]
fn a_pipe_wall_along_a_cubic_stays_its_radius_from_the_path() {
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host
        .run(
            r#"let m = model();
let rail = m.profile(xy()).path(0.0, 0.0)
    .spline_to([[12.0, -16.0], [32.0, -16.0], [44.0, 0.0]], #{ degree: 3 }).close().sketch("rail").extrude(blind(1.0));
rail.pipe(edges().nurbs().at_height(0.0), 20.0, #{ label: "pipe" });
m"#,
            "sweep_refine",
        )
        .unwrap_or_else(|e| panic!("{e}"));
    let mut model = run.model.expect("model");
    let b = model.to_part().expect("the pipe").brep().clone();
    let path = NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(12.0, -16.0, 0.0),
            Point3::new(32.0, -16.0, 0.0),
            Point3::new(44.0, 0.0, 0.0),
        ],
        weights: vec![1.0; 4],
    };
    let mut worst = 0.0_f64;
    let mut walls = 0;
    for face in b.faces.values() {
        let Some(SurfaceGeometry::Nurbs(n)) = face.surface.as_ref().map(|s| &s.geometry) else {
            continue;
        };
        walls += 1;
        let ((u0, u1), (v0, v1)) = SurfaceGeometry::Nurbs(n.clone()).parameter_range();
        for i in 0..=40 {
            for j in 0..=12 {
                let p = n.point_at(
                    u0 + (u1 - u0) * f64::from(j) / 12.0,
                    v0 + (v1 - v0) * f64::from(i) / 40.0,
                );
                let t = closest_point_in_window(&path, &p, 0.0, 1.0);
                worst = worst.max((path.point_at(t).distance_to(&p) - 20.0).abs());
            }
        }
    }
    assert!(walls > 0, "the pipe has a B-spline wall");
    assert!(worst < 1e-5, "the wall strays {worst:.3e} from radius 20");
}
