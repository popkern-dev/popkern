//! A wedge a millionth steep: a cube cut by a plane meeting its top at
//! 1e-6 rad (PK-1131). The wedge is 5e-6 thick at its far end and
//! tapers to an edge where the plane meets the top; its side faces are
//! needles 5 long and at most 5e-6 tall, its front face a strip 5e-6
//! tall. At the tolerance's resolution the wedge is a valid body: its
//! vertices lie more than the tolerance apart, its thin edge within it
//! of both faces. It declined three ways: the plane read as the
//! cube's top (planes within the angular tolerance through a common
//! line, though 5e-6 apart across the face), the front vertical edges'
//! crossings 5e-6 below the corners read as the corners (an end test in
//! parameter units: 5e-6 on an edge 10 long is half a millionth of its
//! parameter), and the graph's on-segment test ignored the same
//! crossings for the same reason.

use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn block(lo: (f64, f64, f64), hi: (f64, f64, f64), tol: Tolerance) -> BRep {
    let rect = |z: f64| {
        vec![
            Point3::new(lo.0, lo.1, z),
            Point3::new(hi.0, lo.1, z),
            Point3::new(hi.0, hi.1, z),
            Point3::new(lo.0, hi.1, z),
        ]
    };
    BRep::polygon_frustum(&rect(lo.2), &rect(hi.2), tol).expect("block")
}

#[test]
fn a_wedge_a_millionth_steep_serves_its_volume() {
    let tol = Tolerance::default();
    let turned =
        Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.7, 0.3), 0.6)
            .expect("axis");
    let tilt = Transform3::rotation_around_axis(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        1e-6f64.atan(),
    )
    .expect("axis");
    // Thickness 1e-6 |y| over y in [-5, 0], 10 wide: 10 * 1e-6 * 12.5.
    let want = 1.25e-4;
    let mut problems = Vec::new();
    for (name, place) in [("identity", Transform3::identity()), ("turned", turned)] {
        let cube = block((-5.0, -5.0, -10.0), (5.0, 5.0, 0.0), tol)
            .transform(&place)
            .expect("place");
        let above = block((-20.0, -20.0, 0.0), (20.0, 20.0, 20.0), tol)
            .transform(&tilt)
            .and_then(|b| b.transform(&place))
            .expect("place");
        match popkern_bool::boolean(
            popkern_bool::BoolOperation::Intersection,
            cube,
            above,
            tol,
            None,
        ) {
            Ok(r) => {
                let v = popkern_measure::body_volume(&r, &tol).expect("volume");
                if let Some(first) = r.validate_all().first() {
                    problems.push(format!("{name}: {first}"));
                } else if (v - want).abs() > 1e-9 + 1e-6 * want {
                    problems.push(format!("{name}: volume {v} for {want}"));
                } else if r.faces.len() != 5 {
                    problems.push(format!("{name}: {} faces", r.faces.len()));
                }
            }
            Err(e) => problems.push(format!("{name}: declined {e}")),
        }
    }
    assert!(problems.is_empty(), "{problems:#?}");
}
