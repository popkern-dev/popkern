//! Drafted faces that meet: the tangency foundations' draft consumer
//! (design §8.4, F1).
//!
//! A rib 10 long, 1 thick and 10 high, both long sides drafted about
//! the base. Drafted apart or together short of the top, the rib
//! serves its closed-form volume; drafted until the sides meet before
//! the top, or exactly at its top edge (where the top face shrinks to
//! a line), the draft declines naming the consumed top face rather
//! than serving a body whose top crosses itself.

use popkern_core::geometry::surface::SurfaceEval as _;
use popkern_core::geometry::{Point3, PortableMath, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_draft::{
    AngleVariation, DraftError, DraftOptions, DraftSpec, NeutralReference, draft_faces,
};

fn rib(place: &Transform3, tol: Tolerance) -> BRep {
    let rect = |z: f64| {
        vec![
            Point3::new(0.0, 0.0, z),
            Point3::new(10.0, 0.0, z),
            Point3::new(10.0, 1.0, z),
            Point3::new(0.0, 1.0, z),
        ]
    };
    BRep::polygon_frustum(&rect(0.0), &rect(10.0), tol)
        .and_then(|b| b.transform(place))
        .expect("rib")
}

#[test]
fn drafted_faces_serve_apart_and_decline_where_they_meet() {
    let tol = Tolerance::default();
    let turned =
        Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.7, 0.3), 0.6)
            .expect("axis");
    // The sides meet at the top edge at tan a = 0.5 / 10.
    let meeting = 0.05f64.patan();
    let mut problems = Vec::new();
    for (name, place) in [("identity", Transform3::identity()), ("turned", turned)] {
        let body = rib(&place, tol);
        let on_side = |p: Point3| {
            let p = place.apply_point(&p);
            body.faces
                .iter()
                .find(|(_, f)| {
                    f.surface.as_ref().is_some_and(|s| {
                        s.geometry
                            .closest_point(&p)
                            .is_some_and(|(u, v)| s.geometry.point_at(u, v).distance_to(&p) < 1e-9)
                    })
                })
                .map(|(h, _)| h)
                .expect("side face")
        };
        let sides = [
            on_side(Point3::new(5.0, 0.0, 5.0)),
            on_side(Point3::new(5.0, 1.0, 5.0)),
        ];
        // Width at height z is 1 − 2 z tan a (a > 0 tapers).
        let volume = |a: f64| 10.0 * 10.0 * (1.0 - 10.0 * a.ptan());
        for a in [
            2f64.to_radians(),
            -2f64.to_radians(),
            -5f64.to_radians(),
            5f64.to_radians(),
            meeting,
        ] {
            let specs: Vec<DraftSpec> = sides
                .iter()
                .map(|&face| DraftSpec {
                    face,
                    angle: a,
                    pull_direction: place.apply_direction(&Vec3::new(0.0, 0.0, 1.0)),
                    neutral: NeutralReference::Plane {
                        origin: place.apply_point(&Point3::new(0.0, 0.0, 0.0)),
                        normal: place.apply_direction(&Vec3::new(0.0, 0.0, 1.0)),
                    },
                    angle_variation: AngleVariation::Constant,
                })
                .collect();
            let served = draft_faces(&body, &specs, &DraftOptions::default(), tol, None);
            let meets = 10.0 * a.ptan() >= 0.5 - 1e-12;
            match (served, meets) {
                (Err(DraftError::FaceConsumed { .. }), true) => {}
                (Ok(r), false) => {
                    let v = popkern_measure::body_volume(&r.brep, &tol).expect("volume");
                    if let Some(first) = r.brep.validate_all().first() {
                        problems.push(format!("{name} a={a}: {first}"));
                    } else if (v - volume(a)).abs() > 1e-9 * volume(a) {
                        problems.push(format!("{name} a={a}: volume {v} for {}", volume(a)));
                    }
                }
                (other, _) => problems.push(format!(
                    "{name} a={a}: {:?}",
                    other.map(|r| r.brep.faces.len())
                )),
            }
        }
    }
    assert!(problems.is_empty(), "{problems:#?}");
}
