//! A shell whose skin touches itself where a wall closes: the tangency
//! foundations' shell consumer (design §8.4, F3).
//!
//! A box [-10, 10]² × [0, 10] with a through hole of radius 3 along x
//! at height 5 leaves walls 2 thick under and over the hole. Shelled by
//! 0.9 it serves its closed-form volume; shelled by 1 the hole's offset
//! cylinder is tangent to the floor's and the top's offset planes along
//! two lines, the cavity split in two with the material meeting itself
//! there. That body touches itself with no shared topology; it was
//! served as a plain solid, validate_all clean. It now declines with
//! the contact (`ShellError::SelfContact`), at identity and turned.

use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::{ShellDirection, ShellError, ShellOptions, shell};
use std::f64::consts::{FRAC_PI_2, PI};

fn holed_box(place: &Transform3, tol: Tolerance) -> BRep {
    let rect = |z: f64| {
        vec![
            Point3::new(-10.0, -10.0, z),
            Point3::new(10.0, -10.0, z),
            Point3::new(10.0, 10.0, z),
            Point3::new(-10.0, 10.0, z),
        ]
    };
    let bx = BRep::polygon_frustum(&rect(0.0), &rect(10.0), tol).expect("box");
    let rod = BRep::cylinder(3.0, 30.0, tol)
        .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, -15.0)))
        .and_then(|c| c.transform(&Transform3::rotation_y(FRAC_PI_2)))
        .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, 5.0)))
        .expect("rod");
    popkern_bool::boolean(popkern_bool::BoolOperation::Subtraction, bx, rod, tol, None)
        .and_then(|b| b.transform(place).map_err(Into::into))
        .expect("hole")
}

#[test]
fn a_shell_whose_wall_closes_exactly_declines_with_the_contact() {
    let tol = Tolerance::default();
    let turned =
        Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.7, 0.3), 0.6)
            .expect("axis");
    // The body less the cavity: the box shrunk by t less the hole grown
    // by t, through the cavity's length 20 − 2t.
    let body = 4000.0 - 9.0 * PI * 20.0;
    let shelled = |t: f64| {
        body - (20.0 - 2.0 * t).powi(2) * (10.0 - 2.0 * t)
            + PI * (3.0 + t).powi(2) * (20.0 - 2.0 * t)
    };
    let mut problems = Vec::new();
    for (name, place) in [("identity", Transform3::identity()), ("turned", turned)] {
        let input = holed_box(&place, tol);
        let run = |t: f64| {
            shell(
                &input,
                t,
                ShellDirection::Inward,
                &[],
                &[],
                &ShellOptions::default(),
                tol,
                None,
            )
        };
        match run(0.9) {
            Ok(r) => {
                let v = popkern_measure::body_volume(&r.brep, &tol).expect("volume");
                if (v - shelled(0.9)).abs() > 1e-9 * shelled(0.9) {
                    problems.push(format!("{name} 0.9: volume {v} for {}", shelled(0.9)));
                }
            }
            Err(e) => problems.push(format!("{name} 0.9: declined {e}")),
        }
        match run(1.0) {
            Err(ShellError::SelfContact { location, .. }) => {
                let p = place.inverse().expect("inverse").apply_point(&location);
                let on_line =
                    p.y.abs() < 1e-3 && ((p.z - 1.0).abs() < 1e-5 || (p.z - 9.0).abs() < 1e-5);
                if !on_line {
                    problems.push(format!("{name} 1.0: contact at {p:?}"));
                }
            }
            other => problems.push(format!(
                "{name} 1.0: {:?}",
                other.map(|r| r.brep.faces.len())
            )),
        }
    }
    assert!(problems.is_empty(), "{problems:#?}");
}
