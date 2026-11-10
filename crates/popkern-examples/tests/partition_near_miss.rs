//! The three-body tangent-records fixture moved off its touch (PK-1130):
//! a box [-10, 10]² × [-5, 0] and two balls of radius 3 on its top face,
//! their centres 6 + g apart, for g of ±half a tolerance, ±2 and ±10
//! tolerances, at identity and turned, under both policies. Within the
//! band the outcome is the touch's (2000 + 36 pi); past it a valid body
//! within |g| times the operands' area of it, or a typed decline that
//! is not a validation failure.
//!
//! At half a tolerance of overlap the union served the box ALONE
//! (volume 2000, validate clean): sphere x sphere read the balls' lens,
//! half a tolerance deep, as a crossing (its circle 1.2e-3 wide) where
//! the surfaces stay within the tolerance of each other throughout, and
//! the record, cut across a seam it grazed, lost both domes. A touch is
//! now judged by depth.

use std::f64::consts::PI;

use popkern_bool::{BoolError, BoolOptions, NonManifoldPolicy, SetExpr, nary_boolean};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// (placement, g) runs that fail today, each with its task.
const KNOWN: &[(&str, f64, &str)] = &[("turned", -2e-6, "PK-1149")];

#[test]
fn the_tangent_records_fixture_moved_off_its_touch() {
    let tol = Tolerance::default();
    let rect = |z: f64| {
        vec![
            Point3::new(-10.0, -10.0, z),
            Point3::new(10.0, -10.0, z),
            Point3::new(10.0, 10.0, z),
            Point3::new(-10.0, 10.0, z),
        ]
    };
    let turned = Transform3::rotation_around_axis(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.7, 0.3),
        0.645_771_823_237_9,
    )
    .expect("axis");
    let want = 2000.0 + 36.0 * PI;
    // The operands' area: the box and the two spheres.
    let area = 1200.0 + 2.0 * 36.0 * PI;
    let mut problems = Vec::new();
    for (name, place) in [("identity", Transform3::identity()), ("turned", turned)] {
        for g in [0.5e-6, -0.5e-6, 2e-6, -2e-6, 1e-5, -1e-5] {
            let known = KNOWN.iter().any(|&(n, kg, _)| n == name && kg == g);
            let bx = BRep::polygon_frustum(&rect(-5.0), &rect(0.0), tol)
                .and_then(|b| b.transform(&place))
                .expect("box");
            let ball = |x: f64| {
                BRep::sphere(3.0, tol)
                    .and_then(|b| b.transform(&Transform3::translation(x, 0.0, 0.0)))
                    .and_then(|b| b.transform(&place))
                    .expect("ball")
            };
            let (a, b) = (ball(-3.0 - g / 2.0), ball(3.0 + g / 2.0));
            for policy in [NonManifoldPolicy::Decline, NonManifoldPolicy::General] {
                let verdict = match nary_boolean(
                    &[&bx, &a, &b],
                    &SetExpr::Any,
                    tol,
                    &BoolOptions {
                        non_manifold: policy,
                    },
                ) {
                    Ok(body) if !body.validate_all().is_empty() => {
                        Err(format!("invalid: {:?}", body.validate_all()))
                    }
                    Ok(body) => {
                        let v = popkern_measure::body_volume(&body, &tol).expect("volume");
                        let bound = 1e-6 * want + g.abs() * area;
                        if (v - want).abs() <= bound {
                            Ok(())
                        } else {
                            Err(format!("volume {v} for {want}"))
                        }
                    }
                    Err(e @ BoolError::ValidationFailed(_)) => Err(format!("{e:?}")),
                    Err(_) => Ok(()),
                };
                match (verdict, known) {
                    (Err(why), false) => problems.push(format!("{name} g={g} {policy:?}: {why}")),
                    (Ok(()), true) => problems.push(format!(
                        "{name} g={g} {policy:?} now passes: drop it from KNOWN"
                    )),
                    _ => {}
                }
            }
        }
    }
    assert!(problems.is_empty(), "{problems:#?}");
}
