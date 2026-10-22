//! PK-913: a mixed-sense pair (an outward disc, an inward prism boss)
//! assembles on the cellular path itself, with no rescue retry, for
//! every operation in both orders and both PK-911 crossing depths,
//! scoped as the serve scopes it and unscoped; a same-sense pair is
//! left alone.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::BoolOperation;

fn boss(z0: f64, tol: Tolerance) -> BRep {
    let pts = [
        Point3::new(-2.0, -2.0, z0),
        Point3::new(2.0, -2.0, z0),
        Point3::new(2.0, 2.0, z0),
        Point3::new(-2.0, 2.0, z0),
    ];
    BRep::prism(&pts, Vec3::new(0.0, 0.0, 1.0), 2.0, tol).expect("boss")
}

fn merge(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    scoped: bool,
) -> Result<BRep, String> {
    let hits = crate::intersect::find_intersections(a, b, &tol).map_err(|e| format!("{e:?}"))?;
    let scopes = crate::scope::serving_scopes(a, b, &hits, op, &tol);
    crate::cellular::boolean::cellular_graph_merge(a, b, &hits, op, &tol, scoped.then_some(&scopes))
        .map(|(body, _)| body)
        .map_err(|d| format!("{d:?}").chars().take(300).collect())
}

#[test]
fn a_mixed_sense_pair_assembles_on_the_cellular_path() {
    let tol = Tolerance::default();
    let disc = BRep::cylinder(30.0, 1.5, tol).expect("disc");
    let mut failures = Vec::new();
    for z0 in [1.4, 0.5] {
        let boss = boss(z0, tol);
        for op in [
            BoolOperation::Union,
            BoolOperation::Intersection,
            BoolOperation::Subtraction,
        ] {
            for (order, a, b) in [("disc, boss", &disc, &boss), ("boss, disc", &boss, &disc)] {
                for scoped in [false, true] {
                    let label = format!("{op:?} ({order}) z0 {z0} scoped {scoped}");
                    match merge(op, a, b, tol, scoped) {
                        Ok(body) if body.validate().is_ok() => {}
                        Ok(body) => failures.push(format!("{label}: {:?}", body.validate().err())),
                        Err(e) => failures.push(format!("{label}: {e}")),
                    }
                }
            }
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}
