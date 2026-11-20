//! Two exact circles close together are two records, not one loop and its
//! duplicate.

use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// A ball 1e-5 wider than a cylinder on its axis crosses its wall in two
/// exact circles 8e-3 apart. The loop dedupe compared samples against the
/// kept loop's eight-chord polyline within a slop gate of 5% of the loop's
/// extent, and dropped the second circle: the operands' common and
/// differences served inside out (the W1.3 corpus's YS-K cells). A loop
/// is a duplicate only where its samples lie on the kept loop's curves
/// within the lateral band.
#[test]
fn two_close_crossing_circles_both_survive() {
    let tol = Tolerance::default();
    let (rc, delta) = (0.8, 1e-5);
    let cylinder = BRep::cylinder(rc, 3.0, tol)
        .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, -1.5)))
        .expect("cylinder");
    let ball = BRep::sphere(rc + delta, tol).expect("ball");
    let ix = crate::intersect::find_intersections(&cylinder, &ball, &tol).expect("intersect");
    let circles = ix.records.iter().filter(|r| r.curve.is_closed).count();
    assert_eq!(circles, 2, "{} records", ix.records.len());
}
