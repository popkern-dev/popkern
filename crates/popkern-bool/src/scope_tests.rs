//! The 64-slot disc's last cut (PK-813): three faces active, the 252
//! other slot walls far, and the unscoped result honours every premise;
//! a scope that leaves out the faces the tool touches is refused.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{active_scope, scope_violation};
use crate::error::BoolOperation;

fn slot(k: usize, tol: Tolerance) -> BRep {
    let a = std::f64::consts::TAU * k as f64 / 64.0;
    let (ca, sa, w) = (a.cos(), a.sin(), 0.6);
    let pts = [
        Point3::new(22.0 * ca + w * sa, 22.0 * sa - w * ca, -1.0),
        Point3::new(27.0 * ca + w * sa, 27.0 * sa - w * ca, -1.0),
        Point3::new(27.0 * ca - w * sa, 27.0 * sa + w * ca, -1.0),
        Point3::new(22.0 * ca - w * sa, 22.0 * sa + w * ca, -1.0),
    ];
    BRep::prism(&pts, Vec3::new(0.0, 0.0, 1.0), 4.0, tol).expect("slot prism")
}

#[test]
fn the_last_slot_cut_touches_three_faces_and_carries_the_rest() {
    let tol = Tolerance::default();
    let mut body = BRep::cylinder(30.0, 1.5, tol).expect("disc");
    for k in 0..63 {
        body = crate::boolean(BoolOperation::Subtraction, body, slot(k, tol), tol, None)
            .expect("slot");
    }
    let tool = slot(63, tol);
    let hits = crate::intersect::find_intersections(&body, &tool, &tol).expect("intersections");
    let scope = active_scope(&body, &hits.aabbs_a, &hits.aabbs_b, 10.0 * tol.linear());
    assert_eq!(
        (scope.active.len(), scope.far.len(), scope.frontier.len()),
        (3, 252, 504),
        "top, bottom and outer wall active; the other slots' walls far"
    );
    // Stage 2b: the top and bottom faces' 63 earlier slot holes each lie
    // away from the tool and are carried like far faces.
    assert_eq!(scope.far_holes.values().map(Vec::len).sum::<usize>(), 126);
    let result =
        crate::boolean(BoolOperation::Subtraction, body.clone(), tool, tol, None).expect("cut");
    assert_eq!(
        scope_violation(BoolOperation::Subtraction, &body, &hits, &scope, &result),
        None
    );
    // Negative control: a pad that empties the tool box calls the faces
    // it touches far.
    let blind = active_scope(&body, &hits.aabbs_a, &hits.aabbs_b, -3.0);
    assert!(
        scope_violation(BoolOperation::Subtraction, &body, &hits, &blind, &result).is_some(),
        "a scope missing the touched faces must be refused"
    );
}

/// Stage 2 (PK-814): the scoped merge — graph from the three active
/// faces, the 252 far faces carried — assembles the same body as the
/// unscoped merge: counts, and every face's surface kind with its loop
/// vertex positions bit for bit. A scope that calls a touched face far
/// cannot stitch its frontier and declines (the caller reruns unscoped).
#[test]
fn the_scoped_merge_assembles_the_unscoped_body() {
    let tol = Tolerance::default();
    let mut body = BRep::cylinder(30.0, 1.5, tol).expect("disc");
    for k in 0..63 {
        body = crate::boolean(BoolOperation::Subtraction, body, slot(k, tol), tol, None)
            .expect("slot");
    }
    let tool = slot(63, tol);
    let hits = crate::intersect::find_intersections(&body, &tool, &tol).expect("intersections");
    let scope = active_scope(&body, &hits.aabbs_a, &hits.aabbs_b, 10.0 * tol.linear());
    let merge = |s: Option<&super::ActiveScope>| {
        crate::cellular::boolean::cellular_graph_merge(
            &body,
            &tool,
            &hits,
            BoolOperation::Subtraction,
            &tol,
            s,
        )
    };
    let (unscoped, _) = merge(None).expect("unscoped merge");
    let (scoped, carried) = merge(Some(&scope)).expect("scoped merge");
    let counts = |b: &BRep| (b.faces.len(), b.edges.len(), b.vertices.len());
    assert_eq!(counts(&scoped), counts(&unscoped));
    assert!(scoped.validate().is_ok());
    let keys = |b: &BRep| {
        let mut k: Vec<_> = b.faces.keys().map(|f| super::face_key(b, f)).collect();
        k.sort();
        k
    };
    assert_eq!(keys(&scoped), keys(&unscoped), "same faces, bit for bit");
    assert!(
        carried.edges.len() >= scope.frontier.len(),
        "the far walls' edges are carried: {}",
        carried.edges.len()
    );
    // Negative control: the blind scope's far faces include the ones the
    // tool cuts; their frontier edges never come through unsplit.
    let blind = active_scope(&body, &hits.aabbs_a, &hits.aabbs_b, -3.0);
    let declined = merge(Some(&blind))
        .map(|(b, _)| b.validate().is_err() || counts(&b) != counts(&unscoped))
        .unwrap_or(true);
    assert!(declined, "a blind scope must not assemble the cut body");
}
