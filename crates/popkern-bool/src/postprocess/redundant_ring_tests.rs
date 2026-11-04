use crate::BoolOperation;
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn cube(x: f64, y: f64, z: f64, s: f64) -> BRep {
    let tol = Tolerance::default();
    BRep::unit_cube(tol)
        .and_then(|c| c.transform(&Transform3::scale(s, s, s)))
        .and_then(|c| c.transform(&Transform3::translation(x, y, z)))
        .expect("cube")
}

/// A 2-block with a 1-cube pocket open at its top, refilled flush.
fn refilled() -> BRep {
    let tol = Tolerance::default();
    let plug = || cube(0.5, 0.5, 1.0, 1.0);
    let pocketed = crate::boolean(
        BoolOperation::Subtraction,
        cube(0.0, 0.0, 0.0, 2.0),
        plug(),
        tol,
        None,
    )
    .expect("pocket");
    assert_eq!(pocketed.faces.len(), 11);
    crate::boolean(BoolOperation::Union, pocketed, plug(), tol, None).expect("refill")
}

#[test]
fn a_flush_refilled_pocket_serves_a_plain_block() {
    let body = refilled();
    assert_eq!(
        (body.vertices.len(), body.edges.len(), body.faces.len()),
        (8, 12, 6),
        "the pocket's rim dissolved with its plug face"
    );
    assert!(body.faces.values().all(|f| f.inner_loops.is_empty()));
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
}

#[test]
fn a_plug_on_another_surface_is_not_dissolved() {
    // The pocket refilled only half way: the plug's top lies 0.5 below
    // the block's, on another plane; nothing merges there.
    let tol = Tolerance::default();
    let plug = cube(0.5, 0.5, 1.0, 1.0);
    let pocketed = crate::boolean(
        BoolOperation::Subtraction,
        cube(0.0, 0.0, 0.0, 2.0),
        plug,
        tol,
        None,
    )
    .expect("pocket");
    let half = BRep::unit_cube(tol)
        .and_then(|c| c.transform(&Transform3::scale(1.0, 1.0, 1.0)))
        .and_then(|c| c.transform(&Transform3::translation(0.5, 0.5, 0.5)))
        .expect("half plug");
    let body = crate::boolean(BoolOperation::Union, pocketed, half, tol, None).expect("union");
    // The block's top keeps the opening (inner loop) and the pocket's
    // floor now sits at z = 1.5 — still a pocket, 0.5 deep.
    assert!(body.faces.values().any(|f| !f.inner_loops.is_empty()));
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
}
