//! PK-276 token-unification proof: ONE `CancellationToken` instance
//! cancels operations across crates.
//!
//! Before PK-276 the kernel carried three structurally identical but
//! distinct token types (`popkern_core::healing`, `popkern_bool`,
//! `popkern_blend`), so a UI Cancel button driving a feature rebuild
//! that spans a boolean and a fillet needed one token per crate. All
//! three module paths now re-export `popkern_core::cancellation::
//! CancellationToken`; this test passing to both entry points from a
//! single binding is the compile-time proof, and the assertions are
//! the runtime half.

use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
fn one_token_cancels_both_a_fillet_and_a_boolean() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit cube builds");

    // A single token — the same *instance* — handed to both crates.
    let token = popkern_core::cancellation::CancellationToken::new();
    token.cancel();

    let specs: Vec<popkern_blend::FilletSpec> = cube
        .edges
        .keys()
        .take(1)
        .map(|e| popkern_blend::FilletSpec {
            edge: e,
            radius: popkern_blend::FilletRadius::Constant(0.1),
        })
        .collect();
    let fillet_result = popkern_blend::fillet_edges_cancellable(
        &cube,
        &specs,
        &popkern_blend::FilletOptions::default(),
        tol,
        None,
        &token,
    );
    assert!(
        matches!(fillet_result, Err(popkern_blend::BlendError::Cancelled)),
        "shared token must cancel the fillet, got {fillet_result:?}"
    );

    let other = BRep::unit_cube(tol)
        .expect("unit cube builds")
        .transform(&Transform3::translation(0.5, 0.3, 0.2))
        .expect("translation applies");
    let bool_result = popkern_bool::boolean_cancellable(
        popkern_bool::BoolOperation::Union,
        &cube,
        &other,
        tol,
        &token,
    );
    assert!(
        matches!(bool_result, Err(popkern_bool::BoolError::Cancelled)),
        "shared token must cancel the boolean, got {bool_result:?}"
    );

    // The healing path names the same type: the token round-trips
    // through the healing module path unchanged.
    let via_healing: popkern_core::healing::CancellationToken = token.clone();
    assert!(via_healing.is_cancelled());

    // Untouched inputs: cancelled ops read only clones.
    assert!(cube.validate().is_ok());
}
