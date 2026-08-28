//! Cellular→graph-merge fallback honesty (PK-351).
//!
//! The cellular selection is the boolean pipeline's primary route;
//! graph-merge is its explicit fallback. These pins hold two facts:
//!
//! 1. On a representative operand corpus the cellular route serves
//!    every operation — the fallback registry stays EMPTY. If this
//!    ever fails, the recorded reasons say exactly which operation
//!    fell back and why (the reason used to be discarded with
//!    `.ok()`, unlogged on all builds).
//! 2. The recording mechanism itself works: forcing the kill-switch
//!    records `CellularDisabled` and the operation still serves
//!    through graph-merge.
//!
//! Process-global registry: safe under nextest (process per test).

use popkern_bool::merge_fallback::{MergeFallbackReason, take_merge_fallback_events};
use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn shifted_cube(dx: f64, dy: f64, dz: f64, scale: f64) -> BRep {
    BRep::unit_cube(tol())
        .expect("cube builds")
        .transform(&Transform3::scale_uniform(scale))
        .expect("scales")
        .transform(&Transform3::translation(dx, dy, dz))
        .expect("translates")
}

#[test]
fn cellular_serves_the_corpus_without_fallback() {
    let _ = take_merge_fallback_events(); // drain any prior state
    let cube = shifted_cube(0.0, 0.0, 0.0, 1.0);
    let offset_cube = shifted_cube(0.5, 0.5, 0.5, 1.0);
    let cylinder = BRep::cylinder(0.4, 2.0, tol()).expect("cylinder builds");
    // Sphere at the cube's CENTER — the generic curved-surface case.
    // (A sphere centered exactly ON a cube vertex is the recorded
    // degeneracy pinned separately below.)
    let sphere = BRep::sphere(0.7, tol())
        .expect("sphere builds")
        .transform(&Transform3::translation(0.5, 0.5, 0.5))
        .expect("translates");

    let cases: Vec<(&str, BoolOperation, &BRep, &BRep)> = vec![
        ("cube-union", BoolOperation::Union, &cube, &offset_cube),
        (
            "cube-subtract",
            BoolOperation::Subtraction,
            &cube,
            &offset_cube,
        ),
        (
            "cube-intersect",
            BoolOperation::Intersection,
            &cube,
            &offset_cube,
        ),
        (
            "cube-cylinder",
            BoolOperation::Subtraction,
            &cube,
            &cylinder,
        ),
        ("cube-sphere", BoolOperation::Union, &cube, &sphere),
    ];
    for (label, op, a, b) in cases {
        let result = boolean(op, a.clone(), b.clone(), tol(), None)
            .unwrap_or_else(|e| panic!("{label}: {e}"));
        assert!(result.validate().is_ok(), "{label}: result must validate");
    }

    let events = take_merge_fallback_events();
    assert!(
        events.is_empty(),
        "cellular selection must serve the corpus; fallbacks fired: {events:?}"
    );
}

/// A sphere centered exactly ON a cube vertex: the intersection curve
/// runs through the corner vertex itself. Cellular assembly leaves a
/// shell gap (recorded as the fallback reason), and graph-merge's
/// strand gate then declines transactionally rather than serving a
/// body with missing regions. This pin holds BOTH honesty contracts:
/// the fallback is recorded, and the decline is a typed error — never
/// a silent route change or a silently-wrong body.
#[test]
fn corner_centered_sphere_union_records_fallback_and_declines() {
    let _ = take_merge_fallback_events();
    let cube = shifted_cube(0.0, 0.0, 0.0, 1.0);
    let sphere = BRep::sphere(0.7, tol()).expect("sphere builds");
    let result = boolean(BoolOperation::Union, cube, sphere, tol(), None);
    let events = take_merge_fallback_events();
    match result {
        Err(_) => {
            assert!(
                events
                    .iter()
                    .any(|e| matches!(e.reason, MergeFallbackReason::CellularInvalid(_))),
                "the cellular decline must be recorded before graph-merge declines: {events:?}"
            );
        }
        // If a future campaign makes this configuration serve, it must
        // do so through the primary route with the registry silent.
        Ok(body) => {
            assert!(body.validate().is_ok(), "served body must validate");
            assert!(
                events.is_empty(),
                "a served result must not hide fallback events: {events:?}"
            );
        }
    }
}

#[test]
fn kill_switch_fallback_is_recorded_and_still_serves() {
    // Process-scoped: safe under nextest's process-per-test model.
    // SAFETY-adjacent note: set_var is fine here — single-threaded
    // test entry, before any boolean runs.
    unsafe { std::env::set_var("POPKERN_CELLULAR", "0") };
    let _ = take_merge_fallback_events();

    let cube = shifted_cube(0.0, 0.0, 0.0, 1.0);
    let other = shifted_cube(0.4, 0.4, 0.4, 1.0);
    let result =
        boolean(BoolOperation::Union, cube, other, tol(), None).expect("graph-merge serves");
    assert!(result.validate().is_ok(), "graph-merge result validates");

    let events = take_merge_fallback_events();
    assert!(
        !events.is_empty(),
        "the forced fallback must be recorded, not silent"
    );
    assert!(
        events
            .iter()
            .all(|e| e.reason == MergeFallbackReason::CellularDisabled),
        "reason must be CellularDisabled: {events:?}"
    );
}
