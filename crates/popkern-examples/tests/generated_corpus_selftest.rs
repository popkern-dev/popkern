//! The generated corpus checks itself (PK-1227, design §9): its operands
//! are valid, its contact placement is exact against closed forms, and
//! its oracles fail a scene whose expectation is wrong (negative
//! controls: a corpus that cannot fail measures nothing).

#[path = "generated_corpus/mod.rs"]
mod generated_corpus;

use generated_corpus::judge::judge;
use generated_corpus::operand::{KINDS, Kind, operand};
use generated_corpus::scene::{Cell, Config, Rung, draw, place};

#[test]
fn every_operand_kind_builds_valid_operands_whose_two_volumes_agree() {
    for kind in KINDS {
        for variant in 0..3u8 {
            let op = operand(kind, variant);
            // An inscribed tessellation loses up to area × chord.
            assert!(
                (op.volume - op.mesh_volume).abs() <= 1.5 * op.area * op.chord + 1e-9 * op.volume,
                "{kind:?}/{variant}: B-rep {} mesh {}",
                op.volume,
                op.mesh_volume
            );
            assert!(
                !op.relevant.is_empty(),
                "{kind:?}/{variant}: no relevant face"
            );
        }
    }
}

#[test]
fn an_external_sphere_touch_is_placed_at_the_closed_form_distance() {
    let cell = Cell {
        a: Kind::S,
        b: Kind::S,
        config: Config::E,
    };
    for index in 0..6 {
        let scene = draw(cell, 7, index, Rung::FREE).expect("a sphere pair always touches");
        let placed = place(&scene);
        let centre = |b: &popkern_core::topology::BRep| {
            b.faces
                .values()
                .find_map(|f| match f.surface.as_ref().map(|s| &s.geometry) {
                    Some(popkern_core::geometry::surface::SurfaceGeometry::Sphere {
                        center,
                        radius,
                        ..
                    }) => Some((*center, *radius)),
                    _ => None,
                })
                .expect("sphere")
        };
        let ((ca, ra), (cb, rb)) = (centre(&placed.a), centre(&placed.b));
        let gap = ca.distance_to(&cb) - ra - rb;
        assert!(
            gap.abs() <= 1e-12 * (ra + rb),
            "scene {index}: centres {gap} off tangency"
        );
    }
}

#[test]
fn the_oracles_fail_a_scene_whose_expectation_is_wrong() {
    // A real touch judged as if apart: the union's decline must fail.
    let cell = Cell {
        a: Kind::S,
        b: Kind::Y,
        config: Config::E,
    };
    let scene = draw(cell, 3, 0, Rung::FREE).expect("draw");
    assert!(
        judge(&scene).failures.is_empty(),
        "the true expectation passes"
    );
    let wrong = generated_corpus::scene::Scene {
        offset: 5.0,
        ..scene.clone()
    };
    assert!(
        !judge(&wrong).failures.is_empty(),
        "claiming apart must fail"
    );
    // The same geometry claimed as an internal touch must fail too.
    let flipped = generated_corpus::scene::Scene {
        cell: Cell {
            config: Config::I,
            ..cell
        },
        ..scene
    };
    assert!(
        !judge(&flipped).failures.is_empty(),
        "a flipped classification must fail"
    );
}

#[test]
fn a_containment_claimed_of_a_disjoint_pair_fails() {
    let cell = Cell {
        a: Kind::P,
        b: Kind::S,
        config: Config::D,
    };
    let scene = draw(cell, 5, 0, Rung::FREE).expect("draw");
    assert!(
        judge(&scene).failures.is_empty(),
        "disjoint passes as disjoint"
    );
    let wrong = generated_corpus::scene::Scene {
        cell: Cell {
            config: Config::N,
            ..cell
        },
        ..scene
    };
    assert!(
        !judge(&wrong).failures.is_empty(),
        "containment claimed of a disjoint pair must fail"
    );
}

/// A kernel panic is the scene's failure, not the pair's: a panicking
/// judge yields a `panic` failure carrying the message (PK-1260; before,
/// a mesher panic in one scene lost the NURBS ball pair's whole report).
#[test]
fn a_panicking_scene_is_recorded_as_a_panic_failure() {
    let v = generated_corpus::judge::caught(|| panic!("index out of bounds"));
    assert_eq!(v.failures.len(), 1);
    assert_eq!(v.failures[0].oracle, "panic");
    assert!(v.failures[0].detail.contains("index out of bounds"));
}
