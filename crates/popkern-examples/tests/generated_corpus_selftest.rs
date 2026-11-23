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

/// A tangent crossing's alignment flag is checked against the placed
/// geometry (negative control): the stored full tier held a cylinder pair
/// flagged opposed over an internal touch, drawn before a generator fix,
/// and the judge expected the union to pinch where the differences do.
/// Flipped, every drawn V scene fails its independent check.
#[test]
fn a_tangent_crossing_with_its_alignment_flipped_fails_its_check() {
    let mut checked = 0;
    for (a, b) in [(Kind::Y, Kind::Y), (Kind::S, Kind::Y), (Kind::Y, Kind::C)] {
        let cell = Cell {
            a,
            b,
            config: Config::V,
        };
        for index in 0..40 {
            let Ok(scene) = draw(cell, 7, index, Rung::FREE) else {
                continue;
            };
            assert!(
                generated_corpus::scene_check::check(&scene).is_ok(),
                "{a:?}{b:?}/{index}: a drawn scene meets its check"
            );
            let mut flipped = scene.clone();
            flipped.aligned = !flipped.aligned;
            assert!(
                generated_corpus::scene_check::check(&flipped).is_err(),
                "{a:?}{b:?}/{index}: aligned {} flipped still passes",
                scene.aligned
            );
            checked += 1;
            if checked % 3 == 0 {
                break;
            }
        }
    }
    assert!(checked >= 6, "only {checked} V scenes drawn");
}

/// A sphere steered onto its seam holding the other inside, touching it
/// (radii 1.15 and 0.9, centres 0.25 apart; the baseline's scene 696 of
/// SS-E, drawn before the steered normal took its hint from the nearest
/// triangle) is no external touch: its B lies outside A away from the
/// contact, but A lies inside B. Passed the check before A's side was read.
#[test]
fn an_external_touch_holding_the_other_inside_fails_its_check() {
    let scene: generated_corpus::scene::Scene =
        serde_json::from_str(r#"{"cell":{"a":"S","b":"S","config":"E"},"seed":1,"index":11099,"variant_a":0,"variant_b":1,"b_scale":1.0,"place_b":[[-0.8510899893929686,-0.32897063503569396,0.40917496409151033,0.22873863270652106],[0.4968120093744296,-0.2526118894752733,0.8302801097443403,-0.08094072980193162],[-0.16977531415496783,0.9099261258915166,0.3784320125513132,0.06022654038446695],[0.0,0.0,0.0,1.0]],"scene":[[-0.4821250291088469,-0.870957902134397,0.09480395042635947,2.0103910287641273],[0.54474156302414,-0.21326448499848893,0.8110332230878993,2.7348451249351724],[-0.6861574788784354,0.44266306842561354,0.577267114974697,-1.312568888918559],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":"Free","steer_b":{"Seam":0.0}},"aligned":false,"shift":0.0,"contact":[2.133060804763025,2.048285710246563,-0.7437206700070339],"curve":[],"curve_closed":false,"third":null}"#).expect("a scene descriptor");
    assert!(matches!(
        generated_corpus::scene_check::check(&scene),
        Err(generated_corpus::scene::Reject::WrongSide)
    ));
}
