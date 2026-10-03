//! PK-630: interference and clearance between part occurrences —
//! exact overlap volumes, certified clearances, contact, containment.

use popkern_assembly::{Assembly, InterferenceOptions, PairStatus, Parent, Placement};
use popkern_script::Plane;
use popkern_script::prelude::*;

fn block(w: f64, h: f64, t: f64) -> Model {
    let mut m = Model::new();
    m.profile(Plane::XY)
        .rect(w, h)
        .sketch("block")
        .expect("sketch")
        .extrude(&mut m, End::Blind(t.into()))
        .expect("extrude");
    m
}

fn rod(r: f64, t: f64) -> Model {
    let mut m = Model::new();
    m.profile(Plane::XY)
        .circle(r)
        .sketch("rod")
        .expect("sketch")
        .extrude(&mut m, End::Blind(t.into()))
        .expect("extrude");
    m
}

fn pair(a: Model, b: Model, at: Placement) -> Assembly {
    let mut asm = Assembly::new("pair");
    let da = asm.add_part("a", a).expect("a");
    let db = asm.add_part("b", b).expect("b");
    asm.place(Parent::Root, da, Some("a"), Placement::identity())
        .expect("place a");
    asm.place(Parent::Root, db, Some("b"), at).expect("place b");
    asm
}

fn check(asm: &mut Assembly, clearance: f64) -> Vec<PairStatus> {
    let report = asm
        .interference(&InterferenceOptions {
            clearance,
            ..InterferenceOptions::default()
        })
        .expect("report");
    report.pairs.into_iter().map(|p| p.status).collect()
}

#[test]
fn overlapping_blocks_report_their_exact_overlap_volume() {
    let mut asm = pair(
        block(1.0, 1.0, 1.0),
        block(1.0, 1.0, 1.0),
        Placement::translation(0.5, 0.0, 0.0),
    );
    match check(&mut asm, 0.0).as_slice() {
        [
            PairStatus::Interfering {
                volume: Some(v),
                decline: None,
            },
        ] => {
            assert!((v - 0.5).abs() < 1e-9, "{v}")
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_block_inside_another_is_interference_by_containment() {
    let mut asm = pair(
        block(2.0, 2.0, 2.0),
        block(0.5, 0.5, 0.5),
        Placement::translation(0.0, 0.0, 0.5),
    );
    match check(&mut asm, 0.0).as_slice() {
        [
            PairStatus::Interfering {
                volume: Some(v), ..
            },
        ] => assert!((v - 0.125).abs() < 1e-9, "{v}"),
        other => panic!("{other:?}"),
    }
}

#[test]
fn separated_blocks_report_a_certified_clearance() {
    let mut asm = pair(
        block(1.0, 1.0, 1.0),
        block(1.0, 1.0, 1.0),
        Placement::translation(1.25, 0.0, 0.0),
    );
    match check(&mut asm, 1.0).as_slice() {
        [
            PairStatus::Clear {
                distance,
                lower_bound,
                refined,
            },
        ] => {
            assert!((distance - 0.25).abs() < 1e-12, "{distance}");
            assert!(
                (lower_bound - 0.25).abs() < 1e-12,
                "planar faces mesh exactly: {lower_bound}"
            );
            assert!(refined);
        }
        other => panic!("{other:?}"),
    }
    // Beyond the requested clearance: not reported.
    assert!(check(&mut asm, 0.1).is_empty());
}

#[test]
fn blocks_sharing_a_face_touch() {
    let mut asm = pair(
        block(1.0, 1.0, 1.0),
        block(1.0, 1.0, 1.0),
        Placement::translation(1.0, 0.0, 0.0),
    );
    assert!(
        matches!(
            check(&mut asm, 0.0).as_slice(),
            [PairStatus::Touching { .. }]
        ),
        "{:?}",
        check(&mut asm, 0.0)
    );
}

#[test]
fn a_rod_through_a_block_reports_the_analytic_overlap() {
    // A block [0, 2] x [-1, 1] x [0, 1] cuts a radius-0.5 rod of height
    // 1 centred on the origin in half: pi r^2 h / 2.
    let mut asm = pair(
        rod(0.5, 1.0),
        block(2.0, 2.0, 1.0),
        Placement::translation(1.0, 0.0, 0.0),
    );
    let want = std::f64::consts::PI * 0.25 / 2.0;
    match check(&mut asm, 0.0).as_slice() {
        [
            PairStatus::Interfering {
                volume: Some(v), ..
            },
        ] => {
            assert!((v - want).abs() < 2e-3 * want, "{v} vs {want}")
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_pin_in_a_bore_reports_its_radial_clearance() {
    // A plate with a radius-0.075 bore, a radius-0.07 pin through it:
    // 0.005 radial clearance, curved faces on both sides.
    let mut plate = Model::new();
    let body = plate
        .profile(Plane::XY)
        .rect(1.0, 1.0)
        .sketch("plate")
        .expect("sketch")
        .extrude(&mut plate, End::Blind(0.2.into()))
        .expect("extrude");
    let bore = plate
        .tool(|m| {
            m.profile(Plane::XY)
                .circle(0.075)
                .sketch("bore")?
                .extrude(m, End::Blind(1.0.into()))
        })
        .expect("bore");
    body.subtract(&mut plate, bore).expect("subtract");
    let mut asm = pair(
        plate,
        rod(0.07, 0.6),
        Placement::translation(0.0, 0.0, -0.2),
    );
    match check(&mut asm, 0.05).as_slice() {
        [
            PairStatus::Clear {
                distance,
                lower_bound,
                ..
            },
        ] => {
            assert!((distance - 0.005).abs() < 1e-6, "{distance}");
            assert!(
                *lower_bound > 0.0 && *lower_bound <= *distance,
                "{lower_bound}"
            );
        }
        other => panic!("{other:?}"),
    }
}
