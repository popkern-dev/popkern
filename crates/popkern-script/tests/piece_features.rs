//! A body's pieces taken apart as recorded features (PK-1061, design:
//! operations on the cell complex, 5.4): `keep_pieces` keeps the solids
//! a geometric rule selects; `separate_pieces` moves them into a new
//! body, the source keeping the rest. Both name nothing new (every
//! kept entity keeps its name), save, reopen and replay with moved
//! tools.
use popkern_measure::body_volume;
use popkern_script::Plane;
use popkern_script::prelude::*;

const MIDDLE: f64 = 0.4 * 1.6 * 0.4;
const SIDE: f64 = 0.9 * 1.6 * 0.4;
/// The three pieces: the plate less the two 0.1 wide slots.
const FULL: f64 = 2.0 * SIDE + MIDDLE;

/// A 2.4 x 1.6 x 0.4 plate less two 0.1 wide slots through it, either
/// side of a 0.4 wide middle centred at parameter `x`: one body of three
/// pieces apart (x < -0.3, the middle, x > 0.3 at x = 0). Pieces apart
/// are lumps of one body; a split's touching pieces are bodies of their
/// own (PK-1120).
fn three_pieces(m: &mut Model) -> Body {
    m.param("x", 0.0, Unit::Meters).expect("param");
    let mut p = m
        .profile(Plane::XY)
        .rect(2.4, 1.6)
        .sketch("plate")
        .expect("sketch")
        .extrude(m, End::Blind(0.4.into()))
        .expect("extrude");
    for at in ["x - 250mm", "x + 250mm"] {
        let slot = m
            .tool(|m| {
                m.profile(Plane::xy_at(-0.5))
                    .rect(0.1, 3.0)
                    .at(popkern_script::scalar::expr(at), 0.0)
                    .sketch("slot")?
                    .extrude(m, End::Blind(2.0.into()))
            })
            .expect("slot");
        p = p.subtract(m, slot).expect("cut");
    }
    p
}

/// Solids, volume and lowest x of `body`; every face named, valid.
fn measured(m: &mut Model, body: Body) -> (usize, f64, f64) {
    let brep = m.brep_of(body).expect("brep");
    assert!(brep.validate_all().is_empty(), "{:?}", brep.validate_all());
    assert!(
        brep.faces.values().all(|f| f.name.is_some()),
        "every face is named"
    );
    let tol = popkern_core::tolerance::Tolerance::default();
    let lo = brep
        .vertices
        .values()
        .map(|v| v.position.x)
        .fold(f64::MAX, f64::min);
    (
        brep.solid_regions().len(),
        body_volume(&brep, &tol).expect("volume"),
        lo,
    )
}

#[test]
fn keep_pieces_keeps_the_largest_and_a_fillet_after_it_resolves() {
    let mut m = Model::new();
    let p = three_pieces(&mut m);
    let p = p.keep_pieces(&mut m, PieceRule::Largest(2)).expect("keep");
    let (solids, volume, _) = measured(&mut m, p);
    assert_eq!(solids, 2);
    assert!((volume - 2.0 * SIDE).abs() < 1e-9, "{volume}");
    let p = p
        .keep_pieces(&mut m, PieceRule::Containing([1.0, 0.0, 0.2]))
        .expect("keep the right one");
    p.fillet(&mut m, edges().straight().direction(Dir::Z), 0.05)
        .expect("fillet by name on the kept piece");
    let (solids, volume, lo) = measured(&mut m, p);
    assert_eq!(solids, 1);
    assert!(volume < SIDE && volume > SIDE - 0.01, "{volume}");
    assert!((lo - 0.3).abs() < 1e-9, "{lo}");
    // The tool moves: the pieces are cut again and the rules select
    // again; the fillet's edge names still resolve.
    m.set_param("x", 0.3).expect("move");
    m.replay().expect("replay with the tool moved");
    let (solids, _, lo) = measured(&mut m, p);
    assert_eq!(solids, 1);
    assert!((lo - 0.6).abs() < 1e-9, "the kept piece moved: {lo}");
}

#[test]
fn separate_pieces_moves_the_selected_pieces_into_a_new_body() {
    let mut m = Model::new();
    let p = three_pieces(&mut m);
    let right = p
        .separate_pieces(&mut m, PieceRule::Containing([1.0, 0.0, 0.2]))
        .expect("separate");
    assert_ne!(right, p);
    assert!(p.is_alive(&m) && right.is_alive(&m));
    let (solids, volume, lo) = measured(&mut m, right);
    assert_eq!(solids, 1);
    assert!((volume - SIDE).abs() < 1e-9 && (lo - 0.3).abs() < 1e-9);
    let (solids, volume, _) = measured(&mut m, p);
    assert_eq!(solids, 2);
    assert!((volume - (FULL - SIDE)).abs() < 1e-9, "{volume}");
    // A later feature on the source sees only what it kept.
    let p = p.keep_pieces(&mut m, PieceRule::Largest(1)).expect("keep");
    let (_, volume, lo) = measured(&mut m, p);
    assert!((volume - SIDE).abs() < 1e-9 && (lo + 1.2).abs() < 1e-9);

    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("piece_features");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join("separate.pkn");
    m.save(&path).expect("save");
    let mut again = Model::open(&path).expect("open");
    again.replay().expect("replay");
    assert_eq!(measured(&mut again, right), measured(&mut m, right));
    assert_eq!(measured(&mut again, p), measured(&mut m, p));
    // The tool moves under the separated body.
    m.set_param("x", -0.3).expect("move");
    m.replay().expect("replay with the tool moved");
    let (_, volume, lo) = measured(&mut m, right);
    assert!(lo.abs() < 1e-9, "{lo}");
    assert!((volume - 1.2 * 1.6 * 0.4).abs() < 1e-9, "{volume}");
}

#[test]
fn piece_rules_that_cannot_answer_fail_with_their_reason() {
    let mut m = Model::new();
    let p = three_pieces(&mut m);
    assert!(matches!(
        p.keep_pieces(&mut m, PieceRule::InsideTools(SetExpr::Body(0))),
        Err(ScriptError::Argument { .. })
    ));
    let nothing = p
        .keep_pieces(&mut m, PieceRule::Containing([9.0, 0.0, 0.0]))
        .expect_err("no piece there");
    assert!(
        format!("{nothing:?}").contains("PieceRuleMatchedNothing"),
        "{nothing:?}"
    );
    let every = p
        .separate_pieces(&mut m, PieceRule::All)
        .expect_err("the source would be empty");
    assert!(format!("{every:?}").contains("every piece"), "{every:?}");
    // Each failure left the document as it was.
    let (solids, volume, _) = measured(&mut m, p);
    assert_eq!(solids, 3);
    assert!((volume - FULL).abs() < 1e-9);
}
