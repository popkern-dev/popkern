//! The partition's operations as recorded features (PK-1052, design:
//! operations on the cell complex, 5.4): a split keeping the pieces a
//! rule selects, the N-ary boolean and symmetric difference, the
//! imprint. Each consumes its tools, names its result, saves and
//! replays.
use popkern_measure::body_volume;
use popkern_part::prelude::*;
use popkern_script::Plane;
use popkern_script::prelude::*;

/// A 2.4 x 1.6 x 0.4 plate on z = 0.
fn plate(m: &mut Model) -> Body {
    m.profile(Plane::XY)
        .rect(2.4, 1.6)
        .sketch("plate")
        .expect("sketch")
        .extrude(m, End::Blind(0.4.into()))
        .expect("extrude")
}

/// A slab `width` wide in x, centred at `x`, through the plate.
fn slab(m: &mut Model, x: f64, width: f64) -> Body {
    m.tool(|m| {
        m.profile(Plane::xy_at(-0.5))
            .rect(width, 3.0)
            .at(x, 0.0)
            .sketch("slab")?
            .extrude(m, End::Blind(2.0.into()))
    })
    .expect("slab")
}

fn measured(m: &mut Model) -> (usize, f64) {
    let part = m.to_part().expect("part");
    let tol = popkern_core::tolerance::Tolerance::default();
    let brep = part.brep();
    assert!(brep.validate().is_ok(), "{:?}", brep.validate());
    assert!(
        brep.faces.values().all(|f| f.name.is_some()),
        "every face is named"
    );
    (
        brep.solid_regions().len(),
        body_volume(brep, &tol).expect("volume"),
    )
}

/// Each of `bodies`: one solid, valid with nothing touching inside it,
/// every face named; their count and total volume.
fn measured_bodies(m: &mut Model, bodies: &[Body]) -> (usize, f64) {
    let tol = popkern_core::tolerance::Tolerance::default();
    let mut volume = 0.0;
    for &b in bodies {
        let brep = m.brep_of(b).expect("brep");
        assert!(brep.validate_all().is_empty(), "{:?}", brep.validate_all());
        assert_eq!(brep.solid_regions().len(), 1, "one piece per body");
        assert!(brep.faces.values().all(|f| f.name.is_some()));
        volume += body_volume(&brep, &tol).expect("volume");
    }
    (bodies.len(), volume)
}

fn split_with(keep: PieceRule) -> Result<(usize, f64), ScriptError> {
    let mut m = Model::new();
    let p = plate(&mut m);
    let s = slab(&mut m, 0.0, 0.4);
    let bodies = p.split(&mut m, &[s], keep)?;
    assert!(!s.is_alive(&m), "the tool is consumed");
    assert_eq!(bodies[0], p, "the plate continues as one piece");
    Ok(measured_bodies(&mut m, &bodies))
}

#[test]
fn a_split_keeps_the_pieces_its_rule_selects() {
    let full = 2.4 * 1.6 * 0.4;
    let middle = 0.4 * 1.6 * 0.4;
    let side = (full - middle) / 2.0;
    let cases = [
        (PieceRule::All, 3, full),
        (PieceRule::InsideTools(SetExpr::Body(0)), 1, middle),
        (
            PieceRule::InsideTools(SetExpr::Minus(
                Box::new(SetExpr::Any),
                Box::new(SetExpr::Body(0)),
            )),
            0,
            0.0,
        ),
        (PieceRule::Containing([1.0, 0.0, 0.2]), 1, side),
        (PieceRule::Largest(2), 2, 2.0 * side),
    ];
    for (rule, solids, volume) in cases {
        let at = format!("{rule:?}");
        if solids == 0 {
            // Any minus body 0 over one tool keeps nothing: no piece is
            // inside a tool other than the slab.
            let failure = format!("{:?}", split_with(rule).expect_err(&at));
            assert!(
                failure.contains("PieceRuleMatchedNothing") && failure.contains("pieces: 3"),
                "{at}: {failure}"
            );
            continue;
        }
        let (got_solids, got_volume) = split_with(rule).unwrap_or_else(|e| panic!("{at}: {e}"));
        assert_eq!(got_solids, solids, "{at}");
        assert!((got_volume - volume).abs() < 1e-9, "{at}: {got_volume}");
    }
}

#[test]
fn the_symmetric_difference_and_nary_boolean_record_and_replay() {
    let full = 2.4 * 1.6 * 0.4;
    let mut m = Model::new();
    let p = plate(&mut m);
    // Two slabs overlapping the plate and each other.
    let a = slab(&mut m, -0.3, 0.8);
    let b = slab(&mut m, 0.3, 0.8);
    m.label_body("a", a).expect("label");
    let kept = p
        .nary_boolean(
            &mut m,
            &[a, b],
            SetExpr::AtLeast(2),
            NonManifoldPolicy::Decline,
        )
        .expect("at least two");
    assert!(!a.is_alive(&m) && !b.is_alive(&m));
    assert_eq!(kept, p);
    // In at least two of the three: the plate within either slab (x in
    // [-0.7, 0.7]) and the slabs' own overlap (x in [-0.1, 0.1], all
    // 3 x 2 of it), the triple overlap counted once.
    let (_, v) = measured(&mut m);
    let want = 1.4 * 1.6 * 0.4 + 0.2 * 3.0 * 2.0 - 0.2 * 1.6 * 0.4;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");

    let mut m = Model::new();
    let p = plate(&mut m);
    let s = slab(&mut m, 0.0, 0.4);
    let xor = p.symmetric_difference(&mut m, s).expect("xor");
    // Its pieces touch along the crossing curves: each its own body.
    let bodies = xor.served_bodies(&m);
    assert!(bodies.len() > 1, "{bodies:?}");
    // The plate minus the slab's crossing, plus the slab's part outside
    // the plate (0.4 x 3 x 2 less the 0.4 x 1.6 x 0.4 inside).
    let (_, v) = measured_bodies(&mut m, &bodies);
    let inside = 0.4 * 1.6 * 0.4;
    assert!(
        (v - (full - inside + 0.4 * 3.0 * 2.0 - inside)).abs() < 1e-9,
        "{v}"
    );

    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("partition_features");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join("xor.pkn");
    m.save(&path).expect("save");
    let mut again = Model::open(&path).expect("open");
    again.replay().expect("replay");
    let (_, w) = measured_bodies(&mut again, &bodies);
    assert!((w - v).abs() < 1e-12, "{w} vs {v}");
}

/// A split serves each kept piece as its own body (PK-1120): when the
/// tool moves, each body keeps its piece (found by a face name only it
/// carries, or in order when its faces were renamed), a reopened
/// document serves the same pieces as the same bodies, and a feature on
/// a piece's body acts on that piece alone.
#[test]
fn split_pieces_are_bodies_that_keep_their_pieces() {
    let mut m = Model::new();
    m.param("x", 0.0, Unit::Meters).expect("param");
    let p = plate(&mut m);
    let s = m
        .tool(|m| {
            m.profile(Plane::xy_at(-0.5))
                .rect(0.4, 3.0)
                .at(popkern_script::scalar::p("x"), 0.0)
                .sketch("slab")?
                .extrude(m, End::Blind(2.0.into()))
        })
        .expect("slab");
    let bodies = p.split(&mut m, &[s], PieceRule::All).expect("split");
    assert_eq!(bodies.len(), 3);
    let full = 2.4 * 1.6 * 0.4;
    let (_, v) = measured_bodies(&mut m, &bodies);
    assert!((v - full).abs() < 1e-9, "{v}");
    // Each body's piece by where it lies in x.
    let spans = |m: &mut Model| -> Vec<(f64, f64)> {
        bodies
            .iter()
            .map(|&b| {
                let brep = m.brep_of(b).expect("brep");
                let xs = brep.vertices.values().map(|v| v.position.x);
                let lo = xs.clone().fold(f64::MAX, f64::min);
                (lo, xs.fold(f64::MIN, f64::max))
            })
            .collect()
    };
    let before = spans(&mut m);
    let right = (0..3).find(|&i| before[i].1 > 1.1).expect("right piece");
    // A feature on one piece's body moves that piece alone.
    bodies[right]
        .translate(&mut m, 0.0, 0.0, 1.0)
        .expect("lift the right piece");
    let lifted = |m: &mut Model, b: Body| {
        let brep = m.brep_of(b).expect("brep");
        brep.vertices
            .values()
            .map(|v| v.position.z)
            .fold(f64::MAX, f64::min)
    };
    for (i, &b) in bodies.iter().enumerate() {
        let want = if i == right { 1.0 } else { 0.0 };
        assert!((lifted(&mut m, b) - want).abs() < 1e-9, "body {i}");
    }
    m.set_param("x", 0.5).expect("move the tool");
    m.replay().expect("replay with the tool moved");
    let after = spans(&mut m);
    for (i, (b, a)) in before.iter().zip(&after).enumerate() {
        // Every piece moves at most by the tool's travel, and keeps the
        // side it was on.
        assert!((a.0 - b.0).abs() <= 0.5 + 1e-9, "body {i}: {b:?} -> {a:?}");
        assert!((a.1 - b.1).abs() <= 0.5 + 1e-9, "body {i}: {b:?} -> {a:?}");
    }
    assert!((after[right].0 - 0.7).abs() < 1e-9, "{:?}", after[right]);
    assert!(
        (lifted(&mut m, bodies[right]) - 1.0).abs() < 1e-9,
        "still lifted"
    );
    let (n, v) = measured_bodies(&mut m, &bodies);
    assert_eq!(n, 3);
    assert!((v - full).abs() < 1e-9, "{v}");

    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("partition_features");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join("split_bodies.pkn");
    m.save(&path).expect("save");
    let mut again = Model::open(&path).expect("open");
    again.replay().expect("replay");
    let reopened: Vec<_> = bodies
        .iter()
        .map(|&b| {
            let brep = again.brep_of(b).expect("brep");
            brep.vertices
                .values()
                .map(|v| v.position.x)
                .fold(f64::MAX, f64::min)
        })
        .collect();
    let now: Vec<_> = after.iter().map(|s| s.0).collect();
    assert_eq!(reopened, now, "the same pieces as the same bodies");
}

#[test]
fn an_imprint_splits_faces_and_keeps_the_volume() {
    let mut m = Model::new();
    let p = plate(&mut m);
    let faces = m.to_part().expect("part").face_count();
    let s = slab(&mut m, 0.0, 0.4);
    p.imprint(&mut m, s).expect("imprint");
    assert!(!s.is_alive(&m));
    let part = m.to_part().expect("part");
    assert!(part.face_count() > faces, "{} faces", part.face_count());
    let (_, v) = measured(&mut m);
    assert!((v - 2.4 * 1.6 * 0.4).abs() < 1e-9, "{v}");
}

#[test]
fn partition_features_reject_bad_tools() {
    let mut m = Model::new();
    let p = plate(&mut m);
    let s = slab(&mut m, 0.0, 0.4);
    assert!(matches!(
        p.split(&mut m, &[], PieceRule::All),
        Err(ScriptError::Argument { .. })
    ));
    assert!(matches!(
        p.split(&mut m, &[p], PieceRule::All),
        Err(ScriptError::Argument { .. })
    ));
    assert!(matches!(
        p.split(&mut m, &[s, s], PieceRule::All),
        Err(ScriptError::Argument { .. })
    ));
    assert!(matches!(
        p.nary_boolean(&mut m, &[s], SetExpr::Body(2), NonManifoldPolicy::Decline),
        Err(ScriptError::Argument { .. })
    ));
}

/// A fillet recorded by stable names on a split's kept piece resolves
/// again when the cutting tool moves: the piece's faces are named by
/// the entities that cut them, never by handles or piece order
/// (design 8, 10: moved-tool replay).
#[test]
fn a_moved_tool_replays_the_split_and_the_fillet_after_it() {
    let mut m = Model::new();
    m.param("x", 0.0, Unit::Meters).expect("param");
    let p = plate(&mut m);
    let s = m
        .tool(|m| {
            m.profile(Plane::xy_at(-0.5))
                .rect(0.4, 3.0)
                .at(popkern_script::scalar::p("x"), 0.0)
                .sketch("slab")?
                .extrude(m, End::Blind(2.0.into()))
        })
        .expect("slab");
    let p = p
        .split(&mut m, &[s], PieceRule::InsideTools(SetExpr::Body(0)))
        .expect("split")[0];
    p.fillet(&mut m, edges().straight().direction(Dir::Z), 0.05)
        .expect("fillet");
    let (_, before) = measured(&mut m);
    let lo = |m: &mut Model| {
        m.to_part()
            .expect("part")
            .brep()
            .vertices
            .values()
            .map(|v| v.position.x)
            .fold(f64::MAX, f64::min)
    };
    assert!((lo(&mut m) + 0.2).abs() < 1e-9);
    m.set_param("x", 0.5).expect("move the tool");
    m.replay().expect("replay with the tool moved");
    let (solids, after) = measured(&mut m);
    assert_eq!(solids, 1);
    assert!(
        (lo(&mut m) - 0.3).abs() < 1e-9,
        "the kept piece moved with the tool"
    );
    assert!((after - before).abs() < 1e-9, "{after} vs {before}");
}

/// Each face of an imprinted plate (design 5.4, PK-1062): its class
/// relative to the tool is queried by stable name, survives a reopened
/// document and follows the tool when it moves.
#[test]
fn imprint_classes_are_queried_by_stable_name() {
    let mut m = Model::new();
    m.param("x", 0.0, Unit::Meters).expect("param");
    let p = plate(&mut m);
    let s = m
        .tool(|m| {
            m.profile(Plane::xy_at(-0.5))
                .rect(0.4, 3.0)
                .at(popkern_script::scalar::p("x"), 0.0)
                .sketch("slab")?
                .extrude(m, End::Blind(2.0.into()))
        })
        .expect("slab");
    p.imprint(&mut m, s).expect("imprint");
    // Each face's class by name, and where its centre lies in x.
    let classes = |m: &mut Model| {
        let brep = m.brep_of(p).expect("brep");
        let mut out = Vec::new();
        for (handle, face) in &brep.faces {
            let name = brep
                .name_pool
                .get(face.name.expect("named"))
                .expect("pooled")
                .clone();
            let class = m.imprint_class(&name).expect("query");
            let xs: Vec<f64> = brep
                .face_outer_loop_positions(handle)
                .iter()
                .map(|q| q.x)
                .collect();
            let mid = xs.iter().sum::<f64>() / xs.len() as f64;
            out.push((mid, class));
        }
        out
    };
    let check = |at: f64, got: &[(f64, Option<FaceClass>)]| {
        assert!(got.iter().all(|(_, c)| c.is_some()), "{got:?}");
        for (mid, class) in got {
            let inside = (mid - at).abs() < 0.2 - 1e-9;
            if inside {
                assert_eq!(*class, Some(FaceClass::Inside), "x {mid}: {got:?}");
            }
        }
        assert!(got.iter().any(|(_, c)| *c == Some(FaceClass::Inside)));
        assert!(got.iter().any(|(_, c)| *c == Some(FaceClass::Outside)));
    };
    let first = classes(&mut m);
    check(0.0, &first);
    // A name no imprint classified.
    use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
    let other = EntityName {
        operation: OperationId(999),
        role: EntityRole {
            category: RoleCategory::Generated,
            index: 0,
        },
        ancestors: Vec::new(),
    };
    assert_eq!(m.imprint_class(&other).expect("query"), None);

    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("partition_features");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join("imprint.pkn");
    m.save(&path).expect("save");
    let mut again = Model::open(&path).expect("open");
    again.replay().expect("replay");
    let mut reopened = classes(&mut again);
    let mut before = first.clone();
    let key = |v: &mut Vec<(f64, Option<FaceClass>)>| {
        v.sort_by(|a, b| {
            a.0.total_cmp(&b.0)
                .then(format!("{:?}", a.1).cmp(&format!("{:?}", b.1)))
        })
    };
    key(&mut reopened);
    key(&mut before);
    assert_eq!(reopened, before);

    m.set_param("x", 0.6).expect("move the tool");
    m.replay().expect("replay");
    check(0.6, &classes(&mut m));
}

/// A fillet recorded on a split piece's edges follows the piece when
/// the tool moves either way and back (PK-1121): the faces and edges the
/// split cut are named by what they are pieces of and by their
/// neighbours (naming scheme 5), so every replay serves the fillet
/// with no rebind and no body matched by order.
#[test]
fn a_fillet_on_a_split_piece_follows_the_moving_tool() {
    let mut m = Model::new();
    m.param("x", 0.0, Unit::Meters).expect("param");
    let p = plate(&mut m);
    let s = m
        .tool(|m| {
            m.profile(Plane::xy_at(-0.5))
                .rect(0.4, 3.0)
                .at(popkern_script::scalar::p("x"), 0.0)
                .sketch("slab")?
                .extrude(m, End::Blind(2.0.into()))
        })
        .expect("slab");
    let bodies = p.split(&mut m, &[s], PieceRule::All).expect("split");
    let right = bodies[2];
    right
        .fillet(&mut m, edges().straight().direction(Dir::Z), 0.05)
        .expect("fillet on the right piece");
    let tol = popkern_core::tolerance::Tolerance::default();
    for x in [0.5, -0.3, 0.0] {
        m.set_param("x", x).expect("move the tool");
        let report = m.replay().expect("replay");
        for row in &report.features {
            assert!(
                matches!(
                    row.outcome,
                    popkern_history::replay_report::ReplayOutcome::Served
                ),
                "x = {x}: {} {:?}",
                row.label,
                row.outcome
            );
            assert!(
                row.fallbacks.is_empty(),
                "x = {x}: {} {:?}",
                row.label,
                row.fallbacks
            );
        }
        let brep = m.brep_of(right).expect("brep");
        assert!(brep.validate_all().is_empty(), "{:?}", brep.validate_all());
        let lo = brep
            .vertices
            .values()
            .map(|v| v.position.x)
            .fold(f64::MAX, f64::min);
        assert!(
            (lo - (x + 0.2)).abs() < 1e-9,
            "x = {x}: the right piece starts at {lo}"
        );
        // The piece less its four filleted vertical edges.
        let block = (1.2 - (x + 0.2)) * 1.6 * 0.4;
        let volume = body_volume(&brep, &tol).expect("volume");
        assert!(
            volume < block - 1e-4 && volume > block - 0.01,
            "x = {x}: {volume} vs {block}"
        );
    }
}
