//! Phase 2e (PK-385): tool models are bodies of the same tree.
//! Booleans take bodies by handle or label, the receiver is the body
//! the op is called on (not always the first body), a consumed tool is
//! no longer alive, and body labels survive save → open.
use popkern_history::SnapshotOrigin;
use popkern_part::prelude::*;
use popkern_script::Plane;
use popkern_script::prelude::*;

fn plate(m: &mut Model) -> Body {
    m.profile(Plane::XY)
        .rect(2.4, 1.6)
        .sketch("plate_sketch")
        .expect("sketch")
        .extrude(m, End::Blind(0.4.into()))
        .expect("extrude")
}

fn signature(part: &Part) -> (usize, usize, i64) {
    let mesh = part
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut v6 = 0.0;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        v6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    (
        part.face_count(),
        part.edge_count(),
        (v6 / 6.0 * 1e9).round() as i64,
    )
}

#[test]
fn boss_plate_tool_bodies_by_label_match_the_direct_build() {
    let direct = Profile::on(popkern_part::Plane::XY)
        .rect(2.4, 1.6)
        .extrude(0.4)
        .expect("plate")
        .union(
            &Profile::on(popkern_part::Plane::xy_at(0.2))
                .circle(0.35)
                .extrude(0.8)
                .expect("boss"),
        )
        .expect("union")
        .subtract(
            &Profile::on(popkern_part::Plane::xy_at(-0.5))
                .circle(0.15)
                .extrude(2.0)
                .expect("bore"),
        )
        .expect("subtract");

    let mut m = Model::new();
    let plate = plate(&mut m);
    let boss = m
        .tool(|m| {
            m.profile(Plane::xy_at(0.2))
                .circle(0.35)
                .sketch("boss_sketch")?
                .extrude(m, End::Blind(0.8.into()))
        })
        .expect("boss");
    m.label_body("boss", boss).expect("label");
    let bore = m
        .tool(|m| {
            m.profile(Plane::xy_at(-0.5))
                .circle(0.15)
                .sketch("bore_sketch")?
                .extrude(m, End::Blind(2.0.into()))
        })
        .expect("bore");
    m.label_body("bore", bore).expect("label");
    assert_eq!(m.body_label(boss).as_deref(), Some("boss"));
    let boss_by_label = m.body("boss").expect("lookup");
    assert_eq!(boss_by_label, boss);
    let plate = plate.union(&mut m, boss_by_label).expect("union");
    let bore_by_label = m.body("bore").expect("lookup");
    let plate = plate.subtract(&mut m, bore_by_label).expect("subtract");
    assert_eq!(m.active().expect("active"), plate);
    assert_eq!(signature(&m.to_part().expect("part")), signature(&direct));

    // Consumed tools are not alive; the receiver is.
    assert!(!boss.is_alive(&m));
    assert!(!bore.is_alive(&m));
    assert!(plate.is_alive(&m));
    let alive: Vec<Body> = m.alive_bodies();
    assert_eq!(alive, vec![plate]);
    assert!(matches!(
        boss.fillet(&mut m, edges(), 0.05),
        Err(ScriptError::BodyNotAlive(_))
    ));

    // Labels survive save → open, including labels of consumed bodies
    // (they name history, not live geometry).
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("body_bookkeeping");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join("boss_plate.pkn");
    m.save(&path).expect("save");
    let mut again = Model::open(&path).expect("open");
    again.replay().expect("replay");
    assert_eq!(again.body("boss").expect("label survives"), boss);
    assert_eq!(again.body("bore").expect("label survives"), bore);
    assert_eq!(
        signature(&again.to_part().expect("part")),
        signature(&direct)
    );
}

#[test]
fn the_receiver_is_the_body_the_op_is_called_on() {
    // Three bodies; union the second INTO the third, leaving the first alone.
    let mut m = Model::new();
    let a = plate(&mut m);
    let b = m
        .tool(|m| {
            m.profile(Plane::xy_at(3.0))
                .rect(1.0, 1.0)
                .sketch("b")?
                .extrude(m, End::Blind(1.0.into()))
        })
        .expect("b");
    let c = m
        .tool(|m| {
            m.profile(Plane::xy_at(3.5))
                .rect(1.0, 1.0)
                .sketch("c")?
                .extrude(m, End::Blind(1.0.into()))
        })
        .expect("c");
    let merged = c.union(&mut m, b).expect("union into c");
    assert_eq!(merged, c, "the receiver keeps its identity");
    assert!(a.is_alive(&m) && c.is_alive(&m) && !b.is_alive(&m));
    let mut alive = m.alive_bodies();
    alive.sort_by_key(|b| b.id.0);
    assert_eq!(alive, vec![a, c]);
    m.set_active(c).expect("activate");
    let part = m.to_part().expect("part");
    assert_eq!(part.face_count(), 6, "c ∪ b is one 1×1×1.5 box");
    m.set_active(a).expect("activate");
    assert_eq!(m.to_part().expect("part").face_count(), 6);
    assert!(
        m.set_active(b).is_err(),
        "a consumed body cannot be activated"
    );
}

#[test]
fn raw_bodies_are_bodies() {
    let mut m = Model::new();
    let plate = plate(&mut m);
    let cube =
        popkern_core::topology::BRep::unit_cube(popkern_core::tolerance::Tolerance::default())
            .expect("cube");
    let tool = m
        .raw_body(
            cube,
            SnapshotOrigin::Import {
                source: "unit_cube".into(),
            },
        )
        .expect("raw body");
    m.label_body("cutter", tool).expect("label");
    let cut = plate.subtract(&mut m, tool).expect("subtract");
    assert_eq!(cut, plate);
    let part = m.to_part().expect("part");
    assert!(
        part.face_count() > 6,
        "a corner notch adds faces: {}",
        part.face_count()
    );
    assert!(!tool.is_alive(&m));
    assert_eq!(m.body("cutter").expect("label"), tool);
}
