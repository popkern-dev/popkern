//! Phase 2d (PK-384): migration verbs are the canonical code form of
//! history-tree edits. Every verb is checked by REPLAY EQUIVALENCE —
//! the migrated document replays to the same mesh as a document
//! recorded directly in the target shape — and by transactionality.
use popkern_part::prelude::*;
use popkern_script::Plane;
use popkern_script::prelude::*;

fn signature(part: &Part) -> (usize, usize, Vec<[i64; 3]>) {
    let mesh = part
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut verts: Vec<[i64; 3]> = mesh
        .vertices
        .iter()
        .map(|v| {
            let p = v.position;
            [
                (p.x * 1e9).round() as i64,
                (p.y * 1e9).round() as i64,
                (p.z * 1e9).round() as i64,
            ]
        })
        .collect();
    verts.sort();
    (part.face_count(), part.edge_count(), verts)
}

/// The boss plate, with the fillet's radius/continuity as arguments.
fn boss_plate(m: &mut Model, radius: f64, continuity: ContinuityKind) -> Body {
    let plate = m
        .profile(Plane::XY)
        .rect(2.4, 1.6)
        .sketch("plate_sketch")
        .expect("sketch")
        .extrude(m, End::Blind(0.4.into()))
        .expect("extrude");
    let boss = m
        .tool(|m| {
            m.profile(Plane::xy_at(0.2))
                .circle(0.35)
                .sketch("boss_sketch")?
                .extrude(m, End::Blind(0.8.into()))
        })
        .expect("boss");
    let plate = plate.union(m, boss).expect("union");
    plate
        .fillet_continuity(
            m,
            edges()
                .between(SurfaceKind::Cylinder, SurfaceKind::Plane)
                .at_height(0.4)
                .or(edges().straight().at_height(0.4))
                .or(edges().straight().direction(Dir::Z)),
            radius,
            continuity,
            "rounds",
        )
        .expect("fillet")
}

fn saved_boss_plate(name: &str) -> std::path::PathBuf {
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("migration_verbs");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join(format!("{name}.pkn"));
    let mut m = Model::new();
    boss_plate(&mut m, 0.08, ContinuityKind::G1);
    m.save(&path).expect("save");
    path
}

#[test]
fn modify_fillet_radius_and_continuity_matches_a_direct_recording() {
    let path = saved_boss_plate("modify");
    let mut m = Model::open(&path).expect("open");
    let rounds = m.feature("rounds").expect("label");
    assert_eq!(
        m.describe(&rounds).expect("describe"),
        "rounds — fillet 0.08 on 9 edges: edges().between(cylinder, plane).at_height(0.4).or(edges().straight().at_height(0.4)).or(edges().straight().direction(z))"
    );
    // (G2 on this selection declines in the kernel at every radius —
    // PK-389; the continuity edit is exercised on a cube below.)
    m.modify(&rounds, Edit::new().radius(0.1)).expect("modify");
    let report = m.replay().expect("replay");
    assert_eq!(report.failures().count(), 0);
    m.assert_watertight().expect("watertight");
    assert!(
        m.describe(&rounds)
            .expect("describe")
            .starts_with("rounds — fillet 0.1 on 9 edges:")
    );

    let mut direct = Model::new();
    boss_plate(&mut direct, 0.1, ContinuityKind::G1);
    assert_eq!(
        signature(&m.to_part().expect("part")),
        signature(&direct.to_part().expect("part")),
        "migrated document replays like a direct recording"
    );
    m.save(&path).expect("save");
    let mut again = Model::open(&path).expect("reopen");
    again.replay().expect("replay");
    assert_eq!(
        signature(&again.to_part().expect("part")),
        signature(&direct.to_part().expect("part"))
    );
}

#[test]
fn modify_is_transactional_and_typed() {
    let path = saved_boss_plate("modify_txn");
    let mut m = Model::open(&path).expect("open");
    let rounds = m.feature("rounds").expect("label");
    let before = signature(&m.to_part().expect("part"));
    // A radius the kernel cannot serve: the edit is rolled back.
    let err = m
        .modify(&rounds, Edit::new().radius(0.9))
        .expect_err("r = 0.9 cannot serve");
    assert!(matches!(err, ScriptError::FeatureFailed { .. }), "{err}");
    assert_eq!(signature(&m.to_part().expect("part")), before);
    assert!(
        m.params()
            .iter()
            .any(|p| p.name == "rounds.radius" && (p.value - 0.08).abs() < 1e-12),
        "radius restored: {:?}",
        m.params()
    );
    // A role the feature does not have is an argument error, not a kernel run.
    let extrude = m.feature("extrude-1").expect("label");
    assert!(matches!(
        m.modify(&extrude, Edit::new().continuity(ContinuityKind::G2)),
        Err(ScriptError::Argument { .. })
    ));
    // Re-selecting edges resolves against the body the feature sees.
    m.modify(
        &rounds,
        Edit::new().edges(edges().straight().direction(Dir::Z)),
    )
    .expect("re-select");
    assert!(m.describe(&rounds).expect("d").contains("on 4 edges"));
}

#[test]
fn insert_reorder_suppress_delete_are_replay_equivalent() {
    let path = saved_boss_plate("insert");
    let mut m = Model::open(&path).expect("open");
    let rounds = m.feature("rounds").expect("label");
    // Insert a boss-top chamfer before the fillet (its edge is not one
    // the fillet touches; a bottom-rim chamfer would make the vertical
    // fillets run out against bevels — a registered decline class).
    let chamfer = m
        .insert_before(&rounds, |m, body| {
            body.chamfer(m, edges().circular().at_height(1.0), 0.05)
        })
        .expect("insert");
    assert_eq!(chamfer.operation, "chamfer");
    let labels: Vec<String> = m.features().iter().map(|f| f.label.clone()).collect();
    assert_eq!(labels[labels.len() - 2], chamfer.label);
    assert_eq!(labels[labels.len() - 1], "rounds");
    let mut direct = Model::new();
    {
        let plate = direct
            .profile(Plane::XY)
            .rect(2.4, 1.6)
            .sketch("plate_sketch")
            .expect("s")
            .extrude(&mut direct, End::Blind(0.4.into()))
            .expect("e");
        let boss = direct
            .tool(|m| {
                m.profile(Plane::xy_at(0.2))
                    .circle(0.35)
                    .sketch("boss_sketch")?
                    .extrude(m, End::Blind(0.8.into()))
            })
            .expect("boss");
        let plate = plate.union(&mut direct, boss).expect("union");
        let plate = plate
            .chamfer(&mut direct, edges().circular().at_height(1.0), 0.05)
            .expect("chamfer");
        plate
            .fillet_continuity(
                &mut direct,
                edges()
                    .between(SurfaceKind::Cylinder, SurfaceKind::Plane)
                    .at_height(0.4)
                    .or(edges().straight().at_height(0.4))
                    .or(edges().straight().direction(Dir::Z)),
                0.08,
                ContinuityKind::G1,
                "rounds",
            )
            .expect("fillet");
    }
    assert_eq!(
        signature(&m.to_part().expect("part")),
        signature(&direct.to_part().expect("part")),
        "insert_before ≡ recording in that order"
    );
    // Suppress the chamfer: equals the original document.
    m.suppress(&chamfer).expect("suppress");
    let mut original = Model::open(&path).expect("open");
    assert_eq!(
        signature(&m.to_part().expect("part")),
        signature(&original.to_part().expect("part")),
        "suppressed ≡ absent"
    );
    m.unsuppress(&chamfer).expect("unsuppress");
    // Reorder the chamfer after the fillet (a fillet then a chamfer).
    let last = m.features().len() - 1;
    m.reorder(&chamfer, last).expect("reorder");
    let labels: Vec<String> = m.features().iter().map(|f| f.label.clone()).collect();
    assert_eq!(labels[last], chamfer.label);
    assert_eq!(m.report().failures().count(), 0);
    // Delete it: back to the original.
    let deleted = m.delete(&chamfer).expect("delete");
    assert_eq!(deleted.len(), 1);
    assert_eq!(
        signature(&m.to_part().expect("part")),
        signature(&original.to_part().expect("part")),
        "deleted ≡ original"
    );
    assert!(m.feature(&chamfer.label).is_err());
    // Dependents: the fillet depends on the union, which depends on both extrudes.
    let union = m.feature("boolean-1").expect("union");
    let deps: Vec<String> = m
        .dependents(&union)
        .iter()
        .map(|f| f.label.clone())
        .collect();
    assert_eq!(deps, vec!["rounds"]);
}

#[test]
fn undo_redo_and_rollback_round_trip() {
    let path = saved_boss_plate("undo");
    let mut m = Model::open(&path).expect("open");
    let before = signature(&m.to_part().expect("part"));
    let rounds = m.feature("rounds").expect("label");
    m.suppress(&rounds).expect("suppress");
    let suppressed = signature(&m.to_part().expect("part"));
    assert_ne!(before, suppressed);
    assert!(m.undo().expect("undo"));
    assert_eq!(signature(&m.to_part().expect("part")), before);
    assert!(m.redo().expect("redo"));
    assert_eq!(signature(&m.to_part().expect("part")), suppressed);
    m.undo().expect("undo");
    // Rollback to just before the fillet: the union body is served.
    let report = m.rollback_to(&rounds).expect("rollback");
    assert_eq!(report.failures().count(), 0);
    assert_eq!(m.to_part().expect("part").face_count(), 8);
    m.rollback_to_end().expect("end");
    assert_eq!(signature(&m.to_part().expect("part")), before);
}

#[test]
fn parameters_and_configurations_surface() {
    let mut m = Model::new();
    m.param("wall", 0.4, Unit::Meters).expect("param");
    let plate = m
        .profile(Plane::XY)
        .rect(2.4, 1.6)
        .sketch("plate")
        .expect("s")
        .extrude(&mut m, End::Blind(popkern_script::scalar::p("wall")))
        .expect("e");
    let _ = plate;
    let z_max = |m: &mut Model| {
        m.to_part()
            .expect("part")
            .brep()
            .vertices
            .values()
            .map(|v| v.position.z)
            .fold(f64::MIN, f64::max)
    };
    assert!((z_max(&mut m) - 0.4).abs() < 1e-12);
    m.set_param("wall", 0.6).expect("edit");
    m.replay().expect("replay");
    assert!((z_max(&mut m) - 0.6).abs() < 1e-12);
    let names: Vec<String> = m.params().iter().map(|p| p.name.clone()).collect();
    assert!(names.contains(&"wall".to_string()));
    let distance = m
        .params()
        .into_iter()
        .find(|p| p.name == "extrude-1.distance")
        .expect("feature param");
    assert_eq!(distance.expression.as_deref(), Some("wall"));
}

#[test]
fn modify_continuity_matches_a_direct_recording() {
    let cube = |m: &mut Model, continuity: ContinuityKind| {
        m.profile(Plane::XY)
            .rect(1.0, 1.0)
            .sketch("cube")
            .expect("s")
            .extrude(m, End::Blind(1.0.into()))
            .expect("e")
            .fillet_continuity(
                m,
                edges().straight().direction(Dir::Z),
                0.1,
                continuity,
                "posts",
            )
            .expect("fillet")
    };
    let mut m = Model::new();
    cube(&mut m, ContinuityKind::G1);
    let posts = m.feature("posts").expect("label");
    m.modify(&posts, Edit::new().continuity(ContinuityKind::G2))
        .expect("modify continuity");
    assert!(m.describe(&posts).expect("d").contains("(G2)"));
    let evidence = m.evidence(&posts).expect("evidence");
    let blend = evidence.blend.as_ref().expect("blend evidence");
    assert!(!blend.downgraded(), "G2 served, not downgraded");
    let mut direct = Model::new();
    cube(&mut direct, ContinuityKind::G2);
    assert_eq!(
        signature(&m.to_part().expect("part")),
        signature(&direct.to_part().expect("part"))
    );
}

#[test]
fn configurations_override_parameters_and_suppress_features() {
    let mut m = Model::new();
    m.param("size", 1.0, Unit::Meters).expect("param");
    m.param_expr("height", "size * 0.75", Unit::Meters)
        .expect("expr");
    let cube = m
        .profile(Plane::XY)
        .rect(
            popkern_script::scalar::p("size"),
            popkern_script::scalar::p("size"),
        )
        .sketch("base")
        .expect("s")
        .extrude(&mut m, End::Blind(popkern_script::scalar::p("height")))
        .expect("e");
    cube.fillet(&mut m, edges().straight().direction(Dir::Z), 0.1)
        .expect("fillet");
    let rounds = m.feature("fillet-1").expect("label");
    let z_max = |m: &mut Model| {
        m.to_part()
            .expect("part")
            .brep()
            .vertices
            .values()
            .map(|v| v.position.z)
            .fold(f64::MIN, f64::max)
    };
    let z = z_max(&mut m);
    assert!((z - 0.75).abs() < 1e-12, "z {z}; params {:?}", m.params());
    m.add_config("tall", "double size").expect("config");
    m.set_override("tall", "size", 2.0).expect("override");
    m.set_suppression("tall", &rounds, true).expect("suppress");
    m.set_active_config(Some("tall")).expect("activate");
    assert_eq!(m.active_config().as_deref(), Some("tall"));
    assert!((z_max(&mut m) - 1.5).abs() < 1e-12, "height follows size");
    assert_eq!(
        m.to_part().expect("part").face_count(),
        6,
        "fillet suppressed in tall"
    );
    m.set_active_config(None).expect("base");
    let z = z_max(&mut m);
    assert!(
        (z - 0.75).abs() < 1e-12,
        "base z {z}; params {:?}",
        m.params()
    );
    assert_eq!(m.to_part().expect("part").face_count(), 10);
    assert_eq!(m.configs(), vec!["tall".to_string()]);
}
