//! Phase-2 acceptance (design §10, spec §12): the boss-plate showcase
//! authored through the recorder is mesh-identical to the direct
//! `popkern-part` build, survives save → open → replay, and every
//! recording call is transactional.
use popkern_part::Plane as PartPlane;
use popkern_part::prelude::*;
use popkern_script::Plane;
use popkern_script::prelude::*;

fn direct() -> Part {
    let plate = Profile::on(PartPlane::XY)
        .rect(2.4, 1.6)
        .extrude(0.4)
        .expect("plate");
    let boss = Profile::on(PartPlane::xy_at(0.2))
        .circle(0.35)
        .extrude(0.8)
        .expect("boss");
    plate
        .union(&boss)
        .expect("union")
        .fillet(0.08, |e| {
            e.between(SurfaceKind::Cylinder, SurfaceKind::Plane)
                .at_height(0.4)
                .or(|e| e.straight().at_height(0.4))
                .or(|e| e.straight().direction(Dir::Z))
        })
        .expect("fillet")
}

fn recorded(m: &mut Model) -> Body {
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
            0.08,
            ContinuityKind::G1,
            "rounds",
        )
        .expect("fillet")
}

fn mesh_signature(part: &Part) -> (usize, usize, Vec<[i64; 3]>) {
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
    (mesh.vertices.len(), mesh.triangles.len(), verts)
}

#[test]
fn boss_plate_recorded_is_mesh_identical_to_the_direct_build() {
    let reference = direct();
    let mut m = Model::new();
    let started = std::time::Instant::now();
    let body = recorded(&mut m);
    let record_time = started.elapsed();
    assert_eq!(m.active().expect("active").id, body.id);
    let part = m.to_part().expect("part");
    assert_eq!(part.face_count(), reference.face_count(), "face count");
    assert_eq!(part.edge_count(), reference.edge_count(), "edge count");
    assert_eq!(
        mesh_signature(&part),
        mesh_signature(&reference),
        "mesh identical"
    );
    m.assert_watertight().expect("watertight");
    // Fallback firings are honest facts, not failures (spec §11); the
    // gate is no failed feature and no downgraded blend.
    assert_eq!(m.report().failures().count(), 0, "{:?}", m.report());
    assert_eq!(m.report().downgraded_blends().count(), 0);
    let rounds = m.feature("rounds").expect("label");
    let evidence = m.evidence(&rounds).expect("evidence");
    assert!(
        evidence.blend.is_some(),
        "fillet carries continuity evidence"
    );
    assert!(
        record_time.as_secs_f64() <= 2.0,
        "boss plate must record in <= 2 s (spec §9), took {record_time:?}"
    );
    eprintln!("[pk382] boss plate recorded in {record_time:?}");
}

#[test]
fn every_literal_became_an_editable_feature_parameter() {
    let mut m = Model::new();
    recorded(&mut m);
    let params = m.tree().parameters();
    let names: Vec<String> = params
        .all_params()
        .filter_map(|(id, _)| params.qualified_name(id))
        .collect();
    for expected in [
        "plate_sketch.width",
        "plate_sketch.height",
        "extrude-1.distance",
        "boss_sketch.radius",
        "extrude-2.distance",
        "rounds.radius",
    ] {
        assert!(
            names.iter().any(|n| n == expected),
            "missing {expected} in {names:?}"
        );
    }
    assert!(
        !names.iter().any(|n| n.starts_with("__pending")),
        "no pending names survive registration: {names:?}"
    );
}

#[test]
fn save_open_replay_round_trips() {
    let mut m = Model::new();
    recorded(&mut m);
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("boss_plate_parity");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join("boss_plate.pkn");
    m.save(&path).expect("save");
    let mut again = Model::open(&path).expect("open");
    let report = again.replay().expect("replay");
    assert_eq!(report.failures().count(), 0);
    assert_eq!(report.downgraded_blends().count(), 0);
    let a = mesh_signature(&m.to_part().expect("part"));
    let b = mesh_signature(&again.to_part().expect("part"));
    assert_eq!(a, b, "reloaded document replays to the same mesh");
    assert_eq!(again.features().len(), m.features().len());
    assert!(again.feature("rounds").is_ok());

    // Campaign sentence: reload → parameter edit → replay → gates.
    again.set_param("extrude-1.distance", 0.5).expect("edit");
    let report = again.replay().expect("replay after edit");
    assert_eq!(report.failures().count(), 0);
    let taller = again.to_part().expect("part");
    assert_eq!(taller.face_count(), m.to_part().expect("part").face_count());
    again.assert_watertight().expect("watertight after edit");
    let z_max = taller
        .brep()
        .vertices
        .values()
        .map(|v| v.position.z)
        .fold(f64::MIN, f64::max);
    assert!(
        (z_max - 1.0).abs() < 1e-9,
        "boss still tops out at 1.0, got {z_max}"
    );
    let plate_top = taller
        .brep()
        .vertices
        .values()
        .filter(|v| (v.position.x.abs() - 1.2).abs() < 1e-9)
        .map(|v| v.position.z)
        .fold(f64::MIN, f64::max);
    // The plate's outer corners top out at the fillet rail: 0.5 - 0.08.
    assert!(
        (plate_top - 0.42).abs() < 1e-9,
        "plate top moved to 0.5 (rail at 0.42), got {plate_top}"
    );
}

#[test]
fn a_declined_feature_leaves_the_document_untouched() {
    let mut m = Model::new();
    let cube = m
        .profile(Plane::XY)
        .rect(1.0, 1.0)
        .sketch("cube")
        .expect("sketch")
        .extrude(&mut m, End::Blind(1.0.into()))
        .expect("extrude");
    let before = m.features().len();
    let params_before = m.tree().parameters().all_params().count();
    let err = cube
        .fillet(&mut m, edges(), 0.9)
        .expect_err("r = 0.9 on every unit-cube edge cannot construct");
    match &err {
        ScriptError::FeatureFailed { operation, .. } => assert_eq!(operation, "fillet"),
        other => panic!("expected FeatureFailed, got {other}"),
    }
    assert_eq!(m.features().len(), before, "failed feature removed");
    assert_eq!(
        m.tree().parameters().all_params().count(),
        params_before,
        "its parameters removed"
    );
    assert_eq!(m.to_part().expect("part").face_count(), 6);
}

#[test]
fn empty_selection_is_a_structured_error() {
    let mut m = Model::new();
    let cube = m
        .profile(Plane::XY)
        .rect(1.0, 1.0)
        .sketch("cube")
        .expect("sketch")
        .extrude(&mut m, End::Blind(1.0.into()))
        .expect("extrude");
    let err = cube
        .fillet(&mut m, edges().circular(), 0.1)
        .expect_err("no circular edges on a cube");
    match err {
        ScriptError::EmptySelection {
            selection,
            available,
            ..
        } => {
            assert_eq!(selection, "edges().circular()");
            assert_eq!(available, 12);
        }
        other => panic!("{other}"),
    }
}
