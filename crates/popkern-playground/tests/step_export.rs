//! PK-627: STEP export is the whole model. A multi-body document (the
//! ball bearing: two races and eight balls) leaves as one solid per
//! body, re-imports with every face, and passes the OCC oracle where it
//! is installed; a single-body document exports exactly the bytes its
//! one body exported before.

use std::path::{Path, PathBuf};
use std::process::Command;

use popkern_playground::{HostOptions, Session};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn served(file: &str) -> Session {
    let source = std::fs::read_to_string(root().join("web/playground/examples").join(file))
        .expect("gallery script");
    let mut session = Session::new(HostOptions::default());
    let envelope = session.run(&source, file);
    assert!(
        envelope.outcome.is_served(),
        "{file}: {:?}",
        envelope.outcome
    );
    session
}

fn count(text: &str, entity: &str) -> usize {
    text.matches(&format!("={entity}(")).count()
}

#[test]
fn a_multi_body_model_exports_every_body_as_its_own_solid() {
    let mut session = served("ball_bearing.rhai");
    let model = session.model_mut().expect("model");
    let bodies = model.alive_bodies().len();
    assert_eq!(bodies, 10, "two races and eight balls");
    let bytes = session.export_step("bearing").expect("export");
    let text = String::from_utf8(bytes.clone()).expect("Part 21 is ASCII");
    assert_eq!(
        count(&text, "MANIFOLD_SOLID_BREP"),
        bodies,
        "one solid per body"
    );
    assert_eq!(count(&text, "CLOSED_SHELL"), bodies);

    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("pk627");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join("bearing.step");
    std::fs::write(&path, &bytes).expect("write");
    // Every face and a solid per body come back through the importer.
    let model = session.model_mut().expect("model");
    model
        .assert_step_roundtrip(dir.join("bearing_gate.step"))
        .expect("round-trip gate");

    let oracle = Path::new(&std::env::var("HOME").unwrap_or_default()).join("venvs/ocp/bin/python");
    let check = root().join("tools/occ-oracle/occ_check.py");
    if oracle.exists() && check.exists() {
        let output = Command::new(&oracle)
            .arg(&check)
            .arg(&path)
            .output()
            .expect("oracle");
        let report = String::from_utf8_lossy(&output.stdout);
        assert!(report.contains("\"passed\": true"), "OCC oracle: {report}");
    } else {
        eprintln!("[pk627] OCC oracle absent; checked by entity count and re-import only");
    }
}

#[test]
fn a_single_body_model_exports_exactly_its_body() {
    let mut session = served("boss_plate.rhai");
    let model = session.model_mut().expect("model");
    assert_eq!(model.alive_bodies().len(), 1);
    let body_bytes = model
        .to_part()
        .expect("part")
        .step_bytes("boss")
        .expect("part bytes");
    let model_bytes = session.export_step("boss").expect("export");
    assert_eq!(
        model_bytes, body_bytes,
        "byte-identical to the body's own export"
    );
}
