// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use serde_json::{Value, json};

use super::{Buffer, Protocol};

const PLATE: &str = "let m = model();\nlet p = m.profile(xy()).rect(1.0, 1.0).sketch(\"plate\").extrude(blind(0.2));\nlet b = m.tool(|m| m.profile(xy()).circle(0.2).sketch(\"pin\").extrude(blind(0.5)));\np.union(b);\nm\n";

fn init() -> Protocol {
    let mut p = Protocol::new();
    let reply = p.handle(&json!({ "id": 1, "type": "init", "options": {} }));
    assert_eq!(reply.json["ok"], true);
    assert_eq!(reply.json["protocol_version"], crate::PROTOCOL_VERSION);
    p
}

#[test]
fn requests_before_init_and_unknown_types_have_the_one_error_shape() {
    let mut p = Protocol::new();
    let early = p.handle(&json!({ "id": 7, "type": "run", "source": "" }));
    assert_eq!(early.json["outcome"]["outcome"], "gate_failure");
    assert_eq!(early.json["id"], 7);
    let mut p = init();
    let bogus = p.handle(&json!({ "id": 8, "type": "frobnicate" }));
    assert_eq!(bogus.json["ok"], false);
    assert_eq!(bogus.json["outcome"]["step"], "worker");
}

#[test]
fn a_run_answers_the_envelope_with_every_body_buffer_attached() {
    let mut p = init();
    let reply = p.handle(&json!({ "id": 2, "type": "run", "source": PLATE, "name": "plate.rhai" }));
    assert_eq!(reply.json["ok"], true);
    assert_eq!(reply.json["envelope"]["outcome"], "served");
    let bodies = reply.json["envelope"]["bodies"]
        .as_array()
        .expect("bodies")
        .len();
    assert_eq!(
        reply.attachments.len(),
        bodies * 6,
        "positions, normals, indices, face_ids, gaussian, mean"
    );
    assert_eq!(reply.attachments[0].path, "buffers.0.positions");
    assert!(matches!(reply.attachments[2].buffer, Buffer::U32(_)));
    let stats = p.handle(&json!({ "id": 3, "type": "stats" }));
    assert_eq!(stats.json["body_count"], bodies);
}

#[test]
fn timeline_and_step_meshes_answer_like_the_worker() {
    let mut p = init();
    p.handle(&json!({ "id": 2, "type": "run", "source": PLATE, "name": "plate.rhai" }));
    let timeline = p.handle(&json!({ "id": 3, "type": "timeline" }));
    let steps = timeline.json["parts"][0]["steps"]
        .as_array()
        .expect("steps")
        .len();
    assert_eq!(steps, 5);
    let first =
        p.handle(&json!({ "id": 4, "type": "step_meshes", "part": "", "index": 2, "known": [] }));
    let key = first.json["rows"][0]["key"].clone();
    assert!(
        first
            .attachments
            .iter()
            .any(|a| a.path == "rows.0.positions")
    );
    let second = p
        .handle(&json!({ "id": 5, "type": "step_meshes", "part": "", "index": 4, "known": [key] }));
    let rows = second.json["rows"].as_array().expect("rows");
    let reused = rows
        .iter()
        .position(|r| r["key"] == first.json["rows"][0]["key"])
        .expect("plate kept");
    assert_eq!(rows[reused]["buffers"], false);
    assert!(
        !second
            .attachments
            .iter()
            .any(|a| a.path.starts_with(&format!("rows.{reused}.")))
    );
}

#[test]
fn exports_attach_bytes_and_failures_answer_outcomes() {
    let mut p = init();
    let none = p.handle(&json!({ "id": 2, "type": "export_pkn" }));
    assert_eq!(none.json["ok"], false);
    p.handle(&json!({ "id": 3, "type": "run", "source": PLATE }));
    let pkn = p.handle(&json!({ "id": 4, "type": "export_pkn" }));
    assert_eq!(pkn.json["ok"], true);
    assert_eq!(pkn.attachments[0].path, "bytes");
    let bytes = pkn.attachments[0].buffer.bytes();
    let doc: Value = serde_json::from_slice(&bytes).expect("pkn json");
    assert!(doc.get("tree").is_some());
}
