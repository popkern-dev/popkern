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

#[test]
fn an_assembly_timeline_and_its_steps_answer_over_the_protocol() {
    let mut p = init();
    let source = "let a = assembly(\"pair\");\nlet b = a.part(\"b\", |m| { m.profile(xy()).rect(1.0, 1.0).sketch(\"s\").extrude(blind(0.1)); });\na.place(b, \"one\");\na.place(b, \"two\").translate(2.0, 0.0, 0.0);\na\n";
    let run = p.handle(&json!({ "id": 2, "type": "run", "source": source, "name": "pair.rhai" }));
    assert_eq!(run.json["envelope"]["outcome"], "served");
    let timeline = p.handle(&json!({ "id": 3, "type": "timeline" }));
    let steps = timeline.json["assembly"]["steps"]
        .as_array()
        .expect("steps")
        .clone();
    assert_eq!(steps.len(), 2);
    assert_eq!(steps[1]["kind"], "place");
    assert_eq!(timeline.json["assembly"]["recorded"], true);
    let one = p.handle(&json!({ "id": 4, "type": "assembly_step", "index": 1 }));
    assert_eq!(one.json["instances"].as_array().map(Vec::len), Some(1));
    let both = p.handle(&json!({ "id": 5, "type": "assembly_step", "index": 2 }));
    assert_eq!(
        both.json["instances"],
        run.json["envelope"]["assembly"]["instances"]
    );
    // A part document has no assembly timeline and refuses steps.
    p.handle(&json!({ "id": 6, "type": "run", "source": PLATE, "name": "plate.rhai" }));
    let part = p.handle(&json!({ "id": 7, "type": "timeline" }));
    assert_eq!(part.json["assembly"], Value::Null);
    let refused = p.handle(&json!({ "id": 8, "type": "assembly_step", "index": 1 }));
    assert_eq!(refused.json["ok"], false);
}

#[test]
fn a_pose_of_joints_answers_one_placement_envelope() {
    let mut p = init();
    let source = "let a = assembly(\"chain\");\n\
        let b = a.part(\"b\", |m| { m.profile(xy()).rect(1.0, 0.2).sketch(\"s\").extrude(blind(0.1)); });\n\
        let base = a.place(b, \"base\");\n\
        let arm = a.place(b, \"arm\");\n\
        let tip = a.place(b, \"tip\");\n\
        a.ground(base);\n\
        a.revolute(\"j1\", base.axis([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), arm.axis([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]),\n\
            #{ zero: [base.datum(\"x\"), arm.datum(\"x\")], angle: 0.0, limits: [-360.0, 360.0] });\n\
        a.revolute(\"j2\", arm.axis([1.0, 0.0, 0.0], [0.0, 0.0, 1.0]), tip.axis([1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),\n\
            #{ zero: [arm.datum(\"x\"), tip.datum(\"x\")], angle: 0.0, limits: [-360.0, 360.0] });\n\
        a\n";
    let run = p.handle(&json!({ "id": 2, "type": "run", "source": source, "name": "chain.rhai" }));
    assert_eq!(run.json["envelope"]["outcome"], "served", "{}", run.json);
    let posed =
        p.handle(&json!({ "id": 3, "type": "set_joints", "values": { "j1": 30.0, "j2": 45.0 } }));
    let envelope = &posed.json["envelope"];
    assert_eq!(envelope["outcome"], "served", "{}", posed.json);
    assert_eq!(envelope["bodies_unchanged"], true);
    let value = |name: &str| {
        envelope["params"]
            .as_array()
            .and_then(|rows| rows.iter().find(|r| r["name"] == name))
            .and_then(|r| r["value"].as_f64())
            .unwrap_or(f64::NAN)
    };
    // Degrees held in SI (radians) and shown back.
    assert!((value("j1") - 30.0).abs() < 1e-9 && (value("j2") - 45.0).abs() < 1e-9);
    // The same pose one joint at a time places every occurrence alike.
    p.handle(&json!({ "id": 4, "type": "set_joints", "values": { "j1": 0.0, "j2": 0.0 } }));
    p.handle(&json!({ "id": 5, "type": "set_param", "name": "j1", "value": 30.0 }));
    let single = p.handle(&json!({ "id": 6, "type": "set_param", "name": "j2", "value": 45.0 }));
    let worlds = |v: &Value| -> Vec<f64> {
        v["assembly"]["instances"]
            .as_array()
            .map(|a| {
                a.iter()
                    .flat_map(|i| i["world"].as_array().cloned().unwrap_or_default())
                    .filter_map(|x| x.as_f64())
                    .collect()
            })
            .unwrap_or_default()
    };
    let (a, b) = (worlds(envelope), worlds(&single.json["envelope"]));
    assert!(!a.is_empty() && a.len() == b.len());
    assert!(
        a.iter().zip(&b).all(|(x, y)| (x - y).abs() < 1e-5),
        "{a:?} vs {b:?}"
    );
    // A value outside a joint's limits refuses the whole pose.
    let refused =
        p.handle(&json!({ "id": 7, "type": "set_joints", "values": { "j1": 10.0, "j2": 400.0 } }));
    assert_ne!(
        refused.json["envelope"]["outcome"], "served",
        "{}",
        refused.json
    );
    let after = p.handle(&json!({ "id": 8, "type": "set_param", "name": "j1", "value": 30.0 }));
    let j2 = after.json["envelope"]["params"]
        .as_array()
        .and_then(|rows| rows.iter().find(|r| r["name"] == "j2"))
        .and_then(|r| r["value"].as_f64())
        .unwrap_or(f64::NAN);
    assert!(
        (j2 - 45.0).abs() < 1e-9,
        "a refused pose kept part of itself: j2 {j2}"
    );
}
