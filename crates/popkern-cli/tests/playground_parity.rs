//! PK-410: the browser engine and `popkern run` are two hosts of one
//! class partition — for every §10 class the PK-403 parity fixture
//! reaches, the playground's envelope outcome equals the command's
//! summary outcome field for field (class, position, payload,
//! remediation text).

use std::path::{Path, PathBuf};

use popkern_cli::{execute, parse_args};
use popkern_playground::{HostOptions, Session};
use serde_json::Value;

const CUBE: &str = r#"let m = model();
let cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
"#;

/// The ten scenarios of `popkern-rhai/tests/error_parity.rs`.
const SCENARIOS: &[(&str, &str, &str)] = &[
    (
        "oversized_fillet",
        "FeatureFailed",
        "cube.fillet(edges(), 0.9);",
    ),
    (
        "negative_radius",
        "FeatureFailed",
        "cube.fillet(edges(), -0.1);",
    ),
    (
        "empty_selection",
        "EmptySelection",
        "cube.fillet(edges().circular(), 0.1);",
    ),
    (
        "cardinality",
        "SelectionCardinality",
        "m.label_face(\"top\", faces().kind(\"plane\"));",
    ),
    (
        "unknown_parameter",
        "ArgumentError",
        "m.set_param(\"nope\", 1.0);",
    ),
    ("label_not_found", "LabelError", "m.feature(\"nope\");"),
    (
        "duplicate_label",
        "LabelError",
        "m.rename(m.feature(\"cube\"), \"extrude-1\");",
    ),
    (
        "consumed_tool",
        "TimelineError",
        "let t = m.tool(|m| m.profile(xy_at(0.5)).circle(0.2).sketch(\"t\").extrude(blind(1.0)));\ncube = cube.subtract(t);\ncube.subtract(t);",
    ),
    (
        "missing_document",
        "IoFailed",
        "open(\"/nonexistent/dir/model.pkn\");",
    ),
    (
        "unknown_expression_parameter",
        "ArgumentError",
        "m.param_expr(\"h\", \"wall * 2\");",
    ),
];

fn write(name: &str, tail: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("playground_parity");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join(format!("{name}.rhai"));
    std::fs::write(&path, format!("{CUBE}{tail}\n")).expect("write");
    path
}

fn strip_report(mut value: Value) -> Value {
    if let Some(payload) = value.get_mut("payload").and_then(Value::as_object_mut) {
        payload.remove("report");
    }
    value
}

#[test]
fn every_reachable_class_lands_in_the_same_partition_with_the_same_payload() {
    for (name, kind, tail) in SCENARIOS {
        let path = write(name, tail);
        let file = path
            .file_name()
            .expect("name")
            .to_string_lossy()
            .into_owned();
        let args = vec!["run".to_string(), path.to_string_lossy().into_owned()];
        let cli = parse_args(&args).expect("args");
        let summary = execute(&cli);
        let mut session = Session::new(HostOptions::default());
        let envelope = session.run(&std::fs::read_to_string(&path).expect("read"), &file);
        let command = strip_report(serde_json::to_value(&summary.outcome).expect("json"));
        let page = strip_report(serde_json::to_value(&envelope.outcome).expect("json"));
        assert_eq!(command, page, "[{name}] outcome parity");
        assert_eq!(summary.outcome.exit_code(), envelope.outcome.exit_code());
        let class = command["kind"]
            .as_str()
            .or_else(|| command["step"].as_str())
            .unwrap_or("");
        assert_eq!(class, *kind, "[{name}] class");
        assert!(
            command["rendered"]
                .as_str()
                .or_else(|| command["detail"].as_str())
                .is_some_and(|r| r.contains(&file)),
            "[{name}] positioned at the script: {command}"
        );
    }
}
