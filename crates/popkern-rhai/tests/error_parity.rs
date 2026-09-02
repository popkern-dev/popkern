//! PK-403: the Rhai error map is the same §10 payload the Rust skin
//! projects (`ScriptError::payload`), field for field and float for
//! float, for every class the recorder fixtures reach; the raising
//! call's line and column are stamped; `InvalidParameter.parameter`
//! is the surface name with the kernel path beside it; uncaught
//! errors render with their remediation.

use popkern_rhai::{ScriptHost, remediation_lines};
use popkern_script::prelude::*;
use popkern_script::{Plane, ScriptError};
use serde_json::Value;

const CUBE: &str = r#"let m = model();
let cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
"#;

fn cube(m: &mut Model) -> Body {
    m.profile(Plane::XY)
        .rect(1.0, 1.0)
        .sketch("cube")
        .expect("sketch")
        .extrude(m, End::Blind(1.0.into()))
        .expect("extrude")
}

fn rust_error(scenario: impl FnOnce(&mut Model, Body) -> Result<(), ScriptError>) -> Value {
    let mut m = Model::new();
    let body = cube(&mut m);
    scenario(&mut m, body)
        .expect_err("scenario must fail")
        .payload()
}

fn rhai_error(tail: &str) -> popkern_rhai::ScriptRunError {
    let mut host = ScriptHost::new();
    host.run(&format!("{CUBE}{tail}\n"), "parity.rhai")
        .expect_err("scenario must fail")
}

/// A scenario: the Rust recorder call and the Rhai spelling of it.
struct Scenario {
    name: &'static str,
    kind: &'static str,
    rust: fn(&mut Model, Body) -> Result<(), ScriptError>,
    rhai: &'static str,
}

fn scenarios() -> Vec<Scenario> {
    vec![
        Scenario {
            name: "oversized fillet",
            kind: "FeatureFailed",
            rust: |m, b| b.fillet(m, edges(), 0.9).map(|_| ()),
            rhai: "cube.fillet(edges(), 0.9);",
        },
        Scenario {
            name: "negative radius",
            kind: "FeatureFailed",
            rust: |m, b| b.fillet(m, edges(), -0.1).map(|_| ()),
            rhai: "cube.fillet(edges(), -0.1);",
        },
        Scenario {
            name: "empty selection",
            kind: "EmptySelection",
            rust: |m, b| b.fillet(m, edges().circular(), 0.1).map(|_| ()),
            rhai: "cube.fillet(edges().circular(), 0.1);",
        },
        Scenario {
            name: "cardinality",
            kind: "SelectionCardinality",
            rust: |m, _| {
                m.label_face("top", &faces().kind(SurfaceKind::Plane))
                    .map(|_| ())
            },
            rhai: "m.label_face(\"top\", faces().kind(\"plane\"));",
        },
        Scenario {
            name: "unknown parameter",
            kind: "ArgumentError",
            rust: |m, _| m.set_param("nope", 1.0),
            rhai: "m.set_param(\"nope\", 1.0);",
        },
        Scenario {
            name: "label not found",
            kind: "LabelError",
            rust: |m, _| m.feature("nope").map(|_| ()),
            rhai: "m.feature(\"nope\");",
        },
        Scenario {
            name: "duplicate label",
            kind: "LabelError",
            rust: |m, _| {
                let f = m.feature("cube")?;
                m.rename(&f, "extrude-1")
            },
            rhai: "m.rename(m.feature(\"cube\"), \"extrude-1\");",
        },
        Scenario {
            name: "consumed tool",
            kind: "TimelineError",
            rust: |m, b| {
                let tool = m.tool(|m| {
                    m.profile(Plane::xy_at(0.5))
                        .circle(0.2)
                        .sketch("t")?
                        .extrude(m, End::Blind(1.0.into()))
                })?;
                let b = b.subtract(m, tool)?;
                b.subtract(m, tool).map(|_| ())
            },
            rhai: "let t = m.tool(|m| m.profile(xy_at(0.5)).circle(0.2).sketch(\"t\").extrude(blind(1.0)));\ncube = cube.subtract(t);\ncube.subtract(t);",
        },
        Scenario {
            name: "missing document",
            kind: "IoFailed",
            rust: |_, _| Model::open("/nonexistent/dir/model.pkn").map(|_| ()),
            rhai: "open(\"/nonexistent/dir/model.pkn\");",
        },
        Scenario {
            name: "unknown expression parameter",
            kind: "ArgumentError",
            rust: |m, _| m.param_expr("h", "wall * 2", Unit::Meters),
            rhai: "m.param_expr(\"h\", \"wall * 2\");",
        },
    ]
}

#[test]
fn every_reachable_class_crosses_field_for_field_and_float_for_float() {
    for s in scenarios() {
        let expected = rust_error(s.rust);
        let err = rhai_error(s.rhai);
        assert_eq!(err.kind, s.kind, "[{}] {err}", s.name);
        assert_eq!(expected["kind"], s.kind, "[{}] rust kind", s.name);
        // The partial replay report carries run-specific timings; the
        // classes and fields are compared with it removed.
        let strip = |mut v: Value| {
            if let Value::Object(map) = &mut v {
                map.remove("report");
            }
            v
        };
        assert_eq!(
            strip(*err.payload.clone()),
            strip(expected),
            "[{}] payload parity",
            s.name
        );
        assert!(err.line.is_some(), "[{}] positioned", s.name);
    }
}

#[test]
fn oversized_radius_reports_max_radius_and_remediation_at_the_exact_position() {
    let err = rhai_error("\n  cube.fillet(edges(), 0.9);");
    assert_eq!(err.kind, "FeatureFailed", "{err}");
    // CUBE is two lines, the tail starts with a blank line, the call
    // is indented two columns: line 4, column 8 — Rhai positions a
    // method call at the method name (`cube.` occupies columns 3–7).
    assert_eq!((err.line, err.column), (Some(4), Some(8)), "{err}");
    let p = &err.payload;
    assert_eq!(p["feature"]["label"], "fillet-1");
    assert_eq!(p["feature"]["operation"], "fillet");
    assert_eq!(p["error"]["kind"], "BlendDeclined");
    let failures = p["error"]["failures"].as_array().expect("failures[]");
    assert!(!failures.is_empty());
    for f in failures {
        assert_eq!(f["reason_kind"], "RadiusTooLarge");
        assert!(f["max_radius"].as_f64().is_some_and(|m| m > 0.0 && m < 0.9));
        assert_eq!(f["remediation"]["action"], "ReduceRadius");
        assert!(f["remediation"]["value"].as_f64().is_some_and(|v| v < 0.9));
        assert!(f["edge_name"].is_object(), "stable edge name: {f}");
    }
    let rendered = err.render();
    assert!(
        rendered.starts_with("FeatureFailed at parity.rhai:4:8: "),
        "{rendered}"
    );
    assert!(rendered.contains("ReduceRadius"), "{rendered}");
    assert_eq!(remediation_lines(&err.payload).len(), 1 + failures.len());
}

#[test]
fn invalid_parameter_carries_the_surface_name_and_the_kernel_path() {
    let err = rhai_error("cube.fillet(edges(), -0.1);");
    let e = &err.payload["error"];
    assert_eq!(e["kind"], "BlendDeclined", "{}", err.payload);
    assert_eq!(e["blend_kind"], "InvalidParameter", "{}", err.payload);
    assert_eq!(e["parameter"], "radius");
    assert_eq!(e["kernel_path"], "specs[0].radius");
    assert!(e["reason"].is_string());
    // Rust projection agrees (the translation lives in popkern-script).
    let rust = rust_error(|m, b| b.fillet(m, edges(), -0.1).map(|_| ()));
    assert_eq!(rust["error"]["parameter"], "radius");
    assert_eq!(rust["error"]["kernel_path"], "specs[0].radius");
    assert_eq!(
        popkern_script::surface_parameter("options.profile.rho"),
        "profile"
    );
    assert_eq!(popkern_script::surface_parameter("specs[2].edge"), "edges");
    assert_eq!(popkern_script::surface_parameter("dist_b"), "distance_2");
}

#[test]
fn skin_level_argument_errors_run_no_kernel_work() {
    // A misspelled enum leaves the document exactly as it was: no
    // feature, no parameter, no replay.
    let mut host = ScriptHost::new();
    let source = format!(
        "{CUBE}let before = m.features().len();\nlet params = m.params().len();\ntry {{ cube.fillet(edges(), 0.1, #{{ continuity: \"g3\" }}); }} catch (e) {{ print(e.kind + \" \" + e.argument); }}\nprint(m.features().len() == before);\nprint(m.params().len() == params);\nm"
    );
    let outcome = host
        .run(&source, "enum.rhai")
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(outcome.printed[0], "ArgumentError continuity");
    assert_eq!(outcome.printed[1], "true");
    assert_eq!(outcome.printed[2], "true", "no parameter was bound");
}

#[test]
fn caught_map_matches_the_uncaught_payload() {
    let uncaught = rhai_error("cube.fillet(edges(), 0.9);");
    let mut host = ScriptHost::new();
    let source = format!(
        "{CUBE}let e = ();\ntry {{ cube.fillet(edges(), 0.9); }} catch (err) {{ e = err; }}\ne"
    );
    let outcome = host
        .run(&source, "catch.rhai")
        .unwrap_or_else(|e| panic!("{e}"));
    let _ = outcome;
    // Re-run returning the map as JSON text through the script itself.
    let source = format!(
        "{CUBE}let e = ();\ntry {{ cube.fillet(edges(), 0.9); }} catch (err) {{ e = err; }}\nprint(e.error.failures[0].max_radius);\nprint(e.error.total);\nprint(e.kind);"
    );
    let outcome = host
        .run(&source, "catch.rhai")
        .unwrap_or_else(|e| panic!("{e}"));
    let expected = uncaught.payload["error"]["failures"][0]["max_radius"]
        .as_f64()
        .expect("max_radius");
    assert_eq!(outcome.printed[0].parse::<f64>().expect("float"), expected);
    assert_eq!(outcome.printed[1], "12");
    assert_eq!(outcome.printed[2], "FeatureFailed");
}
