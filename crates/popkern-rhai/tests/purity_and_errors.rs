//! Script purity (design invariant 9) and the error map: no clock,
//! RNG, environment, file, or `eval` built-in is reachable; a
//! recorder decline reaches the script as a §10 map with `kind` and
//! fields, positioned at the raising call, and `try`/`catch` receives
//! the same map.

use popkern_rhai::{FORBIDDEN_FUNCTIONS, ScriptHost, build_engine, registered_function_names};

const CUBE: &str = r#"
let m = model();
let cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
"#;

#[test]
fn no_impure_function_is_registered() {
    let engine = build_engine(|_| {});
    let names = registered_function_names(&engine);
    assert!(
        names.iter().any(|n| n == "fillet"),
        "dialect registered: {names:?}"
    );
    for forbidden in FORBIDDEN_FUNCTIONS {
        assert!(
            !names.iter().any(|n| n == forbidden),
            "'{forbidden}' must not be callable; registered: {names:?}"
        );
    }
    for suspicious in [
        "rand",
        "random",
        "env",
        "open_file",
        "read_file",
        "exit",
        "spawn",
    ] {
        assert!(
            !names.iter().any(|n| n == suspicious),
            "'{suspicious}' must not be callable"
        );
    }
}

#[test]
fn calling_a_forbidden_built_in_fails_at_its_line() {
    let mut host = ScriptHost::new();
    // Unregistered functions fail at the call (skin-level, before any
    // kernel work); `eval` is a keyword and fails at parse.
    for (source, line, kind) in [
        (
            "let m = model();\nlet t = timestamp();\n",
            2,
            "ArgumentError",
        ),
        ("let m = model();\nsleep(1);\n", 2, "ArgumentError"),
        ("let m = model();\nexit();\n", 2, "ArgumentError"),
        (
            "let m = model();\n\nlet x = eval(\"1 + 1\");\n",
            3,
            "SyntaxError",
        ),
    ] {
        let err = host.run(source, "impure.rhai").expect_err("must fail");
        assert_eq!(err.kind, kind, "{err}");
        assert_eq!(err.line, Some(line), "{err}");
    }
}

#[test]
fn a_declined_fillet_is_a_blend_declined_map_at_the_fillet_line() {
    let mut host = ScriptHost::new();
    let source = format!("{CUBE}\ncube.fillet(edges(), 0.9);\n");
    let err = host
        .run(&source, "cube.rhai")
        .expect_err("r = 0.9 declines");
    assert_eq!(err.kind, "FeatureFailed", "{err}");
    assert_eq!(err.line, Some(5), "raised at the fillet call: {err}");
    let payload = &err.payload;
    assert_eq!(payload["feature"]["operation"], "fillet");
    assert_eq!(payload["error"]["kind"], "BlendDeclined");
    assert_eq!(payload["error"]["blend_kind"], "PartialFailure");
    assert_eq!(payload["error"]["total"], 12);
    let first = &payload["error"]["failures"][0];
    assert_eq!(first["reason_kind"], "RadiusTooLarge");
    assert!(
        first["remediation"]["value"]
            .as_f64()
            .is_some_and(|v| v > 0.0 && v < 0.9)
    );
    assert!(
        payload["message"]
            .as_str()
            .is_some_and(|m| m.contains("of 12 edges failed"))
    );
}

#[test]
fn try_catch_receives_the_same_map_and_the_document_is_untouched() {
    let mut host = ScriptHost::new();
    let source = format!(
        r#"{CUBE}
let caught = ();
try {{ cube.fillet(edges(), 0.9); }} catch (e) {{ caught = e; }}
print(caught.kind);
print(caught.error.kind);
print(caught.error.failures.len());
print(m.features().len());
m
"#
    );
    let outcome = host
        .run(&source, "catch.rhai")
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(outcome.printed[0], "FeatureFailed");
    assert_eq!(outcome.printed[1], "BlendDeclined");
    assert!(outcome.printed[2].parse::<usize>().is_ok_and(|n| n >= 1));
    assert_eq!(outcome.printed[3], "2", "failed feature removed again");
    let mut model = outcome.model.expect("model");
    assert_eq!(model.to_part().expect("part").face_count(), 6);
}

#[test]
fn misspelled_enums_and_unknown_options_fail_before_the_kernel_runs() {
    let mut host = ScriptHost::new();
    for (tail, argument) in [
        (
            "cube.fillet(edges(), 0.1, #{ continuity: \"g3\" });",
            "continuity",
        ),
        ("cube.fillet(edges(), 0.1, #{ colour: \"red\" });", "colour"),
        (
            "cube.shell(faces().facing(\"z\"), 0.1, #{ direction: \"sideways\" });",
            "direction",
        ),
        (
            "cube.fillet(edges().between(\"cylinder\", \"cube\"), 0.1);",
            "b",
        ),
    ] {
        let source = format!("{CUBE}\n{tail}\nm");
        let err = host.run(&source, "bad.rhai").expect_err(tail);
        assert_eq!(err.kind, "ArgumentError", "{tail}: {err}");
        assert_eq!(err.payload["argument"], argument, "{tail}: {}", err.payload);
        assert_eq!(err.line, Some(5), "{tail}: {err}");
    }
}

#[test]
fn empty_selection_and_syntax_errors_are_classified() {
    let mut host = ScriptHost::new();
    let err = host
        .run(
            &format!("{CUBE}\ncube.fillet(edges().circular(), 0.1);"),
            "empty.rhai",
        )
        .expect_err("no circular edges on a cube");
    assert_eq!(err.kind, "EmptySelection");
    assert_eq!(err.payload["selection"], "edges().circular()");
    assert_eq!(err.payload["available"], 12);
    let err = host
        .run("let m = model(\n", "syntax.rhai")
        .expect_err("unbalanced");
    assert_eq!(err.kind, "SyntaxError");
    assert!(err.line.is_some());
}
