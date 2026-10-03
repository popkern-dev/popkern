//! PK-629/PK-630 Rhai surface: `assembly()`, parts recorded by a
//! closure, sub-assemblies, chained rigid placements, part families
//! (`a.part(name, defaults, |m, p| …)` + `.variant(…)`), errors raised in
//! a part closure, and show-as-code round-tripping the tree.

use popkern_assembly::{Assembly, Document};
use popkern_rhai::{ScriptHost, document_to_rhai};

const BEARING: &str = r#"
let a = assembly("balls");
let ball = a.part("ball", |m| {
    print("building ball");
    m.profile(xz())
        .path(0.0, 0.05).arc_to([0.15, 0.2], [0.0, 0.35])
        .line_to(0.0, 0.05).close()
        .sketch("ball")
        .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
});
let pi = 3.141592653589793;
for i in 0..4 {
    let t = i * 2.0 * pi / 4;
    a.place(ball, "ball_" + i).translate(0.7 * cos(t), 0.7 * sin(t), 0.0);
}
let pair = a.subassembly("pair", |s| {
    s.place(ball, "low");
    s.place(ball, "high").translate(0.0, 0.0, 1.0);
});
a.place(pair).rotate(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), 90.0).translate(2.0, 0.0, 0.0);
a.place(pair, "hidden").hide();
a
"#;

fn run(source: &str) -> (Assembly, Vec<String>) {
    let mut host = ScriptHost::new();
    let mut outcome = host.run(source, "asm.rhai").expect("run");
    let printed = outcome.printed.clone();
    match outcome.take_document() {
        Some(Document::Assembly(a)) => (*a, printed),
        _ => panic!("the script returns an assembly"),
    }
}

#[test]
fn a_script_builds_each_part_once_and_places_instances() {
    let (mut asm, printed) = run(BEARING);
    assert_eq!(printed, ["building ball"], "the closure ran once");
    assert_eq!(asm.definitions().len(), 2);
    let rows: Vec<(String, bool)> = asm
        .tree()
        .iter()
        .map(|r| (r.path.to_string(), r.visible))
        .collect();
    assert_eq!(
        rows,
        [
            ("ball_0".into(), true),
            ("ball_1".into(), true),
            ("ball_2".into(), true),
            ("ball_3".into(), true),
            ("pair_1".into(), true),
            ("pair_1/low".into(), true),
            ("pair_1/high".into(), true),
            ("hidden".into(), false),
            ("hidden/low".into(), false),
            ("hidden/high".into(), false),
        ]
    );
    // pair_1 turned 90° about z THEN moved +2 in x: its "high" ball
    // (local z +1) sits at (2, 0, 1).
    let high = asm
        .tree()
        .into_iter()
        .find(|r| r.path.to_string() == "pair_1/high")
        .expect("high");
    let o = high.world.origin();
    assert!(
        (o.x - 2.0).abs() < 1e-12 && o.y.abs() < 1e-12 && (o.z - 1.0).abs() < 1e-12,
        "{o:?}"
    );
    assert_eq!(
        asm.flatten().expect("flatten").len(),
        8,
        "eight part occurrences"
    );
}

#[test]
fn errors_inside_a_part_name_the_part_and_the_closure_line() {
    let source = "let a = assembly();\nlet p = a.part(\"plate\", |m| {\n    m.profile(xy()).rect(1.0, 1.0).sketch(\"s\").extrude(blind(1.0))\n        .fillet(edges().circular(), 0.1);\n});\na";
    let mut host = ScriptHost::new();
    let error = host.run(source, "bad.rhai").expect_err("refused");
    assert_eq!(error.line, Some(4), "{error}");
    assert_eq!(error.payload["part"], "plate", "{}", error.payload);
    assert!(
        error.message.starts_with("in part 'plate': "),
        "{}",
        error.message
    );

    let duplicate = "let a = assembly();\nlet p = a.part(\"x\", |m| { m.profile(xy()).rect(1.0, 1.0).sketch(\"s\").extrude(blind(1.0)); });\na.place(p, \"one\");\na.place(p, \"one\");\na";
    let error = ScriptHost::new()
        .run(duplicate, "dup.rhai")
        .expect_err("refused");
    assert_eq!(error.kind, "AssemblyError");
    assert_eq!(error.payload["reason"], "duplicate_instance_name");
    assert_eq!(error.line, Some(4));

    let edge_axis = "let a = assembly();\nlet p = a.part(\"x\", |m| { m.profile(xy()).rect(1.0, 1.0).sketch(\"s\").extrude(blind(1.0)); });\na.place(p).rotate(axis(edges().straight()), 10.0);\na";
    let error = ScriptHost::new()
        .run(edge_axis, "axis.rhai")
        .expect_err("refused");
    assert_eq!(error.kind, "ArgumentError", "{error}");
}

#[test]
fn a_family_builds_one_variant_per_parameter_set() {
    let source = r#"
let a = assembly("kit");
let block = a.part("block", #{ w: 1.0, h: 2, tag: "a" }, |m, p| {
    print("build " + p.w + " " + p.h + " " + p.tag);
    m.profile(xy()).rect(p.w, p.h).sketch("s").extrude(blind(0.5));
});
a.place(block);
a.place(block.variant(#{ w: 3.0 }));
a.place(block.variant(#{ w: 3.0 }));
a.place(block.variant(#{ tag: "b" }));
a
"#;
    let (asm, printed) = run(source);
    assert_eq!(
        printed,
        ["build 1.0 2.0 a", "build 3.0 2.0 a", "build 1.0 2.0 b"]
    );
    let names: Vec<&str> = asm.definitions().iter().map(|d| d.name.as_str()).collect();
    assert_eq!(names, ["block", "block[w=3]", "block[tag=b]"]);
    assert_eq!(asm.tree().len(), 4);

    for (with, reason) in [
        (
            "#{ width: 3.0 }",
            "has no parameter 'width' (parameters: h, tag, w)",
        ),
        (
            "#{ w: \"wide\" }",
            "parameter 'w' is a number, got a string",
        ),
    ] {
        let bad = source.replace(
            "block.variant(#{ tag: \"b\" })",
            &format!("block.variant({with})"),
        );
        let error = ScriptHost::new()
            .run(&bad, "bad.rhai")
            .expect_err("refused");
        assert!(error.message.contains(reason), "{}", error.message);
    }
}

#[test]
fn show_as_code_rebuilds_the_same_tree() {
    let (mut asm, _) = run(BEARING);
    let code = document_to_rhai(&mut Document::Assembly(Box::new(Assembly::new("tmp"))))
        .expect("empty emits");
    assert!(code.starts_with("let a = assembly(\"tmp\");"));
    let mut doc = Document::Assembly(Box::new(std::mem::replace(&mut asm, Assembly::new(""))));
    let code = document_to_rhai(&mut doc).expect("emit");
    assert!(code.contains("a.part(\"ball\", |m| {"), "{code}");
    assert!(code.contains(".frame(["), "{code}");
    assert!(code.contains(".hide()"), "{code}");
    let Document::Assembly(mut original) = doc else {
        unreachable!()
    };
    let (mut again, _) = run(&code);
    let (a, b) = (original.tree(), again.tree());
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(&b) {
        assert_eq!(x.path, y.path);
        assert_eq!(x.visible, y.visible);
        assert!(
            x.world.transform().approx_eq(y.world.transform(), 1e-12),
            "{}",
            x.path
        );
    }
    assert_eq!(
        original.flatten().expect("a").len(),
        again.flatten().expect("b").len()
    );
    // And the emitted code is a fixed point.
    let code_again = document_to_rhai(&mut Document::Assembly(Box::new(again))).expect("emit");
    assert_eq!(code_again, code);
}

#[test]
fn an_error_inside_a_tool_closure_or_a_script_fn_keeps_its_class_and_line() {
    let tool = "let m = model();\nlet b = m.tool(|m| {\n    m.profile(xy()).rect(1.0, 1.0).sketch(\"s\").extrude(blind(1.0))\n        .fillet(edges().circular(), 0.1)\n});\nm";
    let error = ScriptHost::new()
        .run(tool, "tool.rhai")
        .expect_err("refused");
    assert_eq!(error.kind, "EmptySelection", "{error}");
    assert_eq!(error.line, Some(4), "{error}");
    let in_fn = "fn plate(m) {\n    m.profile(xy()).rect(1.0, 1.0).sketch(\"s\").extrude(blind(1.0)).fillet(edges().circular(), 0.1)\n}\nlet m = model();\nplate(m);\nm";
    let error = ScriptHost::new()
        .run(in_fn, "fn.rhai")
        .expect_err("refused");
    assert_eq!(error.kind, "EmptySelection", "{error}");
    assert_eq!(error.line, Some(2), "{error}");
}

#[test]
fn show_as_code_keeps_assembly_parameters_and_drives() {
    let source = r#"
let a = assembly("plates");
let thin = a.part("thin", |m| {
    m.param("t", 0.1, "m");
    m.profile(xy()).rect(1.0, 1.0).sketch("thin").extrude(blind(p("t")));
});
a.place(thin);
a.param("t", 0.25, "m");
a.drive("t", [thin]);
a
"#;
    let (asm, _) = run(source);
    let mut doc = Document::Assembly(Box::new(asm));
    let code = document_to_rhai(&mut doc).expect("emit");
    assert!(code.contains("a.param(\"t\", 0.25, \"m\");"), "{code}");
    assert!(code.contains("a.drive(\"t\", d0, \"t\");"), "{code}");
    let (again, _) = run(&code);
    let Document::Assembly(original) = doc else {
        unreachable!()
    };
    assert_eq!(again.params(), original.params());
    assert_eq!(again.drives(), original.drives());
}
