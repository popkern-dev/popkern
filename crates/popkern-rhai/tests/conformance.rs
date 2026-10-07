//! PK-405: the front-end conformance checklist is generated from the
//! scripting-surface spec, checked in, and audits both skins. A spec
//! row without an implementation fails here; a stale checked-in
//! checklist fails here (regenerate with `POPKERN_UPDATE_CHECKLIST=1`).

use std::path::{Path, PathBuf};

use popkern_rhai::{
    Checklist, ScriptHost, audit_error_classes, audit_rhai, audit_rust, build_engine,
    checklist_from_spec,
};

const SPEC: &str = "docs/architecture/scripting-surface.md";
const CHECKLIST: &str = "docs/reference/scripting-surface-checklist.json";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("workspace root")
}

fn spec() -> String {
    std::fs::read_to_string(root().join(SPEC)).expect("spec")
}

fn generated() -> Checklist {
    checklist_from_spec(&spec(), SPEC)
}

fn sources(paths: &[&str]) -> String {
    paths
        .iter()
        .map(|p| std::fs::read_to_string(root().join(p)).expect(p))
        .collect::<Vec<_>>()
        .join("\n")
}

fn crate_sources(dir: &str) -> String {
    let mut out = String::new();
    for entry in std::fs::read_dir(root().join(dir)).expect("src dir") {
        let path = entry.expect("entry").path();
        if path.extension().is_some_and(|e| e == "rs") {
            out.push_str(&std::fs::read_to_string(&path).expect("source"));
            out.push('\n');
        }
    }
    out
}

#[test]
fn checked_in_checklist_matches_the_spec() {
    let checklist = generated();
    assert!(
        checklist.ops.len() > 60,
        "ops parsed: {}",
        checklist.ops.len()
    );
    assert!(
        checklist.predicates.len() >= 15,
        "{:?}",
        checklist.predicates
    );
    assert!(
        checklist.error_classes.len() >= 15,
        "{:?}",
        checklist.error_classes
    );
    let text = serde_json::to_string_pretty(&checklist).expect("json") + "\n";
    let path = root().join(CHECKLIST);
    if std::env::var_os("POPKERN_UPDATE_CHECKLIST").is_some() {
        std::fs::write(&path, &text).expect("write checklist");
    }
    let checked_in = std::fs::read_to_string(&path).unwrap_or_default();
    assert_eq!(
        checked_in, text,
        "{CHECKLIST} is stale: regenerate with POPKERN_UPDATE_CHECKLIST=1"
    );
    // Idempotent: generating twice yields the same artifact.
    assert_eq!(generated(), checklist);
}

#[test]
fn rhai_skin_conforms() {
    let engine = build_engine(|_| {});
    let findings = audit_rhai(&engine, &generated());
    assert!(
        findings.is_empty(),
        "Rhai conformance findings:\n{}",
        findings
            .iter()
            .map(|f| format!("  {}: {}", f.subject, f.detail))
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn rust_skin_conforms() {
    let findings = audit_rust(&generated(), &crate_sources("crates/popkern-script/src"));
    assert!(
        findings.is_empty(),
        "Rust conformance findings:\n{}",
        findings
            .iter()
            .map(|f| format!("  {}: {}", f.subject, f.detail))
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn every_error_class_is_produced_with_its_fields() {
    let text = sources(&[
        "crates/popkern-script/src/error_payload.rs",
        "crates/popkern-history/src/failure.rs",
        "crates/popkern-history/src/failure_blend.rs",
    ]);
    let findings = audit_error_classes(&generated(), &text);
    assert!(
        findings.is_empty(),
        "error-class findings:\n{}",
        findings
            .iter()
            .map(|f| format!("  {}: {}", f.subject, f.detail))
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn a_spec_row_without_an_implementation_fails_the_audit() {
    let mut text = spec();
    let anchor = "| `chamfer` | stable |";
    let at = text.find(anchor).expect("chamfer row");
    text.insert_str(
        at,
        "| `frobnicate` | stable | `frobnicate(edges, amount)` | `Feature::Frobnicate` | |\n",
    );
    let checklist = checklist_from_spec(&text, SPEC);
    assert!(checklist.ops.iter().any(|o| o.name == "frobnicate"));
    let engine = build_engine(|_| {});
    let findings = audit_rhai(&engine, &checklist);
    assert!(
        findings.iter().any(|f| f.subject == "frobnicate"),
        "{findings:?}"
    );
    assert!(
        audit_rust(&checklist, &crate_sources("crates/popkern-script/src"))
            .iter()
            .any(|f| f.subject == "frobnicate")
    );
}

/// Every stable op runs in one script: the checklist's stable set is
/// the set this script exercises.
#[test]
fn every_stable_op_runs_in_the_exercise_script() {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("conformance");
    std::fs::create_dir_all(&dir).expect("dir");
    let pkn = dir.join("exercise.pkn");
    let step = dir.join("exercise.step");
    let script = format!(
        r#"let m = model(#{{ unit: "m" }});
m.param("wall", 0.4, "m");
m.param_expr("lip", "wall * 0.5", "m");
let plate = m.profile(xy()).rect(2.4, 1.6).sketch("plate_sketch").extrude(blind(p("wall")));
let boss = m.tool(|m| m.profile(xy_at(0.2)).circle(0.35).at(0.0, 0.0).sketch("boss_sketch").extrude(blind(0.8)));
plate = plate.union(boss);
plate = plate.fillet(edges().between("cylinder", "plane").at_height(0.4).or(edges().straight().at_height(0.4)).or(edges().straight().direction("z")), 0.08, #{{ label: "rounds" }});
let post = m.tool(|m| m.profile(yz()).polygon([[-0.1, -0.1], [0.1, -0.1], [0.1, 0.1], [-0.1, 0.1]]).sketch("post_sketch").extrude(blind(0.3)));
post = post.translate(1.6, 0.0, 0.2);
post = post.rotate(axis_at([1.6, 0.0, 0.2], [0.0, 0.0, 1.0]), 15.0);
let lug = m.tool(|m| m.profile(xz()).path(-0.2, 0.0).line_to(0.2, 0.0).arc_to([0.15, 0.2], [0.0, 0.3]).close().sketch("lug_sketch").extrude(blind(0.1)));
let cup = m.tool(|m| m.profile(xy_at(2.0)).rect(0.6, 0.6).sketch("cup_sketch").extrude(blind(0.6)));
cup = cup.chamfer(edges().straight().at_height(2.6), 0.05);
cup = cup.shell(faces().facing("z"), 0.05);
let disc = m.tool(|m| m.profile(xz()).rect(0.2, 0.1).at(0.6, 3.0).sketch("disc_sketch").revolve(axis_at([0.0, 0.0, 3.0], [0.0, 0.0, 1.0]), full_circle()));
let cutter = m.tool(|m| m.profile(xy_at(0.3)).circle(0.05).at(0.9, 0.5).sketch("cutter_sketch").extrude(blind(0.3)));
plate = plate.subtract(cutter);
let bar = m.tool(|m| m.profile(xy_at(4.0)).rect(0.3, 0.3).sketch("bar_sketch").extrude(blind(0.3)));
let blob = m.tool(|m| m.profile(xy_at(4.1)).circle(0.2).sketch("blob_sketch").extrude(blind(0.3)));
bar = bar.intersect(blob);
let spring = m.tool(|m| m.profile(xy_at(5.0)).circle(0.05).sketch("spring_wire").helical_sweep(axis_at([0.0, 0.0, 5.0], [0.0, 0.0, 1.0]), 0.4, 0.2, 3.0, #{{ start_angle: 30.0, label: "spring_coil" }}));
let s0 = m.tool(|m| m.profile(xy_at(6.0)).circle(0.3).sketch("loft_base"));
let s1 = m.tool(|m| m.profile(xy_at(6.3)).circle(0.2).sketch("loft_tip"));
let horn = m.loft([s0, s1], #{{ label: "horn" }});
let rail = m.tool(|m| m.profile(xz()).rect(0.2, 0.1).at(0.7, 7.05).sketch("rail_section").revolve(axis_at([0.0, 0.0, 7.0], [0.0, 0.0, 1.0]), mid_plane(45.0), #{{ label: "rail_ring" }}));
let tube = rail.pipe(edges().at_height(7.1).chord_length_between(0.55, 0.7), 0.04, #{{ label: "rail_tube" }});
let ribbon = m.tool(|m| m.profile(xy_at(8.0)).rect(0.06, 0.06).sketch("ribbon_section"));
let swept = ribbon.sweep(rail, edges().at_height(7.0).chord_length_between(0.55, 0.7), #{{ twist: 45.0, label: "ribbon_sweep" }});
print(spring.is_alive() && horn.is_alive() && tube.is_alive() && swept.is_alive());
m.set_active(plate);
m.label_edge("rim", edges().circular().at_height(1.0));
m.label_face("top", faces().facing("z").at_height(0.4));
m.rename_label("rim", "boss_rim");
m.unlabel("boss_rim");
print(m.labels().len());
let rounds = m.feature("rounds");
print(rounds.describe());
print(type_of(rounds.evidence()));
rounds.modify(#{{ radius: 0.1 }});
rounds.suppress();
rounds.unsuppress();
m.rename(rounds, "rounds");
m.set_param("wall", 0.4);
m.set_expr("lip", "wall * 0.25");
m.clear_expr("lip");
m.set_bounds("wall", 0.1, 1.0);
m.clear_bounds("wall");
m.rename_param("lip", "rim_lip");
print(m.params().len() > 0);
print(m.features().len());
print(m.bodies().len() > 0);
m.label_body("main", plate);
print(m.body("main").id == plate.id);
let stud = m.threaded_rod("M6", 4.0, #{{ tip_chamfer: false, label: "stud" }});
print(stud.is_alive());
let h = m.hinge("butt", 0.6, #{{ only: "leaf_a", knuckles: 3, label: "ex_hinge" }});
let hb = m.hinge_component("ex_hinge", "leaf_b");
let flex = m.living_hinge("straight", 0.2, 0.4, #{{ label: "flex" }});
print(h.leaf_a.is_alive() && hb.is_alive() && flex.is_alive());
m.set_active(plate);
m.replay();
print(m.report().failures.len());
m.assert_watertight();
m.to_step({step:?});
m.assert_step_roundtrip({step:?});
m.save({pkn:?});
let again = open({pkn:?});
print(again.features().len());
m
"#
    );
    let mut host = ScriptHost::new();
    let outcome = host
        .run(&script, "exercise.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    assert!(outcome.model.is_some());
    let stable: Vec<String> = generated()
        .ops
        .iter()
        .filter(|o| o.mark == "stable")
        .map(|o| o.name.clone())
        .collect();
    let bound = |op: &str| -> String {
        popkern_rhai::conformance_audit::RHAI_BINDINGS
            .iter()
            .find(|b| b.op == op)
            .map(|b| b.function.to_string())
            .unwrap_or_else(|| op.to_string())
    };
    let missing: Vec<&String> = stable
        .iter()
        .filter(|op| !script.contains(&format!("{}(", bound(op))))
        .collect();
    assert!(missing.is_empty(), "stable ops not exercised: {missing:?}");
}
