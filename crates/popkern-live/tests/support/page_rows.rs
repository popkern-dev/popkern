//! The page harness rows for the time limit and the pre-built gallery
//! examples (split from `page.rs` for its line budget).

use serde_json::Value;

/// The time limit: presets up to "no limit", default 24 s (20 s kernel
/// budget + the page's grace), "no limit" lifts both budgets. Pre-built
/// examples, when the bundle was baked (`tools/bake-live.sh`):
/// shown without the kernel, certified, and rebuilt live on demand.
pub fn assert_time_limit_and_bakes(out: &Value) {
    // The time limit: presets up to "no limit", default 24 s (20 s
    // kernel budget + the page's grace), "no limit" lifts both budgets.
    let limit = &out["time_limit"];
    assert_eq!(limit["initial"], "24", "{limit}");
    assert!(
        limit["options"]
            .as_array()
            .is_some_and(|o| o.len() >= 5 && o.last().is_some_and(|l| l == "no limit")),
        "{limit}"
    );
    assert!(limit["none"]["host_timeout"].is_null(), "{limit}");
    assert!(limit["none"]["host_ops"].is_null(), "{limit}");
    assert_eq!(limit["none"]["stored"], "none");
    assert_eq!(limit["restored_timeout"], 20000.0);
    // Pre-built examples, when the bundle was baked
    // (tools/bake-live.sh): shown without the kernel, certified,
    // and rebuilt live on demand.
    if out["bakes_present"] == true {
        let prebuilt = &out["prebuilt"];
        assert!(
            prebuilt["status"]
                .as_str()
                .is_some_and(|s| s.starts_with("served (pre-built)")),
            "{prebuilt}"
        );
        assert_eq!(prebuilt["live"], false);
        assert!(
            prebuilt["wall_ms"].as_u64().is_some_and(|ms| ms < 5000),
            "a pre-built example shows in seconds, not a kernel run: {prebuilt}"
        );
        assert!(prebuilt["bodies_on_screen"].as_u64().unwrap_or(0) >= 1);
        assert!(
            prebuilt["report_text"]
                .as_str()
                .is_some_and(|t| t.contains("certified") && !t.contains("certifying")),
            "{prebuilt}"
        );
        let returning = &out["returning"];
        assert!(
            returning["status"]
                .as_str()
                .is_some_and(|s| s.starts_with("served (pre-built)")),
            "a returning visitor with an edited draft still opens pre-built: {returning}"
        );
        assert_eq!(returning["live"], false, "{returning}");
        assert_eq!(returning["notice_visible"], true, "{returning}");
        assert_eq!(
            returning["playground_draft_key_cleared"], true,
            "{returning}"
        );
        assert_eq!(returning["editor_is_default"], true, "{returning}");
        assert_eq!(returning["editor_restored"], true, "{returning}");
        assert_eq!(returning["notice_after"], false, "{returning}");
        assert_eq!(returning["scratch_open"], true, "{returning}");
        // PK-662: a read-only example, its editable copy, a conflict.
        assert_eq!(prebuilt["read_only"], true, "{prebuilt}");
        let copy = &out["copy_to_edit"];
        assert_eq!(copy["copy_editable"], true, "{copy}");
        assert_eq!(copy["copy_saved"], true, "{copy}");
        assert_eq!(copy["example_unchanged"], true, "{copy}");
        let project = &out["project_example"];
        assert_eq!(project["outcome"], "served", "{project}");
        assert_eq!(project["features"], 6, "{project}");
        assert_eq!(
            project["versions"],
            serde_json::json!(["001_plate", "002_boss", "003_rounds"]),
            "{project}"
        );
        assert_eq!(project["read_only"], true, "{project}");
        let conflict = &out["conflict"];
        assert_eq!(conflict["bar"], true, "{conflict}");
        assert_eq!(
            conflict["theirs_kept"], true,
            "nothing overwritten: {conflict}"
        );
        assert_eq!(conflict["loaded_theirs"], true, "{conflict}");
        assert_eq!(conflict["bar_after"], false, "{conflict}");
        let cheap = &out["prebuilt_cheap"];
        assert_eq!(cheap["made_live"], true, "{cheap}");
        assert!(
            cheap["status_after_live"]
                .as_str()
                .is_some_and(|s| s.starts_with("served:")),
            "{cheap}"
        );
    } else {
        eprintln!("[bakes] SKIP the pre-built rows: no web/live/bakes (run tools/bake-live.sh)");
    }
}

/// An assembly in the page (PK-629): the tree lists every occurrence,
/// each part definition is ONE GPU geometry however many occurrences
/// draw it, hiding a sub-assembly hides what it holds, opening a part
/// scopes the inspector, a `part:name` parameter edit moves every
/// instance, picks spell the occurrence, STEP carries the product
/// structure, and visibility is saved in the `.pkn`.
pub fn assert_assembly(out: &Value) {
    let a = &out["assembly"];
    assert_eq!(a["outcome"], "served", "{a}");
    assert_eq!(a["section_visible"], true, "{a}");
    assert_eq!(a["summary"], "pins: 1 part(s), 6 occurrence(s)", "{a}");
    assert_eq!(
        a["rows"],
        serde_json::json!([
            "pin_0",
            "pin_1",
            "pin_2",
            "pin_3",
            "pair_1",
            "pair_1/low",
            "pair_1/high"
        ]),
        "{a}"
    );
    assert_eq!(a["meshes"], 6, "one mesh per occurrence: {a}");
    assert_eq!(
        a["geometries"], 1,
        "one shared geometry for the one part: {a}"
    );
    assert_eq!(a["drawn_after_hide"], 4, "{a}");
    assert_eq!(
        a["hidden_rows"],
        serde_json::json!(["pair_1", "pair_1/low", "pair_1/high"]),
        "{a}"
    );
    assert_eq!(a["open_name"], "pin", "{a}");
    let params = a["open_params"].as_array().expect("params");
    assert!(
        params.iter().any(|p| p == "pin:radius")
            && params
                .iter()
                .all(|p| p.as_str().is_some_and(|p| p.starts_with("pin:"))),
        "an opened part shows its own (qualified) parameters only: {a}"
    );
    assert_eq!(
        a["open_report_rows"], 2,
        "the pin's sketch and extrude: {a}"
    );
    assert!(
        a["pick"]
            .as_str()
            .is_some_and(|p| p.starts_with("face pin_2:")),
        "an instance-path pick: {a}"
    );
    assert_eq!(a["param_outcome"], "served", "{a}");
    assert_eq!(a["radius_grew_everywhere"], true, "{a}");
    assert_eq!(
        a["drawn_after_edit"], 4,
        "the hidden sub-assembly stays hidden after a rebuild: {a}"
    );
    assert_eq!(a["step_occurrences"], 7, "4 + pair + 2 in the pair: {a}");
    assert_eq!(a["step_solids"], 1, "the pin's solid is written once: {a}");
    // Clearance (PK-630): three neighbouring pins 0.1 apart, nothing
    // interfering, the rows shown.
    let clear: Vec<&Value> = a["clearance"]
        .as_array()
        .expect("clearance")
        .iter()
        .collect();
    assert_eq!(clear.len(), 3, "{a}");
    for row in &clear {
        assert_eq!(row[2], "clear", "{a}");
        assert!((row[3].as_f64().unwrap_or(0.0) - 0.1).abs() < 1e-6, "{a}");
    }
    assert_eq!(a["clearance_rows"], 3, "{a}");
    // Exploded view: transforms only, inside a 60 fps frame budget.
    assert_eq!(a["explode_visible"], true, "{a}");
    assert_eq!(a["explode_moved"], true, "{a}");
    assert!(
        a["explode_frame_ms"]["p95"]
            .as_f64()
            .is_some_and(|ms| ms < 16.7),
        "explode frame over the 60 fps budget: {a}"
    );
    assert_eq!(a["pkn_version"], 5, "{a}");
    assert_eq!(a["pkn_hidden"], serde_json::json!(["pair_1"]), "{a}");
}

/// Open STEP (PK-631): OpenCASCADE's assembly opens view-only with its
/// tree, names and one geometry per shared part.
pub fn assert_open_step(out: &Value) {
    let s = &out["open_step"];
    assert_eq!(s["outcome"], "served", "{s}");
    assert_eq!(
        s["rows"],
        serde_json::json!([
            "top=plate",
            "leg_front=leg",
            "leg_back=leg",
            "pair=leg_pair",
            "pair/left=leg",
            "pair/right=leg"
        ]),
        "{s}"
    );
    assert_eq!(
        s["geometries"], 2,
        "the leg is drawn four times from one geometry: {s}"
    );
    assert_eq!(s["meshes"], 5, "{s}");
    assert_eq!(s["tab"], "occ_bench_assembly.rhai", "{s}");
    assert_eq!(
        s["build_kept"], true,
        "Build keeps the imported document: {s}"
    );
    assert!(
        s["printed"][0]
            .as_str()
            .is_some_and(|p| p.starts_with("imported occ_bench_assembly: 2 part(s)")),
        "{s}"
    );
}

/// Same-origin asset loads only: `fetch("./relative/path")`.
pub fn network_violations(source: &str) -> Vec<String> {
    let mut found = Vec::new();
    for (n, line) in source.lines().enumerate() {
        let trimmed = line.trim_start();
        if trimmed.starts_with("//") {
            continue;
        }
        for api in [
            "XMLHttpRequest",
            "WebSocket",
            "sendBeacon",
            "EventSource",
            "importScripts",
        ] {
            if line.contains(api) {
                found.push(format!("line {}: {api}", n + 1));
            }
        }
        let mut rest = line;
        while let Some(at) = rest.find("fetch(") {
            let after = &rest[at + "fetch(".len()..];
            let literal = after.trim_start();
            let relative = (literal.starts_with("\"./") || literal.starts_with("`./"))
                && !literal.contains("http");
            if !relative {
                found.push(format!("line {}: fetch with a non-relative target", n + 1));
            }
            rest = after;
        }
    }
    found
}

/// PK-661 timeline and PK-660 project requests through the node worker.
pub fn assert_history_rows(out: &serde_json::Value) {
    // PK-661 timeline through the worker: lines, then key reuse.
    let timeline = &out["timeline"];
    assert_eq!(timeline["ok"], true, "{timeline}");
    assert_eq!(
        timeline["labels"].as_array().map(Vec::len),
        Some(6),
        "{timeline}"
    );
    assert_eq!(
        timeline["lines"][0], 5,
        "plate sketch on line 5: {timeline}"
    );
    assert_eq!(timeline["plate_buffers"], true, "{timeline}");
    let boss_rows = timeline["boss_rows"].as_array().expect("rows");
    let plate = boss_rows
        .iter()
        .find(|r| r["key"] == timeline["plate_key"])
        .expect("the plate keeps its key");
    assert_eq!(plate["buffers"], false, "known key not re-sent: {timeline}");
    assert_eq!(plate["sent"], false);
    assert!(
        boss_rows
            .iter()
            .any(|r| r["buffers"] == true && r["sent"] == true)
    );
    // PK-660 project through the worker.
    let project = &out["project"];
    assert_eq!(project["outcome"], "served", "{project}");
    assert_eq!(
        project["versions"],
        serde_json::json!([["001_plate", "built", 2], ["002_pin", "built", 2]])
    );
}
