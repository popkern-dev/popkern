//! The page harness rows for the time limit and the pre-built gallery
//! examples (split from `page.rs` for its line budget).

use serde_json::Value;

/// The time limit: presets up to "no limit", default 24 s (20 s kernel
/// budget + the page's grace), "no limit" lifts both budgets. Pre-built
/// examples, when the bundle was baked (`tools/bake-playground.sh`):
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
    // (tools/bake-playground.sh): shown without the kernel, certified,
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
        assert_eq!(returning["editor_is_default"], true, "{returning}");
        assert_eq!(returning["editor_restored"], true, "{returning}");
        assert_eq!(returning["notice_after"], false, "{returning}");
        let cheap = &out["prebuilt_cheap"];
        assert_eq!(cheap["made_live"], true, "{cheap}");
        assert!(
            cheap["status_after_live"]
                .as_str()
                .is_some_and(|s| s.starts_with("served:")),
            "{cheap}"
        );
    } else {
        eprintln!(
            "[bakes] SKIP the pre-built rows: no web/playground/bakes (run tools/bake-playground.sh)"
        );
    }
}
