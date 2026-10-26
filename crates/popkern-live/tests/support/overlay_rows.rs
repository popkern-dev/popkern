//! The page harness's overlay rows (split from `page_rows.rs` for its
//! line budget).

use serde_json::Value;

/// The overlay rows: a kernel decline's class, marked line and
/// remediation, and a generative design declined unsolved, re-solved
/// with its progress on the status line, then served.
pub fn assert_overlay_rows(out: &Value) {
    // PK-956: a generative design declines unsolved with re-solve on
    // the overlay; re-solve shows its progress and serves the design.
    assert_eq!(
        out["generative"]["deferred"]["outcome"], "kernel_decline",
        "{out}"
    );
    assert_eq!(out["generative"]["deferred"]["kind"], "NeedsSolve", "{out}");
    assert_eq!(out["generative"]["resolve_visible"], true, "{out}");
    let resolved = &out["generative"]["resolved"];
    assert!(
        resolved["progress"]
            .as_str()
            .is_some_and(|t| t.contains("update")),
        "{resolved}"
    );
    assert_eq!(resolved["live"], true, "{resolved}");
    assert_eq!(resolved["overlay_hidden"], true, "{resolved}");
    assert!(
        resolved["row"].as_str().is_some_and(|t| t.contains("kept")),
        "{resolved}"
    );
    // PK-959: a film evolves the re-solved design from its history:
    // the first record shows less than the final design, the space
    // faint throughout; the body hidden while it evolves, shown after,
    // when the layer has faded; closing takes the layer away.
    let film = &out["generative"]["film"];
    assert_eq!(film["ok"], true, "{film}");
    assert!(film["records"].as_u64().unwrap_or(0) >= 10, "{film}");
    let triangles = |at: &str| film[at]["triangles"].as_u64().unwrap_or(0);
    assert!(triangles("end") > 0, "{film}");
    assert_ne!(triangles("start"), triangles("end"), "{film}");
    for at in ["start", "end"] {
        assert_eq!(film[at]["carved"], true, "{at}: {film}");
        assert_eq!(film[at]["body"], false, "{at}: {film}");
        assert!(
            film[at]["space"]
                .as_f64()
                .is_some_and(|o| o > 0.0 && o < 0.5),
            "{at}: {film}"
        );
    }
    assert_eq!(film["after"]["carved"], false, "{film}");
    assert_eq!(film["after"]["body"], true, "{film}");
    assert_eq!(film["in_scene"], true, "{film}");
    assert_eq!(film["removed"], true, "{film}");
    // A kernel decline: class, position and remediation.
    assert_eq!(out["decline"]["outcome"], "kernel_decline");
    assert_eq!(out["decline"]["kind"], "FeatureFailed");
    assert_eq!(out["decline"]["overlay_class"], "kernel decline");
    assert_eq!(out["decline"]["marked_line"], 3);
    assert!(
        out["decline"]["overlay_text"]
            .as_str()
            .is_some_and(|t| t.contains("ReduceRadius")),
        "remediation shown: {}",
        out["decline"]["overlay_text"]
    );
}
