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
