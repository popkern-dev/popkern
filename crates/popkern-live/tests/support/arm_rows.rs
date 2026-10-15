//! Browser-harness rows for the six-axis arm's page features: the
//! section view (§4.6), the placement-only joint slider (§9.8) and a
//! large model's fit.

use serde_json::Value;

/// Every arm row.
pub fn assert_arm_rows(out: &Value) {
    assert_section(out);
    assert_clipping(out);
    assert_joint(out);
    assert_large_fit(out);
}

/// The section view (six-axis-arm §4.6): every material clipped, one
/// cap per visible body that draws real pixels, a move within a frame,
/// and "off" restoring the unclipped scene.
pub fn assert_section(out: &Value) {
    let s = &out["section"];
    assert_eq!(s["slider_shown"], true, "{s}");
    assert_eq!(s["clipped"], true, "{s}");
    // One cap per body the plane passes through (the other part lies
    // wholly beyond the cut: removed, nothing to cap).
    assert!(s["crossed"].as_u64().unwrap_or(0) >= 1, "{s}");
    assert_eq!(s["caps"], s["crossed"], "one cap per crossed body: {s}");
    assert!(
        s["cap_pixels"].as_u64().unwrap_or(0) > 100,
        "the caps draw: {s}"
    );
    assert!(
        s["move_ms"].as_f64().unwrap_or(f64::MAX) < 16.7,
        "a move fits a frame: {s}"
    );
    assert_eq!(s["cleared"], true, "{s}");
}

/// The clipping view (docs/architecture/clipping-view.md §9).
pub fn assert_clipping(out: &Value) {
    let c = &out["clipping"];
    let t = &c["tilted"];
    assert_eq!(
        t["clipped"], true,
        "a tilted plane clips every material: {c}"
    );
    assert!(t["crossed"].as_u64().unwrap_or(0) >= 1, "{c}");
    assert_eq!(t["caps"], t["crossed"], "one cap per crossed body: {c}");
    assert!(
        t["outline_segments"].as_u64().unwrap_or(0) > 8,
        "the cut is outlined: {c}"
    );
    assert!(
        c["flip_pixels"].as_u64().unwrap_or(0) > 500,
        "flip changes the view: {c}"
    );
    assert_eq!(
        c["any_keeps_probe"], false,
        "\"any\" removes the x>cx, y<cy quadrant: {c}"
    );
    assert_eq!(
        c["all_keeps_probe"], true,
        "\"all\" removes only the corner: {c}"
    );
    let x = &c["excluded"];
    assert_eq!(x["plain_unclipped"], true, "{c}");
    assert_eq!(x["housing_clipped"], true, "{c}");
    assert_eq!(x["caps_before"], x["crossed_before"], "{c}");
    assert_eq!(
        x["caps_after"].as_u64().unwrap_or(99) + 1,
        x["caps_before"].as_u64().unwrap_or(0),
        "the excluded part (crossed by the plane) loses its cap: {c}"
    );
    assert!(x["link"].as_str().unwrap_or("").contains("~x:"), "{c}");
    assert!(
        x["toggle_shown"].as_u64().unwrap_or(0) >= 2,
        "tree rows show the section toggle: {c}"
    );
    let g = &c["gizmo"];
    assert_eq!(g["grabbed"], true, "the arrow takes the pointer: {c}");
    assert_eq!(g["orbited"], false, "a handle drag does not orbit: {c}");
    let (moved, expected) = (
        g["moved"].as_f64().unwrap_or(0.0),
        g["expected"].as_f64().unwrap_or(1.0),
    );
    assert!(
        expected > 0.0 && (moved - expected).abs() < 0.05 * expected,
        "40 px along the arrow slides the plane 40 px worth: {c}"
    );
    let f = &c["facing"];
    assert_eq!(
        f["grabbed"], true,
        "the square takes a facing plane's drag: {c}"
    );
    assert!(
        f["away"].as_f64().unwrap_or(0.0) > 0.0,
        "dragging up pushes it away: {c}"
    );
    let l = &c["link"];
    assert!(l["written"].as_str().unwrap_or("").contains("~all"), "{c}");
    assert!(l["cleared_link"].is_null(), "clear drops #section: {c}");
    assert_eq!(l["restored"], 2, "{c}");
    assert_eq!(l["combine"], "all", "{c}");
    assert!(
        c["outline_ms"].as_f64().unwrap_or(f64::MAX) < 16.7,
        "the outline fits a frame: {c}"
    );
    assert_eq!(c["cleared"], true, "{c}");
}

/// A joint slider (§9.8, PK-725): an assembly parameter that drives no
/// part re-places occurrences without re-sending a mesh — the page's
/// geometries are the same objects, the link and its follower move,
/// the ground stays — and a step of the drag fits the wasm budget.
pub fn assert_joint(out: &Value) {
    let j = &out["joint"];
    assert_eq!(j["outcome"], "served", "{j}");
    assert_eq!(j["placed"], true, "a placement-only envelope: {j}");
    assert_eq!(j["sent_buffers"], 0, "{j}");
    assert_eq!(j["same_geometries"], true, "{j}");
    assert_eq!(j["base_moved"], false, "{j}");
    assert_eq!(j["link_moved"], true, "{j}");
    assert_eq!(j["rotor_moved"], true, "{j}");
    let value = j["slider_value"].as_f64().unwrap_or(f64::NAN);
    assert!((value - 40.0).abs() < 1e-9, "{j}");
    assert_eq!(
        j["row_kept"], true,
        "the panel row is updated in place: {j}"
    );
    assert_eq!(j["row_value"], j["slider_value"], "{j}");
    assert_eq!(j["bodies_kept"], 1, "{j}");
    assert_eq!(j["report_kept"], true, "{j}");
    // The fastest step is the step's cost: headless Chrome's software
    // GL adds frame backpressure that scales with machine load (the
    // median read 20 ms alone, 38 ms beside the whole suite).
    assert!(
        j["frame_ms"]["fastest"].as_f64().unwrap_or(f64::MAX) < 30.0,
        "a joint step within the wasm budget: {j}"
    );
}

/// A large model (two metres in millimetres) fits and draws: the depth
/// range follows the fitted size (a fixed far plane hid the arm). The
/// feature tree renders the arm's fin-root warning by its message.
pub fn assert_large_fit(out: &Value) {
    let f = &out["large_fit"];
    assert_eq!(f["outcome"], "served", "{f}");
    assert!(
        f["drawn_fraction"].as_f64().unwrap_or(0.0) > 0.05,
        "the model draws after a fit: {f}"
    );
    // The feature tree: a structured warning by its message, a label
    // as text.
    assert_eq!(f["report"]["warning"], "sliver retained", "{f}");
    assert_eq!(f["report"]["label"], "<b>fins</b>", "{f}");
    assert_eq!(f["report"]["markup"], 0, "{f}");
}
