//! The page harness's Simulate row (PK-948): the live boss plate held
//! on its underside and pushed on its boss top, solved from the panel.

use serde_json::Value;

pub fn assert_simulate(out: &Value) {
    let s = &out["simulate"];
    assert_eq!(s["status"], "solved", "{s}");
    assert!(s["materials"].as_u64().unwrap_or(0) >= 6, "{s}");
    assert_eq!(s["overlay"], true, "the field is drawn: {s}");
    assert_eq!(s["bodies_hidden"], true, "the bodies give way to it: {s}");
    assert_eq!(s["legend"], true, "{s}");
    assert!(
        s["legend_max"].as_str().is_some_and(|t| t.ends_with("MPa")),
        "{s}"
    );
    let rows: Vec<&str> = s["rows"]
        .as_array()
        .map(|r| r.iter().filter_map(Value::as_str).collect())
        .unwrap_or_default();
    for wanted in ["max deflection", "peak stress", "safety factor", "material"] {
        assert!(
            rows.iter().any(|r| r.starts_with(wanted)),
            "{wanted} missing from {rows:?}"
        );
    }
    assert_eq!(
        s["fields"],
        serde_json::json!([
            "displacement_mm",
            "max_principal_mpa",
            "von_mises_mpa",
            "zz_error_mpa"
        ]),
        "{s}"
    );
    assert_eq!(s["cleared"], true, "a re-show clears the field: {s}");
    eprintln!(
        "[pk948] headless simulate on the boss plate: {:.0} ms",
        s["solve_ms"].as_f64().unwrap_or(-1.0)
    );
}
