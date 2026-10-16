//! Browser-harness rows for the model view's corner and the phone
//! layout (PK-867): the standard views, and films on a phone.

use serde_json::Value;

/// Every view row.
pub fn assert_view_rows(out: &Value) {
    assert_standard_views(out);
    assert_phone(out);
}

/// Each standard view looks along its axis (Z up; the front looks
/// along +Y), and the pad offers all seven.
pub fn assert_standard_views(out: &Value) {
    let v = &out["views"];
    assert_eq!(v["pad_buttons"], 7, "{v}");
    let expect = [
        ("front", [0.0, 1.0, 0.0]),
        ("back", [0.0, -1.0, 0.0]),
        ("left", [1.0, 0.0, 0.0]),
        ("right", [-1.0, 0.0, 0.0]),
        ("top", [0.0, 0.0, -1.0]),
        ("bottom", [0.0, 0.0, 1.0]),
    ];
    for (name, dir) in expect {
        let got: Vec<f64> = v["look"][name]
            .as_array()
            .unwrap_or_else(|| panic!("{name}: {v}"))
            .iter()
            .map(|x| x.as_f64().unwrap_or(f64::NAN))
            .collect();
        for (g, d) in got.iter().zip(dir) {
            assert!(
                (g - d).abs() < 2e-3,
                "{name} looks along {got:?}, not {dir:?}"
            );
        }
    }
    // Isometric: from the front-right-top octant, toward the model.
    let iso: Vec<f64> = v["look"]["iso"]
        .as_array()
        .expect("iso")
        .iter()
        .map(|x| x.as_f64().unwrap_or(f64::NAN))
        .collect();
    let third = 1.0 / 3.0_f64.sqrt();
    for (g, d) in iso.iter().zip([-third, third, -third]) {
        assert!((g - d).abs() < 2e-3, "iso looks along {iso:?}");
    }
}

/// On a phone a project's film is listed on the History sheet; opening
/// it brings that sheet up, shows the viewport's own play button (which
/// plays), keeps the frame's shade inside the model view; closing it
/// lists it again. The status chip clears the corner's tools and the
/// view modes keep their clustered card.
pub fn assert_phone(out: &Value) {
    let p = &out["phone"];
    let open = &p["open"];
    assert_eq!(open["active"], true, "{p}");
    assert_eq!(
        open["sheet"], "history",
        "a film opening shows its controls: {p}"
    );
    assert_eq!(open["bar_shown"], true, "{p}");
    assert_eq!(open["play_view_shown"], true, "{p}");
    assert_eq!(open["played"], true, "the viewport play button plays: {p}");
    assert_eq!(open["list_hidden"], true, "{p}");
    assert_eq!(open["clipped"], true, "{p}");
    assert_eq!(open["frame_in_pane"], true, "{p}");
    assert_eq!(open["hole_in_pane"], true, "{p}");
    assert_eq!(p["closed"]["active"], false, "{p}");
    assert_eq!(p["closed"]["list_shown"], true, "{p}");
    assert_eq!(p["closed"]["list_items"], 1, "{p}");
    assert_eq!(p["status_clear_of_tools"], true, "{p}");
    assert_eq!(p["modes_clustered"], true, "{p}");
    assert_eq!(p["orient_shown"], true, "{p}");
}
