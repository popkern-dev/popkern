//! Browser-harness rows for the model view's corner and the phone
//! layout (PK-867): the standard views, and films on a phone.

use serde_json::Value;

/// Every view row.
pub fn assert_view_rows(out: &Value) {
    assert_standard_views(out);
    assert_phone(out);
    assert_projects_ui(out);
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

/// On a phone (PK-867, PK-868) a project's film is listed on the
/// History sheet. A film opens as a stage over the model view and the
/// sheet stays where it was; the stage has exactly one Close and one
/// Play, a bar inside the model view under the film's frame, one
/// chapter per step, a timeline that seeks where it is clicked, and the
/// corner tools step aside. Preparing shows what it waits on, and a
/// film opened while a build runs says so, waits, then plays. Closing
/// lists it again; the status chip clears the corner's tools and the
/// view modes keep their clustered card.
pub fn assert_phone(out: &Value) {
    let p = &out["phone"];
    let open = &p["open"];
    assert_eq!(open["active"], true, "{p}");
    assert_eq!(
        open["sheet"], "code",
        "the stage leaves the sheet alone: {p}"
    );
    assert_eq!(open["bar_shown"], true, "{p}");
    assert_eq!(open["bar_in_pane"], true, "{p}");
    assert_eq!(open["frame_above_bar"], true, "{p}");
    assert_eq!(open["closes"], 1, "one way out: {p}");
    assert_eq!(open["plays"], 1, "one play button: {p}");
    assert_eq!(open["chapters"], open["steps"], "{p}");
    assert_eq!(open["corner_hidden"], true, "{p}");
    assert_eq!(open["played"], true, "{p}");
    let seek = open["seek_fraction"].as_f64().unwrap_or(-1.0);
    assert!(
        (seek - 0.75).abs() < 0.08,
        "a timeline click seeks there: {p}"
    );
    assert_eq!(open["clipped"], true, "{p}");
    assert_eq!(open["frame_in_pane"], true, "{p}");
    assert_eq!(open["hole_in_pane"], true, "{p}");
    let named = |row: &Value, stage: &str| {
        row["stages"].as_array().is_some_and(|a| {
            a.iter()
                .any(|s| s.as_str().is_some_and(|s| s.starts_with(stage)))
        })
    };
    assert!(
        named(open, "Running the film script"),
        "preparing names its stage: {p}"
    );
    assert_eq!(p["closed"]["active"], false, "{p}");
    assert_eq!(p["closed"]["stage_gone"], true, "{p}");
    assert_eq!(p["closed"]["list_shown"], true, "{p}");
    assert_eq!(p["closed"]["list_items"], 1, "{p}");
    assert_eq!(
        p["during_build"]["active"], true,
        "a film opened mid-build plays: {p}"
    );
    assert!(
        named(&p["during_build"], "Building the model"),
        "it says it waits on the build: {p}"
    );
    assert_eq!(p["status_clear_of_tools"], true, "{p}");
    assert_eq!(p["modes_clustered"], true, "{p}");
    assert_eq!(p["orient_shown"], true, "{p}");
}

/// Themes switch and are remembered, Scratch carries the default
/// model's film, the open Scratch deletes (another project opens), and
/// a film selected the moment its example opens plays on that example
/// (PK-869: it played on the previous example's model, whose targets
/// the new film's script could not name).
pub fn assert_projects_ui(out: &Value) {
    let p = &out["projects_ui"];
    let t = &p["themes"];
    assert_eq!(t["options"], 3, "{p}");
    assert_eq!(t["stored"], "light", "{p}");
    assert_ne!(
        t["light"], t["default"],
        "the light theme changes the ground: {p}"
    );
    assert_eq!(t["back"], true, "{p}");
    assert_eq!(p["scratch_film"], true, "{p}");
    assert_eq!(p["deleted"]["scratch_gone"], true, "{p}");
    assert_eq!(p["deleted"]["reopened"], true, "{p}");
    assert_eq!(p["deleted"]["fresh_film"], true, "{p}");
    let race = p["race"].as_array().expect("race rows");
    assert_eq!(race.len(), 2, "{p}");
    for row in race {
        assert_eq!(row["active"], true, "{p}");
        assert!(row["overlay"].is_null(), "{p}");
    }
    assert_eq!(
        race[1]["kind"], "part",
        "the door hinge's own document: {p}"
    );
    assert_eq!(race[1]["bodies"], 3, "{p}");
    // The door hinge's leaf turns by t = 4 (the section sweep moves no
    // body): a seek moves bodies on screen.
    assert_eq!(
        race[1]["moved"], true,
        "a seek moves the film's bodies: {p}"
    );
}

/// One module owns `localStorage`, behind a guarded accessor; the only
/// other reader is the pre-paint theme script (a classic script cannot
/// import a module), which reads the theme key inside a try and never
/// writes.
pub fn assert_storage_owned<S: AsRef<str>>(own: &[S], read: impl Fn(&str) -> String) {
    let store = read("local-store.js");
    assert!(store.contains("localStorage"));
    assert!(store.contains("function storage(action) {\n  try {"));
    let boot = read("theme-boot.js");
    assert!(
        boot.contains("try {") && !boot.contains("setItem"),
        "{boot}"
    );
    for file in own.iter().map(AsRef::as_ref) {
        if file == "local-store.js" || file == "theme-boot.js" {
            continue;
        }
        assert!(
            !read(file).contains("localStorage"),
            "{file} reaches localStorage outside local-store.js"
        );
    }
}
