// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use super::{Gallery, Manifest, PROJECT_FORMAT};
use crate::error::ProjectError;

fn line_of(text: &str) -> (usize, String) {
    match Manifest::parse(text) {
        Err(ProjectError::ManifestInvalid { line, reason, .. }) => (line, reason),
        other => panic!("expected ManifestInvalid, got {other:?}"),
    }
}

#[test]
fn a_full_manifest_round_trips_through_text() {
    let text = r#"
# A project.
[project]
name = "Bracket \"A\""   # quoted
format = 1
base = "imports/vendor.step"

[gallery]
title = "Bracket"
blurb = "Two # signs: # inside a string stay"
thumbnail = "thumb.png"
"#;
    let manifest = Manifest::parse(text).expect("parses");
    assert_eq!(manifest.name, "Bracket \"A\"");
    assert_eq!(manifest.format, PROJECT_FORMAT);
    assert_eq!(manifest.base.as_deref(), Some("imports/vendor.step"));
    let gallery = manifest.gallery.clone().expect("gallery");
    assert_eq!(
        gallery.blurb.as_deref(),
        Some("Two # signs: # inside a string stay")
    );
    assert_eq!(
        Manifest::parse(&manifest.to_text()).expect("reparses"),
        manifest
    );
    assert_eq!(gallery.demonstrates, Gallery::default().demonstrates);
}

#[test]
fn mistakes_name_their_line() {
    assert_eq!(line_of("[project]\nname = \"x\"\nnmae = \"y\"\n").0, 3);
    assert_eq!(line_of("[project]\nname = x\n").0, 2);
    assert_eq!(line_of("name = \"x\"\n").0, 1);
    assert_eq!(line_of("[project]\nname = \"x\"\n[extras]\n").0, 3);
    assert_eq!(line_of("[project]\nname = \"x\"\nname = \"y\"\n").0, 3);
    assert_eq!(line_of("[project]\nname = \"x\"\nformat = 99\n").0, 3);
    let (line, reason) = line_of("[gallery]\ntitle = \"t\"\n");
    assert_eq!(line, 0);
    assert!(reason.contains("[project]"), "{reason}");
}
