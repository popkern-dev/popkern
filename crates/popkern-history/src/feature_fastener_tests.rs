//! Fastener feature serialization tests.

use crate::feature::Feature;
use crate::feature_fastener::{FastenerKind, ThreadDesc};
use crate::parameter::ParameterId;

fn nut(corner_chamfer: bool, mouth_chamfer: bool) -> Feature {
    Feature::Fastener {
        kind: FastenerKind::HexNut,
        thread: ThreadDesc {
            major_diameter: ParameterId::default(),
            pitch: ParameterId::default(),
            left_handed: false,
            starts: 1,
        },
        length: ParameterId::default(),
        thread_length: None,
        clearance: None,
        tip_chamfer: true,
        corner_chamfer,
        mouth_chamfer,
    }
}

#[test]
fn nut_chamfer_toggles_round_trip() {
    let json = serde_json::to_string(&nut(false, false)).unwrap_or_default();
    let back: Result<Feature, _> = serde_json::from_str(&json);
    let again = back
        .ok()
        .map(|f| serde_json::to_string(&f).unwrap_or_default());
    assert_eq!(again, Some(json));
}

#[test]
fn documents_without_the_toggles_keep_both_chamfers() {
    // A pre-PK-686 writer never emitted the two keys.
    let Ok(mut value) = serde_json::to_value(nut(false, false)) else {
        panic!("serializes");
    };
    let Some(fields) = value.get_mut("Fastener").and_then(|v| v.as_object_mut()) else {
        panic!("externally tagged Fastener variant: {value}");
    };
    fields.remove("corner_chamfer");
    fields.remove("mouth_chamfer");
    let back: Result<Feature, _> = serde_json::from_value(value);
    let as_json = |f: &Feature| serde_json::to_value(f).ok();
    assert_eq!(
        back.ok().as_ref().and_then(as_json),
        as_json(&nut(true, true))
    );
}
