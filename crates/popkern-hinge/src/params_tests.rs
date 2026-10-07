// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use popkern_fastener::ThreadSpec;

use crate::HingeError;
use crate::params::{HingeComponent, HingeHoles, HingeKind, HingeParams, hinge_layout};

#[test]
fn default_butt_hinge_lays_out_five_knuckles() {
    let p = HingeParams::new(60.0);
    let layout = hinge_layout(HingeKind::Butt, &p).expect("layout");
    assert_eq!(layout.segments.len(), 5);
    // Rolled leaf: R = pin/2 + clearance + thickness = 2 + 0 + 2.
    assert!((layout.knuckle_radius - 4.0).abs() < 1e-12);
    assert!((layout.bore_radius - 2.0).abs() < 1e-12);
    assert_eq!(layout.owned(HingeComponent::LeafA).len(), 3);
    assert_eq!(layout.owned(HingeComponent::LeafB).len(), 2);
    let (first, last) = (layout.segments[0], layout.segments[4]);
    assert_eq!(first.0, 0.0);
    assert!((last.1 - 60.0).abs() < 1e-9);
    // Exact fit still keeps the notch off the knuckle cylinder.
    assert!(layout.notch > layout.knuckle_radius);
}

#[test]
fn every_decline_names_its_field_and_limit() {
    let base = HingeParams::new(60.0);
    let decline = |kind, p: &HingeParams| hinge_layout(kind, p).expect_err("declines");
    assert!(matches!(
        decline(
            HingeKind::Butt,
            &HingeParams {
                knuckles: 1,
                ..base.clone()
            }
        ),
        HingeError::InvalidSpec {
            field: "knuckles",
            ..
        }
    ));
    // Print-in-place needs leaf B between two of leaf A's knuckles.
    assert!(matches!(
        decline(
            HingeKind::PrintInPlace,
            &HingeParams {
                knuckles: 2,
                ..HingeParams::print_in_place(60.0)
            }
        ),
        HingeError::InvalidSpec {
            field: "knuckles",
            ..
        }
    ));
    assert!(matches!(
        decline(
            HingeKind::Butt,
            &HingeParams {
                knuckles: 41,
                gap: 1.5,
                ..base.clone()
            }
        ),
        HingeError::KnucklesDoNotFit { knuckles: 41, .. }
    ));
    assert!(matches!(
        decline(
            HingeKind::Butt,
            &HingeParams {
                knuckle_diameter: Some(3.9),
                ..base.clone()
            }
        ),
        HingeError::NoKnuckleWall { .. }
    ));
    assert!(matches!(
        decline(
            HingeKind::Butt,
            &HingeParams {
                leaf_width: 4.0,
                ..base.clone()
            }
        ),
        HingeError::PlateTooNarrow { .. }
    ));
    assert!(matches!(
        decline(HingeKind::PrintInPlace, &HingeParams::new(60.0)),
        HingeError::PrintClearance {
            field: "clearance",
            ..
        }
    ));
    assert!(matches!(
        decline(
            HingeKind::PrintInPlace,
            &HingeParams {
                gap: 0.0,
                ..HingeParams::print_in_place(60.0)
            }
        ),
        HingeError::PrintClearance { field: "gap", .. }
    ));
    assert!(matches!(
        decline(
            HingeKind::Butt,
            &HingeParams {
                angle_deg: 360.0,
                ..base.clone()
            }
        ),
        HingeError::InvalidSpec { field: "angle", .. }
    ));
    let m8 = ThreadSpec::parse("M8").expect("M8");
    assert!(matches!(
        decline(
            HingeKind::Butt,
            &HingeParams {
                holes: Some(HingeHoles {
                    per_leaf: 2,
                    thread: m8,
                    countersink: true
                }),
                ..base.clone()
            }
        ),
        HingeError::HoleDoesNotFit {
            dimension: "plate width",
            ..
        }
    ));
    let text = decline(
        HingeKind::Butt,
        &HingeParams {
            knuckles: 41,
            gap: 1.5,
            ..base
        },
    )
    .to_string();
    assert!(
        text.contains("41 knuckles") && text.contains("fewer knuckles"),
        "{text}"
    );
}

#[test]
fn holes_sit_centred_on_the_plate_beyond_the_notch() {
    let p = HingeParams {
        holes: Some(HingeHoles {
            per_leaf: 3,
            thread: ThreadSpec::parse("M3").expect("M3"),
            countersink: true,
        }),
        ..HingeParams::new(60.0)
    };
    let layout = hinge_layout(HingeKind::Butt, &p).expect("layout");
    assert_eq!(layout.holes.len(), 3);
    let expected = layout.notch + (p.leaf_width - layout.notch) / 2.0;
    for (k, hole) in layout.holes.iter().enumerate() {
        assert!((hole.offset - expected).abs() < 1e-12);
        assert!((hole.z - 20.0 * (k as f64 + 0.5)).abs() < 1e-9);
    }
}

#[test]
fn spellings_round_trip() {
    for kind in [HingeKind::Butt, HingeKind::PrintInPlace] {
        assert_eq!(HingeKind::parse(kind.spelling()), Some(kind));
    }
    for c in [
        HingeComponent::LeafA,
        HingeComponent::LeafB,
        HingeComponent::Pin,
    ] {
        assert_eq!(HingeComponent::parse(c.spelling()), Some(c));
    }
    assert_eq!(HingeKind::parse("print"), None);
}
