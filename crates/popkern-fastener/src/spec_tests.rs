//! ThreadSpec parsing and derived-dimension pins.

use super::FastenerError;
use super::spec::{Handedness, ThreadSpec};
use super::standards::{self, ClearanceFit};

#[test]
fn parses_iso_coarse() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    assert!((spec.major_diameter() - 6.0).abs() < 1e-12);
    assert!((spec.pitch() - 1.0).abs() < 1e-12);
    assert_eq!(spec.handedness(), Handedness::Right);
    assert_eq!(spec.starts(), 1);
    let spec = ThreadSpec::parse("m8").expect("case-insensitive");
    assert!((spec.pitch() - 1.25).abs() < 1e-12);
}

#[test]
fn parses_explicit_pitch_and_unified() {
    let fine = ThreadSpec::parse("M6x0.75").expect("fine");
    assert!((fine.pitch() - 0.75).abs() < 1e-12);
    let quarter = ThreadSpec::parse("1/4-20").expect("unified fractional");
    assert!((quarter.major_diameter() - 6.35).abs() < 1e-9);
    assert!((quarter.pitch() - 1.27).abs() < 1e-9);
    let numbered = ThreadSpec::parse("#8-32").expect("unified numbered");
    assert!((numbered.major_diameter() - (0.060 + 0.013 * 8.0) * 25.4).abs() < 1e-9);
    assert!((numbered.pitch() - 25.4 / 32.0).abs() < 1e-9);
}

#[test]
fn unknown_designations_name_the_grammars() {
    for bad in ["M7", "garbage", "M6y1", "1/0-20", "#-32"] {
        let err = ThreadSpec::parse(bad).expect_err(bad);
        let text = err.to_string();
        assert!(
            text.contains(bad.trim()),
            "error for '{bad}' names the input: {text}"
        );
    }
    // The off-table metric error tells the user the explicit-pitch form.
    let err = ThreadSpec::parse("M7").expect_err("M7 off-table");
    assert!(
        err.to_string().contains("M7x"),
        "suggests M7x<pitch>: {err}"
    );
}

#[test]
fn iso_derived_dimensions() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    // H = (√3/2)·P at 60°; minor = major − (5/8)·2H·... = 6 − 1.0825·P.
    let h = spec.fundamental_height();
    assert!((h - 0.866_025_403_784_438_6).abs() < 1e-12);
    assert!((spec.minor_diameter() - (6.0 - 1.082_531_754_730_548_3)).abs() < 1e-12);
    assert!((spec.lead() - 1.0).abs() < 1e-12);
    let two_start = spec.clone().with_starts(2);
    assert!((two_start.lead() - 2.0).abs() < 1e-12);
}

#[test]
fn validate_rejects_degenerate_specs() {
    assert!(matches!(
        ThreadSpec::custom(-1.0, 1.0).validate(),
        Err(FastenerError::InvalidSpec {
            field: "major_diameter",
            ..
        })
    ));
    assert!(matches!(
        ThreadSpec::custom(6.0, 0.0).validate(),
        Err(FastenerError::InvalidSpec { field: "pitch", .. })
    ));
    assert!(matches!(
        ThreadSpec::custom(6.0, 1.0).with_starts(0).validate(),
        Err(FastenerError::InvalidSpec {
            field: "starts",
            ..
        })
    ));
    // Depth reaching the axis: tiny major with a huge pitch.
    assert!(matches!(
        ThreadSpec::custom(1.0, 2.0).validate(),
        Err(FastenerError::ProfileOverlap { .. })
    ));
}

#[test]
fn standards_tables_and_derivations() {
    assert_eq!(standards::iso_coarse_pitch(10.0), Some(1.5));
    assert_eq!(standards::iso_coarse_pitch(7.0), None);
    let hex = standards::hex_head(6.0);
    assert!((hex.across_flats - 10.0).abs() < 1e-12);
    assert!((hex.height - 4.0).abs() < 1e-12);
    // Off-table sizes derive complete proportions.
    let odd = standards::hex_head(7.0);
    assert!((odd.across_flats - 11.2).abs() < 1e-9);
    let socket = standards::socket_head(6.0);
    assert!((socket.socket_across_flats - 5.0).abs() < 1e-12);
    assert!((standards::hex_nut_thickness(6.0) - 5.2).abs() < 1e-12);
    assert!((standards::clearance_diameter(6.0, ClearanceFit::Normal) - 6.6).abs() < 1e-12);
    assert!((standards::clearance_diameter(7.0, ClearanceFit::Loose) - 8.4).abs() < 1e-9);
    assert!((standards::clearance_diameter(6.0, ClearanceFit::Custom(6.9)) - 6.9).abs() < 1e-12);
}
