use super::{DimensionKind, Unit};

const EPSILON: f64 = 1e-12;

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() < EPSILON
}

// ── Round-trip tests ──────────────────────────────────────────────────────────

#[test]
fn meters_round_trip() {
    let v = 1.234;
    assert!(approx_eq(Unit::Meters.from_si(Unit::Meters.to_si(v)), v));
}

#[test]
fn millimeters_round_trip() {
    let v = 500.0;
    assert!(approx_eq(
        Unit::Millimeters.from_si(Unit::Millimeters.to_si(v)),
        v
    ));
}

#[test]
fn centimeters_round_trip() {
    let v = 42.5;
    assert!(approx_eq(
        Unit::Centimeters.from_si(Unit::Centimeters.to_si(v)),
        v
    ));
}

#[test]
fn inches_round_trip() {
    let v = 12.0;
    assert!(approx_eq(Unit::Inches.from_si(Unit::Inches.to_si(v)), v));
}

#[test]
fn radians_round_trip() {
    let v = std::f64::consts::PI / 4.0;
    assert!(approx_eq(Unit::Radians.from_si(Unit::Radians.to_si(v)), v));
}

#[test]
fn degrees_round_trip() {
    let v = 90.0;
    assert!(approx_eq(Unit::Degrees.from_si(Unit::Degrees.to_si(v)), v));
}

#[test]
fn ratio_round_trip() {
    let v = 0.75;
    assert!(approx_eq(Unit::Ratio.from_si(Unit::Ratio.to_si(v)), v));
}

// ── Known-value conversion tests ─────────────────────────────────────────────

#[test]
fn one_inch_is_0_0254_meters() {
    assert!(approx_eq(Unit::Inches.to_si(1.0), 0.0254));
}

#[test]
fn one_meter_is_1000_mm() {
    assert!(approx_eq(Unit::Millimeters.from_si(1.0), 1000.0));
}

#[test]
fn one_meter_is_100_cm() {
    assert!(approx_eq(Unit::Centimeters.from_si(1.0), 100.0));
}

#[test]
fn ninety_degrees_is_pi_over_2() {
    assert!(approx_eq(
        Unit::Degrees.to_si(90.0),
        std::f64::consts::PI / 2.0
    ));
}

#[test]
fn pi_radians_is_180_degrees() {
    assert!(approx_eq(
        Unit::Degrees.from_si(std::f64::consts::PI),
        180.0
    ));
}

// ── kind() tests ──────────────────────────────────────────────────────────────

#[test]
fn meters_kind_is_length() {
    assert_eq!(Unit::Meters.kind(), DimensionKind::Length);
}

#[test]
fn millimeters_kind_is_length() {
    assert_eq!(Unit::Millimeters.kind(), DimensionKind::Length);
}

#[test]
fn centimeters_kind_is_length() {
    assert_eq!(Unit::Centimeters.kind(), DimensionKind::Length);
}

#[test]
fn inches_kind_is_length() {
    assert_eq!(Unit::Inches.kind(), DimensionKind::Length);
}

#[test]
fn radians_kind_is_angle() {
    assert_eq!(Unit::Radians.kind(), DimensionKind::Angle);
}

#[test]
fn degrees_kind_is_angle() {
    assert_eq!(Unit::Degrees.kind(), DimensionKind::Angle);
}

#[test]
fn ratio_kind_is_dimensionless() {
    assert_eq!(Unit::Ratio.kind(), DimensionKind::Dimensionless);
}

// ── suffix() tests ────────────────────────────────────────────────────────────

#[test]
fn suffix_meters() {
    assert_eq!(Unit::Meters.suffix(), "m");
}

#[test]
fn suffix_millimeters() {
    assert_eq!(Unit::Millimeters.suffix(), "mm");
}

#[test]
fn suffix_centimeters() {
    assert_eq!(Unit::Centimeters.suffix(), "cm");
}

#[test]
fn suffix_inches() {
    assert_eq!(Unit::Inches.suffix(), "in");
}

#[test]
fn suffix_radians() {
    assert_eq!(Unit::Radians.suffix(), "rad");
}

#[test]
fn suffix_degrees() {
    assert_eq!(Unit::Degrees.suffix(), "deg");
}

#[test]
fn suffix_ratio() {
    assert_eq!(Unit::Ratio.suffix(), "");
}
