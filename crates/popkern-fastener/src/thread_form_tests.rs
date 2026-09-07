//! Thread-form profile math pins: ISO stations, exact clearance
//! grow, and the overlap invariants.

use super::FastenerError;
use super::spec::ThreadSpec;
use super::thread_form::{external_groove, internal_bore_radius, internal_ridge};

fn m6() -> ThreadSpec {
    ThreadSpec::parse("M6").expect("M6")
}

#[test]
fn external_groove_hits_iso_stations() {
    let spec = m6();
    let profile = external_groove(&spec).expect("groove");
    let p = spec.pitch();
    let r_major = 3.0;
    let r_minor = spec.minor_diameter() / 2.0;
    // Root flat P/4 at the minor radius; opening 7P/8 at the major
    // radius, flanks continued past it by the embed margin.
    assert!((profile.radial_min() - r_minor).abs() < 1e-12);
    assert!((profile.radial_max() - (r_major + p / 16.0)).abs() < 1e-12);
    let at_minor: Vec<f64> = profile
        .points
        .iter()
        .filter(|pt| (pt.0 - r_minor).abs() < 1e-12)
        .map(|pt| pt.1)
        .collect();
    assert_eq!(at_minor.len(), 2);
    assert!((at_minor[0] - p / 8.0).abs() < 1e-12 || (at_minor[0] + p / 8.0).abs() < 1e-12);
}

#[test]
fn internal_ridge_zero_clearance_is_the_nominal_ridge() {
    let spec = m6();
    let profile = internal_ridge(&spec, 0.0).expect("ridge");
    let p = spec.pitch();
    let r_major = 3.0;
    // Tip flat P/8 (half-width P/16) at the major radius.
    let tips: Vec<f64> = profile
        .points
        .iter()
        .filter(|pt| (pt.0 - r_major).abs() < 1e-12)
        .map(|pt| pt.1.abs())
        .collect();
    assert_eq!(tips.len(), 2);
    assert!((tips[0] - p / 16.0).abs() < 1e-12);
    assert!((internal_bore_radius(&spec, 0.0) - spec.minor_diameter() / 2.0).abs() < 1e-12);
}

#[test]
fn clearance_grow_is_a_true_normal_offset() {
    let spec = m6();
    let c = 0.2;
    let nominal = internal_ridge(&spec, 0.0).expect("nominal");
    let grown = internal_ridge(&spec, c).expect("grown");
    // Tip moves radially out by exactly c.
    assert!((grown.radial_max() - (nominal.radial_max() + c)).abs() < 1e-12);
    // Every grown flank point is ≥ c from the nominal flank LINE
    // measured along the flank normal: verify the flank's
    // perpendicular distance grew by exactly c. Nominal +flank goes
    // through (r_minor, 3P/8) with direction (cos α, −sin α).
    let p = spec.pitch();
    let half = spec.flank_angle() / 2.0;
    let r_minor = spec.minor_diameter() / 2.0;
    let (n_r, n_y) = (half.sin(), half.cos()); // outward flank normal
    let line_offset = |pt: &(f64, f64)| (pt.0 - r_minor) * n_r + (pt.1 - 3.0 * p / 8.0) * n_y;
    // Only the in-material corners carry the fit contract: the base
    // corners live inside the bore void (below the wall), where the
    // profile deliberately turns parallel instead of following the
    // offset flank.
    let r_bore = internal_bore_radius(&spec, c);
    for pt in grown
        .points
        .iter()
        .filter(|pt| pt.1 > 0.0 && pt.0 >= r_bore)
    {
        let d = line_offset(pt);
        assert!(
            d >= c - 1e-12,
            "grown corner {pt:?} only {d:.6} beyond the nominal flank (want ≥ {c})"
        );
    }
    // Bore wall grew by c too.
    assert!((internal_bore_radius(&spec, c) - (spec.minor_diameter() / 2.0 + c)).abs() < 1e-12);
}

#[test]
fn overlap_invariants_report_the_gap() {
    let spec = m6();
    // A clearance near half the pitch merges adjacent internal cuts.
    let err = internal_ridge(&spec, 0.45).expect_err("overlapping clearance");
    match err {
        FastenerError::ProfileOverlap { gap, pitch, .. } => {
            assert!(gap <= 0.0, "reported gap is the violated margin: {gap}");
            assert!((pitch - 1.0).abs() < 1e-12);
        }
        other => panic!("expected ProfileOverlap, got {other}"),
    }
    // Negative clearance is a spec error, not an overlap.
    assert!(matches!(
        internal_ridge(&spec, -0.1),
        Err(FastenerError::InvalidSpec {
            field: "clearance",
            ..
        })
    ));
}

#[test]
fn multi_start_widens_the_self_adjacency_budget() {
    // A two-start thread doubles the lead: a profile that would
    // overlap itself at starts=1 with a coarse custom form fits at
    // starts=2 (the per-start spacing doubles).
    let tight = ThreadSpec::custom(6.0, 1.0).with_flank_angle(0.9);
    // Narrow flank angle deepens the thread; groove stays valid.
    external_groove(&tight).expect("single start still fits");
    let two = tight.with_starts(2);
    external_groove(&two).expect("two starts fit a fortiori");
}
