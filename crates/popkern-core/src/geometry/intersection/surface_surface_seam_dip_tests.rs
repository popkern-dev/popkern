//! A seam dip is two opposite jumps with no net winding.

use std::f64::consts::TAU;

use super::dips_across_seam;
use crate::geometry::Point3;
use crate::geometry::intersection::marching::MarchSample;

fn trace(v1: &[f64]) -> Vec<MarchSample> {
    v1.iter()
        .enumerate()
        .map(|(i, &v)| MarchSample {
            u1: 1.0,
            v1: v,
            u2: 0.5,
            v2: TAU - 0.1,
            point: Point3::new(i as f64, 0.0, 0.0),
        })
        .collect()
}

const PERIODS: [Option<f64>; 4] = [Some(TAU), Some(TAU), None, Some(TAU)];

/// A closed trace dipping below the v = 0 seam and back.
#[test]
fn a_trace_dipping_below_the_seam_and_back_dips() {
    let s = trace(&[0.25, 0.1, TAU - 0.05, TAU - 0.2, TAU - 0.1, 0.05, 0.25]);
    assert!(dips_across_seam(&s, PERIODS));
}

/// A trace winding once around the axis is not a dip.
#[test]
fn a_winding_trace_does_not_dip() {
    let s = trace(&[0.3, 0.0001, TAU - 0.5, 4.0, 2.0, 0.6, 0.31]);
    assert!(!dips_across_seam(&s, PERIODS));
}

/// No jump, or a non-periodic axis: no dip.
#[test]
fn a_trace_without_jumps_or_period_does_not_dip() {
    assert!(!dips_across_seam(&trace(&[0.3, 0.2, 0.1, 0.2]), PERIODS));
    assert!(!dips_across_seam(&trace(&[0.1, TAU - 0.1, 0.1]), [None; 4]));
}

/// One sample registered a whole period away and straight back jumps
/// twice without crossing the seam: no dip.
#[test]
fn a_sample_on_another_period_image_does_not_dip() {
    let s = trace(&[0.83, 0.831 + TAU, 0.833]);
    assert!(!dips_across_seam(&s, PERIODS));
}
