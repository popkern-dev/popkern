use super::*;
use crate::merge::walk::chart_embedding::EmbeddingSource;

fn run(uv: Vec<(f64, f64)>) -> ChartEmbedding {
    ChartEmbedding {
        uv,
        source: EmbeddingSource::OwnPcurve,
        seam_like: false,
        refined_first: None,
        refined_last: None,
    }
}

/// Two lobe halves leaving a pole: both start along −v (the meridian
/// direction every run at a pole has) and end at one node, at
/// different azimuths: distinct.
#[test]
fn runs_leaving_a_pole_at_different_azimuths_are_distinct() {
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};
    let a = run(vec![(3.0 * FRAC_PI_4, PI), (2.40, 2.5), (PI, FRAC_PI_2)]);
    let b = run(vec![(FRAC_PI_4, PI), (0.74, 2.5), (0.0, FRAC_PI_2)]);
    assert!(!chart_runs_overlay(&a, &b));
}

/// A record registered twice, its samples a rounding apart: one
/// segment.
#[test]
fn a_run_and_its_copy_overlay() {
    let a = run(vec![(0.0, 0.0), (0.5, 0.25), (1.0, 1.0)]);
    let b = run(vec![
        (0.0, 0.0),
        (0.25, 0.0625),
        (0.5, 0.25 + 1e-9),
        (1.0, 1.0),
    ]);
    assert!(chart_runs_overlay(&a, &b));
}
