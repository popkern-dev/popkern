//! Cut-domain placement (PK-459 stage 2): seam-axis classification and
//! the v fold of a doubly-wrapped face's runs.

use std::f64::consts::TAU;

use super::chart_cut_domain::fold_into_cut_domain;
use super::chart_embedding::{ChartEmbedding, EmbeddingSource};
use super::chart_seam_axes::{ChartAxis, constant_axis};
use super::chart_subdivision::FaceChartSubdivision;

fn run(uv: Vec<(f64, f64)>, seam_like: bool) -> ChartEmbedding {
    ChartEmbedding {
        uv,
        source: EmbeddingSource::OwnPcurve,
        seam_like,
        refined_first: None,
        refined_last: None,
    }
}

#[test]
fn seam_runs_classify_by_their_constant_axis() {
    let eps = TAU * 1.0e-6;
    let meridian = [(0.0, 0.0), (0.0, 3.0), (0.0, TAU)];
    let latitude = [(0.0, 1.0), (3.0, 1.0), (TAU, 1.0)];
    let record = [(0.5, 0.5), (0.6, 0.7)];
    assert_eq!(constant_axis(&meridian, eps, eps), Some(ChartAxis::U));
    assert_eq!(constant_axis(&latitude, eps, eps), Some(ChartAxis::V));
    assert_eq!(constant_axis(&record, eps, eps), None);
}

#[test]
fn v_seam_halves_take_opposite_borders_and_records_fold_by_interior() {
    let mut sub = FaceChartSubdivision {
        sigma: 1.0,
        ..Default::default()
    };
    // The v-seam pair at v = 0: u-increasing and u-decreasing halves,
    // both authored on the low copy.
    sub.embeddings
        .insert(0, run(vec![(0.0, 0.0), (TAU, 0.0)], true));
    sub.embeddings
        .insert(1, run(vec![(TAU, 0.0), (0.0, 0.0)], true));
    // A drill ring's lower half crossing the v-seam: endpoints ON the
    // seam (low copy), interior just below it (v < 0 ≡ just under 2π).
    sub.embeddings
        .insert(2, run(vec![(0.66, 0.0), (0.70, -0.2), (0.74, 0.0)], false));
    // Its upper half, interior above the seam.
    sub.embeddings
        .insert(3, run(vec![(0.74, 0.0), (0.70, 0.2), (0.66, 0.0)], false));
    let domain = fold_into_cut_domain(&mut sub, 0.0, TAU, TAU, 1.0);

    assert_eq!(
        domain.v_seam.iter().copied().collect::<Vec<_>>(),
        vec![0, 1]
    );
    // σ = +1: u-increasing travel bounds the interior from the LOW
    // border, u-decreasing from the HIGH one.
    assert!(sub.embeddings[&0].uv.iter().all(|q| q.1 == 0.0));
    assert!(
        sub.embeddings[&1]
            .uv
            .iter()
            .all(|q| (q.1 - TAU).abs() < 1e-12)
    );
    // The lower half moved whole to the top of the domain: interior
    // inside it, endpoints on the high copy.
    let lower = &sub.embeddings[&2].uv;
    assert!((lower[1].1 - (TAU - 0.2)).abs() < 1e-12, "{lower:?}");
    assert!((lower[0].1 - TAU).abs() < 1e-12 && (lower[2].1 - TAU).abs() < 1e-12);
    // The upper half stays on the low copy.
    assert_eq!(sub.embeddings[&3].uv[1].1, 0.2);

    // σ = −1 mirrors the border rule.
    let mut mirrored = FaceChartSubdivision::default();
    mirrored
        .embeddings
        .insert(0, run(vec![(0.0, 0.0), (TAU, 0.0)], true));
    fold_into_cut_domain(&mut mirrored, 0.0, TAU, TAU, -1.0);
    assert!(
        mirrored.embeddings[&0]
            .uv
            .iter()
            .all(|q| (q.1 - TAU).abs() < 1e-12)
    );
}
