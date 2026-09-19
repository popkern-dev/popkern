//! Projection-seed budget tests (PK-461) — split from `tests.rs`
//! under the 500-line file budget.

use crate::geometry::Point3;
use crate::geometry::nurbs_surface::NurbsSurface;

/// PK-461: the projection-seed subsample is budgeted on the TOTAL
/// pair count, not a blind 48 per axis — a coiled wall with few
/// u-spans and hundreds of v-spans must keep (near) full v-span seed
/// coverage, or the coarse ranking walks the adjacent turn's basin
/// and every projection consumer misreads on-surface points as ~a
/// pitch off (the phantom that declined a valid 25-turn thread
/// serve).
#[test]
fn projection_seeds_spend_pair_budget_on_the_dense_axis() {
    // 1 u-span (degree 1), 300 v-spans (degree 1, uniform knots).
    let num_v = 301;
    let mut knots_v = vec![0.0];
    for i in 0..num_v {
        knots_v.push(i as f64 / (num_v - 1) as f64);
    }
    knots_v.push(1.0);
    let mut control_points = Vec::new();
    let mut weights = Vec::new();
    for iu in 0..2 {
        for iv in 0..num_v {
            control_points.push(Point3::new(iu as f64, iv as f64, 0.0));
            weights.push(1.0);
        }
    }
    let surface = NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v,
        control_points,
        num_u: 2,
        num_v,
        weights,
    };
    let seeds = surface.projection_seeds().expect("many-span chart seeds");
    let mut vs: Vec<f64> = seeds.iter().map(|&(_, v)| v).collect();
    vs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    vs.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
    assert!(
        vs.len() >= 250,
        "v-axis seed coverage crushed to {} distinct values (want every span while the \
         total pair budget allows)",
        vs.len()
    );
    // The total stays within the square budget (+ borders).
    assert!(seeds.len() <= (48 * 48) + 2 * (vs.len() + 4));
}
