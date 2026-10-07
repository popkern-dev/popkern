use super::*;
use std::f64::consts::TAU;

const PERIODS: Periods = (TAU, TAU);

/// A chart polyline.
type Polyline = Vec<(f64, f64)>;

/// One single-run loop per polyline.
fn loops_of(polys: &[&Polyline]) -> Vec<ChartLoop> {
    polys
        .iter()
        .map(|p| ChartLoop::from_runs([p.as_slice()]))
        .collect()
}

/// A full ring at height `v`, sampled at 16 stations, walked toward
/// +u (`forward`) or −u, starting on the period branch `u0`.
fn ring(v: f64, forward: bool, u0: f64) -> Polyline {
    (0..16)
        .map(|k| {
            let t = TAU * k as f64 / 16.0;
            (if forward { u0 + t } else { u0 - t }, v)
        })
        .collect()
}

/// The band v in [1, 2] around the ring direction: material LEFT of
/// travel under σ = +1 means the lower rim walks +u, the upper −u.
fn band() -> (Polyline, Polyline) {
    (ring(1.0, true, 0.0), ring(2.0, false, 0.0))
}

#[test]
fn a_band_around_the_torus_contains_its_own_rows_only() {
    let (lower, upper) = band();
    let loops = loops_of(&[&lower, &upper]);
    assert!(bounds_a_region(&loops, PERIODS));
    let r = interior_reference(&loops, 1.0, PERIODS).expect("reference");
    assert!(r.1 > 1.0 && r.1 < 2.0, "reference {r:?} inside the band");
    for u in [0.3, 2.0, 4.0, 6.0, -5.0, 19.0] {
        assert!(same_side(&loops, r, (u, 1.5), PERIODS), "u {u} in band");
        assert!(!same_side(&loops, r, (u, 3.0), PERIODS), "u {u} above");
        assert!(!same_side(&loops, r, (u, 0.5), PERIODS), "u {u} below");
        // Below the band through the tube seam: v wraps to 2π − 0.2.
        assert!(
            !same_side(&loops, r, (u, -0.2), PERIODS),
            "u {u} across the seam"
        );
    }
}

/// The complement of the band is the same two rims walked the other
/// way: the walk's handedness, not the rims, says which side it is.
#[test]
fn reversed_rims_bound_the_complement() {
    let (mut lower, mut upper) = band();
    lower.reverse();
    upper.reverse();
    let loops = loops_of(&[&lower, &upper]);
    let r = interior_reference(&loops, 1.0, PERIODS).expect("reference");
    assert!(!same_side(&loops, r, (1.0, 1.5), PERIODS));
    assert!(same_side(&loops, r, (1.0, 3.5), PERIODS));
    assert!(same_side(&loops, r, (1.0, 0.5), PERIODS));
    // A left-handed chart (σ = −1) reads the reversed rims as the band.
    let r = interior_reference(&loops, -1.0, PERIODS).expect("reference");
    assert!(same_side(&loops, r, (1.0, 1.5), PERIODS));
    assert!(!same_side(&loops, r, (1.0, 3.5), PERIODS));
}

/// Rims authored on drifting period branches (the PK-685 torus-groove
/// polygon spanned u in [−12.6, 6.3]) are the same curves on the torus:
/// each run sits on its own branch, and the junctions glue them.
#[test]
fn branch_drift_between_runs_changes_nothing() {
    let (lower, upper) = band();
    let lower_runs: Vec<Polyline> = lower
        .chunks(5)
        .enumerate()
        .map(|(k, c)| c.iter().map(|&(u, v)| (u - TAU * k as f64, v)).collect())
        .collect();
    let upper_runs: Vec<Polyline> = vec![upper.iter().map(|&(u, v)| (u - 2.0 * TAU, v)).collect()];
    let loops = loops_of(&[&lower, &upper]);
    let drifted = vec![
        ChartLoop::from_runs(lower_runs.iter().map(|r| r.as_slice())),
        ChartLoop::from_runs(upper_runs.iter().map(|r| r.as_slice())),
    ];
    assert!(bounds_a_region(&drifted, PERIODS));
    let r = interior_reference(&loops, 1.0, PERIODS).expect("reference");
    let rd = interior_reference(&drifted, 1.0, PERIODS).expect("reference");
    for x in [(0.3, 1.5), (4.0, 3.0), (-9.4, 1.2), (2.0, 0.4)] {
        assert_eq!(
            same_side(&loops, r, x, PERIODS),
            same_side(&drifted, rd, x, PERIODS),
            "{x:?}"
        );
    }
}

/// A two-sample run along the tube seam spanning HALF a period (the
/// half-torus band's seam piece, v from π to 2π) is read as lifted, not
/// as its ambiguous short representative.
#[test]
fn a_half_period_endpoint_run_keeps_its_direction() {
    use std::f64::consts::PI;
    let down: Polyline = (0..8).map(|k| (0.0, TAU - PI * k as f64 / 8.0)).collect();
    let rim: Polyline = (0..8).map(|k| (TAU * k as f64 / 8.0, PI)).collect();
    let up: Polyline = vec![(TAU, PI), (TAU, TAU)];
    let top: Polyline = (0..8).map(|k| (TAU - TAU * k as f64 / 8.0, TAU)).collect();
    let band = vec![ChartLoop::from_runs([
        down.as_slice(),
        rim.as_slice(),
        up.as_slice(),
        top.as_slice(),
    ])];
    assert!(bounds_a_region(&band, PERIODS));
    let r = interior_reference(&band, 1.0, PERIODS).expect("reference");
    assert!(same_side(&band, r, (3.0, 4.5), PERIODS));
    assert!(!same_side(&band, r, (3.0, 1.5), PERIODS));
}

/// A contractible disc and the torus with that disc removed.
#[test]
fn a_disc_and_its_complement() {
    let disc: Vec<(f64, f64)> = vec![(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)];
    let loops = loops_of(&[&disc]);
    assert!(bounds_a_region(&loops, PERIODS));
    let r = interior_reference(&loops, 1.0, PERIODS).expect("reference");
    assert!(same_side(&loops, r, (1.5, 1.5), PERIODS));
    assert!(!same_side(&loops, r, (4.0, 4.0), PERIODS));
    assert!(same_side(&loops, r, (1.5 + TAU, 1.5 - TAU), PERIODS));
    let hole: Vec<(f64, f64)> = disc.iter().rev().copied().collect();
    let loops = loops_of(&[&hole]);
    let r = interior_reference(&loops, 1.0, PERIODS).expect("reference");
    assert!(!same_side(&loops, r, (1.5, 1.5), PERIODS));
    assert!(same_side(&loops, r, (4.0, 4.0), PERIODS));
    assert!(same_side(&loops, r, (1.5, 4.0), PERIODS));
}

/// One ring alone bounds nothing on the torus: the decision declines.
#[test]
fn a_lone_non_contractible_ring_bounds_no_region() {
    let lower = ring(1.0, true, 0.0);
    let loops = loops_of(&[&lower]);
    assert!(!bounds_a_region(&loops, PERIODS));
}

/// A path through a sample vertex counts the crossing once.
#[test]
fn a_path_through_a_vertex_counts_once() {
    let (lower, upper) = band();
    let loops = loops_of(&[&lower, &upper]);
    // Stations sit at u = kτ/16; the vertical path at u = τ/16 meets a
    // station of each rim exactly.
    let r = (TAU / 16.0, 1.5);
    assert!(same_side(&loops, r, (TAU / 16.0, 1.6), PERIODS));
    assert!(!same_side(&loops, r, (TAU / 16.0, 2.5), PERIODS));
    assert!(!same_side(&loops, r, (TAU / 16.0, 0.5), PERIODS));
}

#[test]
fn the_widest_material_arc_lies_in_the_region() {
    // The band minus a disc hole straddling u = 3: the row v = 1.5
    // meets the hole's sides, so its material arcs are the two pieces
    // of the ring row either side of the hole, the wider one first.
    let (lower, upper) = band();
    let hole: Vec<(f64, f64)> = vec![(2.5, 1.2), (2.5, 1.8), (3.5, 1.8), (3.5, 1.2)];
    let loops = loops_of(&[&lower, &upper, &hole]);
    let r = interior_reference(&loops, 1.0, PERIODS).expect("reference");
    let u = widest_material_arc(&loops, r, 1.5, PERIODS).expect("arc");
    assert!(!(2.5..=3.5).contains(&u), "sample {u} in the hole");
    assert!(same_side(&loops, r, (u, 1.5), PERIODS));
    // A row outside the band has no material arc.
    assert_eq!(widest_material_arc(&loops, r, 3.0, PERIODS), None);
}

/// The torus-groove shape (PK-685): a band v in [1, 2] whose loop also
/// walks the u = 0 tube seam piece inside the band in and back out.
/// The rims' walk directions are arbitrary (self-loop claims), yet the
/// slit alone places the reference in the band.
#[test]
fn a_seam_slit_places_the_reference_in_the_band_whatever_the_rim_directions() {
    let slit: Vec<(f64, f64)> = (0..=8).map(|k| (0.0, 1.0 + k as f64 / 8.0)).collect();
    for (forward_lower, forward_upper) in [(true, false), (true, true), (false, false)] {
        let lower = ring(1.0, forward_lower, 0.0);
        let upper = ring(2.0, forward_upper, 0.0);
        let back: Vec<(f64, f64)> = slit.iter().rev().map(|&(u, v)| (u + TAU, v)).collect();
        let loops = loops_of(&[&lower, &upper, &slit, &back]);
        let r = slit_reference(&loops, &slit, PERIODS).expect("reference");
        assert!(r.1 > 1.0 && r.1 < 2.0, "reference {r:?}");
        assert!(same_side(&loops, r, (3.0, 1.5), PERIODS));
        assert!(!same_side(&loops, r, (3.0, 3.5), PERIODS));
        assert!(!same_side(&loops, r, (3.0, 0.5), PERIODS));
    }
}
