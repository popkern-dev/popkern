//! Unit tests for the chart-orientation repair (PK-249).

use super::*;

/// A gently domed 3-patch ring in the three-quad layout: patch k's
/// quad has corners `V_{k-1}` at (0,0), the hex mid-vertex `M_k` at
/// (0,1), the ring center `C` at (1,0), and `V_k` at (1,1) — so the
/// `u=1` row is spoke k (C→V_k) and the `v=0` column is spoke k−1
/// reversed (V_{k-1}→C), matching the ring builder's conventions.
/// Grids are bilinear over the quad with a mild dome so charts are
/// healthy by construction. Adjacent patches share the spoke
/// boundary rows exactly (both bilinear maps interpolate the same
/// spoke segment), like any built ring.
fn healthy_ring() -> Vec<RectangularGregoryBiNonic> {
    let rim = |angle_deg: f64| {
        let a = angle_deg.to_radians();
        Point3::new(a.cos(), a.sin(), 0.0)
    };
    let c = Point3::new(0.0, 0.0, 0.25);
    let mut ring = Vec::new();
    for k in 0..3 {
        let v_prev = rim(120.0 * (k as f64) - 120.0);
        let m = rim(120.0 * (k as f64) - 60.0);
        let v = rim(120.0 * (k as f64));
        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
        for (i, row) in grid.iter_mut().enumerate() {
            let u = (i as f64) / 9.0;
            for (j, cell) in row.iter_mut().enumerate() {
                let vpar = (j as f64) / 9.0;
                let x = (1.0 - u) * (1.0 - vpar) * v_prev.x
                    + (1.0 - u) * vpar * m.x
                    + u * (1.0 - vpar) * c.x
                    + u * vpar * v.x;
                let y = (1.0 - u) * (1.0 - vpar) * v_prev.y
                    + (1.0 - u) * vpar * m.y
                    + u * (1.0 - vpar) * c.y
                    + u * vpar * v.y;
                let z = (1.0 - u) * (1.0 - vpar) * v_prev.z
                    + (1.0 - u) * vpar * m.z
                    + u * (1.0 - vpar) * c.z
                    + u * vpar * v.z;
                *cell = Point3::new(x, y, z);
            }
        }
        ring.push(RectangularGregoryBiNonic { grid });
    }
    ring
}

/// Fold one interior row past its rail so the natural normal
/// reverses in a band along the `u = 1` side.
fn banded_ring() -> Vec<RectangularGregoryBiNonic> {
    let mut ring = healthy_ring();
    for j in 3..=6 {
        let boundary = ring[1].grid[9][j];
        let inner = ring[1].grid[8][j];
        // Reflect row 8 across the boundary row — the chart's
        // u-tangent reverses locally while the boundary stays put.
        ring[1].grid[8][j] = Point3::new(
            2.0 * boundary.x - inner.x,
            2.0 * boundary.y - inner.y,
            2.0 * boundary.z - inner.z,
        );
    }
    ring
}

#[test]
fn healthy_ring_counts_zero_bands_and_declines_repair() {
    let ring = healthy_ring();
    assert_eq!(binonic_reversal_band_count(&ring), 0);
    let result = repair_chart_orientation(&ring, None);
    assert!(
        matches!(result, Err(ref reason) if reason == "already_clean"),
        "clean rings must decline as already_clean"
    );
}

#[test]
fn banded_ring_is_detected() {
    let ring = banded_ring();
    assert!(
        binonic_reversal_band_count(&ring) > 0,
        "the reflected row must read as a reversal band"
    );
}

#[test]
fn repair_clears_bands_and_keeps_boundaries_bit_identical() {
    let ring = banded_ring();
    let outcome =
        repair_chart_orientation(&ring, None).expect("the synthetic band must be repairable");
    assert!(outcome.pre_bands > 0);
    assert_eq!(binonic_reversal_band_count(&outcome.ring), 0);
    // The OUTER boundary rows are not degrees of freedom:
    // bit-identical. The spoke rows may move, but only as one
    // shared curve with fixed endpoints.
    for (k, (before, after)) in ring.iter().zip(outcome.ring.iter()).enumerate() {
        for j in 0..10 {
            assert_eq!(before.grid[0][j].x, after.grid[0][j].x, "row0 moved");
            assert_eq!(before.grid[j][9].x, after.grid[j][9].x, "col9 moved");
        }
        assert_eq!(before.grid[9][0].x, after.grid[9][0].x, "C moved");
        assert_eq!(before.grid[9][9].x, after.grid[9][9].x, "V moved");
        let next = &outcome.ring[(k + 1) % 3];
        for j in 0..10 {
            let a = after.grid[9][j];
            let b = next.grid[9 - j][0];
            assert!(
                a.distance_to(&b) < 1e-12,
                "spoke {k} shared-curve identity broken at j {j}"
            );
        }
    }
    // Continuity did not silently collapse: the ring G1 stays
    // within the repair's own gate.
    assert!(
        outcome.post_g1 <= 5.0e-2,
        "repaired ring G1 {:.3e} above the commercial gate",
        outcome.post_g1
    );
}

#[test]
fn interior_relocation_is_rejected_not_accepted() {
    // A ring whose mid-face reference itself is fine but where the
    // repair target is unreachable because the boundary forces the
    // fold: collapse the outer boundary row to a single point so
    // every chart with that boundary is degenerate along it.
    let mut ring = healthy_ring();
    let pinch = ring[1].grid[0][0];
    for j in 0..10 {
        ring[1].grid[0][j] = pinch;
    }
    for j in 3..=6 {
        let boundary = ring[1].grid[9][j];
        let inner = ring[1].grid[8][j];
        ring[1].grid[8][j] = Point3::new(
            2.0 * boundary.x - inner.x,
            2.0 * boundary.y - inner.y,
            2.0 * boundary.z - inner.z,
        );
    }
    // Whatever the solver does, the result must be a decline or a
    // genuinely clean candidate — never an accepted candidate with
    // residual violations.
    match repair_chart_orientation(&ring, None) {
        Ok(outcome) => {
            assert_eq!(binonic_reversal_band_count(&outcome.ring), 0);
        }
        Err(reason) => {
            assert!(!reason.is_empty());
        }
    }
}
