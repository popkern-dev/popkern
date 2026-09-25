//! Unit tests for the chart-orientation repair (PK-249).

use super::*;

/// A domed 3-patch ring in the three-quad layout: patch k's quad has
/// corners `V_{k-1}` at (0,0), the hex mid-vertex `M_k` at (0,1), the
/// ring center `C` at (1,0), and `V_k` at (1,1) — so the `u=1` row is
/// spoke k (C→V_k) and the `v=0` column is spoke k−1 reversed
/// (V_{k-1}→C), matching the ring builder's conventions.
///
/// Every patch is the SAME smooth surface — the dome
/// `z = 0.25 (1 − x² − y²)` over its bilinear quad — so every spoke is
/// exactly G1 (PK-539: the former fixture was bilinear quads under
/// a linear dome, which crease at every spoke; the coplanarity metric
/// read those creases as G1). The lifted map has degree (2, 2), so the
/// degree-9 grid interpolating it at 10 × 10 Chebyshev nodes
/// reproduces it exactly.
fn healthy_ring() -> Vec<RectangularGregoryBiNonic> {
    let rim = |angle_deg: f64| {
        let a = angle_deg.to_radians();
        (a.cos(), a.sin())
    };
    let nodes: Vec<f64> = (0..10)
        .map(|a| 0.5 - 0.5 * ((2 * a + 1) as f64 * std::f64::consts::PI / 20.0).cos())
        .collect();
    let inverse = collocation_inverse(&nodes);
    let mut ring = Vec::new();
    for k in 0..3 {
        let v_prev = rim(120.0 * (k as f64) - 120.0);
        let m = rim(120.0 * (k as f64) - 60.0);
        let v = rim(120.0 * (k as f64));
        let surface = |u: f64, s: f64| {
            let x = (1.0 - u) * (1.0 - s) * v_prev.0 + (1.0 - u) * s * m.0 + u * s * v.0;
            let y = (1.0 - u) * (1.0 - s) * v_prev.1 + (1.0 - u) * s * m.1 + u * s * v.1;
            Point3::new(x, y, 0.25 * (1.0 - x * x - y * y))
        };
        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
        for (i, row) in grid.iter_mut().enumerate() {
            for (jj, cell) in row.iter_mut().enumerate() {
                let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
                for (a, &na) in nodes.iter().enumerate() {
                    for (b, &nb) in nodes.iter().enumerate() {
                        let weight = inverse[i][a] * inverse[jj][b];
                        let p = surface(na, nb);
                        x += weight * p.x;
                        y += weight * p.y;
                        z += weight * p.z;
                    }
                }
                *cell = Point3::new(x, y, z);
            }
        }
        ring.push(RectangularGregoryBiNonic { grid });
    }
    ring
}

/// Inverse of the degree-9 Bernstein collocation matrix at `nodes`
/// (`inv[i][a]`: control point i from node value a), by Gauss–Jordan
/// with partial pivoting.
fn collocation_inverse(nodes: &[f64]) -> [[f64; 10]; 10] {
    let binomial =
        |k: usize| -> f64 { (0..k).fold(1.0, |acc, m| acc * (9 - m) as f64 / (m + 1) as f64) };
    let mut m = [[0.0_f64; 20]; 10];
    for (a, &t) in nodes.iter().enumerate() {
        for (i, entry) in m[a].iter_mut().take(10).enumerate() {
            *entry = binomial(i) * t.powi(i as i32) * (1.0 - t).powi(9 - i as i32);
        }
        m[a][10 + a] = 1.0;
    }
    for col in 0..10 {
        let pivot = (col..10)
            .max_by(|&x, &y| m[x][col].abs().total_cmp(&m[y][col].abs()))
            .unwrap_or(col);
        m.swap(col, pivot);
        let scale = m[col][col];
        for entry in m[col].iter_mut() {
            *entry /= scale;
        }
        let pivot_row = m[col];
        for (row, r) in m.iter_mut().enumerate() {
            if row != col {
                let factor = r[col];
                for (x, p) in r.iter_mut().zip(pivot_row.iter()) {
                    *x -= factor * p;
                }
            }
        }
    }
    let mut inv = [[0.0_f64; 10]; 10];
    for (row, source) in inv.iter_mut().zip(m.iter()) {
        row.copy_from_slice(&source[10..20]);
    }
    inv
}

/// Fold one interior row past its rail on the OUTER side (`u = 0`,
/// the rim V_{k-1}→M_k) so the natural normal reverses in a band there,
/// away from the spokes.
fn banded_ring() -> Vec<RectangularGregoryBiNonic> {
    let mut ring = healthy_ring();
    for j in 3..=6 {
        let boundary = ring[1].grid[0][j];
        let inner = ring[1].grid[1][j];
        // Reflect row 1 across the boundary row — the chart's
        // u-tangent reverses locally while the boundary stays put.
        ring[1].grid[1][j] = Point3::new(
            2.0 * boundary.x - inner.x,
            2.0 * boundary.y - inner.y,
            2.0 * boundary.z - inner.z,
        );
    }
    ring
}

/// The same fold on the SPOKE side (`u = 1`): patch 1 folds back
/// across its shared spoke with patch 2.
fn spoke_folded_ring() -> Vec<RectangularGregoryBiNonic> {
    let mut ring = healthy_ring();
    for j in 3..=6 {
        let boundary = ring[1].grid[9][j];
        let inner = ring[1].grid[8][j];
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

/// PK-539: a fold across a spoke is not repairable by an interior
/// re-chart (the spoke is a shared curve), and the repair must decline
/// it — the coplanarity gate read the folded result as G1 and the
/// repair "succeeded" with the fold intact.
#[test]
fn a_fold_across_a_spoke_is_declined_not_passed() {
    let as_gregory = |ring: Vec<RectangularGregoryBiNonic>| -> Vec<GregoryPatch> {
        ring.into_iter()
            .map(|p| GregoryPatch::RectangularBiNonic(Box::new(p)))
            .collect()
    };
    let spoke_sine = super::super::three_quad_kp_promote::worst_ring_normal_sine;
    assert!(
        spoke_sine(&as_gregory(healthy_ring())) < 1e-9,
        "the dome ring is G1"
    );
    let folded = spoke_folded_ring();
    assert_eq!(
        spoke_sine(&as_gregory(folded.clone())),
        1.0,
        "the fold reads as one"
    );
    let result = repair_chart_orientation(&folded, None);
    assert!(
        matches!(result, Err(ref reason) if reason.starts_with("g1(")),
        "a spoke fold must be declined on G1, got {:?}",
        result.map(|o| o.post_g1)
    );
}
