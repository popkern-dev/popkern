//! The print filter against Langelaar 2017 §3.1's pattern tests (solid
//! features on the bed print as drawn, unsupported ones vanish, a
//! ledge stops at 45°) and its adjoint against central differences.

use super::*;

/// A vertical slab of `nx` × `nz` unit cells, one cell deep in y, with
/// entry (i, k) at index i + nx k.
fn slab(nx: usize, nz: usize) -> Vec<DesignEntry> {
    let mut out = Vec::new();
    for k in 0..nz {
        for i in 0..nx {
            let corner = [i as f64, 0.0, k as f64];
            out.push(DesignEntry {
                solid: 0,
                entry: out.len(),
                centre: corner.map(|c| c + 0.5),
                volume: 1.0,
                corner,
            });
        }
    }
    out
}

const UP: BuildDirection = BuildDirection {
    axis: 2,
    negative: false,
};

/// Blueprint of the pattern: a solid column at i = 1 from the bed to
/// the top, a one-layer ledge from it at k = 5 out to i = 6, and a
/// floating 2 × 2 block at i = 8–9, k = 3–4.
fn pattern(nx: usize, nz: usize) -> Vec<f64> {
    let mut x = vec![0.0; nx * nz];
    for k in 0..nz {
        x[1 + nx * k] = 1.0;
    }
    for i in 2..=6 {
        x[i + nx * 5] = 1.0;
    }
    for k in 3..=4 {
        for i in 8..=9 {
            x[i + nx * k] = 1.0;
        }
    }
    x
}

#[test]
fn the_pattern_prints_as_the_paper_shows() {
    let (nx, nz) = (11, 8);
    let filter = PrintFilter::new(&slab(nx, nz), 1.0, UP);
    let x = pattern(nx, nz);
    let printed = filter.apply(&x, &vec![0.0; x.len()]);
    let at = |i: usize, k: usize| printed.densities[i + nx * k];
    // The column prints solid (the smooth min's error under 1%).
    for k in 0..nz {
        assert!((at(1, k) - 1.0).abs() < 0.01, "column {k}: {}", at(1, k));
    }
    // The ledge's first cell leans 45° off the column; the rest of it
    // overhangs and is not printed.
    assert!(at(2, 5) > 0.99, "{}", at(2, 5));
    for i in 3..=6 {
        assert!(at(i, 5) < 0.01, "ledge {i}: {}", at(i, 5));
    }
    // The floating block is not printed.
    for k in 3..=4 {
        for i in 8..=9 {
            assert!(at(i, k) < 0.01, "block {i} {k}: {}", at(i, k));
        }
    }
    // Empty space stays empty.
    assert!(at(5, 2) < 1e-3);
}

/// Built the other way up (the slab's top on the bed), the column still
/// prints and the ledge's underside now faces up: the cells over the
/// ledge's free end hang from nothing and the ledge itself prints only
/// near the column.
#[test]
fn the_build_direction_flips_the_bed() {
    let (nx, nz) = (11, 8);
    let down = BuildDirection {
        axis: 2,
        negative: true,
    };
    let filter = PrintFilter::new(&slab(nx, nz), 1.0, down);
    let mut x = pattern(nx, nz);
    // A strut on the ledge's tip, up to the top layer: built downwards
    // from the bed it hangs the ledge's tip from it.
    for k in 6..nz {
        x[6 + nx * k] = 1.0;
    }
    let printed = filter.apply(&x, &vec![0.0; x.len()]);
    let at = |i: usize, k: usize| printed.densities[i + nx * k];
    for k in 0..nz {
        assert!((at(1, k) - 1.0).abs() < 0.01);
    }
    assert!(at(6, nz - 1) > 0.99);
    // The tip is held by the strut above (now below), and its neighbour
    // by the tip's 45°.
    assert!(
        at(6, 5) > 0.99 && at(5, 5) > 0.99,
        "{} {}",
        at(6, 5),
        at(5, 5)
    );
    assert!(at(4, 5) < 0.01, "{}", at(4, 5));
}

/// Kept material prints regardless and holds up what stands on it.
#[test]
fn kept_material_prints_and_supports() {
    let (nx, nz) = (11, 8);
    let filter = PrintFilter::new(&slab(nx, nz), 1.0, UP);
    let mut x = pattern(nx, nz);
    // A post on the floating block.
    for k in 5..nz {
        x[8 + nx * k] = 1.0;
    }
    let mut floors = vec![0.0; x.len()];
    for k in 3..=4 {
        for i in 8..=9 {
            floors[i + nx * k] = 1.0;
        }
    }
    let printed = filter.apply(&x, &floors);
    let at = |i: usize, k: usize| printed.densities[i + nx * k];
    assert_eq!(at(8, 3), 1.0);
    for k in 5..nz {
        assert!(at(8, k) > 0.99, "post {k}: {}", at(8, k));
    }
}

/// The adjoint (eqs 13–14) against central differences of
/// f = Σ w_i ξ_i on a grey blueprint, every entry checked.
#[test]
fn the_adjoint_matches_central_differences() {
    let (nx, nz) = (6, 5);
    let entries = slab(nx, nz);
    let filter = PrintFilter::new(&entries, 1.0, UP);
    let n = entries.len();
    let x: Vec<f64> = (0..n)
        .map(|i| 0.15 + 0.7 * ((7 * i + 3) % 11) as f64 / 10.0)
        .collect();
    let w: Vec<f64> = (0..n).map(|i| 1.0 + ((5 * i) % 7) as f64).collect();
    let mut floors = vec![0.0; n];
    floors[2 + nx * 2] = 0.95;
    let f = |x: &[f64]| -> f64 {
        filter
            .apply(x, &floors)
            .densities
            .iter()
            .zip(&w)
            .map(|(a, b)| a * b)
            .sum()
    };
    let printed = filter.apply(&x, &floors);
    assert!(printed.held().iter().any(|&h| h));
    let grad = filter.back(&x, &printed, &w);
    let delta = 1e-6;
    for j in 0..n {
        let mut up = x.clone();
        up[j] += delta;
        let mut down = x.clone();
        down[j] -= delta;
        let fd = (f(&up) - f(&down)) / (2.0 * delta);
        assert!(
            (fd - grad[j]).abs() <= 1e-6 * (1.0 + fd.abs()),
            "entry {j}: {fd} against {}",
            grad[j]
        );
    }
}

/// A nearly empty blueprint (densities about 1e-8, as void reaches
/// under a sharp projection): Σ ξᴾ is subnormal, and the adjoint stays
/// finite and exact. (Raising that sum to 1/Q − 1 overflowed to
/// infinite sensitivities, which the update turned to NaN.)
#[test]
fn the_adjoint_is_finite_over_a_nearly_empty_support() {
    let (nx, nz) = (5, 4);
    let entries = slab(nx, nz);
    let filter = PrintFilter::new(&entries, 1.0, UP);
    let n = entries.len();
    let x: Vec<f64> = (0..n)
        .map(|i| 1e-8 * (1.0 + (i % 7) as f64 * 0.03))
        .collect();
    let floors = vec![0.0; n];
    let w: Vec<f64> = (0..n).map(|i| 1.0 + (i % 3) as f64).collect();
    let printed = filter.apply(&x, &floors);
    let grad = filter.back(&x, &printed, &w);
    assert!(grad.iter().all(|g| g.is_finite()), "{grad:?}");
    let f = |x: &[f64]| -> f64 {
        filter
            .apply(x, &floors)
            .densities
            .iter()
            .zip(&w)
            .map(|(a, b)| a * b)
            .sum()
    };
    let delta = 1e-11;
    for j in 0..n {
        let mut up = x.clone();
        up[j] += delta;
        let mut down = x.clone();
        down[j] -= delta;
        let fd = (f(&up) - f(&down)) / (2.0 * delta);
        assert!(
            (fd - grad[j]).abs() <= 1e-5 * (1.0 + fd.abs()),
            "entry {j}: {fd} against {}",
            grad[j]
        );
    }
}
