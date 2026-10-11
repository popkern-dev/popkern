// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]
use super::rectangular_binonic::RectangularGregoryBiNonic;
use super::rectangular_binonic_dual_twist::RectangularGregoryBiNonicDualTwist;
use crate::geometry::Point3;

/// A curved plain bi-nonic fixture: paraboloid-sampled grid.
fn curved_plain() -> RectangularGregoryBiNonic {
    let mut patch = RectangularGregoryBiNonic::zeros();
    for i in 0..10 {
        for j in 0..10 {
            let x = i as f64 / 9.0;
            let y = j as f64 / 9.0;
            patch.grid[i][j] = Point3::new(x, y, 0.3 * x * x + 0.2 * y * y - 0.25 * x * y);
        }
    }
    patch
}

#[test]
fn coincident_pairs_reproduce_the_plain_patch() {
    let plain = curved_plain();
    let dual = RectangularGregoryBiNonicDualTwist::from_plain(&plain);
    for iu in 0..=8 {
        for iv in 0..=8 {
            let (u, v) = (iu as f64 / 8.0, iv as f64 / 8.0);
            let a = plain.point_at(u, v);
            let b = dual.point_at(u, v);
            assert!(
                a.distance_to(&b) < 1.0e-12,
                "coincident-pair dual twist must equal the plain patch at ({u}, {v})"
            );
            let (pa, pb) = (plain.partial_v(u, v), dual.partial_v(u, v));
            let d = pa - pb;
            assert!(
                d.dot(&d).sqrt() < 1.0e-10,
                "partial_v must match the plain patch at ({u}, {v})"
            );
        }
    }
}

#[test]
fn boundary_cross_derivative_sees_only_its_own_candidate() {
    let plain = curved_plain();
    let mut dual = RectangularGregoryBiNonicDualTwist::from_plain(&plain);
    // Perturb only the P_v candidate of corner (0,0): the
    // u = 0 boundary's cross-derivative must be unchanged
    // (it reads P_u), while the v = 0 boundary's must move.
    let perturbed = Point3::new(
        dual.inner_cps[0].1.x,
        dual.inner_cps[0].1.y,
        dual.inner_cps[0].1.z + 0.05,
    );
    dual.inner_cps[0].1 = perturbed;
    let baseline = RectangularGregoryBiNonicDualTwist::from_plain(&plain);
    for sample in 1..8 {
        let t = sample as f64 / 8.0;
        let unchanged = dual.partial_u(0.0, t) - baseline.partial_u(0.0, t);
        assert!(
            unchanged.dot(&unchanged).sqrt() < 1.0e-12,
            "u=0 cross-derivative must ignore the P_v candidate at v={t}"
        );
    }
    let moved = dual.partial_v(0.5, 0.0) - baseline.partial_v(0.5, 0.0);
    assert!(
        moved.dot(&moved).sqrt() > 1.0e-4,
        "v=0 cross-derivative must respond to the P_v candidate"
    );
}

#[test]
fn analytic_partials_match_finite_differences_with_split_pairs() {
    let plain = curved_plain();
    let mut dual = RectangularGregoryBiNonicDualTwist::from_plain(&plain);
    // Split all four pairs to exercise the rational terms.
    for k in 0..4 {
        let (p_u, p_v) = dual.inner_cps[k];
        dual.inner_cps[k] = (
            Point3::new(p_u.x + 0.01, p_u.y - 0.02, p_u.z + 0.03),
            Point3::new(p_v.x - 0.015, p_v.y + 0.01, p_v.z - 0.02),
        );
    }
    let h = 1.0e-6;
    for &(u, v) in &[(0.2, 0.3), (0.5, 0.5), (0.85, 0.15), (0.1, 0.9)] {
        let fd_u = (dual.point_at(u + h, v) - dual.point_at(u - h, v)) * (0.5 / h);
        let an_u = dual.partial_u(u, v);
        let du = fd_u - an_u;
        assert!(
            du.dot(&du).sqrt() < 1.0e-6,
            "partial_u disagrees with finite differences at ({u}, {v})"
        );
        let fd_v = (dual.point_at(u, v + h) - dual.point_at(u, v - h)) * (0.5 / h);
        let an_v = dual.partial_v(u, v);
        let dv = fd_v - an_v;
        assert!(
            dv.dot(&dv).sqrt() < 1.0e-6,
            "partial_v disagrees with finite differences at ({u}, {v})"
        );
    }
}

#[test]
fn boundary_rows_are_plain_controls() {
    let plain = curved_plain();
    let mut dual = RectangularGregoryBiNonicDualTwist::from_plain(&plain);
    for k in 0..4 {
        let (p_u, p_v) = dual.inner_cps[k];
        dual.inner_cps[k] = (
            Point3::new(p_u.x, p_u.y, p_u.z + 0.1),
            Point3::new(p_v.x, p_v.y, p_v.z - 0.1),
        );
    }
    for sample in 0..=8 {
        let t = sample as f64 / 8.0;
        for (a, b) in [
            (dual.point_at(t, 0.0), plain.point_at(t, 0.0)),
            (dual.point_at(t, 1.0), plain.point_at(t, 1.0)),
            (dual.point_at(0.0, t), plain.point_at(0.0, t)),
            (dual.point_at(1.0, t), plain.point_at(1.0, t)),
        ] {
            assert!(
                a.distance_to(&b) < 1.0e-12,
                "boundary positions must be unaffected by pair splits (t={t})"
            );
        }
    }
}
