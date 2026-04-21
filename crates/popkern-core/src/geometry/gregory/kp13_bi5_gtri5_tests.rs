//! Tests for the KP13 bi-5 auxiliary map `g△`.

use super::super::bernstein_tri5::flat_index;
use super::*;

const EPS: f64 = 1e-10;

fn approx_point(a: Point3, b: Point3) -> bool {
    (a.x - b.x).abs() < EPS && (a.y - b.y).abs() < EPS && (a.z - b.z).abs() < EPS
}

// ---- paper label → flat-index convention ----
//
// Paper triple `(a, b, c)` with `a` in the `o`-slot, `b` in the
// `m_k`-slot, `c` in the `m_{k+1}`-slot. Our basis uses
// `B^5_{i, j, k}(u, v, w = 1-u-v)`: at `(u, v) = (0, 0)` (w = 1),
// only `B^5_{0, 0, 5}` is non-zero → CP flat-index 0 is `o`.
// Hence paper `(a, b, c)` → our `(i, j, k) = (b, c, a)` →
// `flat_index(5, b, c)`.

#[test]
fn fig4_label_flat_indices_match() {
    // Paper 500 = star at o. → flat_index(5, 0, 0) = 0
    assert_eq!(flat_index(5, 0, 0), IDX_500);
    // Paper 050 = m_k. → flat_index(5, 5, 0) = 20
    assert_eq!(flat_index(5, 5, 0), IDX_050);
    // Paper 005 = m_{k+1}. → flat_index(5, 0, 5) = 5
    assert_eq!(flat_index(5, 0, 5), IDX_005);
    // Paper 230. → flat_index(5, 3, 0) = 15
    assert_eq!(flat_index(5, 3, 0), IDX_230);
    // Paper 122. → flat_index(5, 2, 2) = 13
    assert_eq!(flat_index(5, 2, 2), IDX_122);
    // Paper 032. → flat_index(5, 3, 2) = 17
    assert_eq!(flat_index(5, 3, 2), IDX_032);
    // Paper 023. → flat_index(5, 2, 3) = 14
    assert_eq!(flat_index(5, 2, 3), IDX_023);
}

#[test]
fn free_bb_indices_are_distinct() {
    let mut ixs = free_bb_indices();
    ixs.sort();
    for w in ixs.windows(2) {
        assert_ne!(w[0], w[1], "duplicate index in free_bb_indices");
    }
    assert_eq!(ixs.len(), 7);
}

// ---- sector_corner ----

#[test]
fn sector_corner_n6() {
    // c_0 = (1, 0), c_3 = (-1, 0) for n = 6.
    let (x0, y0) = sector_corner(6, 0);
    assert!((x0 - 1.0).abs() < EPS);
    assert!(y0.abs() < EPS);
    let (x3, y3) = sector_corner(6, 3);
    assert!((x3 + 1.0).abs() < EPS);
    assert!(y3.abs() < EPS);
}

// ---- Gtri5Sector eval ----

#[test]
fn zero_sector_evaluates_to_origin() {
    let s = Gtri5Sector::zero();
    assert!(approx_point(
        s.evaluate(0.3, 0.3),
        Point3::new(0.0, 0.0, 0.0)
    ));
    assert!(approx_point(
        s.partial_u(0.2, 0.4),
        Point3::new(0.0, 0.0, 0.0)
    ));
    assert!(approx_point(
        s.second_partial_uv(0.1, 0.5),
        Point3::new(0.0, 0.0, 0.0)
    ));
}

#[test]
fn sector_evaluates_corners_to_bb_coefs() {
    let mut s = Gtri5Sector::zero();
    // At (u, v) = (0, 0) → w = 1 → wedge vertex `o` → paper 500.
    let p_star = Point3::new(1.0, 2.0, 3.0);
    // At (u, v) = (1, 0) → wedge vertex `m_k` → paper 050.
    let p_mk = Point3::new(4.0, 5.0, 6.0);
    // At (u, v) = (0, 1) → wedge vertex `m_{k+1}` → paper 005.
    let p_mk1 = Point3::new(7.0, 8.0, 9.0);
    s.cps[IDX_500] = p_star;
    s.cps[IDX_050] = p_mk;
    s.cps[IDX_005] = p_mk1;
    assert!(approx_point(s.evaluate(0.0, 0.0), p_star));
    assert!(approx_point(s.evaluate(1.0, 0.0), p_mk));
    assert!(approx_point(s.evaluate(0.0, 1.0), p_mk1));
}

#[test]
fn sector_set_cp_roundtrip() {
    // Paper "230" = paper triple (2, 3, 0) = our (i, j, k) = (3, 0, 2)
    // → set_cp(3, 0, _).
    let mut s = Gtri5Sector::zero();
    let p = Point3::new(0.1, 0.2, 0.3);
    s.set_cp(3, 0, p);
    assert!(approx_point(s.cp(3, 0), p));
    assert_eq!(s.cps[IDX_230], p);
}

#[test]
#[should_panic]
fn cp_panics_out_of_simplex() {
    let s = Gtri5Sector::zero();
    let _ = s.cp(3, 3); // i + j = 6 > 5
}

// ---- Gtri5 constructors ----

#[test]
fn new_allocates_n_sectors() {
    let g = Gtri5::new(5);
    assert_eq!(g.n, 5);
    assert_eq!(g.sectors.len(), 5);
    for s in &g.sectors {
        for cp in &s.cps {
            assert!(approx_point(*cp, Point3::new(0.0, 0.0, 0.0)));
        }
    }
}

#[test]
#[should_panic]
fn new_below_three_panics() {
    let _ = Gtri5::new(2);
}

#[test]
fn star_point_is_shared_across_sectors() {
    let star = Point3::new(1.5, -0.25, 2.0);
    let g = Gtri5::new_with_star_point(6, star);
    assert!(approx_point(g.star_point(), star));
    for s in &g.sectors {
        assert!(approx_point(s.cps[IDX_500], star));
    }
}

#[test]
fn set_star_point_propagates() {
    let mut g = Gtri5::new(8);
    let star = Point3::new(0.0, 0.0, 1.0);
    g.set_star_point(star);
    for k in 0..8 {
        assert!(approx_point(g.sector(k).cps[IDX_500], star));
    }
}

// ---- Catmull-Clark limit convenience ----

#[test]
fn new_from_cap_net_pins_star_to_cc_limit() {
    // For a trivial cap-net where the center and all outer nodes
    // equal the same point, the CC limit is that point.
    let c = Point3::new(1.5, -0.5, 2.0);
    for n in [5, 6, 7, 8] {
        let c5s = vec![c; n];
        let c6s = vec![c; n];
        let g = Gtri5::new_from_cap_net(n, c, &c5s, &c6s);
        for k in 0..n {
            assert!(
                approx_point(g.sector(k).cps[IDX_500], c),
                "n = {n}, sector {k}: star != {c:?}"
            );
        }
    }
}

#[test]
fn new_from_cap_net_matches_standalone_formula() {
    // CC limit = (1 - n·ν₅ - n·ν₆)·c7 + Σ(ν₅·c5 + ν₆·c6). Verify
    // new_from_cap_net reproduces the standalone helper's output.
    let n = 6;
    let c7 = Point3::new(0.0, 0.0, 1.0);
    let c5s: Vec<Point3> = (0..n).map(|k| Point3::new(k as f64, 0.0, 0.0)).collect();
    let c6s: Vec<Point3> = (0..n).map(|k| Point3::new(0.0, k as f64, 0.0)).collect();
    let g = Gtri5::new_from_cap_net(n, c7, &c5s, &c6s);
    let expected = super::super::kp_polar_shared::catmull_clark_star_point(n, c7, &c5s, &c6s);
    assert!(approx_point(g.star_point(), expected));
}

// ---- Mid-edge invariant ----

#[test]
fn set_mid_edge_updates_both_neighbors() {
    let mut g = Gtri5::new(5);
    let p = Point3::new(3.0, 4.0, 0.0);
    g.set_mid_edge(2, p);
    // Sector 2's IDX_050 holds m_2.
    assert!(approx_point(g.sector(2).cps[IDX_050], p));
    // Sector 1's IDX_005 also holds m_2.
    assert!(approx_point(g.sector(1).cps[IDX_005], p));
    assert!(approx_point(g.mid_edge(2), p));
}

#[test]
fn set_mid_edge_wraps() {
    let mut g = Gtri5::new(5);
    let p = Point3::new(7.0, 0.0, 0.0);
    g.set_mid_edge(0, p);
    // Sector 0 IDX_050 and sector 4 (= (0 - 1) mod 5) IDX_005.
    assert!(approx_point(g.sector(0).cps[IDX_050], p));
    assert!(approx_point(g.sector(4).cps[IDX_005], p));
}

#[test]
fn resync_mid_edges_restores_invariant() {
    let mut g = Gtri5::new(4);
    // Break invariant via sector_mut.
    g.sector_mut(0).cps[IDX_050] = Point3::new(1.0, 0.0, 0.0);
    // (sector 3's IDX_005 was left at origin)
    g.resync_mid_edges();
    // Now all paired pairs agree.
    g.validate().unwrap();
}

// ---- validate() ----

#[test]
fn validate_zero_grid_passes() {
    let g = Gtri5::new(7);
    g.validate().unwrap();
}

#[test]
fn validate_with_star_pinning_passes() {
    let g = Gtri5::new_with_star_point(6, Point3::new(0.5, 0.5, 0.5));
    g.validate().unwrap();
}

#[test]
fn validate_detects_star_mismatch() {
    let mut g = Gtri5::new_with_star_point(5, Point3::new(1.0, 0.0, 0.0));
    g.sectors[2].cps[IDX_500] = Point3::new(1.0, 0.1, 0.0);
    let err = g.validate().unwrap_err();
    assert!(err.contains("star-point"));
}

#[test]
fn validate_detects_mid_edge_mismatch() {
    let mut g = Gtri5::new(4);
    g.sectors[1].cps[IDX_050] = Point3::new(1.0, 0.0, 0.0);
    // Sector 0's IDX_005 is still zero → should be equal.
    let err = g.validate().unwrap_err();
    assert!(err.contains("mid-edge"));
}

// ---- evaluate_extended (v_k domain extension) ----

#[test]
fn evaluate_extended_agrees_with_evaluate_inside_simplex() {
    // Inside the simplex, extended should match the checked variant.
    let mut s = Gtri5Sector::zero();
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let k = 5 - i - j;
            s.set_cp(i, j, Point3::new(i as f64, j as f64, k as f64));
        }
    }
    for (u, v) in &[(0.1, 0.1), (0.3, 0.5), (0.5, 0.3), (0.0, 0.0), (1.0, 0.0)] {
        let a = s.evaluate(*u, *v);
        let b = s.evaluate_extended(*u, *v);
        assert!(
            approx_point(a, b),
            "disagreement at ({u}, {v}): {a:?} vs {b:?}"
        );
    }
}

#[test]
fn evaluate_extended_linear_matches_analytic_outside_simplex() {
    // Linear polynomial f(u, v) = 2u + 3v + 1 on the simplex — all
    // BB-coefs are f(i/5, j/5). Then evaluate_extended at (u, v)
    // outside the simplex (e.g. u + v > 1 or u < 0) should yield
    // the analytic value.
    let mut s = Gtri5Sector::zero();
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let u = i as f64 / 5.0;
            let v = j as f64 / 5.0;
            s.set_cp(i, j, Point3::new(2.0 * u + 3.0 * v + 1.0, 0.0, 0.0));
        }
    }
    // Outside simplex: (1.5, 0.3) → u + v = 1.8 > 1.
    let out = s.evaluate_extended(1.5, 0.3);
    let expected = 2.0 * 1.5 + 3.0 * 0.3 + 1.0;
    assert!(
        (out.x - expected).abs() < 1e-10,
        "extended eval at (1.5, 0.3): got {}, expected {}",
        out.x,
        expected
    );
    // Negative barycentric: (-0.2, 0.5).
    let out = s.evaluate_extended(-0.2, 0.5);
    let expected = 2.0 * (-0.2) + 3.0 * 0.5 + 1.0;
    assert!(
        (out.x - expected).abs() < 1e-10,
        "extended eval at (-0.2, 0.5): got {}, expected {}",
        out.x,
        expected
    );
}

// ---- Partial-derivative FD check ----

#[test]
fn partial_u_matches_fd() {
    // Set a non-trivial test sector with identifiable coefficients.
    let mut s = Gtri5Sector::zero();
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let k = 5 - i - j;
            s.set_cp(i, j, Point3::new(i as f64, j as f64, k as f64));
        }
    }
    let h = 1e-6;
    let (u, v) = (0.3, 0.4);
    let fd_x = (s.evaluate(u + h, v).x - s.evaluate(u - h, v).x) / (2.0 * h);
    let analytic = s.partial_u(u, v);
    assert!(
        (fd_x - analytic.x).abs() < 1e-5,
        "∂u FD mismatch: fd = {}, analytic = {}",
        fd_x,
        analytic.x
    );
}

#[test]
fn partial_v_matches_fd() {
    let mut s = Gtri5Sector::zero();
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let k = 5 - i - j;
            s.set_cp(
                i,
                j,
                Point3::new((i + 1) as f64, (j + 2) as f64, (k + 3) as f64),
            );
        }
    }
    let h = 1e-6;
    let (u, v) = (0.2, 0.3);
    let fd_y = (s.evaluate(u, v + h).y - s.evaluate(u, v - h).y) / (2.0 * h);
    let analytic = s.partial_v(u, v);
    assert!(
        (fd_y - analytic.y).abs() < 1e-5,
        "∂v FD mismatch: fd = {}, analytic = {}",
        fd_y,
        analytic.y
    );
}
