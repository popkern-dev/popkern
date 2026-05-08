#![cfg(test)]
use super::tests_support::*;
use super::*;

// ---- Validation ----

#[test]
fn zero_guide_validates() {
    for n in [3, 4, 5, 6, 10] {
        let g = Bi5Guide::zero(n);
        g.validate().unwrap();
    }
}

#[test]
#[should_panic]
fn zero_guide_below_three_panics() {
    let _ = Bi5Guide::zero(2);
}

#[test]
fn validate_detects_sector_count_mismatch() {
    let mut g = Bi5Guide::zero(5);
    g.sectors.pop();
    assert!(g.validate().is_err());
}

#[test]
fn validate_detects_non_finite_cp() {
    let mut g = Bi5Guide::zero(5);
    g.sectors[2].cps[3][4].x = f64::NAN;
    assert!(g.validate().is_err());
}

// ---- Evaluation ----

#[test]
fn sector_eval_at_corners_equals_corner_cps() {
    let mut sector = Bi5GuideSector::zero();
    sector.cps[0][0] = Point3::new(1.0, 2.0, 3.0);
    sector.cps[5][0] = Point3::new(4.0, 5.0, 6.0);
    sector.cps[0][5] = Point3::new(7.0, 8.0, 9.0);
    sector.cps[5][5] = Point3::new(10.0, 11.0, 12.0);
    assert!(approx_p(sector.eval(0.0, 0.0), sector.cps[0][0], 1e-12));
    assert!(approx_p(sector.eval(1.0, 0.0), sector.cps[5][0], 1e-12));
    assert!(approx_p(sector.eval(0.0, 1.0), sector.cps[0][5], 1e-12));
    assert!(approx_p(sector.eval(1.0, 1.0), sector.cps[5][5], 1e-12));
}

#[test]
fn eval_is_continuous_across_sector_boundary() {
    // Sample the shared edge at v = 0 for sector k+1 vs u = 1 for
    // sector k (both trace out the m_{k+1}-side of the boundary).
    let params = Bi5FreeParams::zero(6);
    let g_tri = apply_free_params(&params, Point3::new(0.0, 0.0, 0.0)).unwrap();
    let g_bi = Bi5Guide::from_gtri5(&g_tri);
    // Sector k's v = 1 edge vs sector (k + 1)'s u = 0 edge (both =
    // the shared m_{k+1} edge; note different parameter directions
    // along the edge in the two sectors, so we compare endpoints).
    for k in 0..6 {
        let k_next = (k + 1) % 6;
        let a = g_bi.sector(k).eval(0.0, 1.0); // m_{k+1} via sector k
        let b = g_bi.sector(k_next).eval(1.0, 0.0); // m_{k+1} via sector k+1
        assert!(
            approx_p(a, b, 1e-10),
            "C⁰ junction at k={k}, k_next={k_next}: {a:?} vs {b:?}"
        );
    }
}
