use super::tests_support::*;
use super::*;

// ---- Bi5Guide::from_gtri5: end-to-end ----

#[test]
fn from_gtri5_preserves_star_per_sector() {
    // All sectors share the same star point (I.11 invariant). After
    // degree-raise, each sector's bi-5 CP (0, 0) should equal the
    // shared star.
    let params = Bi5FreeParams::zero(6);
    let star = Point3::new(0.5, -0.2, 1.3);
    let g_tri = apply_free_params(&params, star).unwrap();
    let g_bi = Bi5Guide::from_gtri5(&g_tri);
    for k in 0..6 {
        let bi_star = g_bi.sector(k).cps[0][0];
        assert!(
            approx_p(bi_star, star, 1e-10),
            "sector {k} bi-5 star {bi_star:?} ≠ tri-5 star {star:?}"
        );
    }
}

#[test]
fn from_gtri5_c0_junction_between_sectors() {
    // Adjacent sectors share a boundary edge. Sector k's (0, 5)
    // corner = m_{k+1} = sector (k+1)'s (5, 0) corner = m_{k+1}.
    // The degree-raised guide inherits this shared edge.
    let params = Bi5FreeParams::zero(7);
    let g_tri = apply_free_params(&params, Point3::new(0.0, 0.0, 0.0)).unwrap();
    let g_bi = Bi5Guide::from_gtri5(&g_tri);
    for k in 0..7 {
        let k_next = (k + 1) % 7;
        let mk1_via_k = g_bi.sector(k).cps[0][5];
        let mk1_via_next = g_bi.sector(k_next).cps[5][0];
        assert!(
            approx_p(mk1_via_k, mk1_via_next, 1e-10),
            "sector {k}→{k_next} mid-edge mismatch: {mk1_via_k:?} vs {mk1_via_next:?}"
        );
    }
}

#[test]
fn from_gtri5_works_for_all_valences() {
    // Including low-valence n ∈ {3, 4} (via KP09a § 6 pseudo-inverse
    // in the circulant upstream).
    for n in [3, 4, 5, 6, 7, 8] {
        let params = Bi5FreeParams::zero(n);
        let g_tri = apply_free_params(&params, Point3::new(0.1, 0.2, 0.3)).unwrap();
        let g_bi = Bi5Guide::from_gtri5(&g_tri);
        g_bi.validate().unwrap_or_else(|e| panic!("n={n}: {e}"));
        assert_eq!(g_bi.n, n);
        assert_eq!(g_bi.sectors.len(), n);
    }
}
