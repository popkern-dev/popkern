//! Tests for KP13 eq (6) post-LSQ m_k exact match.

use super::super::kp13_bi5_boundary::TensorBorder;
use super::super::kp13_bi5_free_params::{Bi5FreeParams, apply_free_params};
use super::super::kp13_bi5_gtri5::Gtri5;
use super::super::kp13_bi5_lsq::{composite_jet_bilinear, solve_g_tri_lsq};
use super::*;

const EPS: f64 = 1e-10;

fn approx_p(a: Point3, b: Point3) -> bool {
    (a.x - b.x).abs() < EPS && (a.y - b.y).abs() < EPS && (a.z - b.z).abs() < EPS
}

// ---- Forward-substitution solve ----

#[test]
fn zero_target_gives_zero_cps() {
    let cps = solve_mk_corner_cps_per_axis(&TaylorJet2::zero());
    for cp in &cps {
        assert!(approx_p(*cp, Point3::new(0.0, 0.0, 0.0)));
    }
}

#[test]
fn solve_then_read_round_trips_value_and_du() {
    // Value-only jet: ∂u, ∂v, ... = 0, value = (1, 0, 0).
    // Expected CPs: flat 20 = (1, 0, 0), flat 18 = (1, 0, 0) (forced
    // by ∂u = 5·(20 − 18) = 0), flat 19 = flat 18 + 0 = (1, 0, 0),
    // etc. All 6 CPs = (1, 0, 0).
    let mut target = TaylorJet2::zero();
    target.value = Point3::new(1.0, 0.0, 0.0);
    let cps = solve_mk_corner_cps_per_axis(&target);
    for cp in &cps {
        assert!(approx_p(*cp, Point3::new(1.0, 0.0, 0.0)));
    }
}

#[test]
fn solve_reproduces_jet_via_tri5_forward_eval() {
    // Write the solved CPs into a Gtri5 sector, then evaluate the
    // composite jet — it should match the original target.
    let target = TaylorJet2 {
        value: Point3::new(1.5, -0.5, 0.25),
        d_u: Point3::new(0.3, 0.2, 0.0),
        d_v: Point3::new(-0.1, 0.4, 0.1),
        d_uu: Point3::new(0.5, 0.0, 0.0),
        d_uv: Point3::new(0.0, 0.1, -0.2),
        d_vv: Point3::new(0.2, -0.3, 0.0),
    };
    let cps = solve_mk_corner_cps_per_axis(&target);
    let mut g = Gtri5::new_with_star_point(5, Point3::new(0.0, 0.0, 0.0));
    for (slot, &flat) in MK_CORNER_FLATS.iter().enumerate() {
        g.sector_mut(0).cps[flat] = cps[slot];
    }
    // Re-derive the jet via composite_jet_bilinear at Mk corner.
    let recovered = composite_jet_bilinear(&g, 0, ChiCorner::Mk);
    assert!(approx_p(recovered.value, target.value));
    assert!(approx_p(recovered.d_u, target.d_u));
    assert!(approx_p(recovered.d_v, target.d_v));
    assert!(approx_p(recovered.d_uu, target.d_uu));
    assert!(approx_p(recovered.d_uv, target.d_uv));
    assert!(approx_p(recovered.d_vv, target.d_vv));
}

// ---- apply_eq6_raw ----

#[test]
fn apply_eq6_raw_makes_jets_match_target_for_consistent_border() {
    // For a consistent tensor-border (same CPs across all sectors),
    // adjacent sectors' jets at the shared m_k vertex agree, so
    // after apply_eq6_raw every sector's composite jet at Mk matches
    // its target exactly.
    let n = 6;
    let mut target = TensorBorder::zero(n);
    let off = 0.25;
    for patch in target.patches.iter_mut() {
        for row in patch.cps.iter_mut() {
            for cp in row.iter_mut() {
                cp.x = off;
            }
        }
        patch.cps[2][0].x += 0.05;
    }
    let cc = Point3::new(0.0, 0.0, 0.0);
    let mut g = solve_g_tri_lsq(&target, cc).unwrap();
    apply_eq6_raw(&mut g, &target);
    for k in 0..n {
        let tgt = target.jet(k, ChiCorner::Mk);
        let actual = composite_jet_bilinear(&g, k, ChiCorner::Mk);
        assert!(approx_p(actual.value, tgt.value), "k={k}: value");
        assert!(approx_p(actual.d_u, tgt.d_u), "k={k}: d_u");
        assert!(approx_p(actual.d_v, tgt.d_v), "k={k}: d_v");
        assert!(approx_p(actual.d_uu, tgt.d_uu), "k={k}: d_uu");
        assert!(approx_p(actual.d_uv, tgt.d_uv), "k={k}: d_uv");
        assert!(approx_p(actual.d_vv, tgt.d_vv), "k={k}: d_vv");
    }
    // Structural invariants still hold.
    g.validate().unwrap();
}

#[test]
fn apply_eq6_raw_preserves_validate_with_inconsistent_border() {
    // With inconsistent tensor-border data (different values per
    // sector), the shared-vertex sync enforces last-write-wins.
    // Structural validate() still passes; the per-sector jet may
    // not match exactly for every sector but that's user-visible
    // behaviour of inconsistent inputs.
    let n = 6;
    let mut target = TensorBorder::zero(n);
    for (k, patch) in target.patches.iter_mut().enumerate() {
        let off = 0.1 * k as f64;
        for row in patch.cps.iter_mut() {
            for cp in row.iter_mut() {
                cp.x = off;
            }
        }
    }
    let cc = Point3::new(0.0, 0.0, 0.0);
    let mut g = solve_g_tri_lsq(&target, cc).unwrap();
    apply_eq6_raw(&mut g, &target);
    g.validate().unwrap();
}

// ---- apply_eq6_with_c2_reconcile ----

#[test]
fn apply_eq6_with_reconcile_keeps_mk_jet_exact_for_consistent_border() {
    let n = 7;
    let mut target = TensorBorder::zero(n);
    // Consistent target: same CPs across all sectors so the shared-
    // vertex values agree.
    let off = 0.2;
    for patch in target.patches.iter_mut() {
        for row in patch.cps.iter_mut() {
            for cp in row.iter_mut() {
                cp.x = off;
            }
        }
    }
    let cc = Point3::new(0.3, 0.0, 0.0);
    let mut g = solve_g_tri_lsq(&target, cc).unwrap();
    apply_eq6_with_c2_reconcile(&mut g, &target).unwrap();

    // After reconcile, Mk jet matches exactly for every sector.
    for k in 0..n {
        let tgt = target.jet(k, ChiCorner::Mk);
        let actual = composite_jet_bilinear(&g, k, ChiCorner::Mk);
        assert!(approx_p(actual.value, tgt.value), "k={k}: value");
        assert!(approx_p(actual.d_u, tgt.d_u), "k={k}: d_u");
        assert!(approx_p(actual.d_uu, tgt.d_uu), "k={k}: d_uu");
    }
    // Structural invariants still hold.
    g.validate().unwrap();
}

#[test]
fn apply_eq6_with_reconcile_passes_validate_for_inconsistent_border() {
    // With inconsistent boundary data, exact per-sector match is
    // impossible (shared-vertex conflict). We only require the
    // reconcile to run without panics and leave a valid Gtri5.
    let n = 6;
    let mut target = TensorBorder::zero(n);
    for (k, patch) in target.patches.iter_mut().enumerate() {
        for row in patch.cps.iter_mut() {
            for cp in row.iter_mut() {
                cp.x = 0.1 * k as f64;
            }
        }
    }
    let cc = Point3::new(0.3, 0.0, 0.0);
    let mut g = solve_g_tri_lsq(&target, cc).unwrap();
    apply_eq6_with_c2_reconcile(&mut g, &target).unwrap();
    g.validate().unwrap();
}

#[test]
fn apply_eq6_zero_target_zero_star_keeps_all_zero() {
    let n = 6;
    let target = TensorBorder::zero(n);
    let cc = Point3::new(0.0, 0.0, 0.0);
    let mut g = apply_free_params(&Bi5FreeParams::zero(n), cc).unwrap();
    apply_eq6_raw(&mut g, &target);
    for k in 0..n {
        for cp in g.sector(k).cps.iter() {
            assert!(approx_p(*cp, Point3::new(0.0, 0.0, 0.0)));
        }
    }
}
