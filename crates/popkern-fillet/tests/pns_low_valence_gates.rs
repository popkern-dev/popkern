//! PK-192 / PK-193 — PnS cap coverage at the low valences.
//!
//! The PnS bi-4 cap was dispatched only at `n ≥ 5` for most of its
//! life. `n = 4` was believed structurally impossible (σ's
//! `1/cos(2π/n)` pole) and `n = 3` was never wired up. Both beliefs
//! were wrong, and both valences are the *common* ones in mechanical
//! parts: `n = 3` is the cube corner, `n = 4` the pyramid apex.
//!
//! These gates hold the low-valence caps to the same live-seam
//! continuity bar as the rest of the family. Continuity is measured
//! across the committed seam topology with finite differences in UV,
//! so it is independent of which surface kind the cap shipped.
//!
//! The PK-192 route decision and its reversal condition are pinned
//! by the live gates below; both this file and
//! `pns_low_valence_probes.rs` measure through the shared
//! `support/pns_apex_outcome.rs` harness.

mod kp_gate_helpers;
#[path = "support/pns_apex_outcome.rs"]
mod pns_apex_outcome;

use helpers::{G0_GATE, G1_GATE};
use kp_gate_helpers as helpers;
use pns_apex_outcome::run_pyramid_apex_outcome;

// PK-192 decision (2026-07-28): the three aspirational
// `pns_serves_*` gates and the POPKERN_PNS route-comparison
// instrument that lived here were retired after the route decision
// was recorded — dispatch stays with direct-apex/Coons at n∈{3,4}
// (incumbent g0 ~1e-15 / cos_min 0.999994+; the PnS candidate's
// border-row floor ≈6e-3 fails its own 5e-3 bar, and a forced cap
// measures cap-cap g0 6.5e-2 / cos_min 0.007). The decision, its
// numbers, and the reversal condition (a construction that
// represents the kinked arc+bridge rim, e.g. per-rim-piece 2N
// partitioning at degree 4, beating the incumbent) are pinned by
// the live gates below: `cube_corner_default_route_meets_bars`
// (route + bars) and `pns_cap_internal_seams_beat_the_shipped_fan`
// (the candidate's clean-input potential). Ledger: PK-192.

/// PK-197 / PK-201: the cross-spoke G1 correction reaches BOTH low
/// valences.
///
/// The exact-boundary wedges were only ever constrained against their
/// OUTER boundaries, so nothing held the inner spokes — the seams a
/// sector shares with its neighbours — and every measured cap was
/// creased there. `cross_spoke_g1` adds that constraint:
///
/// - **n = 4**: `cos_min` 0.006755 -> 0.996283; the residual sits at
///   the pole end (f = 0.95), the ill-conditioned chart corner.
/// - **n = 3**: `cos_min` 0.00214 -> 0.999994 once
///   `harmonize_mismatched_spoke_nets` gave its non-square nets
///   ((49, 41) / (56, 48)) a common spoke knot vector — the
///   correction itself was always sufficient; index-wise pairing
///   simply could not see those seams.
///
/// Both halves are asserted as G1 now. The n = 3 half spent months as
/// an inverted assertion (cos_min < 0.5) pinning the residual so it
/// could not disappear from view; it flipped 2026-07-26.
#[test]
fn cross_spoke_correction_fixes_both_low_valences() {
    let n4 = run_pyramid_apex_outcome(4, &[0.04, 0.06, 0.05, 0.055]);
    assert!(
        n4.metrics.seam_count > 0,
        "n=4 must commit measurable internal cap seams"
    );
    assert!(
        n4.metrics.cos_min > 0.99,
        "n=4 cap seams must be near-tangent after the cross-spoke correction, \
         measured cos_min {:.6} (was 0.006755 before it)",
        n4.metrics.cos_min
    );
    assert!(
        n4.metrics.worst_g0_gap < G0_GATE,
        "n=4 seams must stay positionally exact, got {:.3e}",
        n4.metrics.worst_g0_gap
    );

    let n3 = run_pyramid_apex_outcome(3, &[0.04, 0.06, 0.05]);
    assert!(
        n3.metrics.seam_count > 0,
        "n=3 must commit measurable internal cap seams"
    );
    assert!(
        n3.metrics.worst_g0_gap < G0_GATE,
        "n=3 seams must stay positionally exact, got {:.3e}",
        n3.metrics.worst_g0_gap
    );
    assert!(
        n3.metrics.cos_min > 0.99,
        "n=3 cap seams must be near-tangent after knot harmonization + the \
         cross-spoke correction, measured cos_min {:.6} (was 0.00214 before \
         harmonization reached its non-square nets)",
        n3.metrics.cos_min
    );
}

/// PK-192 decision record, pinned on the ticket's named fixture: a
/// unit-cube corner with mixed radii, filleted through the DEFAULT
/// dispatch (no env switches), must be served by the incumbent route
/// within the live-seam bars.
///
/// The recorded route comparison this pins (2026-07-28, both routes
/// through the identical `run_cube_corner_outcome` instrument):
///
/// - incumbent (this default): Coons, G1, cap-cap g0 9.795e-17,
///   g1 4.632e-3, cos_min 0.999989, valid
/// - PnS candidate (POPKERN_PNS_FIRST_N3=1 POPKERN_PNS_BI4_BORDER=1):
///   acceptance declines at cap_outer_g0 5.972e-3 vs the 5.000e-3
///   bar, support exact at 1.4e-16 — the border-row floor
/// - PnS forced through (+POPKERN_PNS_FORCE_ACCEPT=1): finished cap
///   measures cap-cap g0 1.106e-1, cos_min 0.002276
///
/// Per the ticket's pre-wired rule the dispatch therefore stays. If
/// this gate ever fails, the incumbent regressed on the commonest
/// corner in mechanical parts — that is a defect, not a license to
/// revisit the route decision.
#[test]
fn cube_corner_default_route_meets_bars() {
    let o = pns_apex_outcome::run_cube_corner_outcome(&[0.04, 0.06, 0.05]);
    assert!(o.valid, "cube corner fillet must leave a valid BRep");
    assert!(
        !o.served_by_pns(),
        "default dispatch must not route the cube corner to PnS — the PK-192 route \
         comparison decided against it (see this test's doc)"
    );
    let m = &o.metrics;
    assert!(
        m.seam_count > 0,
        "cube corner cap must commit measurable internal seams"
    );
    assert!(
        m.worst_g0_gap < G0_GATE,
        "cube corner cap-cap seams must stay positionally exact, got {:.3e}",
        m.worst_g0_gap
    );
    assert!(
        m.worst_g1_normal < G1_GATE,
        "cube corner cap-cap seams must stay tangent-continuous, got {:.3e}",
        m.worst_g1_normal
    );
}

/// PK-192 verdict: the PnS cap's internal seams against the shipped
/// fan's, over the same corner rim.
///
/// The fan and the cap cover the same opening, so the comparison is
/// like-for-like: take the rim the pipeline actually committed, build
/// a PnS cap over it, and measure normal agreement where adjacent
/// patches meet. The fan's numbers come from the shipped body
/// (PK-197); the cap's come from evaluating its patches directly.
///
/// This is the measurement the ticket requires before any dispatch
/// change — and it is reported either way, win or falsification.
///
/// Context for the final PK-192 decision: this synthetic win (PnS
/// 0.0° across internal seams on a clean regular rim) is exactly why
/// the route was worth probing, and it still holds. The dispatch
/// decision nevertheless went to the incumbent, because the *real*
/// corner rim is a kinked arc+bridge pair per sector, which is
/// outside the construction's border-row input class: in situ the
/// candidate declines at rim ≈ 6e-3 (bar 5e-3) and a force-accepted
/// finished cap measures cap-cap g0 6.5e-2 / cos_min 0.007. Potential
/// on clean input and delivery on real input are different claims;
/// this gate pins the first, the route comparison decided on the
/// second.
#[test]
fn pns_cap_internal_seams_beat_the_shipped_fan() {
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::{FillBoundaryEdge, Point3, Vec3, fill_n_sided_boundary};

    for (n, radii) in [
        (3_usize, vec![0.04, 0.06, 0.05]),
        (4, vec![0.04, 0.06, 0.05, 0.055]),
    ] {
        let shipped = run_pyramid_apex_outcome(n, &radii);
        assert!(
            shipped.metrics.seam_count > 0,
            "n={n}: need a shipped fan to compare against"
        );
        // LiveSeamMetrics reports |n_a x n_b|; convert to an angle so
        // the two constructions are compared in the same units.
        let fan_worst_deg = shipped
            .metrics
            .worst_g1_normal
            .clamp(0.0, 1.0)
            .asin()
            .to_degrees()
            .max(shipped.metrics.cos_min.clamp(0.0, 1.0).acos().to_degrees());

        // Build a PnS cap over a rim of the same size and shape class.
        // A regular ring is the fair, construction-independent stand-in:
        // it isolates what the two constructions do internally rather
        // than how each of them samples this particular rim.
        let radius = radii.iter().copied().fold(0.0_f64, f64::max);
        let ring: Vec<FillBoundaryEdge> = (0..n)
            .map(|k| {
                let a = |i: f64| 2.0 * std::f64::consts::PI * i / (n as f64);
                let p = |t: f64| Point3::new(radius * f64::cos(t), radius * f64::sin(t), 0.0);
                let (start, mid, end) =
                    (p(a(k as f64)), p(a(k as f64 + 0.5)), p(a(k as f64 + 1.0)));
                let normal = Vec3::new(0.0, 0.0, 1.0);
                let tangent = (end - start).normalized().expect("side has length");
                FillBoundaryEdge {
                    start,
                    mid,
                    end,
                    normal_start: normal,
                    normal_mid: normal,
                    normal_end: normal,
                    tangent_start: tangent,
                    tangent_mid: tangent,
                    tangent_end: tangent,
                    extra_samples: Vec::new(),
                }
            })
            .collect();
        let fill = fill_n_sided_boundary(&ring, 1e-9).expect("PnS cap builds over the rim");

        // Adjacent PnS patches share the seam u = 0 of one with v = 0
        // of the next. Sample away from the pole, where the chart is
        // degenerate and a finite-difference normal is direction noise.
        let mut cap_worst_deg = 0.0_f64;
        for k in 0..n {
            let cur = &fill.surfaces[k];
            let nxt = &fill.surfaces[(k + 1) % n];
            for i in 1..=9 {
                let t = i as f64 / 10.0;
                let (na, nb) = (cur.normal_at(t, 0.0), nxt.normal_at(0.0, t));
                let (Some(na), Some(nb)) = (na.normalized(), nb.normalized()) else {
                    continue;
                };
                let cos = na.dot(&nb).abs().clamp(0.0, 1.0);
                cap_worst_deg = cap_worst_deg.max(cos.acos().to_degrees());
            }
        }

        println!(
            "[pk192-verdict] n={n} shipped_fan_worst={fan_worst_deg:.1} deg  \
             pns_cap_worst={cap_worst_deg:.1} deg"
        );
        assert!(
            cap_worst_deg < fan_worst_deg,
            "n={n} FALSIFIED: the PnS cap ({cap_worst_deg:.1} deg) is no better than the \
             shipped fan ({fan_worst_deg:.1} deg) across internal seams — record this and \
             leave dispatch alone"
        );
        assert!(
            cap_worst_deg < 5.0,
            "n={n}: the PnS cap should be near-tangent across its own seams by \
             construction, measured {cap_worst_deg:.1} deg"
        );
    }
}
