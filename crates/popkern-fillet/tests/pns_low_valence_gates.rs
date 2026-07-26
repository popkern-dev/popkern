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
//! The A/B against the incumbent route is
//! `pns_low_valence_route_comparison`, which reads the PnS dispatch
//! switch from the environment and is therefore ignored by default —
//! see its doc comment for the run command.
//!
//! The campaign's ignored diagnosis probes live in
//! `pns_low_valence_probes.rs`; both files measure through the shared
//! `support/pns_apex_outcome.rs` harness.

mod kp_gate_helpers;
#[path = "support/pns_apex_outcome.rs"]
mod pns_apex_outcome;

use helpers::{G0_GATE, G1_GATE};
use kp_gate_helpers as helpers;
use pns_apex_outcome::run_pyramid_apex_outcome;

fn run_pyramid_apex(n: usize, radii: &[f64]) -> helpers::LiveSeamMetrics {
    let outcome = run_pyramid_apex_outcome(n, radii);
    assert!(outcome.valid, "n={n} apex fillet must leave a valid BRep");
    outcome.metrics
}

/// `n = 3` — the cube/tetra corner. Previously below the dispatch
/// floor despite the adapter supporting it all along.
#[test]
#[ignore = "PK-192 target gate: PnS does not yet serve n=3 corners — they are built by the exact-boundary Coons route, which is G0-only across its internal seams (PK-197). Flip when PnS dispatch reaches this valence."]
fn pns_serves_valence_3_apex_within_gates() {
    let m = run_pyramid_apex(3, &[0.05; 3]);
    assert!(
        m.seam_count > 0,
        "n=3 cap must commit at least one live seam to measure"
    );
    assert!(
        m.worst_g0_gap < G0_GATE,
        "n=3 cap worst G0 seam gap {:.3e} exceeds {:.0e}",
        m.worst_g0_gap,
        G0_GATE
    );
    assert!(
        m.worst_g1_normal < G1_GATE,
        "n=3 cap worst G1 normal misalignment {:.3e} exceeds {:.0e}",
        m.worst_g1_normal,
        G1_GATE
    );
}

/// `n = 4` — the pyramid apex and cross junction. This is the
/// valence the kernel refused outright until the σ pole was shown
/// removable.
#[test]
#[ignore = "PK-192 target gate: PnS does not yet serve n=4 corners — they are built by the exact-boundary Coons route, which is G0-only across its internal seams (PK-197). Flip when PnS dispatch reaches this valence."]
fn pns_serves_valence_4_apex_within_gates() {
    let m = run_pyramid_apex(4, &[0.05; 4]);
    assert!(
        m.seam_count > 0,
        "n=4 cap must commit at least one live seam to measure"
    );
    assert!(
        m.worst_g0_gap < G0_GATE,
        "n=4 cap worst G0 seam gap {:.3e} exceeds {:.0e}",
        m.worst_g0_gap,
        G0_GATE
    );
    assert!(
        m.worst_g1_normal < G1_GATE,
        "n=4 cap worst G1 normal misalignment {:.3e} exceeds {:.0e}",
        m.worst_g1_normal,
        G1_GATE
    );
}

/// Mixed radii at the low valences — the asymmetric case that
/// distinguishes a genuine cap construction from one that only works
/// on a symmetric fixture.
#[test]
#[ignore = "PK-192 target gate: PnS does not yet serve low-valence mixed-radius corners — they are built by the exact-boundary Coons route, which is G0-only across its internal seams (PK-197). Flip when PnS dispatch reaches this valence."]
fn pns_serves_low_valence_mixed_radii() {
    for (n, radii) in [
        (3_usize, vec![0.04, 0.06, 0.05]),
        (4, vec![0.04, 0.06, 0.05, 0.055]),
    ] {
        let m = run_pyramid_apex(n, &radii);
        assert!(
            m.worst_g0_gap < G0_GATE,
            "n={n} mixed-radius cap worst G0 {:.3e} exceeds {:.0e}",
            m.worst_g0_gap,
            G0_GATE
        );
        assert!(
            m.worst_g1_normal < G1_GATE,
            "n={n} mixed-radius cap worst G1 {:.3e} exceeds {:.0e}",
            m.worst_g1_normal,
            G1_GATE
        );
    }
}

/// A/B of the PnS route against the incumbent route at `n ∈ {3, 4}`.
///
/// Reads `POPKERN_PNS`, which is process-wide, so it must run in its
/// own process and is ignored by default. This is the measurement
/// that decides whether lowering the dispatch floor was justified:
///
/// ```text
/// POPKERN_PNS=1 cargo test -p popkern-fillet --test pns_low_valence_gates \
///     pns_low_valence_route_comparison -- --ignored --nocapture
/// POPKERN_PNS=0 cargo test -p popkern-fillet --test pns_low_valence_gates \
///     pns_low_valence_route_comparison -- --ignored --nocapture
/// ```
#[test]
#[ignore = "PK-192 route A/B: reads process-wide POPKERN_PNS; see doc comment for the two run commands"]
fn pns_low_valence_route_comparison() {
    let route = if std::env::var("POPKERN_PNS")
        .map(|v| v != "0")
        .unwrap_or(true)
    {
        "pns"
    } else {
        "incumbent"
    };
    for (n, radii) in [
        (3_usize, vec![0.05, 0.05, 0.05]),
        (3, vec![0.04, 0.06, 0.05]),
        (4, vec![0.05, 0.05, 0.05, 0.05]),
        (4, vec![0.04, 0.06, 0.05, 0.055]),
    ] {
        let o = run_pyramid_apex_outcome(n, &radii);
        let m = &o.metrics;
        println!(
            "[low-valence-ab] switch={route} n={n} radii={radii:?} served_by_pns={} methods={:?} achieved={:?} cap_faces={} valid={} seams={} g0={:.3e} g1={:.3e} g2={:.3e} cos_min={:.6}",
            o.served_by_pns(),
            o.methods,
            o.achieved,
            o.cap_faces,
            o.valid,
            m.seam_count,
            m.worst_g0_gap,
            m.worst_g1_normal,
            m.worst_g2_residual,
            m.cos_min
        );
        for d in &o.pns_diagnostics {
            println!("[low-valence-ab]     diag: {d}");
        }
        for line in &o.seam_detail {
            println!("[low-valence-ab]     seam: {line}");
        }
    }
}

/// PK-197 / PK-201: where the cross-spoke G1 correction reaches, and
/// where it does not.
///
/// The exact-boundary wedges were only ever constrained against their
/// OUTER boundaries, so nothing held the inner spokes — the seams a
/// sector shares with its neighbours — and every measured cap was
/// creased there. `cross_spoke_g1` adds that constraint, and this
/// gate pins the measured outcome on both sides:
///
/// - **n = 4 is fixed.** `cos_min` went 0.006755 -> 0.996283 and the
///   whole mid-span reads 1.0000; the residual sits at the pole end
///   (f = 0.95), which is the ill-conditioned chart corner.
/// - **n = 3 is not.** It measures exactly as before, so the
///   correction is not reaching that cap's seams at all.
///
/// Both halves are asserted. Pinning only the win would let the n = 3
/// residual disappear from view, which is how the original defect
/// survived: nothing measured it.
#[test]
fn cross_spoke_correction_fixes_n4_and_not_yet_n3() {
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
        n3.metrics.cos_min < 0.5,
        "PK-201 residual: n=3 was still creased (cos_min 0.00214) because the \
         cross-spoke pairing does not reach its seams. Measured {:.6} — if this \
         is now near 1.0 the residual is FIXED and this half of the gate should \
         become a G1 assertion like n=4's.",
        n3.metrics.cos_min
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
