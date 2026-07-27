//! PK-203 pipeline probes for the low-valence PnS caps, split from
//! `pns_low_valence_probes.rs` for the 500-line structural gate.
//! Same harness, same env-gated run commands in each doc comment.

mod kp_gate_helpers;
#[path = "support/pns_apex_outcome.rs"]
mod pns_apex_outcome;

use kp_gate_helpers as helpers;
use pns_apex_outcome::run_pyramid_apex_outcome;

/// PK-203 decisive probe: what does a FINISHED PnS cap measure?
///
/// Every quantity the acceptance gate refuses on improves by orders
/// of magnitude between dispatch time and pipeline end (support
/// 110,000x, cap rim 3x from one completion pass) — but the refused
/// cap is rolled back before those passes run, so no finished PnS cap
/// has ever been measured. This forces the continuity refusal open
/// (`POPKERN_PNS_FORCE_ACCEPT`, topology refusals still decline),
/// lets the whole repair pipeline run, and measures the final body:
/// rim/spoke fidelity on the cap faces, live cap-cap seam continuity,
/// and cap/support seam continuity over the cap+neighbour face set.
///
/// Decision rule, wired before the numbers existed: if the finished
/// cap meets the 5.000e-3 bar, PK-203's late re-judged acceptance is
/// the entire remaining fix for PK-192 and no border/spoke geometry
/// mechanism is needed. If it stays over, the residual is genuinely
/// the cap's own and the spoke-curve authorship mechanism is next.
///
/// ```text
/// POPKERN_PNS_FIRST_N3=1 POPKERN_PNS_FORCE_ACCEPT=1 \
///     cargo test -p popkern-fillet --test pns_low_valence_probes \
///     probe_pns_force_accept_pipeline_end -- --ignored --nocapture
/// ```
#[test]
#[ignore = "PK-203 force-accept probe: reads process-wide POPKERN_PNS_FIRST_N3 + POPKERN_PNS_FORCE_ACCEPT; see doc comment"]
fn probe_pns_force_accept_pipeline_end() {
    let armed = std::env::var("POPKERN_PNS_FIRST_N3")
        .map(|v| v != "0")
        .unwrap_or(false)
        && std::env::var("POPKERN_PNS_FORCE_ACCEPT")
            .map(|v| v == "1")
            .unwrap_or(false);
    if !armed {
        println!("[pk203-force] flags not set — nothing measured. Run per the doc comment.");
        return;
    }
    let o = run_pyramid_apex_outcome(3, &[0.04, 0.06, 0.05]);
    println!(
        "[pk203-force] n=3 served_by_pns={} valid={} cap_faces={} methods={:?}",
        o.served_by_pns(),
        o.valid,
        o.cap_faces,
        o.methods
    );
    for d in &o.pns_diagnostics {
        println!("[pk203-force]     diag: {d}");
    }
    assert!(
        o.served_by_pns(),
        "force-accept armed but the PnS cap did not ship — the probe measured \
         the incumbent; diagnostics above say why"
    );
    let brep = &o.brep;
    let split = popkern_fillet::pcurve_fidelity::measure_cap_boundary_fidelity_split(
        brep,
        &o.cap_face_handles,
    );
    let fmt = |f: Option<popkern_fillet::pcurve_fidelity::PcurveFidelity>| {
        f.map(|f| {
            format!(
                "worst={:.4e} nearest={:.4e} edge={:?} face={:?} f={:.2} measured={}",
                f.worst,
                f.worst_nearest_point,
                f.worst_edge,
                f.worst_face,
                f.worst_fraction,
                f.measured
            )
        })
        .unwrap_or_else(|| "unmeasured".into())
    };
    println!("[pk203-force] FINAL rim  : {}", fmt(split.rim));
    println!("[pk203-force] FINAL spoke: {}", fmt(split.spoke));
    println!(
        "[pk203-force] FINAL cap-cap seams: n={} g0={:.3e} g1={:.3e} cos_min={:.6}",
        o.metrics.seam_count, o.metrics.worst_g0_gap, o.metrics.worst_g1_normal, o.metrics.cos_min
    );
    // Cap/support continuity: measure over the cap faces plus every
    // face across a cap edge, which makes the rim seams measurable by
    // the same live-seam instrument the gates use.
    let mut with_neighbours = o.cap_face_handles.clone();
    for &cap in &o.cap_face_handles {
        let Some(fd) = brep.faces.get(cap) else {
            continue;
        };
        let Some(ld) = brep.loops.get(fd.outer_loop) else {
            continue;
        };
        let start = ld.half_edge;
        let mut he = start;
        for _ in 0..32 {
            let Some(hd) = brep.half_edges.get(he) else {
                break;
            };
            if let Some(f) = brep
                .half_edges
                .get(hd.twin)
                .and_then(|t| brep.loops.get(t.loop_))
                .map(|l| l.face)
                && !with_neighbours.contains(&f)
            {
                with_neighbours.push(f);
            }
            he = hd.next;
            if he == start {
                break;
            }
        }
    }
    let rim_seams = helpers::measure_live_cap_seams(brep, &with_neighbours, helpers::N_SAMPLES);
    println!(
        "[pk203-force] FINAL cap+neighbour seams: n={} g0={:.3e} g1={:.3e} cos_min={:.6}",
        rim_seams.seam_count, rim_seams.worst_g0_gap, rim_seams.worst_g1_normal, rim_seams.cos_min
    );
    println!(
        "[pk203-force] verdict inputs: bar=5.000e-3; finished-cap rim/spoke above decide \
         whether late acceptance alone unblocks N=3"
    );
}
