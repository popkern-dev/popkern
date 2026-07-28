//! PK-192 / PK-203 pipeline probes for the low-valence PnS caps,
//! split from `pns_low_valence_probes.rs` for the 500-line structural
//! gate. Same harness, same env-gated run commands in each doc
//! comment.

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
/// POPKERN_PNS_FIRST_N3=1 POPKERN_PNS_FORCE_ACCEPT=1 POPKERN_PNS_BI4_BORDER=1 \
///     cargo test -p popkern-fillet --test pns_low_valence_pipeline_probes \
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

    // PK-192's named fixture class: the same forced-acceptance
    // measurement on a unit-cube corner. Informational — the cube
    // corner's bridge edges are kef-consumed before the blend face
    // forms, so whether PnS even survives to bid is part of what this
    // records; the diagnostics say why when it does not.
    let c = pns_apex_outcome::run_cube_corner_outcome(&[0.04, 0.06, 0.05]);
    println!(
        "[pk203-force] cube corner served_by_pns={} valid={} cap_faces={} methods={:?} \
         seams={} g0={:.3e} g1={:.3e} cos_min={:.6}",
        c.served_by_pns(),
        c.valid,
        c.cap_faces,
        c.methods,
        c.metrics.seam_count,
        c.metrics.worst_g0_gap,
        c.metrics.worst_g1_normal,
        c.metrics.cos_min
    );
    for d in &c.pns_diagnostics {
        println!("[pk203-force]     cube diag: {d}");
    }
}

/// Rim versus spoke on the PnS cap ITSELF, mid-pipeline.
///
/// This is the measurement `probe_cap_rim_versus_spoke_deviation`
/// (in `pns_low_valence_probes.rs`) could not make: that probe runs
/// on the final body, where the PnS cap was rolled back, so it reads
/// the incumbent. The split here is computed INSIDE the dispatch
/// (`cap_rim_spoke_split` diagnostic, emitted at the same moment the
/// acceptance gate judges the cap) and surfaces through the
/// structured diagnostics channel.
///
/// Reads process-wide `POPKERN_PNS_FIRST_N3`, so it must run in its
/// own process and is ignored by default:
///
/// ```text
/// POPKERN_PNS_FIRST_N3=1 POPKERN_PNS_BI4_BORDER=1 cargo test -p popkern-fillet \
///     --test pns_low_valence_pipeline_probes probe_pns_mid_pipeline_rim_spoke_split \
///     -- --ignored --nocapture
/// ```
#[test]
#[ignore = "PK-192 rim/spoke probe: reads process-wide POPKERN_PNS_FIRST_N3; see doc comment for the run command"]
fn probe_pns_mid_pipeline_rim_spoke_split() {
    // Same read as `pns_promote::pns_first_at_low_valence` — the test
    // only labels its output with it; routing reads the real one.
    let released = std::env::var("POPKERN_PNS_FIRST_N3")
        .map(|v| v != "0")
        .unwrap_or(false);
    // The uniform-radius apex collapses to the mitred/bridge
    // construction and never reaches the dispatch (zero corner
    // patches, zero pns diagnostics) — printed for the record, not
    // asserted. The mixed-radius apex is the fixture every PK-192
    // number was measured on, and there the split must be recorded.
    for (n, radii, must_record) in [
        (3_usize, vec![0.05, 0.05, 0.05], false),
        (3, vec![0.04, 0.06, 0.05], true),
    ] {
        let o = run_pyramid_apex_outcome(n, &radii);
        println!(
            "[pk192-mid-split] n={n} radii={radii:?} released_for_pns={released} served_by_pns={}",
            o.served_by_pns()
        );
        for d in &o.pns_diagnostics {
            println!("[pk192-mid-split]     diag: {d}");
        }
        if released && must_record {
            assert!(
                o.pns_diagnostics
                    .iter()
                    .any(|d| d.contains("cap_rim_spoke_split")),
                "n={n}: with N=3 released for PnS the dispatch must record the \
                 rim/spoke split before its verdict; diagnostics were {:?}",
                o.pns_diagnostics
            );
        } else if !released {
            println!(
                "[pk192-mid-split] POPKERN_PNS_FIRST_N3 not set — N=3 stays on the \
                 direct-apex route and the PnS dispatch never runs. Set the env var \
                 per the doc comment; this run measured nothing."
            );
        }
    }

    // Valence-scaling control, attempted and retired (2026-07-28).
    // The PK-192 structural verdict asked for one n >= 5 measurement
    // of the same border-row split, on the theory that the kinked
    // arc+bridge fit error shrinks with sector angle. The control is
    // unreachable: with production options the setback/KP-polar route
    // claims the n = 5 mixed apex upstream (zero pns diagnostics);
    // with the polar ladder disabled the legacy Gregory ladder claims
    // it instead (methods=[Gregory], cap_faces=20, still zero pns
    // diagnostics). The mitred arc+bridge boundary class therefore
    // never reaches the PnS dispatch at n >= 5 AT ALL — the n >= 5
    // PnS parity gates (pk097_pns_*) ride setback boundaries, whose
    // per-sector rims are smooth. The ~6e-3 border-row floor is a
    // low-valence-only exposure by construction, and no scaling
    // number exists to measure.
    let o5 =
        pns_apex_outcome::run_pyramid_apex_outcome_no_polar(5, &[0.04, 0.06, 0.05, 0.055, 0.045]);
    println!(
        "[pk192-mid-split] n=5 mitred-class control: served_by_pns={} methods={:?} \
         cap_faces={} pns_diags={} — PnS never bids on mitred n>=5 corners (see comment)",
        o5.served_by_pns(),
        o5.methods,
        o5.cap_faces,
        o5.pns_diagnostics.len()
    );
}
