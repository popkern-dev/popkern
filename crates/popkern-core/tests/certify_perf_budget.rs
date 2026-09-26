//! PK-562: blend-evidence seam certification over real G2 fillet
//! output — the pentagonal pyramid apex (bi-nonic setback cap) and the
//! boss plate at G2 — with the exact options the evidence uses and the
//! certificates the inline evidence produced when the fixture was
//! captured (2026-09-26, before any certifier optimisation).
//!
//! Every optimisation of the certifier must keep each seam's status
//! and gate decisions; bit-identical certificates are reported so a
//! numerics change is visible, not silent.
//!
//! The corpus run is a release-only perf budget in the PK-272 pattern
//! (it self-skips in debug; `tools/quality-check.sh --perf` runs it):
//!
//! ```text
//! cargo test --release -p popkern-core --test certify_perf_budget -- --nocapture
//! ```
//!
//! A three-seam sample checks the decisions in every build.
//!
//! # Budgets
//!
//! `measured × 5 headroom`, release, single-threaded, 16-core linux
//! dev reference machine, 2026-09-26:
//!
//! | fixture | seams | intervals | captured (before PK-562) | after PK-562 | budget |
//! |---|---|---|---|---|---|
//! | pyramid5_g2 | 60 | 2878 | 16.3 s | 2.4 s | 12 s |
//! | boss_plate_g2 | 138 | 8890 | 57.2 s | 7.6 s | 38 s |
//!
//! PK-562 mechanisms, each measured on this benchmark: the fraction
//! probe's intervals reused by the full-budget rerun (bit-identical,
//! −26/−33%); the midpoint state reused for its single branch
//! (bit-identical, −30%); allocation-free restriction (−2%); iso-u
//! boundary lines evaluated once per seam (bit-identical, −24/−29%);
//! iso-v lines evaluated v-first and cached (every status and
//! decision kept, interval counts unchanged, −62%). A sign-specialised
//! de Casteljau product was measured 15% slower and dropped.

use std::time::Instant;

use popkern_core::geometry::certify::{SeamCertificate, SeamCertifyOptions};
use popkern_core::topology::BRep;
use popkern_core::topology::EdgeHandle;
use popkern_core::topology::edge_seam_certificate::certify_edge_seam;
use serde::Deserialize;

#[derive(Deserialize)]
struct Seam {
    edge: EdgeHandle,
    certificate: Option<SeamCertificate>,
}

#[derive(Deserialize)]
struct Fixture {
    brep: BRep,
    options: SeamCertifyOptions,
    seams: Vec<Seam>,
}

fn load(name: &str) -> Fixture {
    let path = format!(
        "{}/tests/fixtures/certify/{name}.fixture.json",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("{path}: {e}"))
}

/// The per-quantity decision a certificate makes against its gate:
/// `Some(true)` certified pass, `Some(false)` certified fail, `None`
/// undecided (or not covering the seam).
fn decisions(c: &SeamCertificate, gates: [f64; 3]) -> [Option<bool>; 3] {
    if !c.bounds_cover_seam() {
        return [None; 3];
    }
    let upper = [c.g0_upper, c.g1_upper, c.g2_upper];
    std::array::from_fn(|k| {
        if upper[k] <= gates[k] {
            Some(true)
        } else if c.lower[k] > gates[k] {
            Some(false)
        } else {
            None
        }
    })
}

struct Outcome {
    seams: usize,
    intervals: usize,
    bit_identical: usize,
    seconds: f64,
}

fn certify_all(fixture: &Fixture, take: usize) -> Outcome {
    let gates = fixture
        .options
        .decision_gates
        .expect("evidence options carry decision gates");
    let mut out = Outcome {
        seams: 0,
        intervals: 0,
        bit_identical: 0,
        seconds: 0.0,
    };
    for seam in fixture.seams.iter().take(take) {
        let baseline = seam
            .certificate
            .as_ref()
            .expect("gated seams are certified");
        let started = Instant::now();
        let fresh = certify_edge_seam(&fixture.brep, seam.edge, &fixture.options);
        out.seconds += started.elapsed().as_secs_f64();
        assert_eq!(
            fresh.status, baseline.status,
            "edge {:?}: certificate status changed",
            seam.edge
        );
        assert_eq!(
            decisions(&fresh, gates),
            decisions(baseline, gates),
            "edge {:?}: a gate decision changed",
            seam.edge
        );
        out.seams += 1;
        out.intervals += fresh.intervals;
        out.bit_identical += usize::from(&fresh == baseline);
    }
    out
}

fn report(name: &str, o: &Outcome) {
    println!(
        "[pk562] {name}: {} seams, {} intervals, {:.3} s, {:.3} ms/interval, {}/{} bit-identical",
        o.seams,
        o.intervals,
        o.seconds,
        1e3 * o.seconds / o.intervals.max(1) as f64,
        o.bit_identical,
        o.seams
    );
}

#[test]
fn evidence_certification_corpus_keeps_every_decision() {
    if cfg!(debug_assertions) {
        eprintln!("PK-562 certify perf budget: skipped in debug — run with --release");
        return;
    }
    for (name, budget_s) in [("pyramid5_g2", 12.0), ("boss_plate_g2", 38.0)] {
        let fixture = load(name);
        let outcome = certify_all(&fixture, usize::MAX);
        report(name, &outcome);
        assert_eq!(outcome.seams, fixture.seams.len());
        assert!(
            outcome.seconds <= budget_s,
            "{name}: certification took {:.2} s, budget {budget_s} s",
            outcome.seconds
        );
    }
}

#[test]
fn evidence_certification_sample_keeps_every_decision() {
    let fixture = load("pyramid5_g2");
    let outcome = certify_all(&fixture, 3);
    assert_eq!(outcome.seams, 3);
}
