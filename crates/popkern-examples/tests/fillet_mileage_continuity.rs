//! PK-188 — corpus-wide continuity evidence with a recorded baseline.
//!
//! `fillet_mileage_corpus` gates watertightness; this suite gates the
//! MEASURED continuity of every fillet operation across the same corpus
//! parts: per-seam G0/G1/G2 residuals (transverse-aware — see
//! `popkern_blend::blend_continuity_evidence`), achieved level
//! classified against the commercial gates, and requested-vs-achieved
//! so a silent downgrade is a test failure, not a stderr line.
//!
//! The evidence is structured output asserted here (no print-trace):
//! the printed table is the ledger evidence line, the assertions are
//! the gate.

use popkern_blend::ContinuityLevel;
use popkern_part::prelude::*;

#[path = "fillet_mileage_corpus_support/parts.rs"]
mod parts;
use parts::*;

type CorpusBuild = fn() -> Result<Part, PartError>;

fn corpus() -> Vec<(&'static str, CorpusBuild)> {
    vec![
        ("mounting_bracket", build_mounting_bracket),
        ("ribbed_plate", build_ribbed_plate),
        ("bossed_housing", build_bossed_housing),
        ("turned_bushing", build_turned_bushing),
        ("pocketed_tray", build_pocketed_tray),
        ("slotted_block", build_slotted_block),
        ("countersunk_plate", build_countersunk_plate),
        ("oblique_countersunk_plate", build_oblique_countersunk_plate),
        ("drafted_cup", build_drafted_cup),
        ("cross_drilled_boss", build_cross_drilled_boss),
        ("coaxial_bushing", build_coaxial_bushing),
        ("tilted_hole_plate", build_tilted_hole_plate),
        ("offset_bracket", build_offset_bracket),
        ("rim_chamfered_disc", build_rim_chamfered_disc),
        ("star_socket_plate", build_star_socket_plate),
        ("mixed_corner_block", build_mixed_corner_block),
        ("chamfered_peg", build_chamfered_peg),
        ("skewed_dress_bar", build_skewed_dress_bar),
        ("counterbored_plate", build_counterbored_plate),
        ("angled_rib_plate", build_angled_rib_plate),
        ("double_boss_plate", build_double_boss_plate),
        ("notched_flange", build_notched_flange),
        ("pentagon_tower", build_pentagon_tower),
    ]
}

/// Aggregated continuity record for one corpus part.
struct PartContinuity {
    name: &'static str,
    /// `Err` when the part itself declined to build (the pass-rate
    /// baseline in the sibling suite owns that story).
    outcome: Result<PartEvidence, String>,
}

#[derive(Debug, Default)]
struct PartEvidence {
    /// Fillet operations that produced continuity evidence.
    fillet_reports: usize,
    /// Weakest achieved level across those operations.
    min_achieved: Option<ContinuityLevel>,
    /// Reports where achieved < requested.
    downgrades: usize,
    /// Total gated (intended-smooth) seams.
    gated_seams: usize,
    /// Total intended-sharp seams (census-pinned — a rail that
    /// degrades past the intent bound moves between these counts).
    sharp_seams: usize,
    /// Seams that could not be measured (must stay 0).
    unmeasured_seams: usize,
    worst_g0: f64,
    worst_g1: f64,
    worst_g2: f64,
    /// Gated seams whose certificate decides every gate (PK-499).
    certified_seams: usize,
    /// Gated seams classified by the sampled fallback.
    uncertified_seams: usize,
}

fn run_corpus_continuity() -> Vec<PartContinuity> {
    corpus()
        .into_iter()
        .map(|(name, build)| {
            let outcome = match build() {
                Err(e) => Err(format!("{e}")),
                Ok(part) => {
                    let mut ev = PartEvidence::default();
                    for report in part.blend_reports() {
                        let Some(continuity) = report.continuity.as_ref() else {
                            continue; // chamfer: designed crease, no claim
                        };
                        ev.fillet_reports += 1;
                        ev.min_achieved = Some(match ev.min_achieved {
                            Some(m) if m <= continuity.achieved => m,
                            _ => continuity.achieved,
                        });
                        if continuity.downgraded() {
                            ev.downgrades += 1;
                        }
                        ev.gated_seams += continuity.seams.iter().filter(|s| s.gated).count();
                        ev.sharp_seams += continuity.sharp_seams;
                        ev.unmeasured_seams += continuity.unmeasured_seams;
                        ev.worst_g0 = ev.worst_g0.max(continuity.worst_g0_gap);
                        ev.worst_g1 = ev.worst_g1.max(continuity.worst_g1_normal);
                        ev.worst_g2 = ev.worst_g2.max(continuity.worst_g2_residual);
                        ev.certified_seams += continuity.certified_seams;
                        ev.uncertified_seams += continuity.uncertified_seams;
                    }
                    Ok(ev)
                }
            };
            PartContinuity { name, outcome }
        })
        .collect()
}

/// Per-part baseline: weakest allowed achieved level, minimum gated-
/// seam census, residual caps, and the expected downgrade count.
///
/// Pinned from the measured 2026-07-27 baseline run; caps carry ~2x
/// headroom over the measured worst so sampling noise cannot flake,
/// while a level regression, a census drop, or an order-of-magnitude
/// residual regression always trips.
///
/// The honest state this pins (the point of PK-188 — recorded, not
/// hidden):
/// - `ribbed_plate` delivers a genuinely G1 straight-seam fillet
///   (g1 = 4.7e-9).
/// - The circular-seam parts (`bossed_housing`, `turned_bushing`,
///   `double_boss_plate`) read genuinely G1 (g0 ~ 2e-11..6e-11,
///   g1 ~ 2e-10) since PK-210: their boss-ring rail pcurves were
///   23-sample cubic interpolations leaving 2.5e-4..9.2e-4 of
///   sagitta (above the 1e-5 trim gate → achieved G0, silent
///   downgrades); an Arc in a Plane face now gets its exact
///   angle-parameterized elliptical pcurve.
/// - `offset_bracket` reads genuinely G1 since PK-211 (g0 2.8e-11 —
///   its 1.76-displaced strip-end pcurves and phantom 0.43 rail
///   error are gone); `mixed_corner_block` reads genuinely G1 since
///   PK-338 (g0 1.1e-11 — its 1.5e-2 trim deviation was the pcurve
///   endpoint surgery snapping to marginal vertex projections).
struct PartBaseline {
    name: &'static str,
    achieved_floor: ContinuityLevel,
    min_gated_seams: usize,
    g0_cap: f64,
    g1_cap: f64,
    expected_downgrades: usize,
    /// Gated seams whose PK-499 certificate decides every gate.
    min_certified_seams: usize,
    /// Gated seams left to the sampled fallback: none. The
    /// closed boss/bushing ring rails decide every gate since PK-542
    /// (phase-shifted correspondence, merged affine polyline pcurves,
    /// fraction continuation when witnesses coincide). Seams
    /// into a sphere's chart pole certify chart-free since PK-543
    /// (mixed_corner_block 7 -> 9, pentagon_tower 20 -> 40 certified).
    /// Seams on a dual-jet corner piece's boundary certify through its
    /// polynomial boundary surrogate since PK-541 (pentagon_tower
    /// 40 -> 60).
    max_uncertified_seams: usize,
}

const BASELINE: &[PartBaseline] = &[
    PartBaseline {
        name: "ribbed_plate",
        achieved_floor: ContinuityLevel::G1,
        min_gated_seams: 16,
        g0_cap: 1.0e-10,
        g1_cap: 1.0e-8,
        expected_downgrades: 0,
        min_certified_seams: 16,
        max_uncertified_seams: 0,
    },
    PartBaseline {
        name: "bossed_housing",
        achieved_floor: ContinuityLevel::G1,
        min_gated_seams: 2,
        g0_cap: 1.0e-10,
        g1_cap: 1.0e-9,
        expected_downgrades: 0,
        min_certified_seams: 4,
        max_uncertified_seams: 0,
    },
    PartBaseline {
        name: "turned_bushing",
        achieved_floor: ContinuityLevel::G1,
        min_gated_seams: 6,
        g0_cap: 2.0e-10,
        g1_cap: 1.0e-9,
        expected_downgrades: 0,
        min_certified_seams: 8,
        max_uncertified_seams: 0,
    },
    PartBaseline {
        // PK-211 closed: the strip-end ruling pcurves were fit
        // against transitional mid-carve geometry and survived every
        // guard (g0 read 1.76, g1 a phantom 0.43 measured at the
        // displaced image); the end-of-pipeline displaced-pcurve
        // backstop recomputes them on the final geometry — measured
        // g0 2.8e-11, g1 exact zero.
        name: "offset_bracket",
        achieved_floor: ContinuityLevel::G1,
        min_gated_seams: 8,
        g0_cap: 1.0e-10,
        g1_cap: 1.0e-9,
        expected_downgrades: 0,
        min_certified_seams: 8,
        max_uncertified_seams: 0,
    },
    PartBaseline {
        // Genuinely G1 since PK-338 (measured g0 1.1e-11, g1 6.6e-10,
        // 0 downgrades): its long-standing 1.5e-2 trim deviation was
        // the pcurve endpoint surgery dragging control points to
        // marginal vertex projections — the model-tolerance-scaled
        // 3D bypass in `fix_pcurve_endpoints` retired the class.
        name: "mixed_corner_block",
        achieved_floor: ContinuityLevel::G1,
        min_gated_seams: 9,
        g0_cap: 1.0e-10,
        g1_cap: 1.0e-8,
        expected_downgrades: 0,
        min_certified_seams: 9,
        max_uncertified_seams: 0,
    },
    PartBaseline {
        // BLUNT (18°-half-angle) apex under a G2 request through the
        // public API: the ladder dispatches, the fan ACCEPTS under
        // the honest transverse postcondition, the part serves with
        // its accepted inexactness recorded as entity tolerances
        // (g0 9.5e-4 sits inside the recorded band — the evidence's
        // tolerance-widened count keeps it visible), and it meshes
        // watertight (PK-208 closed). The recorded downgrade is the
        // rim-arc tangency (PK-212: 0.119 → 6.7e-2 via the tangent-
        // plane border prolongation + rim cross-jet reconciliation;
        // the residual is the pinned mid-arc bump, gated kernel-side
        // in kp_polar_commercial_gates::rim_tangency). PK-399: the
        // setback corners are now exact trim crossings, so the
        // off-strip corner point no longer splits each strip trim at
        // an extra vertex (measured on the blunt pentagon: faces 31
        // unchanged, edges 80 → 70, vertices 51 → 41 — the short
        // cap↔strip seam stubs next to every corner are gone, no
        // seam reclassified): the census is 60 (was 65), the rim
        // tangency 3.8e-3 (was 6.7e-2), achieved G1.
        // PK-575 (2026-09-27): the G2 request now MEASURES G2. Rim
        // cluster fields and spoke-junction dual-jet corners took the
        // blunt fan's rim transverse curvature from 0.30 to 7e-4.
        // Measured: achieved G2, g0 1.6e-8, g1 7.3e-6, g2 1.1e-3, all 60
        // gated seams certified.
        name: "pentagon_tower",
        achieved_floor: ContinuityLevel::G2,
        min_gated_seams: 60,
        g0_cap: 1.0e-7,
        g1_cap: 1.0e-5,
        expected_downgrades: 0,
        min_certified_seams: 60,
        max_uncertified_seams: 0,
    },
    PartBaseline {
        name: "double_boss_plate",
        achieved_floor: ContinuityLevel::G1,
        min_gated_seams: 4,
        g0_cap: 1.0e-10,
        g1_cap: 1.0e-9,
        expected_downgrades: 0,
        min_certified_seams: 8,
        max_uncertified_seams: 0,
    },
];

#[test]
fn corpus_continuity_evidence_holds_baseline() {
    let records = run_corpus_continuity();
    for r in &records {
        match &r.outcome {
            Ok(ev) if ev.fillet_reports == 0 => {
                eprintln!("[continuity] {}: no fillet evidence (chamfer-only)", r.name);
            }
            Ok(ev) => {
                eprintln!(
                    "[continuity] {}: ops={} achieved={:?} downgrades={} gated={} sharp={} \
                     unmeasured={} g0={:.3e} g1={:.3e} g2={:.3e} certified={} uncertified={}",
                    r.name,
                    ev.fillet_reports,
                    ev.min_achieved,
                    ev.downgrades,
                    ev.gated_seams,
                    ev.sharp_seams,
                    ev.unmeasured_seams,
                    ev.worst_g0,
                    ev.worst_g1,
                    ev.worst_g2,
                    ev.certified_seams,
                    ev.uncertified_seams,
                );
            }
            Err(e) => {
                eprintln!(
                    "[continuity] {}: DECLINED ({})",
                    r.name,
                    &e[..e.len().min(120)]
                );
            }
        }
    }

    for b in BASELINE {
        let name = b.name;
        let record = records
            .iter()
            .find(|r| r.name == name)
            .unwrap_or_else(|| panic!("baseline names unknown part {name}"));
        let ev = match &record.outcome {
            Ok(ev) => ev,
            Err(e) => panic!("[{name}] baselined part declined to build: {e}"),
        };
        let achieved = ev
            .min_achieved
            .unwrap_or_else(|| panic!("[{name}] baselined part produced no fillet evidence"));
        assert!(
            achieved >= b.achieved_floor,
            "[{name}] achieved continuity regressed: {achieved:?} < baseline {:?}",
            b.achieved_floor
        );
        assert!(
            ev.gated_seams >= b.min_gated_seams,
            "[{name}] gated-seam census dropped: {} < baseline {} — either a rail \
             degraded past the smooth-intent bound or seams vanished",
            ev.gated_seams,
            b.min_gated_seams
        );
        assert!(
            ev.worst_g0 <= b.g0_cap,
            "[{name}] worst G0 trim deviation regressed: {:.3e} > cap {:.1e}",
            ev.worst_g0,
            b.g0_cap
        );
        assert!(
            ev.worst_g1 <= b.g1_cap,
            "[{name}] worst G1 residual regressed: {:.3e} > cap {:.1e}",
            ev.worst_g1,
            b.g1_cap
        );
        assert!(
            ev.certified_seams >= b.min_certified_seams,
            "[{name}] certified-seam census dropped: {} < baseline {}",
            ev.certified_seams,
            b.min_certified_seams
        );
        assert!(
            ev.uncertified_seams <= b.max_uncertified_seams,
            "[{name}] sampled-fallback seams grew: {} > baseline {}",
            ev.uncertified_seams,
            b.max_uncertified_seams
        );
        assert_eq!(
            ev.unmeasured_seams, 0,
            "[{name}] a blend seam became unmeasurable"
        );
        assert_eq!(
            ev.downgrades, b.expected_downgrades,
            "[{name}] downgrade count moved (baseline {}): a fillet either silently \
             downgraded below its request or silently un-downgraded (re-pin with evidence)",
            b.expected_downgrades
        );
    }
    let evidencing = records
        .iter()
        .filter(|r| matches!(&r.outcome, Ok(ev) if ev.fillet_reports > 0))
        .count();
    assert!(
        evidencing >= BASELINE.len(),
        "parts producing fillet evidence dropped below the baselined census: {evidencing} < {}",
        BASELINE.len()
    );
}
