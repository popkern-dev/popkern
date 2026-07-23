//! PK-144 — fillet mileage corpus: real-part-class models under the
//! full dress-up feature set, with a recorded pass-rate baseline.
//!
//! Each corpus part is a build recipe (profile → booleans → fillets/
//! chamfers) at laptop-stand-class complexity. The runner records a
//! per-part verdict (built + watertight vs declined/failed) and the
//! suite gates on the PASS-RATE BASELINE — a decline is a recorded
//! data point feeding the campaign ledger, not an automatic red, but
//! the rate must never silently drop.

use popkern_part::prelude::*;

struct CorpusVerdict {
    name: &'static str,
    outcome: Result<usize, String>,
}

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
    ]
}

fn run_corpus() -> Vec<CorpusVerdict> {
    corpus()
        .into_iter()
        .map(|(name, build)| {
            let outcome = match build() {
                Err(e) => Err(format!("{e}")),
                Ok(part) if part.brep().faces.is_empty() => {
                    // An empty body meshes to zero edges and every
                    // downstream gate passes VACUOUSLY (measured: the
                    // notched flange's subtract emptied the body and
                    // "passed watertight" — PK-173's real story).
                    Err("EMPTY BODY (vacuous)".to_string())
                }
                Ok(part) => match part.open_edge_count() {
                    Err(e) => Err(format!("open_edge_count: {e}")),
                    Ok(opens) => Ok(opens),
                },
            };
            CorpusVerdict { name, outcome }
        })
        .collect()
}

/// The corpus gate: per-part verdicts printed as the ledger evidence
/// line, and the pass rate held at the recorded baseline.
#[test]
fn corpus_pass_rate_holds_baseline() {
    let verdicts = run_corpus();
    let mut passed = 0usize;
    for v in &verdicts {
        match &v.outcome {
            Ok(0) => {
                passed += 1;
                eprintln!("[corpus] {}: PASS (watertight)", v.name);
            }
            Ok(opens) => {
                eprintln!("[corpus] {}: OPEN ({opens} open edges)", v.name);
            }
            Err(e) => {
                eprintln!(
                    "[corpus] {}: DECLINED/FAILED ({})",
                    v.name,
                    &e[..e.len().min(160)]
                );
            }
        }
    }
    let total = verdicts.len();
    eprintln!("[corpus] pass rate: {passed}/{total}");
    // Baseline history: 3/6 (2026-07-18, PK-144 stage 1); 13/22
    // (2026-07-22, PK-161 — corpus grown to 22 parts); 14/22
    // (2026-07-22, PK-170 — corner-breaching oblique bore served:
    // seam-strip component placement + material-side fold sampling);
    // 15/22 (2026-07-22, PK-170 cross-drill leg — analytic
    // cylinder×cylinder, LSQ fit knots, wrapping-wall chart roles);
    // 16/22 (2026-07-22, PK-171/PK-160 — corner sphere promotion
    // gates made radius-relative; the stale-outlier fit score is
    // r(√2−1) by construction and the 0.5-absolute gate silently
    // stopped promoting exact spheres above r ≈ 1.2); 17/22
    // (2026-07-22, PK-171 — notched_flange chamfer over the
    // interrupted rim: arc_to_nurbs split at parameter-of-angle
    // (rational quadratic parameter ≠ angle inside a 90° span), rim
    // extension candidates ranked by vertex distance, exact NURBS
    // trim slicing + vertex-onto-trim snapping, terminal Cone conic
    // stubs, second edge_reconcile pass after terminal settle);
    // 18/22 (2026-07-23, PK-171 — offset_bracket: obstructed
    // strip-end surgery (interior tangency foot + terminal-arc
    // crossing + ceiling face via split/mef/kve) and seam-slit
    // dissolve after coplanar RedundantEdge face merges in
    // sliver_absorb).
    // Failing
    // classes, each routed to its campaign item:
    //   mounting_bracket — flush-overlap L union (coincident-boolean
    //     class, PK-106 family);
    //   pocketed_tray — pocket-floor concave fillet legs (PK-140);
    //   slotted_block — flush through-slot, 60 opens (coincident);
    //   angled_rib_plate — strip interference decline on the skewed
    //     rib junction (legitimate transactional decline; recorded).
    // Raise when campaigns close classes; NEVER lower — a drop is a
    // regression.
    const BASELINE: usize = 18;
    assert!(
        passed >= BASELINE,
        "corpus pass rate {passed}/{total} fell below the recorded baseline {BASELINE}/{total}"
    );
}

/// The OCC-oracle acceptance leg (PK-161): every corpus part that
/// builds watertight must also pass BRepCheck + free-bounds via the
/// OCP oracle (tools/occ-oracle/occ_check.py). Skips with a notice
/// when the OCP venv is absent (CI without the oracle stays green);
/// with the oracle present, a watertight part failing OCC is a red.
#[test]
fn corpus_occ_oracle_leg() {
    let oracle_python = std::path::Path::new(&std::env::var("HOME").unwrap_or_default())
        .join("venvs/ocp/bin/python");
    let oracle_script = std::path::Path::new("../../tools/occ-oracle/occ_check.py");
    if !oracle_python.exists() || !oracle_script.exists() {
        eprintln!("[corpus-occ] OCP oracle unavailable — leg skipped");
        return;
    }
    let out_dir = std::path::Path::new("../../examples/output/fixtures/corpus");
    std::fs::create_dir_all(out_dir).expect("corpus fixture dir");
    let mut checked = 0usize;
    let mut failures: Vec<String> = Vec::new();
    for (name, build) in corpus() {
        let Ok(part) = build() else { continue };
        if part.brep().faces.is_empty() {
            continue;
        }
        let Ok(0) = part.open_edge_count() else {
            continue;
        };
        let path = out_dir.join(format!("{name}.step"));
        if let Err(e) = part.export_step(&path) {
            failures.push(format!("{name}: export failed: {e}"));
            continue;
        }
        let output = std::process::Command::new(&oracle_python)
            .arg(oracle_script)
            .arg(&path)
            .output()
            .expect("oracle invocation");
        let report = String::from_utf8_lossy(&output.stdout);
        let passed = report.contains("\"passed\": true");
        checked += 1;
        eprintln!(
            "[corpus-occ] {name}: {}",
            if passed { "OCC PASS" } else { "OCC FAIL" }
        );
        if !passed {
            failures.push(format!(
                "{name}: OCC rejected: {}",
                report.chars().take(400).collect::<String>()
            ));
        }
    }
    let occ_passed = checked - failures.len();
    eprintln!("[corpus-occ] {occ_passed}/{checked} watertight part(s) OCC-clean");
    // OCC baseline (2026-07-22, PK-161): 10/14; 14/17 (2026-07-22,
    // PK-170/PK-171 — bore, cross-drill, mixed-corner and
    // notched-flange classes serve and export OCC-clean); 15/18
    // (2026-07-23, PK-171 — offset_bracket obstructed-end class);
    // 18/18 (2026-07-23, PK-172 — torus-band EDGE_CURVE senses: the
    // periodic chain normalization now trusts stored pcurve
    // directions over the assumed chart winding, closing the
    // SelfIntersectingWire + UnorientableShape class on every
    // curved-seam Torus face). No recorded failing classes.
    // Raise when classes close; NEVER lower.
    const OCC_BASELINE: usize = 18;
    for f in &failures {
        eprintln!("[corpus-occ] recorded: {}", &f[..f.len().min(200)]);
    }
    assert!(
        occ_passed >= OCC_BASELINE,
        "OCC pass count {occ_passed}/{checked} fell below the baseline {OCC_BASELINE}"
    );
    assert!(
        checked > 0,
        "oracle present but no watertight parts checked"
    );
}
