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

fn build_mounting_bracket() -> Result<Part, PartError> {
    // L-bracket: base plate + upright wall, mounting holes, filleted
    // outer verticals and chamfered hole rims.
    let base = Profile::on(Plane::XY).rect(60.0, 40.0).extrude(6.0)?;
    let wall = Profile::on(Plane::XY)
        .at(-27.0, 0.0)
        .rect(6.0, 40.0)
        .extrude(35.0)?;
    let body = base.union(&wall)?;
    let hole = Profile::on(Plane::XY)
        .at(15.0, 0.0)
        .circle(4.0)
        .extrude(6.0)?;
    let holed = body.subtract(&hole)?;
    holed.fillet(3.0, |e| {
        e.straight()
            .direction(Dir::Z)
            .chord_length_between(5.0, 7.0)
    })
}

fn build_ribbed_plate() -> Result<Part, PartError> {
    // Plate with a stiffening rib; fillet the rib-plate junction
    // (the classic concave line fillet).
    let plate = Profile::on(Plane::XY).rect(80.0, 50.0).extrude(5.0)?;
    let rib = Profile::on(Plane::xy_at(5.0))
        .rect(80.0, 6.0)
        .extrude(12.0)?;
    let body = plate.union(&rib)?;
    body.fillet(2.0, |e| {
        e.straight()
            .direction(Dir::X)
            .chord_length_between(70.0, 90.0)
            .allow_empty()
    })
}

fn build_bossed_housing() -> Result<Part, PartError> {
    // Housing block with a cylindrical boss on top; fillet the boss
    // base seam (curved-seam fillet — the flagship class).
    let block = Profile::on(Plane::XY).rect(50.0, 50.0).extrude(20.0)?;
    let boss = Profile::on(Plane::xy_at(20.0)).circle(10.0).extrude(15.0)?;
    let body = block.union(&boss)?;
    body.fillet(3.0, |e| {
        e.circular()
            .between(SurfaceKind::Cylinder, SurfaceKind::Plane)
            .allow_empty()
    })
}

fn build_turned_bushing() -> Result<Part, PartError> {
    // Turned part: cylinder stack with a rim fillet through the
    // periodic seam vertex (the PK-141 rim class).
    let shaft = Profile::on(Plane::XY).circle(12.0).extrude(30.0)?;
    let flange = Profile::on(Plane::XY).circle(20.0).extrude(6.0)?;
    let body = flange.union(&shaft)?;
    body.fillet(2.0, |e| {
        e.circular()
            .between(SurfaceKind::Cylinder, SurfaceKind::Plane)
            .allow_empty()
    })
}

fn build_pocketed_tray() -> Result<Part, PartError> {
    // Pocketed tray with chamfered outer verticals (the bracket
    // example's class) plus a pocket-floor fillet attempt.
    let plate = Profile::on(Plane::XY).rect(60.0, 36.0).extrude(10.0)?;
    let tool = Profile::on(Plane::xy_at(4.0))
        .rect(40.0, 20.0)
        .extrude(10.0)?;
    let pocketed = plate.subtract(&tool)?;
    let chamfered = pocketed.chamfer(2.5, |e| {
        e.straight()
            .direction(Dir::Z)
            .chord_length_between(8.0, 12.0)
            .allow_empty()
    })?;
    chamfered.fillet(1.5, |e| {
        e.straight()
            .direction(Dir::X)
            .chord_length_between(35.0, 45.0)
            .allow_empty()
    })
}

fn build_slotted_block() -> Result<Part, PartError> {
    // Through-slot interrupting long edges, then chamfer the long
    // top edges (the slot-interrupted class — recorded survivors).
    let block = Profile::on(Plane::XY).rect(70.0, 30.0).extrude(12.0)?;
    let slot = Profile::on(Plane::XY).rect(12.0, 30.0).extrude(12.0)?;
    let slotted = block.subtract(&slot)?;
    slotted.chamfer(1.5, |e| {
        e.straight()
            .direction(Dir::X)
            .chord_length_between(20.0, 80.0)
            .allow_empty()
    })
}

type CorpusBuild = fn() -> Result<Part, PartError>;

fn corpus() -> Vec<(&'static str, CorpusBuild)> {
    vec![
        ("mounting_bracket", build_mounting_bracket),
        ("ribbed_plate", build_ribbed_plate),
        ("bossed_housing", build_bossed_housing),
        ("turned_bushing", build_turned_bushing),
        ("pocketed_tray", build_pocketed_tray),
        ("slotted_block", build_slotted_block),
    ]
}

fn run_corpus() -> Vec<CorpusVerdict> {
    corpus()
        .into_iter()
        .map(|(name, build)| {
            let outcome = match build() {
                Err(e) => Err(format!("{e}")),
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
    // Baseline recorded 2026-07-18 (PK-144 stage 1): 3/6. The three
    // failing classes, each feeding its campaign:
    //   mounting_bracket — flush-overlap L union fails outright
    //     (loop gap = the shared wall span; coincident-boolean class,
    //     PK-106 family);
    //   pocketed_tray — pocket-floor concave fillet, 2 of 4 legs
    //     fail and the operation declines transactionally (PK-140);
    //   slotted_block — flush through-slot subtract ships 60 open
    //     edges (coincident-boolean class).
    // Raise when campaigns close classes; NEVER lower — a drop is a
    // regression.
    const BASELINE: usize = 3;
    assert!(
        passed >= BASELINE,
        "corpus pass rate {passed}/{total} fell below the recorded baseline {BASELINE}/{total}"
    );
}
