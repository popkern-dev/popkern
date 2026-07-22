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

fn build_countersunk_plate() -> Result<Part, PartError> {
    // Perpendicular countersunk bore (the PK-147 class): through-bore
    // plus a coaxial drafted-frustum sink crossing the bore wall.
    let plate = Profile::on(Plane::XY).rect(60.0, 60.0).extrude(30.0)?;
    let bore = Profile::on(Plane::xy_at(-1.0)).circle(10.0).extrude(32.0)?;
    let plate = plate.subtract(&bore)?;
    let sink = Profile::on(Plane::xy_at(24.0))
        .circle(8.0)
        .extrude(8.0)?
        .draft(-30.0_f64.to_radians(), Dir::Z, |f: FaceQuery| {
            f.kind(SurfaceKind::Cylinder)
        })?;
    plate.subtract(&sink)
}

fn build_oblique_countersunk_plate() -> Result<Part, PartError> {
    // The PK-166 stand class: cone countersink through a tilted rest
    // plane over a blind seat.
    let plate = Profile::on(Plane::XY).rect(200.0, 200.0).extrude(140.0)?;
    let shave = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)?
        .rotate_x(18.0_f64.to_radians())?
        .translate(0.0, 0.0, 100.0)?;
    let plate = plate.subtract(&shave)?;
    let seat = Profile::on(Plane::xy_at(70.0)).circle(35.2).extrude(60.0)?;
    let plate = plate.subtract(&seat)?;
    let sink = Profile::on(Plane::xy_at(82.0))
        .circle(34.6)
        .extrude(36.0)?
        .draft(-8.5_f64.to_radians(), Dir::Z, |f: FaceQuery| {
            f.kind(SurfaceKind::Cylinder)
        })?;
    plate.subtract(&sink)
}

fn build_drafted_cup() -> Result<Part, PartError> {
    // Drafted cylinder wall — the exact-frustum class.
    Profile::on(Plane::XY).circle(30.0).extrude(40.0)?.draft(
        6.0_f64.to_radians(),
        Dir::Z,
        |f: FaceQuery| f.kind(SurfaceKind::Cylinder),
    )
}

fn build_cross_drilled_boss() -> Result<Part, PartError> {
    // Radial hole through a cylinder wall (cylinder-x-cylinder rims).
    let body = Profile::on(Plane::XY).circle(16.0).extrude(40.0)?;
    let drill = Profile::on(Plane::XY)
        .circle(6.0)
        .extrude(80.0)?
        .rotate_y(90.0_f64.to_radians())?
        .translate(-40.0, 0.0, 20.0)?;
    body.subtract(&drill)
}

fn build_coaxial_bushing() -> Result<Part, PartError> {
    // Tube: both walls full-period periodic faces.
    let outer = Profile::on(Plane::XY).circle(20.0).extrude(25.0)?;
    let bore = Profile::on(Plane::xy_at(-1.0)).circle(8.0).extrude(27.0)?;
    outer.subtract(&bore)
}

fn build_tilted_hole_plate() -> Result<Part, PartError> {
    // Oblique through-hole: cylinder-x-plane ellipse rims both faces.
    let plate = Profile::on(Plane::XY).rect(60.0, 60.0).extrude(15.0)?;
    let drill = Profile::on(Plane::XY)
        .circle(6.0)
        .extrude(80.0)?
        .rotate_x(25.0_f64.to_radians())?
        .translate(0.0, -15.0, -20.0)?;
    plate.subtract(&drill)
}

fn build_offset_bracket() -> Result<Part, PartError> {
    // L-bracket with the wall OVERLAPPING the base interior (not
    // flush) so the union stays off the coincident class; fillet the
    // outer verticals.
    let base = Profile::on(Plane::XY).rect(60.0, 40.0).extrude(8.0)?;
    let wall = Profile::on(Plane::xy_at(4.0))
        .at(-24.0, 0.0)
        .rect(8.0, 40.0)
        .extrude(35.0)?;
    let body = base.union(&wall)?;
    body.fillet(4.0, |e| {
        e.straight()
            .direction(Dir::Z)
            .chord_length_between(6.0, 10.0)
            .allow_empty()
    })
}

fn build_rim_chamfered_disc() -> Result<Part, PartError> {
    // Large-leg rim chamfer: the conical band class.
    Profile::on(Plane::XY)
        .circle(44.0)
        .extrude(20.0)?
        .chamfer(3.0, |e| e.circular().at_height(20.0))
}

fn build_star_socket_plate() -> Result<Part, PartError> {
    // Hexagram socket through a plate (the star-membrane class).
    let star: Vec<(f64, f64)> = (0..12)
        .map(|i| {
            let a = std::f64::consts::PI / 6.0 * i as f64;
            let r = if i % 2 == 0 { 14.0 } else { 7.0 };
            (r * a.cos(), r * a.sin())
        })
        .collect();
    let plate = Profile::on(Plane::XY).rect(60.0, 60.0).extrude(12.0)?;
    let tool = Profile::on(Plane::xy_at(-1.0))
        .polygon(&star)
        .extrude(14.0)?;
    plate.subtract(&tool)
}

fn build_mixed_corner_block() -> Result<Part, PartError> {
    // Three fillets meeting at one cube corner in a single operation
    // — the V-free DefaultOn corner class via the Part facade.
    let block = Profile::on(Plane::XY).rect(30.0, 30.0).extrude(30.0)?;
    block.fillet(3.0, |e| {
        e.straight().filter(|i| {
            let near_x = i.midpoint.x > 14.9;
            let near_y = i.midpoint.y > 14.9;
            let near_z = i.midpoint.z > 29.9;
            (near_x && near_y) || (near_x && near_z) || (near_y && near_z)
        })
    })
}

fn build_chamfered_peg() -> Result<Part, PartError> {
    // Cylinder peg with chamfer lead-ins on both end rims (PK-135).
    Profile::on(Plane::XY)
        .circle(10.0)
        .extrude(30.0)?
        .chamfer(1.5, |e| e.circular().at_height(30.0))?
        .chamfer(1.5, |e| e.circular().at_height(0.0))
}

fn build_skewed_dress_bar() -> Result<Part, PartError> {
    // Yawed+tilted boolean-cut rest plane with side chamfers — the
    // PK-168 skewed terminal-chord class.
    let bar = Profile::on(Plane::XY).rect(40.0, 60.0).extrude(96.0)?;
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)?
        .rotate_x(18.0_f64.to_radians())?
        .rotate_z(-12.0_f64.to_radians())?
        .translate(0.0, 0.0, 80.0)?;
    let bar = bar.subtract(&upper)?;
    let tilt = 18.0_f64.to_radians();
    let w = -12.0_f64.to_radians();
    let rest_z = move |x: f64, y: f64| 80.0 + (-x * w.sin() + y * w.cos()) * tilt.tan();
    bar.chamfer(1.0, move |e| {
        e.straight().filter(move |i| {
            i.direction.map(|d| d.y.abs() > 0.85).unwrap_or(false)
                && (i.midpoint.z - rest_z(i.midpoint.x, i.midpoint.y)).abs() < 0.5
        })
    })
}

fn build_counterbored_plate() -> Result<Part, PartError> {
    // Bore + larger coaxial counterbore: concentric rims on the top
    // face (the nested-hole walk class).
    let plate = Profile::on(Plane::XY).rect(50.0, 50.0).extrude(20.0)?;
    let bore = Profile::on(Plane::xy_at(-1.0)).circle(5.0).extrude(22.0)?;
    let plate = plate.subtract(&bore)?;
    let cbore = Profile::on(Plane::xy_at(12.0)).circle(9.0).extrude(9.0)?;
    plate.subtract(&cbore)
}

fn build_angled_rib_plate() -> Result<Part, PartError> {
    // Diagonal stiffening rib; concave junction fillet with a
    // non-axis-aligned direction filter.
    let plate = Profile::on(Plane::XY).rect(80.0, 50.0).extrude(5.0)?;
    let rib = Profile::on(Plane::xy_at(5.0))
        .rect(70.0, 6.0)
        .extrude(10.0)?
        .rotate_z(20.0_f64.to_radians())?;
    let body = plate.union(&rib)?;
    body.fillet(1.5, |e| {
        e.straight()
            .chord_length_between(40.0, 80.0)
            .filter(|i| (i.midpoint.z - 5.0).abs() < 0.1)
            .allow_empty()
    })
}

fn build_double_boss_plate() -> Result<Part, PartError> {
    // Two cylindrical bosses, both base seams filleted (multi-seam
    // curved fillets in one operation).
    let plate = Profile::on(Plane::XY).rect(80.0, 40.0).extrude(10.0)?;
    let boss_a = Profile::on(Plane::xy_at(10.0))
        .at(-20.0, 0.0)
        .circle(8.0)
        .extrude(12.0)?;
    let boss_b = Profile::on(Plane::xy_at(10.0))
        .at(20.0, 0.0)
        .circle(6.0)
        .extrude(15.0)?;
    let body = plate.union(&boss_a)?.union(&boss_b)?;
    body.fillet(2.0, |e| {
        e.circular()
            .between(SurfaceKind::Cylinder, SurfaceKind::Plane)
            .allow_empty()
    })
}

fn build_notched_flange() -> Result<Part, PartError> {
    // Turned flange with a keyway notch through the rim, then a rim
    // chamfer over the interrupted circle (chain class).
    let flange = Profile::on(Plane::XY).circle(25.0).extrude(8.0)?;
    let key = Profile::on(Plane::xy_at(-1.0))
        .at(25.0, 0.0)
        .rect(10.0, 8.0)
        .extrude(10.0)?;
    let notched = flange.subtract(&key)?;
    notched.chamfer(1.0, |e| e.circular().at_height(8.0).allow_empty())
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
    // Baseline history: 3/6 (2026-07-18, PK-144 stage 1); 14/22
    // (2026-07-22, PK-161 — corpus grown to 22 parts). Failing
    // classes, each routed to its campaign item:
    //   mounting_bracket — flush-overlap L union (coincident-boolean
    //     class, PK-106 family);
    //   pocketed_tray — pocket-floor concave fillet legs (PK-140);
    //   slotted_block — flush through-slot, 60 opens (coincident);
    //   cross_drilled_boss — Part rotate on a cylinder body breaks
    //     loop closure (PK-169);
    //   tilted_hole_plate — oblique off-center through-hole, 821
    //     opens (PK-170);
    //   offset_bracket — fillet on an overlapping union serves an
    //     edge 1.66 off-surface, check-9 catch (PK-171);
    //   mixed_corner_block — V-free corner class at 30x scale serves
    //     0.21 off-surface, check-9 catch (PK-171);
    //   angled_rib_plate — strip interference decline on the skewed
    //     rib junction (legitimate transactional decline; recorded);
    //   notched_flange — the rim-overlapping key subtract EMPTIES the
    //     body; formerly a vacuous green (empty mesh = zero opens),
    //     now an explicit failure (PK-173).
    // Raise when campaigns close classes; NEVER lower — a drop is a
    // regression.
    const BASELINE: usize = 13;
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
    // OCC baseline (2026-07-22, PK-161): 10/14. Recorded failing
    // classes, routed to their campaign items:
    //   bossed_housing / turned_bushing / double_boss_plate — the
    //     curved-seam TORUS strip exports a wire OCC rejects
    //     (SelfIntersectingWire + UnorientableShape; PK-172);
    // (The notched flange's former 'empty CLOSED_SHELL' entry was a
    // vacuous empty body — excluded above, tracked as PK-173.)
    // Raise when classes close; NEVER lower.
    const OCC_BASELINE: usize = 10;
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
