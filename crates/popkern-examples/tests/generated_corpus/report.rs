//! Running a kind pair's cells, writing the report, and the known-failure
//! ratchet (design §7–§8).

use std::collections::BTreeMap;
use std::path::PathBuf;

use super::cells::{CONFIGS, CellRun, run_cell};
use super::known::KNOWN;
use super::operand::Kind;
use super::scene::{Cell, Config};

/// The run's seed: `CORPUS_SEED` (test configuration), else 1.
pub fn seed() -> u64 {
    std::env::var("CORPUS_SEED")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1)
}

/// Random scenes per cell: `CORPUS_RANDOM`, else `default`.
pub fn random_count(default: usize) -> usize {
    std::env::var("CORPUS_RANDOM")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

fn report_dir() -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR"))
        .join("corpus")
        .join(format!("seed-{}", seed()));
    std::fs::create_dir_all(&dir).expect("report dir");
    dir
}

/// The cells of the unordered pair `(a, b)`: every configuration, the
/// ordered ones both ways.
pub fn pair_cells(a: Kind, b: Kind) -> Vec<Cell> {
    let mut out = Vec::new();
    for config in CONFIGS {
        out.push(Cell { a, b, config });
        if a != b && matches!(config, Config::I | Config::N | Config::K) {
            out.push(Cell { a: b, b: a, config });
        }
    }
    out
}

/// A failure's class key: cell, oracle, operation.
pub fn class_key(cell: &Cell, oracle: &str, op: &str) -> String {
    format!("{}/{}/{}", cell.name(), oracle, op)
}

/// Run a pair, write its report, and return the unexpected failures and
/// the known entries the run showed fixed.
pub fn run_pair(a: Kind, b: Kind, random: usize, smoke: bool) -> (Vec<String>, Vec<String>) {
    let seed = seed();
    let mut runs: Vec<CellRun> = Vec::new();
    let mut table = String::new();
    for cell in pair_cells(a, b) {
        if !super::feasible::FEASIBLE.contains(&cell.name().as_str()) {
            table.push_str(&format!(
                "{:<14} infeasible (frozen table, feasible.rs)\n",
                cell.name()
            ));
            continue;
        }
        let run = run_cell(cell, seed, random, smoke);
        table.push_str(&format!(
            "{:<14} ladder {:>3} random {:>4} passed {:>4} failed {:>3} rejected {:>4}/{:<4} worst {:.1e}\n",
            cell.name(),
            run.ladder_run,
            run.random_run,
            run.passed,
            run.failed.len(),
            run.rejected,
            run.draws,
            run.worst_rel
        ));
        runs.push(run);
    }
    eprintln!("{table}");
    let tier = if smoke { "smoke" } else { "full" };
    let path = report_dir().join(format!("{tier}-{a:?}x{b:?}.json"));
    std::fs::write(&path, serde_json::to_string(&runs).expect("json")).expect("write report");
    // The ratchet.
    let mut seen: BTreeMap<String, usize> = BTreeMap::new();
    let mut unexpected = Vec::new();
    for run in &runs {
        let cell = run.cell.expect("cell");
        if run.draws > 0 && run.rejected * 2 > run.draws {
            unexpected.push(format!(
                "{}: generator rejects {} of {} draws {:?}",
                cell.name(),
                run.rejected,
                run.draws,
                run.reject_reasons
            ));
        }
        for j in &run.failed {
            for f in &j.verdict.failures {
                let key = class_key(
                    &cell,
                    f.oracle,
                    &f.op.map_or("-".to_string(), |o| format!("{o:?}")),
                );
                *seen.entry(key.clone()).or_default() += 1;
                if !KNOWN.iter().any(|(k, _)| *k == key) {
                    // The first failure of each class is minimized.
                    let first = !unexpected
                        .iter()
                        .any(|u: &String| u.starts_with(&format!("{key}:")));
                    let (scene, reduced) = if first {
                        super::minimize::minimize(&j.scene, f.oracle)
                    } else {
                        (j.scene.clone(), false)
                    };
                    unexpected.push(format!(
                        "{key}: {} [{}replay: CORPUS_SCENE='{}']",
                        f.detail,
                        if reduced {
                            "minimized, not placement luck; "
                        } else {
                            ""
                        },
                        serde_json::to_string(&scene).expect("json")
                    ));
                }
            }
        }
    }
    // Only a full run covers every rung a known class may sit on.
    let fixed = KNOWN
        .iter()
        .filter(|_| !smoke)
        .filter(|(k, _)| {
            pair_cells(a, b)
                .iter()
                .any(|c| k.starts_with(&format!("{}/", c.name())))
        })
        .filter(|(k, _)| !seen.contains_key(*k))
        .map(|(k, task)| format!("{k} ({task}) passed this run: remove it from KNOWN"))
        .collect();
    (unexpected, fixed)
}

/// Run the sheet cells of a solid kind (both sheets, every
/// configuration), returning the unexpected failures (the same ratchet:
/// keys `<sheet>x<solid>-<config>/<oracle>/-`).
pub fn run_sheets(solid: Kind, random: usize) -> Vec<String> {
    use super::sheets::{SheetConfig, SheetKind, draw_sheet, judge_sheet};
    let seed = seed();
    let mut unexpected = Vec::new();
    let mut table = String::new();
    for sheet in [SheetKind::Quad, SheetKind::Dome] {
        for config in [SheetConfig::Cross, SheetConfig::Tangent, SheetConfig::Apart] {
            let name = format!("{sheet:?}x{solid:?}-{config:?}");
            let (mut run, mut passed, mut rejected) = (0, 0, 0);
            for k in 0..(4 + random as u64) {
                let Some(scene) =
                    (0..24).find_map(|t| draw_sheet(sheet, solid, config, seed, k * 1009 + t))
                else {
                    rejected += 1;
                    continue;
                };
                run += 1;
                let verdict = judge_sheet(&scene);
                if verdict.failures.is_empty() {
                    passed += 1;
                }
                for f in &verdict.failures {
                    let key = format!("{name}/{}/-", f.oracle);
                    if !KNOWN.iter().any(|(k, _)| *k == key) {
                        unexpected.push(format!(
                            "{key}: {} [scene: {}]",
                            f.detail,
                            serde_json::to_string(&scene).expect("json")
                        ));
                    }
                }
            }
            table.push_str(&format!(
                "{name:<22} run {run:>4} passed {passed:>4} rejected {rejected:>3}\n"
            ));
        }
    }
    eprintln!("{table}");
    unexpected
}

/// Run an explicit list of `(cell, rung)` scenes (the smoke tier),
/// returning the unexpected failures under the ratchet.
pub fn run_listed(list: &[(Cell, super::scene::Rung)]) -> Vec<String> {
    use super::judge::judge;
    use super::scene::draw;
    let seed = seed();
    let mut unexpected = Vec::new();
    for (i, (cell, rung)) in list.iter().enumerate() {
        let Some(scene) = (0..24).find_map(|t| draw(*cell, seed, i as u64 * 1009 + t, *rung).ok())
        else {
            unexpected.push(format!(
                "{}: no scene drawn for {rung:?} (generator defect)",
                cell.name()
            ));
            continue;
        };
        for f in judge(&scene).failures {
            let key = class_key(
                cell,
                f.oracle,
                &f.op.map_or("-".to_string(), |o| format!("{o:?}")),
            );
            if !KNOWN.iter().any(|(k, _)| *k == key) {
                unexpected.push(format!(
                    "{key}: {} [replay: CORPUS_SCENE='{}']",
                    f.detail,
                    serde_json::to_string(&scene).expect("json")
                ));
            }
        }
    }
    unexpected
}
