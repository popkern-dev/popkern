//! Cells and the ladder (design §4–§5), and running a cell.

use popkern_core::tolerance::Tolerance;
use serde::Serialize;

use super::geometry::Steer;
use super::judge::{Verdict, judge_with};
use super::operand::{KINDS, Kind};
use super::scene::{Cell, Config, OFFSETS, Reject, Rung, Scene, draw};

/// The configurations, in table order.
pub const CONFIGS: [Config; 9] = [
    Config::D,
    Config::E,
    Config::I,
    Config::V,
    Config::K,
    Config::X,
    Config::X3,
    Config::N,
    Config::O,
];

/// Whether the configuration depends on which operand is first (B
/// inside A): those run over ordered pairs, the rest over unordered
/// (both differences are judged either way).
fn ordered(config: Config) -> bool {
    matches!(config, Config::I | Config::N | Config::K)
}

/// Every candidate cell; feasibility is decided by trial draws
/// (`feasible`), not by hand rules.
pub fn candidate_cells() -> Vec<Cell> {
    let mut out = Vec::new();
    for config in CONFIGS {
        for (i, &a) in KINDS.iter().enumerate() {
            for (j, &b) in KINDS.iter().enumerate() {
                if !ordered(config) && j < i {
                    continue;
                }
                out.push(Cell { a, b, config });
            }
        }
    }
    out
}

/// The kinds whose relevant chart has a seam and those with poles.
fn has_seam(k: Kind) -> bool {
    !matches!(k, Kind::P | Kind::W | Kind::L | Kind::G)
}

fn has_pole(k: Kind) -> bool {
    matches!(k, Kind::S | Kind::B)
}

fn touching(config: Config) -> bool {
    matches!(config, Config::E | Config::I | Config::V | Config::K)
}

/// The deterministic ladder of a cell (design §5).
pub fn ladder(cell: Cell) -> Vec<Rung> {
    let tol = Tolerance::default().linear();
    let mut out = vec![Rung::FREE];
    if touching(cell.config) {
        for &o in &OFFSETS[1..] {
            out.push(Rung {
                offset: o,
                ..Rung::FREE
            });
        }
    }
    for s in [1e-3, 1e3] {
        out.push(Rung {
            scale: s,
            ..Rung::FREE
        });
    }
    let distances = [0.0, 0.5 * tol, tol, 10.0 * tol, 1e-3];
    if touching(cell.config) || matches!(cell.config, Config::X | Config::X3) {
        for &d in &distances {
            if has_seam(cell.a) {
                out.push(Rung {
                    steer_a: Steer::Seam(d),
                    ..Rung::FREE
                });
            }
            if has_pole(cell.a) {
                out.push(Rung {
                    steer_a: Steer::Pole(d),
                    ..Rung::FREE
                });
            }
            if touching(cell.config) {
                if has_seam(cell.b) {
                    out.push(Rung {
                        steer_b: Steer::Seam(d),
                        ..Rung::FREE
                    });
                }
                if has_pole(cell.b) {
                    out.push(Rung {
                        steer_b: Steer::Pole(d),
                        ..Rung::FREE
                    });
                }
            }
        }
        // Both events at once, and an event at a near miss.
        if has_pole(cell.a) && has_seam(cell.b) && touching(cell.config) {
            out.push(Rung {
                steer_a: Steer::Pole(0.0),
                steer_b: Steer::Seam(0.0),
                ..Rung::FREE
            });
        }
        if has_seam(cell.a) && touching(cell.config) {
            out.push(Rung {
                offset: 0.5,
                steer_a: Steer::Seam(0.0),
                ..Rung::FREE
            });
        }
    }
    out
}

/// The smoke tier's ladder (design §7): the free scene, the half-tolerance
/// near miss, and the contact exactly on a seam and on a pole.
pub fn smoke_ladder(cell: Cell) -> Vec<Rung> {
    let mut out = vec![Rung::FREE];
    if touching(cell.config) {
        out.push(Rung {
            offset: 0.5,
            ..Rung::FREE
        });
        if has_seam(cell.a) {
            out.push(Rung {
                steer_a: Steer::Seam(0.0),
                ..Rung::FREE
            });
        }
        if has_pole(cell.a) {
            out.push(Rung {
                steer_a: Steer::Pole(0.0),
                ..Rung::FREE
            });
        }
    }
    out
}

/// A random scene's rung: offsets and chart events drawn too.
pub fn random_rung(cell: Cell, k: u64) -> Rung {
    let mut rng = super::rng::Rng::new(0x5EED, super::rng::stream_id(&cell.name()), k);
    let offset = if touching(cell.config) {
        OFFSETS[rng.below(OFFSETS.len())]
    } else {
        0.0
    };
    let scale = [1.0, 1.0, 1.0, 1e-3, 1e3][rng.below(5)];
    Rung {
        offset,
        scale,
        ..Rung::FREE
    }
}

/// One judged scene.
#[derive(Debug, Clone, Serialize)]
pub struct Judged {
    pub scene: Scene,
    pub verdict: Verdict,
}

/// A cell's run.
#[derive(Debug, Clone, Serialize, Default)]
pub struct CellRun {
    pub cell: Option<Cell>,
    pub ladder_run: usize,
    pub random_run: usize,
    pub passed: usize,
    pub rejected: usize,
    pub draws: usize,
    pub failed: Vec<Judged>,
    pub worst_rel: f64,
    pub reject_reasons: Vec<String>,
}

/// Draw a rung's scene, retrying rejected draws (counted) up to `tries`.
fn draw_retrying(
    cell: Cell,
    seed: u64,
    base: u64,
    rung: Rung,
    tries: u64,
    run: &mut CellRun,
) -> Option<Scene> {
    for t in 0..tries {
        run.draws += 1;
        let drawn = std::panic::catch_unwind(|| draw(cell, seed, base * 1009 + t, rung));
        let Ok(drawn) = drawn else {
            // A panic placing a scene is a generator or operand-mesh
            // defect: recorded with the rejects, and the rung abandoned.
            run.rejected += 1;
            run.reject_reasons.push("panic while drawing".to_string());
            return None;
        };
        match drawn {
            Ok(scene) => return Some(scene),
            Err(r) => {
                run.rejected += 1;
                if run.reject_reasons.len() < 8 {
                    run.reject_reasons.push(format!("{r:?}"));
                }
                if matches!(r, Reject::NoSample) && !matches!(rung.steer_a, Steer::Free) {
                    // The steer has no event on this operand's faces.
                    return None;
                }
            }
        }
    }
    None
}

/// Run a cell: its ladder (the smoke ladder when `smoke`) and `random`
/// random scenes.
pub fn run_cell(cell: Cell, seed: u64, random: usize, smoke: bool) -> CellRun {
    let mut run = CellRun {
        cell: Some(cell),
        ..CellRun::default()
    };
    let record = |scene: Scene, run: &mut CellRun, ladder: bool| {
        // The OCC oracle on a deterministic tenth of the full tier.
        let occ = !smoke
            && super::rng::stream_id(&format!("{}-{}", cell.name(), scene.index))
                .is_multiple_of(10);
        let verdict = judge_with(&scene, occ);
        if ladder {
            run.ladder_run += 1;
        } else {
            run.random_run += 1;
        }
        run.worst_rel = run.worst_rel.max(verdict.worst_rel);
        if verdict.failures.is_empty() {
            run.passed += 1;
        } else {
            run.failed.push(Judged { scene, verdict });
        }
    };
    let rungs = if smoke {
        smoke_ladder(cell)
    } else {
        ladder(cell)
    };
    for (i, rung) in rungs.into_iter().enumerate() {
        if let Some(scene) = draw_retrying(cell, seed, i as u64, rung, 24, &mut run) {
            record(scene, &mut run, true);
        }
    }
    for k in 0..random as u64 {
        let rung = random_rung(cell, k ^ seed);
        if let Some(scene) = draw_retrying(cell, seed, 10_000 + k, rung, 24, &mut run) {
            record(scene, &mut run, false);
        }
    }
    run
}

/// Whether a cell is feasible: some draw among `tries` is accepted.
pub fn feasible(cell: Cell, tries: u64) -> bool {
    (0..tries).any(|t| draw(cell, 0xF0, t, Rung::FREE).is_ok())
}
