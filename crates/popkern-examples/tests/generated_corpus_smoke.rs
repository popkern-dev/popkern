//! The generated corpus's per-commit tier (PK-1227, design §7): one
//! representative cell per configuration, its free scene and a contact
//! exactly on a seam or a pole, judged by every oracle under the
//! known-failure ratchet. The full tier runs every cell
//! (`--profile foundations`).

#[path = "generated_corpus/mod.rs"]
mod generated_corpus;

use generated_corpus::geometry::Steer;
use generated_corpus::operand::Kind::*;
use generated_corpus::report::run_listed;
use generated_corpus::scene::{Cell, Config::*, Rung};

fn free() -> Rung {
    Rung::FREE
}

fn seam_a() -> Rung {
    Rung {
        steer_a: Steer::Seam(0.0),
        ..Rung::FREE
    }
}

fn pole_a() -> Rung {
    Rung {
        steer_a: Steer::Pole(0.0),
        ..Rung::FREE
    }
}

fn check(list: &[(Cell, Rung)]) {
    let unexpected = run_listed(list);
    assert!(
        unexpected.is_empty(),
        "{} problems:\n{}",
        unexpected.len(),
        unexpected.join("\n")
    );
}

#[test]
fn touches() {
    check(&[
        (
            Cell {
                a: S,
                b: Y,
                config: E,
            },
            free(),
        ),
        (
            Cell {
                a: S,
                b: Y,
                config: E,
            },
            pole_a(),
        ),
        (
            Cell {
                a: P,
                b: S,
                config: I,
            },
            free(),
        ),
        (
            Cell {
                a: Y,
                b: S,
                config: K,
            },
            free(),
        ),
        (
            Cell {
                a: S,
                b: S,
                config: E,
            },
            Rung {
                offset: 0.5,
                ..Rung::FREE
            },
        ),
    ]);
}

#[test]
fn crossings() {
    check(&[
        (
            Cell {
                a: Y,
                b: S,
                config: X,
            },
            seam_a(),
        ),
        (
            Cell {
                a: S,
                b: Y,
                config: X3,
            },
            free(),
        ),
        (
            Cell {
                a: Y,
                b: Y,
                config: V,
            },
            free(),
        ),
    ]);
}

#[test]
fn apart_contained_coincident() {
    check(&[
        (
            Cell {
                a: S,
                b: P,
                config: D,
            },
            free(),
        ),
        (
            Cell {
                a: P,
                b: B,
                config: N,
            },
            free(),
        ),
        (
            Cell {
                a: S,
                b: S,
                config: O,
            },
            free(),
        ),
        (
            Cell {
                a: Y,
                b: Y,
                config: O,
            },
            free(),
        ),
    ]);
}
