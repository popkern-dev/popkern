//! Per-fixture shards of the validation discovery sweep — split
//! from `validation_sweep.rs` for the 500-line structural cap
//! (PK-268: the all-fixtures monolith was a 131–169 s single test;
//! sharding lets the suite schedule parallelize the sweep, and each
//! fixture leg keeps the monolith's collect-all + invariant +
//! improvement assertions).

use super::*;

/// One fixture's leg of the discovery sweep: all 5 variant rungs,
/// rung results printed to stdout, invariants collected and asserted
/// per fixture — the same collect-all semantics the single
/// all-fixtures test had, sharded per fixture so the suite schedule
/// can parallelize the sweep (PK-268: the monolith was a 131–169 s
/// single test; each fixture leg is well under the ~120 s bar).
fn run_sweep_fixture(name: &str, fixture_fn: fn() -> VertexEdges, radii: &[f64]) {
    // cylinder_cap fixture deferred — the cylinder primitive has no
    // valence-3 vertex; a curved-support corner needs its own design.
    let _ = fixture_cylinder_cap_corner;
    let rungs = run_fixture_all_rungs(name, fixture_fn, radii);
    assert_rung_invariants(name, &rungs);
    assert_end_to_end_improvement(name, &rungs);
    print_summary(&[(name.to_string(), rungs)]);
}

#[test]
fn validation_sweep_cube_mixed() {
    run_sweep_fixture("cube_mixed", fixture_cube, &[0.08, 0.10, 0.12]);
}

#[test]
fn validation_sweep_cube_uniform() {
    run_sweep_fixture("cube_uniform", fixture_cube, &[0.10, 0.10, 0.10]);
}

#[test]
fn validation_sweep_cube_extreme() {
    run_sweep_fixture("cube_extreme", fixture_cube, &[0.02, 0.10, 0.20]);
}

#[test]
fn validation_sweep_pyramid4_uniform() {
    run_sweep_fixture("pyramid4_uniform", fixture_square_pyramid, &[0.05; 4]);
}

#[test]
fn validation_sweep_pyramid4_mixed() {
    run_sweep_fixture(
        "pyramid4_mixed",
        fixture_square_pyramid,
        &[0.04, 0.06, 0.05, 0.07],
    );
}

#[test]
fn validation_sweep_pyramid5_uniform() {
    run_sweep_fixture(
        "pyramid5_uniform",
        || fixture_n_gonal_pyramid(5),
        &[0.04; 5],
    );
}

#[test]
fn validation_sweep_pyramid5_mixed() {
    run_sweep_fixture(
        "pyramid5_mixed",
        || fixture_n_gonal_pyramid(5),
        &[0.03, 0.04, 0.045, 0.05, 0.055],
    );
}

#[test]
fn validation_sweep_pyramid6_uniform() {
    run_sweep_fixture(
        "pyramid6_uniform",
        || fixture_n_gonal_pyramid(6),
        &[0.04; 6],
    );
}

#[test]
fn validation_sweep_pyramid6_mixed() {
    run_sweep_fixture(
        "pyramid6_mixed",
        || fixture_n_gonal_pyramid(6),
        &[0.025, 0.035, 0.04, 0.045, 0.05, 0.055],
    );
}
