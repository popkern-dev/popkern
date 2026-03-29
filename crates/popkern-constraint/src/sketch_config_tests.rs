//! Tests for named configuration management.

use crate::constraint::Constraint;
use crate::Sketch;

/// Save a configuration with distance=5, change to distance=10 in a second
/// config, then verify activating the first restores distance=5.
#[test]
fn save_and_activate_restores_distance_value() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(5.0, 0.0);
    sketch.add_constraint(Constraint::Fixed(a));
    let dist = sketch.add_constraint(Constraint::Distance(a, b, 5.0));
    sketch.solve().unwrap();

    // Save "config_5".
    let cfg5 = sketch.create_config("config_5");
    sketch.save_config(cfg5).unwrap();

    // Change distance to 10.
    sketch.set_constraint_value(dist, 10.0).unwrap();
    sketch.solve().unwrap();
    let (bx, _) = sketch.point_position(b).unwrap();
    assert!((bx - 10.0).abs() < 1e-6, "expected b.x=10, got {bx}");

    // Save "config_10".
    let cfg10 = sketch.create_config("config_10");
    sketch.save_config(cfg10).unwrap();

    // Activate config_5 — should restore distance=5 and solve.
    sketch.activate_config(cfg5).unwrap();
    let (bx5, _) = sketch.point_position(b).unwrap();
    assert!(
        (bx5 - 5.0).abs() < 1e-6,
        "after activating cfg5, b.x should be 5.0, got {bx5}"
    );

    // Activate config_10 — should restore distance=10.
    sketch.activate_config(cfg10).unwrap();
    let (bx10, _) = sketch.point_position(b).unwrap();
    assert!(
        (bx10 - 10.0).abs() < 1e-6,
        "after activating cfg10, b.x should be 10.0, got {bx10}"
    );
}

/// Activating a config that references a deleted constraint skips it without
/// panicking.
#[test]
fn activate_config_skips_missing_constraint_handles() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(3.0, 0.0);
    sketch.add_constraint(Constraint::Fixed(a));
    let dist = sketch.add_constraint(Constraint::Distance(a, b, 3.0));
    sketch.solve().unwrap();

    let cfg = sketch.create_config("cfg");
    sketch.save_config(cfg).unwrap();

    // Delete the distance constraint.
    sketch.remove_constraint(dist);

    // Activating should not panic — the stale handle is silently skipped.
    sketch.activate_config(cfg).unwrap();
}

/// Group enabled states are saved and restored by configurations.
#[test]
fn config_saves_and_restores_group_states() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(3.0, 0.0);
    sketch.add_constraint(Constraint::Fixed(a));
    let dist = sketch.add_constraint(Constraint::Distance(a, b, 10.0));
    let group = sketch.create_group("grp");
    sketch.add_to_group(group, dist).unwrap();
    sketch.disable_group(group);
    sketch.solve().unwrap();

    // Save config with group disabled.
    let cfg_disabled = sketch.create_config("disabled");
    sketch.save_config(cfg_disabled).unwrap();

    // Enable the group, save another config.
    sketch.enable_group(group);
    let cfg_enabled = sketch.create_config("enabled");
    sketch.save_config(cfg_enabled).unwrap();

    // Activate the disabled config — group should be disabled again.
    sketch.activate_config(cfg_disabled).unwrap();
    assert_eq!(sketch.is_group_enabled(group), Some(false));
    // Constraint was inactive; b.x should be near 3.0 still.
    let (bx, _) = sketch.point_position(b).unwrap();
    assert!(
        (bx - 3.0).abs() < 0.5,
        "b.x should be near 3 with group disabled, got {bx}"
    );

    // Activate the enabled config — group should be enabled, solver pulls b to 10.
    sketch.activate_config(cfg_enabled).unwrap();
    assert_eq!(sketch.is_group_enabled(group), Some(true));
    let (bx2, _) = sketch.point_position(b).unwrap();
    assert!(
        (bx2 - 10.0).abs() < 1e-6,
        "b.x should be 10 with group enabled, got {bx2}"
    );
}

/// list_configs returns all created configs.
#[test]
fn list_configs_returns_all_configs() {
    let mut sketch = Sketch::new();
    let c1 = sketch.create_config("alpha");
    let c2 = sketch.create_config("beta");

    let list = sketch.list_configs();
    assert_eq!(list.len(), 2);
    let names: Vec<&str> = list.iter().map(|(_, n)| *n).collect();
    assert!(names.contains(&"alpha"), "missing alpha");
    assert!(names.contains(&"beta"), "missing beta");

    sketch.delete_config(c1);
    let list2 = sketch.list_configs();
    assert_eq!(list2.len(), 1);
    assert_eq!(list2[0].0, c2);
    assert_eq!(list2[0].1, "beta");
}

/// Warm-start: activating a saved config with matching param length reuses
/// the saved params vector. The solve should converge in very few iterations
/// since the warm start is already near the solution.
#[test]
fn activate_config_warm_start_converges_quickly() {
    use crate::solver::SolverOptions;

    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(5.0, 0.0);
    sketch.add_constraint(Constraint::Fixed(a));
    sketch.add_constraint(Constraint::Distance(a, b, 5.0));
    sketch.solve().unwrap();

    let cfg = sketch.create_config("warm");
    sketch.save_config(cfg).unwrap();

    // Perturb b slightly.
    let idx = sketch.points[b].param_indices()[0];
    sketch.params[idx] = 4.9;

    // Activate — warm-start from near-solution.
    let opts = SolverOptions {
        max_iterations: 50,
        ..Default::default()
    };
    let result = sketch
        .solve_with_options(&opts)
        .expect("solve should succeed");

    // Warm start means very few iterations needed.
    assert!(
        result.iterations <= 5,
        "expected ≤5 iterations with warm start, got {}",
        result.iterations
    );
}
