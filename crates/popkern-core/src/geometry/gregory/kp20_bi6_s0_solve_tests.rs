use super::*;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn assert_close(a: Point3, b: Point3) {
    assert!(
        a.approx_eq(&b, 1e-12),
        "expected {a:?} to match {b:?}, distance {}",
        a.distance_to(&b)
    );
}

#[test]
fn exact_square_system_recovers_free_points() {
    let samples = vec![
        Kp20S0Sample::new(vec![(0, 1.0)], p(0.0, 0.0, 0.0), p(1.0, 2.0, 3.0)),
        Kp20S0Sample::new(vec![(1, 1.0)], p(0.0, 0.0, 0.0), p(4.0, 5.0, 6.0)),
    ];

    let out = solve_kp20_s0_projection(2, &samples).unwrap();

    assert_close(out.free_points[0], p(1.0, 2.0, 3.0));
    assert_close(out.free_points[1], p(4.0, 5.0, 6.0));
    assert!(out.residual_squared < 1e-24);
}

#[test]
fn overdetermined_system_projects_to_average() {
    let samples = vec![
        Kp20S0Sample::new(vec![(0, 1.0)], p(0.0, 0.0, 0.0), p(1.0, 2.0, 3.0)),
        Kp20S0Sample::new(vec![(0, 1.0)], p(0.0, 0.0, 0.0), p(3.0, 4.0, 5.0)),
    ];

    let out = solve_kp20_s0_projection(1, &samples).unwrap();

    assert_close(out.free_points[0], p(2.0, 3.0, 4.0));
    assert_close(out.projected_samples[0], p(2.0, 3.0, 4.0));
    assert_close(out.projected_samples[1], p(2.0, 3.0, 4.0));
    assert!((out.residual_squared - 6.0).abs() < 1e-12);
}

#[test]
fn constants_are_removed_before_solving_and_restored_after_projection() {
    let samples = vec![Kp20S0Sample::new(
        vec![(0, 2.0)],
        p(1.0, -2.0, 3.0),
        p(5.0, 4.0, 13.0),
    )];

    let out = solve_kp20_s0_projection(1, &samples).unwrap();

    assert_close(out.free_points[0], p(2.0, 3.0, 5.0));
    assert_close(out.projected_samples[0], p(5.0, 4.0, 13.0));
}

#[test]
fn coupled_system_solves_independently_per_axis() {
    let samples = vec![
        Kp20S0Sample::new(
            vec![(0, 1.0), (1, 1.0)],
            p(0.0, 0.0, 0.0),
            p(3.0, 7.0, 11.0),
        ),
        Kp20S0Sample::new(
            vec![(0, 1.0), (1, -1.0)],
            p(0.0, 0.0, 0.0),
            p(-1.0, -1.0, -1.0),
        ),
    ];

    let out = solve_kp20_s0_projection(2, &samples).unwrap();

    assert_close(out.free_points[0], p(1.0, 3.0, 5.0));
    assert_close(out.free_points[1], p(2.0, 4.0, 6.0));
}

#[test]
fn out_of_range_free_index_is_reported() {
    let samples = vec![Kp20S0Sample::new(
        vec![(2, 1.0)],
        p(0.0, 0.0, 0.0),
        p(1.0, 0.0, 0.0),
    )];

    let err = solve_kp20_s0_projection(1, &samples).unwrap_err();

    assert_eq!(
        err,
        Kp20S0SolveError::FreeIndexOutOfRange {
            sample: 0,
            index: 2,
            num_free: 1
        }
    );
}

#[test]
fn singular_normal_matrix_is_reported() {
    let samples = vec![Kp20S0Sample::new(
        vec![(0, 1.0)],
        p(0.0, 0.0, 0.0),
        p(1.0, 0.0, 0.0),
    )];

    let err = solve_kp20_s0_projection(2, &samples).unwrap_err();

    assert_eq!(err, Kp20S0SolveError::SingularOnAxis(0));
}
