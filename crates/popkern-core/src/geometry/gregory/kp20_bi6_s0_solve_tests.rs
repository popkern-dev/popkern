use super::*;
use crate::geometry::gregory::{Kp20Bi6Patch, kp20_g2_strip_weights};

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
    assert!(out.normal_condition_number.is_finite());
    assert!(out.normal_condition_number >= 1.0);
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

#[test]
fn g2_strip_samples_encode_kp20_eq13_linear_rows() {
    let mut p_patch = Kp20Bi6Patch::zero();
    let mut q_patch = Kp20Bi6Patch::zero();
    p_patch.cps[4][0] = p(10.0, 0.0, 0.0);
    p_patch.cps[4][1] = p(11.0, 0.0, 0.0);
    q_patch.cps[4][1] = p(9.0, 0.0, 0.0);
    p_patch.cps[4][2] = p(14.0, 0.0, 0.0);
    q_patch.cps[4][2] = p(6.0, 0.0, 0.0);
    p_patch.cps[5][0] = p(2.0, 0.0, 0.0);
    p_patch.cps[6][0] = p(5.0, 0.0, 0.0);
    p_patch.cps[5][2] = p(11.0, 0.0, 0.0);
    q_patch.cps[5][2] = p(7.0, 0.0, 0.0);
    p_patch.cps[6][2] = p(13.0, 0.0, 0.0);
    q_patch.cps[6][2] = p(17.0, 0.0, 0.0);

    let samples = build_kp20_g2_strip_samples(&p_patch, &q_patch, 2.0);

    assert_eq!(samples.len(), 12);
    let (w1, w2, w3) = kp20_g2_strip_weights(2.0);
    let expected_const = 10.0 - w3 * 2.0 + w3 * 5.0;
    assert!((samples[10].constant.x - expected_const).abs() < 1e-12);
    assert_eq!(
        samples[10].coeffs,
        vec![(0, 0.25), (3, -0.25), (1, w2), (4, -w2), (2, w1), (5, -w1)]
    );
    assert_eq!(
        samples[11].coeffs,
        vec![(0, -0.25), (3, 0.25), (1, -w2), (4, w2), (2, -w1), (5, w1)]
    );
}

#[test]
fn project_g2_strip_pair_writes_eq13_consistent_rows() {
    let mut p_patch = Kp20Bi6Patch::zero();
    let mut q_patch = Kp20Bi6Patch::zero();
    for i in 4..=6 {
        p_patch.cps[i][0] = p(i as f64, 0.0, 0.0);
        q_patch.cps[i][0] = p(100.0 + i as f64, 0.0, 0.0);
        p_patch.cps[i][1] = p(30.0 + i as f64, 0.0, 0.0);
        q_patch.cps[i][1] = p(40.0 + i as f64, 0.0, 0.0);
        p_patch.cps[i][2] = p(10.0 + i as f64, 0.0, 0.0);
        q_patch.cps[i][2] = p(20.0 + i as f64, 0.0, 0.0);
    }

    let out = project_kp20_g2_strip_pair(&mut p_patch, &mut q_patch, 2.0).unwrap();

    assert_eq!(out.free_points.len(), 6);
    let expected_p51 = p_patch.cps[5][0].x + (p_patch.cps[5][2].x - q_patch.cps[5][2].x) / 4.0;
    let expected_q51 = p_patch.cps[5][0].x - (p_patch.cps[5][2].x - q_patch.cps[5][2].x) / 4.0;
    assert!((p_patch.cps[5][1].x - expected_p51).abs() < 1e-10);
    assert!((q_patch.cps[5][1].x - expected_q51).abs() < 1e-10);
}
