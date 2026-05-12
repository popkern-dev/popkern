use super::*;
use crate::geometry::Point3;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn sample_sector(offset: f64) -> [[Point3; 4]; 3] {
    [
        [
            p(offset + 0.0, 0.0, 0.0),
            p(offset + 1.0, 0.0, 0.0),
            p(offset + 2.0, 0.0, 0.0),
            p(offset + 3.0, 0.0, 0.0),
        ],
        [
            p(offset + 0.0, 1.0, 0.0),
            p(offset + 1.0, 1.0, 0.0),
            p(offset + 2.0, 1.0, 0.0),
            p(offset + 3.0, 1.0, 0.0),
        ],
        [
            p(offset + 0.0, 2.0, 0.0),
            p(offset + 1.0, 2.0, 0.0),
            p(offset + 2.0, 2.0, 0.0),
            p(offset + 3.0, 2.0, 0.0),
        ],
    ]
}

fn assert_close(a: Point3, b: Point3) {
    assert!(
        a.approx_eq(&b, 1e-12),
        "expected {a:?} to match {b:?}, distance {}",
        a.distance_to(&b)
    );
}

#[test]
fn split_tensor_border_preserves_boundary_endpoints() {
    let border = TensorBorderBi3 {
        sectors: vec![sample_sector(0.0), sample_sector(10.0), sample_sector(20.0)],
    };

    let split = split_kp20_tensor_border(&border).unwrap();

    assert_eq!(split.len(), 3);
    assert_close(split[0].near_c.rows[0][0], p(0.0, 0.0, 0.0));
    assert_close(split[0].near_c.rows[0][6], p(1.5, 0.0, 0.0));
    assert_close(split[0].near_m.rows[0][0], p(1.5, 0.0, 0.0));
    assert_close(split[0].near_m.rows[0][6], p(3.0, 0.0, 0.0));
}

#[test]
fn split_tensor_border_degree_elevates_linear_rows() {
    let border = TensorBorderBi3 {
        sectors: vec![sample_sector(0.0), sample_sector(10.0), sample_sector(20.0)],
    };

    let split = split_kp20_tensor_border(&border).unwrap();

    for i in 0..7 {
        assert_close(split[0].near_c.rows[0][i], p(0.25 * i as f64, 0.0, 0.0));
        assert_close(
            split[0].near_m.rows[0][i],
            p(1.5 + 0.25 * i as f64, 0.0, 0.0),
        );
        assert_close(split[0].near_c.rows[1][i], p(0.25 * i as f64, 1.0, 0.0));
    }
}

#[test]
fn split_tensor_border_applies_eq14_second_cross_reparameterisation() {
    let border = TensorBorderBi3 {
        sectors: vec![sample_sector(0.0), sample_sector(10.0), sample_sector(20.0)],
    };

    let split = split_kp20_tensor_border(&border).unwrap();
    let c = (core::f64::consts::TAU / 3.0).cos();
    let d1_3 = 2.0 * c / 3.0;
    let delta = d1_3 / 4.0 - c / 12.0;

    assert_close(split[0].near_c.rows[2][0], p(0.0, 2.0, 0.0));
    assert_close(split[0].near_c.rows[2][6], p(1.5, 2.0 + 0.5 * delta, 0.0));
    assert_close(split[0].near_m.rows[2][0], p(1.5, 2.0 + 0.5 * delta, 0.0));
    assert_close(split[0].near_m.rows[2][6], p(3.0, 2.0 + 0.5 * d1_3, 0.0));
}

#[test]
fn initialise_macro_patches_assigns_m_side_boundary_rows() {
    let border = TensorBorderBi3 {
        sectors: vec![sample_sector(0.0), sample_sector(10.0), sample_sector(20.0)],
    };

    let macros = initialise_kp20_macro_patches_from_tensor_border(&border).unwrap();

    assert_eq!(macros.len(), 3);
    assert_close(
        macros[1].piece(Kp20PieceId::NearMidpointLeft).cps[0][0],
        p(10.0, 0.0, 0.0),
    );
    assert_close(
        macros[1].piece(Kp20PieceId::NearMidpointRight).cps[0][6],
        p(13.0, 0.0, 0.0),
    );
    assert_close(
        macros[1].piece(Kp20PieceId::NearOriginLeft).cps[0][0],
        p(0.0, 0.0, 0.0),
    );
}

#[test]
fn c2_upper_rows_follow_kp20_eq12() {
    let mut lower = Kp20Bi6Patch::zero();
    let mut upper = Kp20Bi6Patch::zero();
    for j in 0..7 {
        lower.cps[2][j] = p(2.0, j as f64, 0.0);
        lower.cps[4][j] = p(6.0, j as f64, 0.0);
        lower.cps[6][j] = p(10.0, j as f64, 0.0);
    }

    kp20_propagate_c2_upper_from_lower(&lower, &mut upper);

    for j in 0..7 {
        assert_close(upper.cps[0][j], p(10.0, j as f64, 0.0));
        assert_close(upper.cps[5][j], p(11.0, j as f64, 0.0));
        assert_close(upper.cps[1][j], p(9.0, j as f64, 0.0));
    }
}

#[test]
fn g2_strip_sets_rows_from_kp20_eq13() {
    let mut p_patch = Kp20Bi6Patch::zero();
    let mut q_patch = Kp20Bi6Patch::zero();
    p_patch.cps[4][0] = p(10.0, 0.0, 0.0);
    p_patch.cps[4][2] = p(14.0, 0.0, 0.0);
    q_patch.cps[4][2] = p(6.0, 0.0, 0.0);
    p_patch.cps[5][0] = p(2.0, 0.0, 0.0);
    p_patch.cps[6][0] = p(5.0, 0.0, 0.0);
    p_patch.cps[5][2] = p(11.0, 0.0, 0.0);
    q_patch.cps[5][2] = p(7.0, 0.0, 0.0);
    p_patch.cps[6][2] = p(13.0, 0.0, 0.0);
    q_patch.cps[6][2] = p(17.0, 0.0, 0.0);

    kp20_apply_g2_strip(&mut p_patch, &mut q_patch, 2.0);

    let (w1, w2, w3) = kp20_g2_strip_weights(2.0);
    let expected_correction = -w3 * (2.0 - 5.0) + w2 * (11.0 - 7.0) + w1 * (13.0 - 17.0);
    let expected_quarter = (14.0 - 6.0) / 4.0;
    assert_close(
        p_patch.cps[4][1],
        p(10.0 + expected_quarter + expected_correction, 0.0, 0.0),
    );
    assert_close(
        q_patch.cps[4][1],
        p(10.0 - expected_quarter + expected_correction, 0.0, 0.0),
    );
}
