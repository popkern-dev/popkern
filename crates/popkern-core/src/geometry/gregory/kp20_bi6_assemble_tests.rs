use super::*;
use crate::geometry::Point3;
use crate::geometry::gregory::{TensorBorderBi3, kp20_apply_g2_strip};

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn assert_close(a: Point3, b: Point3) {
    assert!(
        a.approx_eq(&b, 1e-10),
        "expected {a:?} to match {b:?}, distance {}",
        a.distance_to(&b)
    );
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

fn tensor_border() -> TensorBorderBi3 {
    TensorBorderBi3 {
        sectors: vec![sample_sector(0.0), sample_sector(10.0), sample_sector(20.0)],
    }
}

fn sampled_macro_patches(n: usize) -> Vec<Kp20MacroPatch> {
    let mut out = vec![Kp20MacroPatch::zero(); n];
    for (k, macro_patch) in out.iter_mut().enumerate() {
        for id in Kp20PieceId::all() {
            let piece = macro_patch.piece_mut(id);
            for i in 0..7 {
                for j in 0..7 {
                    piece.cps[i][j] = p(
                        1000.0 + 100.0 * k as f64 + 10.0 * id.index() as f64 + i as f64,
                        j as f64,
                        0.0,
                    );
                }
            }
        }
    }
    out
}

#[test]
fn assembler_rejects_valence_mismatch() {
    let input = Kp20Bi6MacroCapInput {
        tensor_border: tensor_border(),
        sampled_macro_patches: sampled_macro_patches(2),
    };

    let err = assemble_kp20_bi6_macro_cap(&input).unwrap_err();

    assert!(matches!(
        err,
        Kp20Bi6AssembleError::ValenceMismatch {
            tensor_border: 3,
            sampled_macro_patches: 2
        }
    ));
}

#[test]
fn assembler_applies_tensor_border_and_internal_c2_rows() {
    let input = Kp20Bi6MacroCapInput {
        tensor_border: tensor_border(),
        sampled_macro_patches: sampled_macro_patches(3),
    };

    let out = assemble_kp20_bi6_macro_cap(&input).unwrap();

    assert_eq!(out.n, 3);
    assert_eq!(out.macro_patches.len(), 3);
    assert_eq!(out.patches_bi6.len(), 12);
    let sector = &out.macro_patches[0];
    assert_close(
        sector.piece(Kp20PieceId::NearMidpointLeft).cps[0][0],
        p(0.0, 0.0, 0.0),
    );
    assert_close(
        sector.piece(Kp20PieceId::NearMidpointRight).cps[0][6],
        p(3.0, 0.0, 0.0),
    );
    assert_close(
        sector.piece(Kp20PieceId::NearOriginLeft).cps[0][3],
        sector.piece(Kp20PieceId::NearMidpointLeft).cps[6][3],
    );
    assert_close(
        sector.piece(Kp20PieceId::NearOriginRight).cps[0][3],
        sector.piece(Kp20PieceId::NearMidpointRight).cps[6][3],
    );
}

#[test]
fn assembler_projects_origin_strip_to_eq13_constraints() {
    let input = Kp20Bi6MacroCapInput {
        tensor_border: tensor_border(),
        sampled_macro_patches: sampled_macro_patches(3),
    };

    let out = assemble_kp20_bi6_macro_cap(&input).unwrap();

    assert_eq!(out.strip_projections.len(), 3);
    let sector = &out.macro_patches[1];
    let left = sector.piece(Kp20PieceId::NearOriginLeft);
    let right = sector.piece(Kp20PieceId::NearOriginRight);

    for i in 5..=6 {
        let expected_left = left.cps[i][0] + (left.cps[i][2] - right.cps[i][2]) * 0.25;
        let expected_right = left.cps[i][0] - (left.cps[i][2] - right.cps[i][2]) * 0.25;
        assert_close(left.cps[i][1], expected_left);
        assert_close(right.cps[i][1], expected_right);
    }

    let mut expected_left = *left;
    let mut expected_right = *right;
    kp20_apply_g2_strip(&mut expected_left, &mut expected_right, -1.0);
    assert_close(left.cps[4][1], expected_left.cps[4][1]);
    assert_close(right.cps[4][1], expected_right.cps[4][1]);
}

#[test]
fn assembler_preserves_planarity_for_planar_inputs() {
    let input = Kp20Bi6MacroCapInput {
        tensor_border: tensor_border(),
        sampled_macro_patches: sampled_macro_patches(3),
    };

    let out = assemble_kp20_bi6_macro_cap(&input).unwrap();

    for patch in out.patches_bi6 {
        for row in patch {
            for cp in row {
                assert!(cp.z.abs() < 1e-12);
            }
        }
    }
}
