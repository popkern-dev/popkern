use super::*;
use crate::geometry::Point3;
use crate::geometry::gregory::{Kp20Bi6CapInput, assemble_kp20_bi6_cap};

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

fn tensor_border() -> TensorBorderBi3 {
    TensorBorderBi3 {
        sectors: vec![sample_sector(0.0), sample_sector(10.0), sample_sector(20.0)],
    }
}

#[test]
fn validation_report_accepts_assembled_planar_cap() {
    let tensor_border = tensor_border();
    let output = assemble_kp20_bi6_cap(&Kp20Bi6CapInput {
        star_point: p(0.0, 0.0, 0.0),
        tensor_border: tensor_border.clone(),
    })
    .unwrap();

    let report = measure_kp20_bi6_macro_cap(&output, &tensor_border).unwrap();

    assert!(report.passes(1e-10), "{report:?}");
    assert!(report.max_abs_z <= 1e-12, "{report:?}");
}

#[test]
fn validation_report_detects_c2_drift() {
    let tensor_border = tensor_border();
    let mut output = assemble_kp20_bi6_cap(&Kp20Bi6CapInput {
        star_point: p(0.0, 0.0, 0.0),
        tensor_border: tensor_border.clone(),
    })
    .unwrap();
    output.macro_patches[0]
        .piece_mut(Kp20PieceId::NearOriginLeft)
        .cps[1][3]
        .x += 0.125;

    let report = measure_kp20_bi6_macro_cap(&output, &tensor_border).unwrap();

    assert!(report.max_internal_c2_first_residual > 0.1, "{report:?}");
    assert!(!report.passes(1e-10));
}

#[test]
fn validation_report_rejects_valence_mismatch() {
    let tensor_border = tensor_border();
    let output = assemble_kp20_bi6_cap(&Kp20Bi6CapInput {
        star_point: p(0.0, 0.0, 0.0),
        tensor_border: tensor_border.clone(),
    })
    .unwrap();
    let shorter = TensorBorderBi3 {
        sectors: tensor_border.sectors[0..2].to_vec(),
    };

    let err = measure_kp20_bi6_macro_cap(&output, &shorter).unwrap_err();

    assert!(matches!(
        err,
        Kp20Bi6ValidationError::ValenceMismatch {
            output: 3,
            tensor_border: 2
        }
    ));
}
