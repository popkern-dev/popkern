use super::*;
use crate::geometry::Point3;
use crate::geometry::gregory::{Kp20RegularCNet, TensorBorderBi3, kp20_apply_g2_strip};

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn assert_close(a: Point3, b: Point3) {
    assert!(
        a.approx_eq(&b, 1e-9),
        "expected {a:?} to match {b:?}, distance {}",
        a.distance_to(&b)
    );
}

fn affine_point(alpha: f64, a: Point3, beta: f64, b: Point3) -> Point3 {
    p(
        alpha * a.x + beta * b.x,
        alpha * a.y + beta * b.y,
        alpha * a.z + beta * b.z,
    )
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

fn tensor_border_with_valence(n: usize) -> TensorBorderBi3 {
    TensorBorderBi3 {
        sectors: (0..n)
            .map(|sector| sample_sector(10.0 * sector as f64))
            .collect(),
    }
}

fn curved_sector(offset: f64, twist: f64) -> [[Point3; 4]; 3] {
    let mut sector = [[p(0.0, 0.0, 0.0); 4]; 3];
    for (row_idx, row) in sector.iter_mut().enumerate() {
        for (col_idx, cp) in row.iter_mut().enumerate() {
            let u = col_idx as f64;
            let v = row_idx as f64;
            *cp = p(
                offset + u + 0.13 * v,
                0.7 * v + twist * u * v,
                twist * (u + 1.0) * (v + 0.5),
            );
        }
    }
    sector
}

fn curved_tensor_border(scale: f64) -> TensorBorderBi3 {
    TensorBorderBi3 {
        sectors: vec![
            curved_sector(scale, 0.03),
            curved_sector(10.0 + scale, -0.02),
            curved_sector(20.0 + scale, 0.05),
        ],
    }
}

fn regular_cnet() -> Kp20RegularCNet {
    let mut sectors = Vec::new();
    for sector in 0..3 {
        let mut grid = [[p(0.0, 0.0, 0.0); 4]; 4];
        for (row, cps) in grid.iter_mut().enumerate() {
            for (col, cp) in cps.iter_mut().enumerate() {
                *cp = p(
                    10.0 * sector as f64 + row as f64,
                    col as f64,
                    0.05 * (row * col) as f64,
                );
            }
        }
        sectors.push(grid);
    }
    Kp20RegularCNet { sectors }
}

fn combine_tensor_borders(
    alpha: f64,
    a: &TensorBorderBi3,
    beta: f64,
    b: &TensorBorderBi3,
) -> TensorBorderBi3 {
    let sectors = a
        .sectors
        .iter()
        .zip(&b.sectors)
        .map(|(a_sector, b_sector)| {
            let mut sector = [[p(0.0, 0.0, 0.0); 4]; 3];
            for row in 0..3 {
                for col in 0..4 {
                    sector[row][col] =
                        affine_point(alpha, a_sector[row][col], beta, b_sector[row][col]);
                }
            }
            sector
        })
        .collect();
    TensorBorderBi3 { sectors }
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
fn assembler_accepts_appendix_a_fig6c_tau_as_paper_reparameterization() {
    let input = Kp20Bi6CapInput {
        star_point: p(0.0, 0.0, 0.0),
        tensor_border: tensor_border_with_valence(5),
    };

    let output = assemble_kp20_bi6_cap_with_appendix_a_fig6c_tau(&input)
        .expect("paper-valid Fig.6c tau assembly should not require triangle containment");

    assert_eq!(output.sector_samples.len(), 5);
    assert_eq!(output.cap.macro_patches.len(), 5);
    for sample in &output.sector_samples {
        assert_eq!(sample.tau.orientation.negative_determinants, 0);
        assert_eq!(sample.tau.orientation.near_zero_determinants, 0);
        assert!(sample.tau.square_domain.stays_inside_square());
    }
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
    // KP20 equation (14) uses the two adjacent sector-local near-c
    // tensor-border halves for the sector-separating near-m strip.
    assert_close(
        sector.piece(Kp20PieceId::NearMidpointLeft).cps[0][0],
        p(20.0, 0.0, 0.0),
    );
    assert_close(
        sector.piece(Kp20PieceId::NearMidpointRight).cps[0][6],
        p(1.5, 0.0, 0.0),
    );
    // After paper §4.4 step (i) overlay and the polynomial-
    // consistency restoration that follows it, NearOrigin no
    // longer carries C²-propagated values from NearMidpoint.
    // The prior assertion (`cps[0][3] == NearMidpoint.cps[6][3]`)
    // documented the kernel's pre-paper-fix behavior; it is no
    // longer a meaningful invariant after the paper §4.4 step (i)
    // wire-up. The other invariants in this test (n, lengths,
    // tensor-border endpoints) remain valid.
}

#[test]
fn assembler_can_use_sampled_tensor_border_source() {
    let sampled_macro_patches = sampled_macro_patches(3);
    let input = Kp20Bi6MacroCapInput {
        tensor_border: tensor_border(),
        sampled_macro_patches: sampled_macro_patches.clone(),
    };

    let out = assemble_kp20_bi6_macro_cap_with_sampled_tensor_border_source(&input).unwrap();

    let sector = &out.macro_patches[0];
    assert_close(
        sector.piece(Kp20PieceId::NearMidpointLeft).cps[0][0],
        sampled_macro_patches[0]
            .piece(Kp20PieceId::NearMidpointLeft)
            .cps[0][0],
    );
    assert_close(
        sector.piece(Kp20PieceId::NearMidpointRight).cps[2][6],
        sampled_macro_patches[0]
            .piece(Kp20PieceId::NearMidpointRight)
            .cps[2][6],
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

#[test]
fn cap_entry_point_builds_guides_and_assembles_macro_cap() {
    let input = Kp20Bi6CapInput {
        star_point: p(0.0, 0.0, 0.0),
        tensor_border: tensor_border(),
    };

    let out = assemble_kp20_bi6_cap(&input).unwrap();

    assert_eq!(out.n, 3);
    assert_eq!(out.macro_patches.len(), 3);
    assert_eq!(out.patches_bi6.len(), 12);
    assert_eq!(out.strip_projections.len(), 3);
    for patch in out.patches_bi6 {
        for row in patch {
            for cp in row {
                assert!(cp.z.abs() < 1e-10);
            }
        }
    }
}

#[test]
fn cap_assembly_is_affine_linear_in_star_and_tensor_border() {
    let tensor_a = curved_tensor_border(0.25);
    let tensor_b = curved_tensor_border(-0.75);
    let star_a = p(0.3, -0.2, 0.4);
    let star_b = p(-0.1, 0.5, -0.3);
    let alpha = 1.7;
    let beta = -0.4;

    let out_a = assemble_kp20_bi6_cap(&Kp20Bi6CapInput {
        star_point: star_a,
        tensor_border: tensor_a.clone(),
    })
    .unwrap();
    let out_b = assemble_kp20_bi6_cap(&Kp20Bi6CapInput {
        star_point: star_b,
        tensor_border: tensor_b.clone(),
    })
    .unwrap();
    let out_combined = assemble_kp20_bi6_cap(&Kp20Bi6CapInput {
        star_point: affine_point(alpha, star_a, beta, star_b),
        tensor_border: combine_tensor_borders(alpha, &tensor_a, beta, &tensor_b),
    })
    .unwrap();

    assert_eq!(out_combined.patches_bi6.len(), out_a.patches_bi6.len());
    for patch_idx in 0..out_combined.patches_bi6.len() {
        for i in 0..7 {
            for j in 0..7 {
                let expected = affine_point(
                    alpha,
                    out_a.patches_bi6[patch_idx][i][j],
                    beta,
                    out_b.patches_bi6[patch_idx][i][j],
                );
                assert_close(out_combined.patches_bi6[patch_idx][i][j], expected);
            }
        }
    }
}

#[test]
fn cap_entry_point_from_cnet_matches_explicit_tensor_border_path() {
    let cnet = regular_cnet();
    let tensor_border = kp20_tensor_border_from_regular_cnet(&cnet).unwrap();
    let star_point = p(0.1, -0.2, 0.3);

    let from_cnet =
        assemble_kp20_bi6_cap_from_cnet(&Kp20Bi6CNetCapInput { star_point, cnet }).unwrap();
    let from_tensor = assemble_kp20_bi6_cap(&Kp20Bi6CapInput {
        star_point,
        tensor_border,
    })
    .unwrap();

    assert_eq!(from_cnet.patches_bi6.len(), from_tensor.patches_bi6.len());
    for patch_idx in 0..from_cnet.patches_bi6.len() {
        for i in 0..7 {
            for j in 0..7 {
                assert_close(
                    from_cnet.patches_bi6[patch_idx][i][j],
                    from_tensor.patches_bi6[patch_idx][i][j],
                );
            }
        }
    }
}
