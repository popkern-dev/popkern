use super::*;
use crate::geometry::gregory::cc_characteristic_tensor_border::cc_characteristic_tensor_border;

fn characteristic_border(n: usize) -> Kp20PaperGuideBorder {
    let chi = cc_characteristic_tensor_border(n).expect("chi");
    Kp20PaperGuideBorder::from_characteristic(&chi)
}

/// A genuinely valid curved fixture: the tensor-border of a REAL
/// Catmull-Clark ring on a node-wise paraboloid-lifted characteristic
/// net (exactly C²-compatible across sectors by construction; the
/// ring's height is not a polynomial function of the plane, so the
/// guide fit has genuine misfit). Chord-type fixtures with corners at
/// the spoke points are NOT valid KP20 input — the paper's border is
/// Hermite data of a C² ring — and were measured to break the G rows
/// structurally.
fn production_border(n: usize) -> Kp20PaperGuideBorder {
    let pieces =
        crate::geometry::gregory::cc_characteristic_tensor_border::cc_lifted_ring_border_pieces(
            n,
            |x, y| 0.4 * (x * x + y * y),
        )
        .expect("lifted ring");
    Kp20PaperGuideBorder {
        sectors: pieces
            .into_iter()
            .map(|[near_c, near_m]| {
                [
                    crate::geometry::gregory::kp20_paper_guide::Kp20PaperGuideBorderPiece {
                        rows: near_c,
                    },
                    crate::geometry::gregory::kp20_paper_guide::Kp20PaperGuideBorderPiece {
                        rows: near_m,
                    },
                ]
            })
            .collect(),
    }
}

fn report_line(label: &str, report: &Kp20PaperCapReport) {
    eprintln!(
        "{label} n={n}: c0={c0:e} g1={g1:e} g2={g2:e} border={b:e} lemma4={l4:e} strip_move={sm:e} guide_fit={gf:e} min_reg={mr:e}",
        n = report.n,
        c0 = report.max_c0_residual,
        g1 = report.max_g1_residual,
        g2 = report.max_g2_residual,
        b = report.max_border_residual,
        l4 = report.lemma4_delta,
        sm = report.max_strip_movement,
        gf = report.guide_fit_residual,
        mr = report.min_regularity,
    );
}

/// Characteristic fixture: the paper's own validation universe. With
/// the identity guide and the WP2 τ, steps (ii)/(iv) are near no-ops
/// and every cap-level gate must close at solver tolerance.
#[test]
fn characteristic_cap_closes_all_gates() {
    for n in [3, 5, 6, 8] {
        let report = assemble_kp20_paper_cap(&characteristic_border(n), Point3::new(0.0, 0.0, 0.0))
            .expect("assemble");
        report_line("characteristic", &report);
        assert!(report.max_c0_residual < 1.0e-8, "n={n} c0");
        assert!(report.max_g1_residual < 1.0e-7, "n={n} g1");
        assert!(report.max_g2_residual < 1.0e-7, "n={n} g2");
        assert!(report.max_border_residual < 1.0e-9, "n={n} border");
        assert!(report.lemma4_delta < 1.0e-8, "n={n} lemma4");
        assert!(report.min_regularity > 1.0e-3, "n={n} regularity");
    }
}

/// Production fixture: a genuinely 3D non-characteristic border. The
/// construction must still deliver G² at solver tolerance — that is
/// the paper's whole claim — regardless of the guide's ~0.14 misfit.
#[test]
fn production_cap_closes_g2_gates() {
    for n in [5, 6] {
        let report = assemble_kp20_paper_cap(&production_border(n), Point3::new(0.0, 0.0, 0.0))
            .expect("assemble");
        report_line("production", &report);
        assert!(report.max_c0_residual < 1.0e-8, "n={n} c0");
        assert!(
            report.max_g1_residual < 1.0e-7,
            "n={n} g1 {:e}",
            report.max_g1_residual,
        );
        assert!(
            report.max_g2_residual < 1.0e-7,
            "n={n} g2 {:e}",
            report.max_g2_residual,
        );
        assert!(report.max_border_residual < 1.0e-9, "n={n} border");
        assert!(report.lemma4_delta < 1.0e-7, "n={n} lemma4");
        assert!(report.min_regularity > 1.0e-3, "n={n} regularity");
    }
}

/// Diagnostic: per-BB-row eq (8)/(9) residual profile along curve 0
/// of the production cap — localizes which end (o = row 0, m = last
/// row) and which half-pair carries the remaining G residual.
#[test]
fn production_cap_g_residual_profile_is_reported() {
    use crate::geometry::gregory::kp20_bi6_cap::Kp20PieceId;
    use crate::geometry::gregory::kp20_bi6_local_constraint::{
        kp20_eq8_jacobian_dense, kp20_eq9_jacobian_dense,
    };
    use crate::geometry::gregory::kp20_bi6_math::{kp20_b_coeffs_near_m, kp20_b_coeffs_near_o};
    let n = 5;
    let report = assemble_kp20_paper_cap(&production_border(n), Point3::new(0.0, 0.0, 0.0))
        .expect("assemble");
    let sectors = &report.sectors;
    let c = (core::f64::consts::TAU / n as f64).cos();
    let p1 = Kp20PieceId::NearOriginLeft;
    let p2 = Kp20PieceId::NearOriginRight;
    let p4 = Kp20PieceId::NearMidpointLeft;
    let prev = n - 1;
    let curr = 0;
    for (label, prev_piece, curr_piece, b) in [
        ("near-o", p1, p1, kp20_b_coeffs_near_o(c)),
        ("near-m", p2, p4, kp20_b_coeffs_near_m(c)),
    ] {
        let eq8: Vec<[f64; 98]> = kp20_eq8_jacobian_dense(b).into_iter().collect();
        let eq9: Vec<[f64; 98]> = kp20_eq9_jacobian_dense(b).into_iter().collect();
        for (eq, jacobians) in [("eq8", &eq8), ("eq9", &eq9)] {
            let mut profile = Vec::new();
            for jac in jacobians {
                let mut acc = [0.0_f64; 3];
                for side in 0..2_usize {
                    for i in 0..7_usize {
                        for j in 0..7_usize {
                            let w = jac[side * 49 + i * 7 + j];
                            if w == 0.0 {
                                continue;
                            }
                            let cp = if side == 0 {
                                sectors[prev].piece(prev_piece).cps[i][j]
                            } else {
                                sectors[curr].piece(curr_piece).cps[j][i]
                            };
                            acc[0] += w * cp.x;
                            acc[1] += w * cp.y;
                            acc[2] += w * cp.z;
                        }
                    }
                }
                profile.push(format!(
                    "{:.1e}",
                    acc[0].abs().max(acc[1].abs()).max(acc[2].abs())
                ));
            }
            eprintln!("{label} {eq}: [{}]", profile.join(", "));
        }
    }
}
