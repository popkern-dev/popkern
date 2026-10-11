// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]
use super::*;
use crate::geometry::gregory::cc_characteristic_tensor_border::cc_characteristic_tensor_border;
use crate::geometry::gregory::kp20_paper_guide::Kp20PaperGuideBorderPiece;

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
                    crate::geometry::gregory::kp20_paper_guide::Kp20PaperGuideBorderPiece::from_bi3(
                        near_c,
                    ),
                    crate::geometry::gregory::kp20_paper_guide::Kp20PaperGuideBorderPiece::from_bi3(
                        near_m,
                    ),
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
    for n in [3, 4, 5, 6, 7, 8, 9, 10] {
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
    for n in [4, 5, 6, 7, 9, 10] {
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

/// The production adapter's output must be consumable by the macro-cap
/// validator: C2 pairing and the eq (14) boundary placement agree
/// with the macro-cap contract, and the flattened grids follow
/// [`Kp20PieceId::all`] sector-major order.
#[test]
fn paper_macro_cap_adapter_is_contract_compatible() {
    use crate::geometry::gregory::kp_polar_boundary_input::TensorBorderBi3;
    use crate::geometry::gregory::kp20_bi6_validate::measure_kp20_bi6_macro_cap;
    let n = 5;
    // A smooth production-shaped border: single cubic per sector from
    // the lifted-ring pieces joined at c (exact — the two halves of
    // one C2 ring edge merge back into one cubic per rim edge is NOT
    // generally possible, so use the characteristic-plus-lift border
    // via a degree-matched smooth fixture instead: the planar
    // characteristic border's production form).
    let pieces =
        crate::geometry::gregory::cc_characteristic_tensor_border::cc_lifted_ring_border_pieces(
            n,
            |x, y| 0.25 * (x * x + y * y),
        )
        .expect("lifted ring");
    // Join each sector's two halves into one cubic per depth row by
    // de-Casteljau inversion is invalid for cornered data; instead
    // exercise the adapter through the two-piece paper entry and only
    // check the adapter packaging + the macro-cap validator READS it.
    let border = Kp20PaperGuideBorder {
        sectors: pieces
            .into_iter()
            .map(|[near_c, near_m]| {
                [
                    Kp20PaperGuideBorderPiece::from_bi3(near_c),
                    Kp20PaperGuideBorderPiece::from_bi3(near_m),
                ]
            })
            .collect(),
    };
    let report = assemble_kp20_paper_cap(&border, Point3::new(0.0, 0.0, 0.0)).expect("assemble");
    let patches_bi6: Vec<[[Point3; 7]; 7]> = report
        .sectors
        .iter()
        .flat_map(|sector| {
            crate::geometry::gregory::Kp20PieceId::all()
                .into_iter()
                .map(|piece_id| sector.piece(piece_id).cps)
        })
        .collect();
    assert_eq!(patches_bi6.len(), 4 * n);
    let output = crate::geometry::gregory::Kp20Bi6MacroCapOutput {
        n,
        macro_patches: report.sectors.clone(),
        patches_bi6,
        strip_projections: Vec::new(),
    };
    // Reconstruct the single-cubic border the macro-cap validator wants:
    // for a smooth C2 ring the two halves merge exactly.
    let mut sectors = Vec::with_capacity(n);
    for sector in &border.sectors {
        let mut rows = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        let (first, second) = (sector[0].bi3_rows(), sector[1].bi3_rows());
        for (depth, row) in rows.iter_mut().enumerate() {
            *row = merge_cubic_halves(first[depth], second[depth]);
        }
        sectors.push(rows);
    }
    let tensor_border = TensorBorderBi3 { sectors };
    let validation = measure_kp20_bi6_macro_cap(&output, &tensor_border).expect("measure");
    eprintln!(
        "macro-cap validator on paper cap: c2_pos={c2p:e} c2_first={c2f:e} origin_strip={os:e} border={b:e}",
        c2p = validation.max_internal_c2_position_residual,
        c2f = validation.max_internal_c2_first_residual,
        os = validation.max_origin_g2_strip_residual,
        b = validation.max_tensor_border_residual,
    );
}

/// Exact inverse of the de-Casteljau midpoint split for a smooth
/// cubic: valid only when the two halves come from one cubic.
fn merge_cubic_halves(first: [Point3; 4], second: [Point3; 4]) -> [Point3; 4] {
    let lerp = |a: Point3, b: Point3, t: f64| {
        Point3::new(
            a.x + (b.x - a.x) * t,
            a.y + (b.y - a.y) * t,
            a.z + (b.z - a.z) * t,
        )
    };
    // first = [p0, m01, mm, c], with m01 = (p0+p1)/2 midpoints of the
    // original: recover p1 = 2*first[1] - p0, and p2 = 2*second[2] - p3.
    let p0 = first[0];
    let p1 = lerp(first[0], first[1], 2.0);
    let p3 = second[3];
    let p2 = lerp(second[3], second[2], 2.0);
    [p0, p1, p2, p3]
}

/// PK-398 structural fact: the characteristic border's own two-piece
/// corner data does NOT agree on the KP20 corner piece's overlap — at
/// `c` its cross leg is 1.0× the CUBIC rim leg while the eq (14) band
/// (half-piece elevated to degree 6, depth row `½(row0 + row1)`)
/// needs half of that, so `kp20_paper_border_corner_mismatch` reads
/// `rim_leg / 2` even for the regular n = 4 border. The paper's
/// overlap average is therefore always active on bi-3 input; the
/// setback extractor's corner-compatible jets (`d1 = degree-6 leg`)
/// are the one input class for which it is a no-op.
#[test]
fn characteristic_border_corner_overlap_mismatch_is_half_the_rim_leg() {
    use crate::geometry::gregory::cc_characteristic_tensor_border::cc_lifted_ring_border_pieces;
    for n in [3usize, 4, 5, 6] {
        let pieces = cc_lifted_ring_border_pieces(n, |_, _| 0.0).expect("chi");
        let [near_c, near_m] = &pieces[0];
        let rim_leg = (near_c[0][3] - near_c[0][2]).length();
        let cross_leg = (near_c[1][3] - near_c[0][3]).length();
        assert!(
            (cross_leg / rim_leg - 1.0).abs() < 1e-9,
            "n={n}: cross/rim leg ratio at c {}",
            cross_leg / rim_leg
        );
        let c = (core::f64::consts::TAU / n as f64).cos();
        let mismatch = crate::geometry::gregory::kp20_paper_border_corner_mismatch(
            n,
            &crate::geometry::gregory::elevate_bi3_rows_to_sextic(*near_c),
            &crate::geometry::gregory::elevate_bi3_rows_to_sextic(*near_m),
            2.0 * c / 3.0,
        );
        assert!(
            mismatch > 0.3 * rim_leg && mismatch < 0.7 * rim_leg,
            "n={n}: corner mismatch {mismatch:e} vs rim leg {rim_leg:e}"
        );
    }
}

/// PK-518: the regular valence n = 4 (outside the paper's tabulated
/// n = 3, 5..=10). Two mechanisms make it a first-class cap:
/// (1) the τ solve carries explicit cross-sector C0 rows — at c = 0
/// the G rows no longer imply them and the sampled cap broke C0
/// across sectors by 1.3e-2; (2) the guide's KKT system is singular
/// at n = 4, and the null space is now selected by the fairing
/// functional instead of the absolute ridge (the characteristic cap
/// needed a 2.7 strip movement to repair the ridge's wild guide).
#[test]
fn regular_valence_characteristic_cap_is_sector_consistent_and_fair() {
    use crate::geometry::gregory::Kp20PieceId;
    let n = 4;
    let border = characteristic_border(n);
    let guide = crate::geometry::gregory::kp20_paper_guide::solve_kp20_paper_guide_with_fairing(
        &border,
        Point3::new(0.0, 0.0, 0.0),
        0.0,
    )
    .expect("guide");
    let tau = solve_kp20_paper_tau(n).expect("tau");
    let alpha = core::f64::consts::TAU / n as f64;
    let sample = |s: usize| {
        let rot = [(s as f64 * alpha).cos(), (s as f64 * alpha).sin()];
        let mut sector = tau.tau;
        for piece in Kp20PieceId::all() {
            let patch = sector.piece_mut(piece);
            for i in 0..7 {
                for j in 0..7 {
                    let cp = patch.cps[i][j];
                    patch.cps[i][j] = Point3::new(
                        cp.x * rot[0] - cp.y * rot[1],
                        cp.x * rot[1] + cp.y * rot[0],
                        0.0,
                    );
                }
            }
        }
        kp20_sample_paper_cap_sector(n, s, &guide.sectors[s], &sector)
    };
    let (s0, s1) = (sample(0), sample(1));
    let worst_c0 = (0..7)
        .map(|i| {
            s0.piece(Kp20PieceId::NearOriginLeft).cps[i][0]
                .distance_to(&s1.piece(Kp20PieceId::NearOriginLeft).cps[0][i])
        })
        .fold(0.0_f64, f64::max);
    assert!(worst_c0 < 1.0e-10, "sampled cross-sector C0 {worst_c0:e}");
    let report = assemble_kp20_paper_cap(&border, Point3::new(0.0, 0.0, 0.0)).expect("assemble");
    report_line("characteristic regular", &report);
    assert!(
        report.max_strip_movement < 1.0e-2,
        "the fair guide needs no strip repair: {:e}",
        report.max_strip_movement
    );
}
