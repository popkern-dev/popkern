use super::*;
use crate::geometry::Point3;
use crate::geometry::gregory::tri6_flat_index;

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

fn affine_guide() -> Kp20Tri6GuideSector {
    let mut guide = Kp20Tri6GuideSector::zero();
    for i in 0..=6 {
        for j in 0..=(6 - i) {
            let u = i as f64 / 6.0;
            let v = j as f64 / 6.0;
            guide.cps[tri6_flat_index(6, i, j)] = p(u, v, 0.0);
        }
    }
    guide
}

fn affine_tau_patch(origin: [f64; 2], du: [f64; 2], dv: [f64; 2]) -> Kp20MacroPatch {
    let mut tau = Kp20MacroPatch::zero();
    for piece_id in Kp20PieceId::all() {
        let piece = tau.piece_mut(piece_id);
        for i in 0..=6 {
            for j in 0..=6 {
                let u = i as f64 / 6.0;
                let v = j as f64 / 6.0;
                piece.cps[i][j] = p(
                    origin[0] + du[0] * u + dv[0] * v,
                    origin[1] + du[1] * u + dv[1] * v,
                    0.0,
                );
            }
        }
    }
    tau
}

#[test]
fn default_piece_maps_cover_four_macro_pieces() {
    let maps = kp20_default_sector_piece_maps();

    assert_eq!(maps[0].0, Kp20PieceId::NearOriginLeft);
    assert_eq!(maps[1].0, Kp20PieceId::NearOriginRight);
    assert_eq!(maps[2].0, Kp20PieceId::NearMidpointRight);
    assert_eq!(maps[3].0, Kp20PieceId::NearMidpointLeft);
}

#[test]
fn default_piece_maps_sample_all_pieces_from_guide() {
    let guide = affine_guide();

    let macro_patch = sample_kp20_macro_patch_from_default_maps(&guide);

    for id in Kp20PieceId::all() {
        let piece = macro_patch.piece(id);
        assert!(piece.cps[0][0].z.abs() < 1e-12);
        assert!(piece.cps[6][6].z.abs() < 1e-12);
    }
    assert_eq!(
        macro_patch.piece(Kp20PieceId::NearOriginLeft).cps[0][0],
        p(0.0, 0.0, 0.0)
    );
    assert_eq!(
        macro_patch.piece(Kp20PieceId::NearMidpointRight).cps[0][0],
        p(0.5, 0.5, 0.0)
    );
}

#[test]
fn orientation_valid_piece_maps_sample_all_pieces_from_guide() {
    let guide = affine_guide();

    let macro_patch =
        sample_kp20_macro_patch_from_piece_maps(&guide, kp20_orientation_valid_sector_piece_maps());

    for id in Kp20PieceId::all() {
        let piece = macro_patch.piece(id);
        assert!(piece.cps[0][0].z.abs() < 1e-12);
        assert!(piece.cps[6][6].z.abs() < 1e-12);
    }
    assert_eq!(
        macro_patch.piece(Kp20PieceId::NearMidpointRight).cps[0][0],
        p(0.5, 0.0, 0.0)
    );
}

#[test]
fn paper_oriented_piece_maps_preserve_default_footprint_with_positive_orientation() {
    let default_report = kp20_default_piece_map_domain_report();
    let paper_report = kp20_piece_map_domain_report(kp20_paper_oriented_sector_piece_maps());

    assert_eq!(default_report.outside_triangle_corners, 0);
    assert_eq!(paper_report.outside_triangle_corners, 0);
    assert_eq!(paper_report.non_positive_orientation_maps, 0);
    assert!(paper_report.orientation_offenders.is_empty());
}

#[test]
fn paper_oriented_piece_maps_sample_all_pieces_from_guide() {
    let guide = affine_guide();

    let macro_patch =
        sample_kp20_macro_patch_from_piece_maps(&guide, kp20_paper_oriented_sector_piece_maps());

    for id in Kp20PieceId::all() {
        let piece = macro_patch.piece(id);
        for row in piece.cps {
            for cp in row {
                assert!(cp.x.is_finite());
                assert!(cp.y.is_finite());
                assert!(cp.z.is_finite());
            }
        }
    }
}

#[test]
fn planar_bi6_tau_validity_report_accepts_bounded_positive_affine_map() {
    const SAMPLES: [f64; 3] = [0.0, 0.5, 1.0];
    let tau = affine_tau_patch([0.05, 0.05], [0.35, 0.0], [0.0, 0.35]);

    let report = kp20_planar_bi6_tau_validity_report(&tau, &SAMPLES, 1.0e-12, 1.0e-12);

    assert!(report.is_paper_valid_reparameterization());
    assert!(report.is_strict_sector_bounded_reparameterization());
    assert_eq!(report.domain.outside_triangle_samples, 0);
    assert_eq!(report.orientation.negative_determinants, 0);
    assert_eq!(report.orientation.near_zero_determinants, 0);
}

#[test]
fn planar_bi6_tau_validity_report_rejects_degenerate_and_unbounded_maps() {
    const SAMPLES: [f64; 3] = [0.0, 0.5, 1.0];
    let degenerate = affine_tau_patch([0.05, 0.05], [0.25, 0.0], [0.5, 0.0]);
    let unbounded = affine_tau_patch([0.75, 0.75], [0.25, 0.0], [0.0, 0.25]);

    let degenerate_report =
        kp20_planar_bi6_tau_validity_report(&degenerate, &SAMPLES, 1.0e-12, 1.0e-12);
    let unbounded_report =
        kp20_planar_bi6_tau_validity_report(&unbounded, &SAMPLES, 1.0e-12, 1.0e-12);

    assert!(!degenerate_report.is_paper_valid_reparameterization());
    assert_eq!(degenerate_report.orientation.negative_determinants, 0);
    assert_eq!(
        degenerate_report.orientation.near_zero_determinants,
        degenerate_report.orientation.sample_count
    );
    assert!(unbounded_report.is_paper_valid_reparameterization());
    assert!(!unbounded_report.is_strict_sector_bounded_reparameterization());
    assert!(unbounded_report.domain.outside_triangle_samples > 0);
}

#[test]
fn valid_planar_bi6_tau_sampler_declines_invalid_reparameterizations() {
    const SAMPLES: [f64; 3] = [0.0, 0.5, 1.0];
    let guide = affine_guide();
    let valid = affine_tau_patch([0.05, 0.05], [0.35, 0.0], [0.0, 0.35]);
    let degenerate = affine_tau_patch([0.05, 0.05], [0.25, 0.0], [0.5, 0.0]);
    let unbounded = affine_tau_patch([0.75, 0.75], [0.25, 0.0], [0.0, 0.25]);

    let sampled = sample_kp20_macro_patch_from_valid_planar_bi6_tau_macro_patch(
        &guide, &valid, &SAMPLES, 1.0e-12, 1.0e-12,
    )
    .expect("bounded positive affine tau should sample");
    assert_close(
        sampled.piece(Kp20PieceId::NearOriginLeft).cps[0][0],
        p(0.05, 0.05, 0.0),
    );

    let degenerate_error = sample_kp20_macro_patch_from_valid_planar_bi6_tau_macro_patch(
        &guide,
        &degenerate,
        &SAMPLES,
        1.0e-12,
        1.0e-12,
    )
    .expect_err("degenerate tau must be rejected before cap sampling");
    assert_eq!(
        degenerate_error.validity.orientation.near_zero_determinants,
        degenerate_error.validity.orientation.sample_count
    );

    let unbounded_error = sample_kp20_macro_patch_from_valid_planar_bi6_tau_macro_patch(
        &guide, &unbounded, &SAMPLES, 1.0e-12, 1.0e-12,
    )
    .expect_err("unbounded tau must be rejected before cap sampling");
    assert!(unbounded_error.validity.domain.outside_triangle_samples > 0);
}

#[test]
fn paper_valid_planar_bi6_tau_sampler_allows_polynomial_guide_extension() {
    const SAMPLES: [f64; 3] = [0.0, 0.5, 1.0];
    let guide = affine_guide();
    let unbounded = affine_tau_patch([0.75, 0.75], [0.25, 0.0], [0.0, 0.25]);
    let degenerate = affine_tau_patch([0.05, 0.05], [0.25, 0.0], [0.5, 0.0]);

    let sampled = sample_kp20_macro_patch_from_paper_valid_planar_bi6_tau_macro_patch(
        &guide, &unbounded, &SAMPLES, 1.0e-12, 1.0e-12,
    )
    .expect("paper-valid tau may sample polynomial guide outside the unit sector");
    assert_close(
        sampled.piece(Kp20PieceId::NearOriginLeft).cps[0][0],
        p(0.75, 0.75, 0.0),
    );

    let degenerate_error = sample_kp20_macro_patch_from_paper_valid_planar_bi6_tau_macro_patch(
        &guide,
        &degenerate,
        &SAMPLES,
        1.0e-12,
        1.0e-12,
    )
    .expect_err("paper-valid tau sampler must still reject singular maps");
    assert_eq!(
        degenerate_error.validity.orientation.near_zero_determinants,
        degenerate_error.validity.orientation.sample_count
    );
}

#[test]
fn c2_stitched_piece_maps_are_inside_positive_and_edge_compatible() {
    let maps = kp20_c2_stitched_sector_piece_maps();
    let domain = kp20_piece_map_domain_report(maps);
    let stitch = kp20_piece_map_stitch_report(maps);

    assert!(domain.stays_inside_triangle());
    assert_eq!(domain.non_positive_orientation_maps, 0);
    assert!(stitch.max_internal_c2_edge_mismatch <= 1.0e-12);
    assert!(stitch.max_origin_strip_edge_mismatch > 1.0e-12);
}

#[test]
fn c2_stitched_piece_maps_sample_all_pieces_from_guide() {
    let guide = affine_guide();

    let macro_patch =
        sample_kp20_macro_patch_from_piece_maps(&guide, kp20_c2_stitched_sector_piece_maps());

    for id in Kp20PieceId::all() {
        let piece = macro_patch.piece(id);
        for row in piece.cps {
            for cp in row {
                assert!(cp.x.is_finite());
                assert!(cp.y.is_finite());
                assert!(cp.z.is_finite());
            }
        }
    }
}

#[test]
fn coupled_stitched_piece_maps_are_inside_positive_and_edge_compatible() {
    let maps = kp20_coupled_stitched_sector_piece_maps();
    let domain = kp20_piece_map_domain_report(maps);
    let stitch = kp20_piece_map_stitch_report(maps);

    assert!(domain.stays_inside_triangle());
    assert_eq!(domain.non_positive_orientation_maps, 0);
    assert!(stitch.max_internal_c2_edge_mismatch <= 1.0e-12);
    assert!(stitch.max_origin_strip_edge_mismatch <= 1.0e-12);
    assert!(!stitch.origin_strip_reversed);
    assert!(stitch.min_piece_center_separation <= 1.0e-12);
}

#[test]
fn nonoverlap_coupled_piece_maps_are_inside_positive_and_separated() {
    let maps = kp20_positive_nonoverlap_c2_sector_piece_maps();
    let domain = kp20_piece_map_domain_report(maps);
    let stitch = kp20_piece_map_stitch_report(maps);

    assert!(domain.stays_inside_triangle());
    assert_eq!(domain.non_positive_orientation_maps, 0);
    assert!(stitch.max_internal_c2_edge_mismatch <= 1.0e-12);
    assert!(stitch.max_origin_strip_edge_mismatch <= 1.0e-12);
    assert!(stitch.origin_strip_reversed);
    assert!(stitch.min_piece_center_separation > 0.24);
}

#[test]
fn paper_quadrant_piece_maps_match_figure8_domain_layout() {
    let maps = kp20_paper_quadrant_sector_piece_maps();
    let domain = kp20_piece_map_domain_report(maps);
    let stitch = kp20_piece_map_stitch_report(maps);

    assert!(!domain.stays_inside_triangle());
    assert_eq!(domain.outside_triangle_corners, 5);
    assert_eq!(domain.non_positive_orientation_maps, 0);
    assert!(stitch.min_piece_center_separation > 0.49);

    assert_eq!(maps[0].0, Kp20PieceId::NearOriginLeft);
    assert_eq!(maps[0].1[0].1.origin, [0.0, 0.0]);
    assert_eq!(maps[1].0, Kp20PieceId::NearOriginRight);
    assert_eq!(maps[1].1[0].1.origin, [0.5, 0.0]);
    assert_eq!(maps[2].0, Kp20PieceId::NearMidpointRight);
    assert_eq!(maps[2].1[0].1.origin, [0.5, 0.5]);
    assert_eq!(maps[3].0, Kp20PieceId::NearMidpointLeft);
    assert_eq!(maps[3].1[0].1.origin, [0.0, 0.5]);
}

#[test]
fn coupled_stitched_piece_maps_sample_all_pieces_from_guide() {
    let guide = affine_guide();

    let macro_patch =
        sample_kp20_macro_patch_from_piece_maps(&guide, kp20_coupled_stitched_sector_piece_maps());

    for id in Kp20PieceId::all() {
        let piece = macro_patch.piece(id);
        for row in piece.cps {
            for cp in row {
                assert!(cp.x.is_finite());
                assert!(cp.y.is_finite());
                assert!(cp.z.is_finite());
            }
        }
    }
}

#[test]
fn nonoverlap_coupled_piece_maps_sample_all_pieces_from_guide() {
    let guide = affine_guide();

    let macro_patch = sample_kp20_macro_patch_from_piece_maps(
        &guide,
        kp20_positive_nonoverlap_c2_sector_piece_maps(),
    );

    for id in Kp20PieceId::all() {
        let piece = macro_patch.piece(id);
        for row in piece.cps {
            for cp in row {
                assert!(cp.x.is_finite());
                assert!(cp.y.is_finite());
                assert!(cp.z.is_finite());
            }
        }
    }
}

#[test]
fn stitch_report_fails_closed_for_missing_semantic_piece() {
    let mut maps = kp20_c2_stitched_sector_piece_maps();
    maps[0].0 = Kp20PieceId::NearMidpointLeft;

    let stitch = kp20_piece_map_stitch_report(maps);

    assert!(stitch.max_internal_c2_edge_mismatch.is_infinite());
    assert!(stitch.max_origin_strip_edge_mismatch.is_infinite());
}

#[test]
fn tau_piece_maps_sample_all_pieces_from_guide() {
    let guide = affine_guide();

    let macro_patch = sample_kp20_macro_patch_from_tau_maps(&guide, 6, 4).unwrap();

    for id in Kp20PieceId::all() {
        let piece = macro_patch.piece(id);
        for row in piece.cps {
            for cp in row {
                assert!(cp.x.is_finite());
                assert!(cp.y.is_finite());
                assert!(cp.z.is_finite());
            }
        }
    }
    assert_eq!(
        macro_patch.piece(Kp20PieceId::NearOriginLeft).cps[0][0],
        p(0.0, 0.0, 0.0)
    );
}

#[test]
fn paper_tau_piece_maps_sample_all_pieces_from_guide() {
    let guide = affine_guide();

    let macro_patch = sample_kp20_macro_patch_from_paper_tau_maps(&guide, 6, 4).unwrap();

    for id in Kp20PieceId::all() {
        let piece = macro_patch.piece(id);
        for row in piece.cps {
            for cp in row {
                assert!(cp.x.is_finite());
                assert!(cp.y.is_finite());
                assert!(cp.z.is_finite());
            }
        }
    }
}

#[test]
fn per_piece_sampler_maps_allow_origin_tau_with_affine_outer_pieces() {
    let guide = affine_guide();
    let samplers = [
        (
            Kp20PieceId::NearOriginLeft,
            Kp20PieceSamplerFamily::Tau(Kp20TauFamily::PaperRefinedApex),
        ),
        (
            Kp20PieceId::NearOriginRight,
            Kp20PieceSamplerFamily::Tau(Kp20TauFamily::PaperRefinedApex),
        ),
        (
            Kp20PieceId::NearMidpointRight,
            Kp20PieceSamplerFamily::AffineGuide,
        ),
        (
            Kp20PieceId::NearMidpointLeft,
            Kp20PieceSamplerFamily::AffineGuide,
        ),
    ];

    let macro_patch = sample_kp20_macro_patch_from_piece_sampler_maps(
        &guide,
        6,
        4,
        samplers,
        kp20_positive_nonoverlap_c2_sector_piece_maps(),
    )
    .unwrap();

    for id in Kp20PieceId::all() {
        let piece = macro_patch.piece(id);
        for row in piece.cps {
            for cp in row {
                assert!(cp.x.is_finite());
                assert!(cp.y.is_finite());
                assert!(cp.z.is_finite());
            }
        }
    }
}

#[test]
fn planar_bi6_tau_sampler_matches_affine_piece_maps_on_affine_guide() {
    let guide = affine_guide();
    let maps = kp20_positive_nonoverlap_c2_sector_piece_maps();
    let expected = sample_kp20_macro_patch_from_piece_maps(&guide, maps);
    let mut tau_macro_patch = Kp20MacroPatch::zero();

    for (piece_id, corner_maps) in maps {
        let base = corner_maps[0].1;
        let piece = tau_macro_patch.piece_mut(piece_id);
        for i in 0..=6 {
            for j in 0..=6 {
                let u = i as f64 / 6.0;
                let v = j as f64 / 6.0;
                piece.cps[i][j] = p(
                    base.origin[0] + u * base.du[0] + v * base.dv[0],
                    base.origin[1] + u * base.du[1] + v * base.dv[1],
                    0.0,
                );
            }
        }
    }

    let actual = sample_kp20_macro_patch_from_planar_bi6_tau_macro_patch(&guide, &tau_macro_patch);

    for piece_id in Kp20PieceId::all() {
        for i in 0..=6 {
            for j in 0..=6 {
                assert_close(
                    actual.piece(piece_id).cps[i][j],
                    expected.piece(piece_id).cps[i][j],
                );
            }
        }
    }
}

#[test]
fn planar_bi6_tau_f3_energy_is_zero_for_affine_maps() {
    let mut tau_macro_patch = Kp20MacroPatch::zero();

    for (piece_id, corner_maps) in kp20_positive_nonoverlap_c2_sector_piece_maps() {
        let base = corner_maps[0].1;
        let piece = tau_macro_patch.piece_mut(piece_id);
        for i in 0..=6 {
            for j in 0..=6 {
                let u = i as f64 / 6.0;
                let v = j as f64 / 6.0;
                piece.cps[i][j] = p(
                    base.origin[0] + u * base.du[0] + v * base.dv[0],
                    base.origin[1] + u * base.du[1] + v * base.dv[1],
                    0.0,
                );
            }
        }
    }

    assert!(kp20_planar_bi6_tau_f3_energy(&tau_macro_patch) <= 1.0e-18);
}

#[test]
fn planar_bi6_tau_orientation_report_accepts_positive_affine_maps() {
    let mut tau_macro_patch = Kp20MacroPatch::zero();

    for (piece_id, corner_maps) in kp20_positive_nonoverlap_c2_sector_piece_maps() {
        let base = corner_maps[0].1;
        let piece = tau_macro_patch.piece_mut(piece_id);
        for i in 0..=6 {
            for j in 0..=6 {
                let u = i as f64 / 6.0;
                let v = j as f64 / 6.0;
                piece.cps[i][j] = p(
                    base.origin[0] + u * base.du[0] + v * base.dv[0],
                    base.origin[1] + u * base.du[1] + v * base.dv[1],
                    0.0,
                );
            }
        }
    }

    let report =
        kp20_planar_bi6_tau_orientation_report(&tau_macro_patch, &[0.0, 0.5, 1.0], 1.0e-12);

    assert_eq!(report.sample_count, 36);
    assert_eq!(report.negative_determinants, 0);
    assert_eq!(report.near_zero_determinants, 0);
    assert!(report.is_orientation_preserving());
    assert!(report.min_determinant > 0.0);
}

#[test]
fn planar_bi6_tau_f3_energy_detects_cubic_curvature() {
    let mut tau_macro_patch = Kp20MacroPatch::zero();

    for piece_id in Kp20PieceId::all() {
        let piece = tau_macro_patch.piece_mut(piece_id);
        for i in 0..=6 {
            for j in 0..=6 {
                let u = i as f64 / 6.0;
                let v = j as f64 / 6.0;
                piece.cps[i][j] = p(u, v + 0.1 * u * u * u, 0.0);
            }
        }
    }

    let energy = kp20_planar_bi6_tau_f3_energy(&tau_macro_patch);

    assert!(energy.is_finite());
    assert!(energy > 1.0e-6);
}

#[test]
fn paper_tau_sampler_matches_tau_value_on_affine_guide_at_piece_corners() {
    let guide = affine_guide();
    let valence = 6;
    let sector = 4;
    let tau = Kp20TauMap::for_sector_local_paper(valence, sector).unwrap();
    let maps = kp20_positive_nonoverlap_c2_sector_piece_maps();

    let macro_patch = sample_kp20_macro_patch_from_tau_family_maps(
        &guide,
        valence,
        sector,
        Kp20TauFamily::Paper,
        maps,
    )
    .unwrap();

    for (piece_id, corner_maps) in maps {
        let piece = macro_patch.piece(piece_id);
        for (corner, map) in corner_maps {
            let value = tau.jet3(map.origin[0], map.origin[1]).d[0][0];
            let expected = p(value[0], value[1], 0.0);
            let (i, j) = match corner {
                Kp20PatchCorner::U0V0 => (0, 0),
                Kp20PatchCorner::U1V0 => (6, 0),
                Kp20PatchCorner::U0V1 => (0, 6),
                Kp20PatchCorner::U1V1 => (6, 6),
            };
            assert_close(piece.cps[i][j], expected);
        }
    }
}

#[test]
fn tau_piece_map_domain_report_separates_kernel_and_paper_tau() {
    let maps = kp20_positive_nonoverlap_c2_sector_piece_maps();

    let kernel = kp20_tau_piece_map_domain_report(maps, 5, 0, Kp20TauFamily::Kernel).unwrap();
    let paper = kp20_tau_piece_map_domain_report(maps, 5, 0, Kp20TauFamily::Paper).unwrap();

    assert_eq!(kernel.tau_family, Kp20TauFamily::Kernel);
    assert_eq!(paper.tau_family, Kp20TauFamily::Paper);
    assert_eq!(kernel.sample_count, 4 * 4 * 6);
    assert_eq!(paper.sample_count, 4 * 4 * 6);
    assert!(kernel.min_u.is_finite());
    assert!(kernel.max_v.is_finite());
    assert!(paper.min_u.is_finite());
    assert!(paper.max_v.is_finite());
    assert!(kernel.max_simplex_excess >= 0.0);
    assert!(paper.max_simplex_excess >= 0.0);
    assert_ne!(
        kernel.outside_triangle_samples,
        usize::MAX,
        "report must populate a real outside count"
    );
    assert_ne!(
        paper.outside_triangle_samples,
        usize::MAX,
        "report must populate a real outside count"
    );
}

#[test]
fn explicit_paper_tau_piece_map_validity_is_measured() {
    let valence = 5;
    let sector = 0;
    let cases = [
        (
            "positive-nonoverlap",
            kp20_positive_nonoverlap_c2_sector_piece_maps(),
        ),
        ("paper-quadrant", kp20_paper_quadrant_sector_piece_maps()),
    ];

    for (label, maps) in cases {
        let validity = kp20_tau_piece_map_validity_report(
            maps,
            valence,
            sector,
            Kp20TauFamily::Paper,
            1.0e-10,
        )
        .unwrap();
        eprintln!(
            "explicit paper tau validity {label}: outside={outside}/{domain_samples} max_excess={max_excess:e} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            max_excess = validity.domain.max_simplex_excess,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            worst = validity.orientation.worst_piece,
            u = validity.orientation.worst_u,
            v = validity.orientation.worst_v,
        );
        assert!(validity.domain.sample_count > 0);
        assert_eq!(
            validity.domain.sample_count,
            validity.orientation.sample_count
        );
        assert!(validity.domain.max_simplex_excess.is_finite());
        assert!(validity.orientation.min_determinant.is_finite());
        assert!(validity.orientation.max_determinant.is_finite());
        assert!(validity.is_paper_valid_reparameterization());
        assert!(!validity.is_strict_sector_bounded_reparameterization());
    }
}

#[test]
fn default_piece_maps_stay_inside_sector_triangle() {
    let report = kp20_default_piece_map_domain_report();

    assert!(report.stays_inside_triangle());
    assert_eq!(report.outside_triangle_corners, 0);
    assert!(report.max_simplex_excess <= 1e-12);
    assert!(report.offending_corners.is_empty());
}

#[test]
fn default_piece_maps_report_non_positive_orientation_scaffold() {
    let report = kp20_default_piece_map_domain_report();

    assert_eq!(report.non_positive_orientation_maps, 3);
    assert!(
        report
            .orientation_offenders
            .iter()
            .any(|offender| offender.starts_with("NearOriginRight:"))
    );
    assert!(
        report
            .orientation_offenders
            .iter()
            .any(|offender| offender.starts_with("NearMidpointRight:"))
    );
    assert!(
        report
            .orientation_offenders
            .iter()
            .any(|offender| offender.starts_with("NearMidpointLeft:"))
    );
}
