use super::*;
use crate::geometry::Point3;
use crate::geometry::gregory::tri6_flat_index;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
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
fn default_piece_maps_stay_inside_sector_triangle() {
    let report = kp20_default_piece_map_domain_report();

    assert!(report.stays_inside_triangle());
    assert_eq!(report.outside_triangle_corners, 0);
    assert!(report.max_simplex_excess <= 1e-12);
    assert!(report.offending_corners.is_empty());
}
