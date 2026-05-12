use super::*;
use crate::geometry::gregory::Kp20PieceId;
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
            let idx = tri6_flat_index(6, i, j);
            let u = i as f64 / 6.0;
            let v = j as f64 / 6.0;
            guide.cps[idx] = p(2.0 * u + 3.0 * v + 1.0, -u + 0.5 * v, 4.0);
        }
    }
    guide
}

#[test]
fn tri6_guide_jet_reproduces_affine_value_and_partials() {
    let guide = affine_guide();

    let jet = guide.jet3(0.2, 0.3);

    assert_close(jet.h[0][0], p(2.3, -0.05, 4.0));
    assert_close(jet.h[1][0], p(2.0, -1.0, 0.0));
    assert_close(jet.h[0][1], p(3.0, 0.5, 0.0));
    assert_close(jet.h[2][0], p(0.0, 0.0, 0.0));
    assert_close(jet.h[1][1], p(0.0, 0.0, 0.0));
    assert_close(jet.h[0][2], p(0.0, 0.0, 0.0));
}

#[test]
fn corner_block_applies_affine_chain_rule() {
    let guide = affine_guide();
    let map = Kp20AffineGuideMap::new([0.1, 0.2], [0.5, 0.0], [0.0, 0.25]);

    let block = kp20_sample_corner_block(&guide, map);

    assert_close(block[0][0], p(1.8, 0.0, 4.0));
    assert_close(block[1][0], p(1.8 + 1.0 / 6.0, -1.0 / 12.0, 4.0));
    assert_close(block[0][1], p(1.8 + 0.75 / 6.0, 0.5 / 24.0, 4.0));
    assert_close(
        block[3][3],
        p(
            1.8 + 0.5 * 2.0 * 0.5 + 0.5 * 3.0 * 0.25,
            -0.25 + 0.0625,
            4.0,
        ),
    );
}

#[test]
fn patch_sampling_covers_full_bi6_grid_and_averages_overlaps() {
    let guide = affine_guide();
    let maps = [
        (
            Kp20PatchCorner::U0V0,
            Kp20AffineGuideMap::new([0.0, 0.0], [0.5, 0.0], [0.0, 0.5]),
        ),
        (
            Kp20PatchCorner::U1V0,
            Kp20AffineGuideMap::new([0.5, 0.0], [-0.5, 0.0], [0.0, 0.5]),
        ),
        (
            Kp20PatchCorner::U0V1,
            Kp20AffineGuideMap::new([0.0, 0.5], [0.5, 0.0], [0.0, -0.5]),
        ),
        (
            Kp20PatchCorner::U1V1,
            Kp20AffineGuideMap::new([0.5, 0.5], [-0.5, 0.0], [0.0, -0.5]),
        ),
    ];

    let patch = kp20_sample_patch_from_guide(&guide, maps);

    for i in 0..7 {
        for j in 0..7 {
            let u = 0.5 * i as f64 / 6.0;
            let v = 0.5 * j as f64 / 6.0;
            assert_close(
                patch.cps[i][j],
                p(2.0 * u + 3.0 * v + 1.0, -u + 0.5 * v, 4.0),
            );
        }
    }
}

#[test]
fn macro_patch_sampling_populates_all_four_quadrants() {
    let guide = affine_guide();
    let maps = [
        (
            Kp20PatchCorner::U0V0,
            Kp20AffineGuideMap::new([0.0, 0.0], [0.25, 0.0], [0.0, 0.25]),
        ),
        (
            Kp20PatchCorner::U1V0,
            Kp20AffineGuideMap::new([0.25, 0.0], [-0.25, 0.0], [0.0, 0.25]),
        ),
        (
            Kp20PatchCorner::U0V1,
            Kp20AffineGuideMap::new([0.0, 0.25], [0.25, 0.0], [0.0, -0.25]),
        ),
        (
            Kp20PatchCorner::U1V1,
            Kp20AffineGuideMap::new([0.25, 0.25], [-0.25, 0.0], [0.0, -0.25]),
        ),
    ];
    let piece_maps = [
        (Kp20PieceId::NearOriginLeft, maps),
        (Kp20PieceId::NearOriginRight, maps),
        (Kp20PieceId::NearMidpointRight, maps),
        (Kp20PieceId::NearMidpointLeft, maps),
    ];

    let macro_patch = kp20_sample_macro_patch_from_guide(&guide, piece_maps);

    for id in Kp20PieceId::all() {
        assert_close(macro_patch.piece(id).cps[0][0], p(1.0, 0.0, 4.0));
        assert_close(macro_patch.piece(id).cps[6][6], p(2.25, -0.125, 4.0));
    }
}
