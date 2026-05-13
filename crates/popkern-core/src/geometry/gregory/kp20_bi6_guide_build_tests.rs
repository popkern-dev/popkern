use super::*;
use crate::geometry::gregory::tri6_flat_index;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn sector(offset: f64) -> [[Point3; 4]; 3] {
    [
        [
            p(offset + 0.0, 0.0, 0.0),
            p(offset + 1.0, 0.0, 0.0),
            p(offset + 2.0, 0.0, 0.0),
            p(offset + 3.0, 0.0, 0.0),
        ],
        [p(0.0, 0.0, 0.0); 4],
        [p(0.0, 0.0, 0.0); 4],
    ]
}

#[test]
fn guide_builder_creates_one_tri6_sector_per_tensor_border_sector() {
    let input = Kp20GuideBuildInput {
        star_point: p(0.0, 0.0, 0.0),
        tensor_border: TensorBorderBi3 {
            sectors: vec![sector(0.0), sector(10.0), sector(20.0)],
        },
    };

    let guide = build_kp20_guide_ring(&input).unwrap();

    assert_eq!(guide.n, 3);
    assert_eq!(guide.sectors.len(), 3);
}

#[test]
fn guide_builder_pins_star_and_outer_boundary() {
    let input = Kp20GuideBuildInput {
        star_point: p(7.0, 8.0, 9.0),
        tensor_border: TensorBorderBi3 {
            sectors: vec![sector(0.0), sector(10.0), sector(20.0)],
        },
    };

    let guide = build_kp20_guide_ring(&input).unwrap();
    let g = guide.sector(0);

    assert_eq!(g.cps[tri6_flat_index(6, 0, 0)], p(7.0, 8.0, 9.0));
    assert_eq!(g.cps[tri6_flat_index(6, 6, 0)], p(0.0, 0.0, 0.0));
    assert_eq!(g.cps[tri6_flat_index(6, 0, 6)], p(3.0, 0.0, 0.0));
}

#[test]
fn guide_builder_preserves_planarity() {
    let input = Kp20GuideBuildInput {
        star_point: p(0.0, 0.0, 0.0),
        tensor_border: TensorBorderBi3 {
            sectors: vec![sector(0.0), sector(10.0), sector(20.0)],
        },
    };

    let guide = build_kp20_guide_ring(&input).unwrap();

    for sector in guide.sectors {
        for cp in sector.cps {
            assert!(cp.z.abs() < 1e-12);
        }
    }
}

#[test]
fn guide_builder_rejects_invalid_valence() {
    let input = Kp20GuideBuildInput {
        star_point: p(0.0, 0.0, 0.0),
        tensor_border: TensorBorderBi3 {
            sectors: vec![sector(0.0), sector(10.0)],
        },
    };

    assert!(build_kp20_guide_ring(&input).is_err());
}

#[test]
fn guide_builder_falls_back_for_regular_valence_four() {
    let input = Kp20GuideBuildInput {
        star_point: p(1.0, 2.0, 3.0),
        tensor_border: TensorBorderBi3 {
            sectors: vec![sector(0.0), sector(10.0), sector(20.0), sector(30.0)],
        },
    };

    let guide = build_kp20_guide_ring(&input).unwrap();

    assert_eq!(guide.n, 4);
    assert_eq!(
        guide.sector(0).cps[tri6_flat_index(6, 0, 0)],
        p(1.0, 2.0, 3.0)
    );
}
