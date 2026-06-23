use super::*;
use crate::geometry::gregory::tri6_flat_index;

const EPS: f64 = 1e-12;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn assert_close(a: Point3, b: Point3) {
    assert!(
        a.distance_to(&b) < EPS,
        "expected {a:?} to be within {EPS} of {b:?}"
    );
}

fn guide_sector() -> Kp20Tri6GuideSector {
    let mut sector = Kp20Tri6GuideSector::zero();
    for i in 0..=6 {
        for j in 0..=(6 - i) {
            let k = 6 - i - j;
            sector.cps[tri6_flat_index(6, i, j)] = p(
                i as f64 + 0.1 * j as f64,
                j as f64 - 0.2 * k as f64,
                k as f64 + 0.3 * i as f64,
            );
        }
    }
    sector
}

fn point_sum(terms: &[(f64, Point3)]) -> Point3 {
    let mut out = p(0.0, 0.0, 0.0);
    for (weight, point) in terms {
        out.x += weight * point.x;
        out.y += weight * point.y;
        out.z += weight * point.z;
    }
    out
}

fn get(sector: &Kp20Tri6GuideSector, i: usize, j: usize, k: usize) -> Point3 {
    assert_eq!(i + j + k, 6);
    sector.cps[tri6_flat_index(6, i, j)]
}

#[test]
fn guide_join_weights_match_kp18_eq1() {
    let (w0, w1, w2) = kp18_guide_join_weights(4).unwrap();
    assert!((w0 + 1.0).abs() < EPS);
    assert!(w1.abs() < EPS);
    assert!((w2 - 2.0).abs() < EPS);
    assert!(kp18_guide_join_weights(2).is_err());
}

#[test]
#[allow(deprecated)]
fn c0_join_identifies_shared_edge() {
    let left = guide_sector();
    let mut right = Kp20Tri6GuideSector::zero();
    kp18_apply_c0_join_degree6(&left, &mut right);

    for i in 0..=6 {
        assert_close(get(&right, 6 - i, i, 0), get(&left, 6 - i, 0, i));
    }
}

#[test]
#[allow(deprecated)]
fn c1_join_matches_kp18_eq3() {
    let left = guide_sector();
    let mut right = Kp20Tri6GuideSector::zero();
    let n = 5;
    let (w0, w1, w2) = kp18_guide_join_weights(n).unwrap();
    kp18_apply_c1_join_degree6(&left, &mut right, n).unwrap();

    for i in 1..=6 {
        let expected = point_sum(&[
            (w0, get(&left, 6 - i, 1, i - 1)),
            (w1, get(&left, 6 - i, 0, i)),
            (w2, get(&left, 6 - i + 1, 0, i - 1)),
        ]);
        assert_close(get(&right, 6 - i, i - 1, 1), expected);
    }
}

#[test]
#[allow(deprecated)]
fn c2_join_matches_kp18_eq4() {
    let left = guide_sector();
    let mut right = Kp20Tri6GuideSector::zero();
    let n = 6;
    let (w0, w1, w2) = kp18_guide_join_weights(n).unwrap();
    kp18_apply_c2_join_degree6(&left, &mut right, n).unwrap();

    for i in 2..=6 {
        let expected = point_sum(&[
            (w0 * w0, get(&left, 6 - i, 2, i - 2)),
            (2.0 * w0 * w1, get(&left, 6 - i, 1, i - 1)),
            (w1 * w1, get(&left, 6 - i, 0, i)),
            (2.0 * w0 * w2, get(&left, 6 - i + 1, 1, i - 2)),
            (2.0 * w1 * w2, get(&left, 6 - i + 1, 0, i - 1)),
            (w2 * w2, get(&left, 6 - i + 2, 0, i - 2)),
        ]);
        assert_close(get(&right, 6 - i, i - 2, 2), expected);
    }
}

#[test]
#[allow(deprecated)]
fn c3_join_matches_kp18_eq5() {
    let left = guide_sector();
    let mut right = Kp20Tri6GuideSector::zero();
    let n = 7;
    let (w0, w1, w2) = kp18_guide_join_weights(n).unwrap();
    kp18_apply_c3_join_degree6(&left, &mut right, n).unwrap();

    for i in 3..=6 {
        let expected = point_sum(&[
            (w0.powi(3), get(&left, 6 - i, 3, i - 3)),
            (3.0 * w0 * w0 * w1, get(&left, 6 - i, 2, i - 2)),
            (3.0 * w0 * w1 * w1, get(&left, 6 - i, 1, i - 1)),
            (w1.powi(3), get(&left, 6 - i, 0, i)),
            (3.0 * w0 * w0 * w2, get(&left, 6 - i + 1, 2, i - 3)),
            (6.0 * w0 * w1 * w2, get(&left, 6 - i + 1, 1, i - 2)),
            (3.0 * w1 * w1 * w2, get(&left, 6 - i + 1, 0, i - 1)),
            (3.0 * w0 * w2 * w2, get(&left, 6 - i + 2, 1, i - 3)),
            (3.0 * w1 * w2 * w2, get(&left, 6 - i + 2, 0, i - 2)),
            (w2.powi(3), get(&left, 6 - i + 3, 0, i - 3)),
        ]);
        assert_close(get(&right, 6 - i, i - 3, 3), expected);
    }
}
