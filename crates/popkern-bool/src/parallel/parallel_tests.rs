//! Tests for parallel intersection dispatch (PK-277).
//!
//! The determinism contract: `find_intersections` output is positionally
//! identical at any thread count, because the task list is a pure
//! function of the two solids and results are collected in task order.

use popkern_core::geometry::CurveEval;
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::intersect::{IntersectionRecord, find_intersections};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn overlapping_cubes() -> (BRep, BRep) {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.3, 0.2))
        .unwrap();
    (a, b)
}

/// Bitwise comparison of two record sequences: same length, same face
/// pairs in the same positions, same tangency flags, and bit-identical
/// 3D curve samples.
fn assert_records_identical(lhs: &[IntersectionRecord], rhs: &[IntersectionRecord]) {
    assert_eq!(lhs.len(), rhs.len(), "record count differs");
    for (i, (l, r)) in lhs.iter().zip(rhs.iter()).enumerate() {
        assert_eq!(l.face_a, r.face_a, "record {i}: face_a differs");
        assert_eq!(l.face_b, r.face_b, "record {i}: face_b differs");
        assert_eq!(l.is_tangent, r.is_tangent, "record {i}: tangency differs");
        let (lt0, lt1) = l.curve.curve_3d.parameter_range();
        let (rt0, rt1) = r.curve.curve_3d.parameter_range();
        assert_eq!(lt0.to_bits(), rt0.to_bits(), "record {i}: t0 differs");
        assert_eq!(lt1.to_bits(), rt1.to_bits(), "record {i}: t1 differs");
        for k in 0..=8 {
            let t = lt0 + (lt1 - lt0) * (k as f64) / 8.0;
            let lp = l.curve.curve_3d.point_at(t);
            let rp = r.curve.curve_3d.point_at(t);
            assert_eq!(
                (lp.x.to_bits(), lp.y.to_bits(), lp.z.to_bits()),
                (rp.x.to_bits(), rp.y.to_bits(), rp.z.to_bits()),
                "record {i}: curve sample {k} differs"
            );
        }
    }
}

#[test]
fn overlapping_cubes_produce_records() {
    let (a, b) = overlapping_cubes();
    let result = find_intersections(&a, &b, &tol()).unwrap();
    assert!(
        !result.records.is_empty(),
        "overlapping cubes should produce intersection curves"
    );
}

#[test]
fn disjoint_cubes_produce_no_records() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(5.0, 0.0, 0.0))
        .unwrap();
    let result = find_intersections(&a, &b, &tol()).unwrap();
    assert!(result.records.is_empty());
}

/// PK-277 acceptance: identical results 1-thread vs N-thread. Runs the
/// same intersection under two in-process rayon pools and compares the
/// record sequences bitwise. With the `parallel` feature disabled this
/// still passes (both runs take the sequential path) but pins nothing
/// new.
#[test]
fn thread_count_does_not_change_records() {
    let (a, b) = overlapping_cubes();

    let single = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(|| find_intersections(&a, &b, &tol()).unwrap());
    let many = rayon::ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .unwrap()
        .install(|| find_intersections(&a, &b, &tol()).unwrap());

    assert_records_identical(&single.records, &many.records);
}
