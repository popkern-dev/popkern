//! The principal-branch landing is a rigid translation (PK-681).

use std::f64::consts::TAU;

use popkern_core::geometry::Point3;
use popkern_core::geometry::nurbs_curve::NurbsCurve;

use super::translate_onto_principal_branch;

/// A degree-1 query through the given u values at v = 0.4.
fn query(us: &[f64]) -> NurbsCurve {
    let n = us.len();
    let mut knots = vec![0.0; 2];
    knots.extend((1..n - 1).map(|i| i as f64 / (n - 1) as f64));
    knots.extend([1.0; 2]);
    NurbsCurve {
        degree: 1,
        knots,
        control_points: us.iter().map(|&u| Point3::new(u, 0.4, 0.0)).collect(),
        weights: vec![1.0; n],
    }
}

fn us(curve: &NurbsCurve) -> Vec<f64> {
    curve.control_points.iter().map(|p| p.x).collect()
}

#[test]
fn a_full_period_circle_on_the_next_branch_keeps_its_seam_end() {
    // The M2 nut mouth-chamfer circle on the crest band: one full
    // period at u in [2pi, 4pi]. A per-point fold sent its closing
    // control point from 2pi to 0 and the circle ran backwards.
    let circle = query(&[TAU, TAU + 1.5, TAU + 3.0, TAU + 4.5, 2.0 * TAU]);
    let landed = us(&translate_onto_principal_branch(&circle, TAU));
    let expected = [0.0, 1.5, 3.0, 4.5, TAU];
    for (got, want) in landed.iter().zip(expected) {
        assert!((got - want).abs() < 1e-12, "{landed:?}");
    }
}

#[test]
fn a_query_on_the_principal_branch_stays_put() {
    let q = query(&[0.2, 1.0, 3.0]);
    assert_eq!(us(&translate_onto_principal_branch(&q, TAU)), us(&q));
    // Mostly on the principal branch with a dip below zero: stays
    // (lowest-point landing moved it a whole period off the boundary
    // and broke the filleted-operand subtraction canary).
    for q in [query(&[-1e-12, 3.0, TAU]), query(&[-0.01, 1.0, 3.0])] {
        assert_eq!(us(&translate_onto_principal_branch(&q, TAU)), us(&q));
    }
}

#[test]
fn a_seam_straddling_query_moves_rigidly() {
    // Straddles 2pi on the second branch: the whole curve moves one
    // period down and stays continuous across the seam.
    let q = query(&[TAU + 5.4, TAU + 6.0, TAU + 7.0]);
    let landed = us(&translate_onto_principal_branch(&q, TAU));
    let expected = [5.4, 6.0, 7.0];
    for (got, want) in landed.iter().zip(expected) {
        assert!((got - want).abs() < 1e-12, "{landed:?}");
    }
    // Negative branches land too (the fold could not reach them).
    let q = query(&[-2.0 * TAU + 1.0, -2.0 * TAU + 2.0]);
    let landed = us(&translate_onto_principal_branch(&q, TAU));
    assert!((landed[0] - 1.0).abs() < 1e-12 && (landed[1] - 2.0).abs() < 1e-12);
}
