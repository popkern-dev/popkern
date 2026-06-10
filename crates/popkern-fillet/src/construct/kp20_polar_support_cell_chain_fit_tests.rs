use popkern_core::geometry::{Point3, Vec3};

use super::super::chain_solver::{ChainCurveConstraints, dbernstein9, fit_constrained_degree9};
use super::eval_row;

fn quadratic_arc_row() -> [Point3; 10] {
    // Degree-9 Bezier reproducing (t, t², 0.3·t) exactly via the
    // Bezier form of the monomials: x linear, y = t² has
    // coefficients i(i−1)/72 at degree 9.
    let mut row = [Point3::new(0.0, 0.0, 0.0); 10];
    for (i, point) in row.iter_mut().enumerate() {
        let t = i as f64 / 9.0;
        let y = (i as f64) * (i as f64 - 1.0) / 72.0;
        *point = Point3::new(t, y, 0.3 * t);
    }
    row
}

fn derivative_at(row: &[Point3; 10], t: f64) -> Vec3 {
    let basis = dbernstein9(t);
    let mut out = Vec3::new(0.0, 0.0, 0.0);
    for (index, weight) in basis.iter().enumerate() {
        out = Vec3::new(
            out.x + row[index].x * weight,
            out.y + row[index].y * weight,
            out.z + row[index].z * weight,
        );
    }
    out
}

#[test]
fn constrained_fit_recovers_polynomial_and_pins_nodes() {
    let source = quadratic_arc_row();
    let mut samples = Vec::new();
    for index in 0..40 {
        let t = index as f64 / 39.0;
        samples.push((t, eval_row(&source, t)));
    }
    let constraints = ChainCurveConstraints {
        positions: vec![
            (0.0, eval_row(&source, 0.0)),
            (0.25, eval_row(&source, 0.25)),
            (0.5, eval_row(&source, 0.5)),
            (1.0, eval_row(&source, 1.0)),
        ],
        tangent_planes: Vec::new(),
    };
    let fitted = fit_constrained_degree9(&samples, &constraints).expect("fit must solve");
    for index in 0..=20 {
        let t = index as f64 / 20.0;
        let expected = eval_row(&source, t);
        let actual = eval_row(&fitted, t);
        assert!(
            expected.distance_to(&actual) < 1.0e-8,
            "fit diverged at t={t}: {expected:?} vs {actual:?}"
        );
    }
    for &(t, position) in &constraints.positions {
        let actual = eval_row(&fitted, t);
        assert!(
            position.distance_to(&actual) < 1.0e-9,
            "node not interpolated at t={t}"
        );
    }
}

#[test]
fn constrained_fit_moves_to_pinned_nodes_off_the_samples() {
    // Samples trace one line; one pinned node sits off the line.
    // The fit must pass through the node exactly while staying
    // near the samples elsewhere.
    let mut samples = Vec::new();
    for index in 0..40 {
        let t = index as f64 / 39.0;
        samples.push((t, Point3::new(t, 0.0, 0.0)));
    }
    let constraints = ChainCurveConstraints {
        positions: vec![(0.5, Point3::new(0.5, 1.0e-3, 0.0))],
        tangent_planes: Vec::new(),
    };
    let fitted = fit_constrained_degree9(&samples, &constraints).expect("fit must solve");
    let at_node = eval_row(&fitted, 0.5);
    assert!(
        constraints.positions[0].1.distance_to(&at_node) < 1.0e-10,
        "pinned node must be interpolated exactly, got {at_node:?}"
    );
    let at_end = eval_row(&fitted, 0.0);
    assert!(
        at_end.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 1.0e-3,
        "ends should stay near the sample line"
    );
}

#[test]
fn tangent_plane_constraint_holds_exactly() {
    // Samples trace a parabola whose end tangent has a y
    // component; constrain the tangent at t=1 to lie in the
    // plane perpendicular to Y. The fitted curve must satisfy
    // C'(1)·Y = 0 exactly while still passing near the samples.
    let source = quadratic_arc_row();
    let mut samples = Vec::new();
    for index in 0..40 {
        let t = index as f64 / 39.0;
        samples.push((t, eval_row(&source, t)));
    }
    let normal = Vec3::new(0.0, 1.0, 0.0);
    let constraints = ChainCurveConstraints {
        positions: vec![(0.0, eval_row(&source, 0.0)), (1.0, eval_row(&source, 1.0))],
        tangent_planes: vec![(1.0, normal)],
    };
    let fitted = fit_constrained_degree9(&samples, &constraints).expect("fit must solve");
    let tangent_dot = derivative_at(&fitted, 1.0).dot(&normal);
    assert!(
        tangent_dot.abs() < 1.0e-9,
        "end tangent must lie in the constrained plane, got {tangent_dot}"
    );
    let source_dot = derivative_at(&source, 1.0).dot(&normal);
    assert!(
        source_dot.abs() > 0.1,
        "fixture's unconstrained tangent must violate the plane"
    );
}

#[test]
fn fit_declines_without_enough_samples() {
    let samples = vec![(0.0, Point3::new(0.0, 0.0, 0.0))];
    assert!(matches!(
        fit_constrained_degree9(&samples, &ChainCurveConstraints::default()),
        Err("sub_chain_sample_count_insufficient")
    ));
}
