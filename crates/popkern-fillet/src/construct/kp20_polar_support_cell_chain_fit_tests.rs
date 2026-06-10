use popkern_core::geometry::Point3;

use super::{eval_row, fit_constrained_degree9};

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

#[test]
fn constrained_fit_recovers_polynomial_and_pins_nodes() {
    let source = quadratic_arc_row();
    let mut samples = Vec::new();
    for index in 0..40 {
        let t = index as f64 / 39.0;
        samples.push((t, eval_row(&source, t)));
    }
    let constraints = vec![
        (0.0, eval_row(&source, 0.0)),
        (0.25, eval_row(&source, 0.25)),
        (0.5, eval_row(&source, 0.5)),
        (1.0, eval_row(&source, 1.0)),
    ];
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
    for &(t, position) in &constraints {
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
    let node = (0.5, Point3::new(0.5, 1.0e-3, 0.0));
    let fitted = fit_constrained_degree9(&samples, &[node]).expect("fit must solve");
    let at_node = eval_row(&fitted, 0.5);
    assert!(
        node.1.distance_to(&at_node) < 1.0e-10,
        "pinned node must be interpolated exactly, got {at_node:?}"
    );
    let at_end = eval_row(&fitted, 0.0);
    assert!(
        at_end.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 1.0e-3,
        "ends should stay near the sample line"
    );
}

#[test]
fn fit_declines_without_enough_samples() {
    let samples = vec![(0.0, Point3::new(0.0, 0.0, 0.0))];
    assert_eq!(
        fit_constrained_degree9(&samples, &[]),
        Err("sub_chain_sample_count_insufficient")
    );
}
