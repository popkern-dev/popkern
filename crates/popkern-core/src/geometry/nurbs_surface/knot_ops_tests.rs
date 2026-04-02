use super::NurbsSurface;
use crate::geometry::Point3;
use crate::geometry::surface::SurfaceEval;

/// Build a biquadratic (degree 2x2) surface with a known curved shape.
/// The surface has a "bump" in the center via the middle control point.
fn biquadratic_bump() -> NurbsSurface {
    NurbsSurface {
        degree_u: 2,
        degree_v: 2,
        knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            // Row u=0
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.5, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            // Row u=0.5
            Point3::new(0.5, 0.0, 0.0),
            Point3::new(0.5, 0.5, 1.0), // bump in Z
            Point3::new(0.5, 1.0, 0.0),
            // Row u=1
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 0.5, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 3,
        num_v: 3,
        weights: vec![1.0; 9],
    }
}

/// Build a bilinear patch (degree 1x1) for simple tests.
fn bilinear_patch() -> NurbsSurface {
    NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0, 1.0, 1.0, 1.0],
    }
}

/// Helper: sample the surface at a grid of (u,v) points and return positions.
fn sample_surface(surface: &NurbsSurface, n: usize) -> Vec<Point3> {
    let mut pts = Vec::with_capacity(n * n);
    for i in 0..n {
        let u = i as f64 / (n - 1).max(1) as f64;
        for j in 0..n {
            let v = j as f64 / (n - 1).max(1) as f64;
            pts.push(surface.point_at(u, v));
        }
    }
    pts
}

/// Insert a knot in u then remove it -- shape should be preserved exactly.
#[test]
fn remove_knot_u_roundtrip_biquadratic() {
    let surface = biquadratic_bump();
    let before = sample_surface(&surface, 5);

    let inserted = surface.insert_knot_u(0.5).unwrap();
    assert_eq!(inserted.num_u, surface.num_u + 1);

    // knots_u after insert: [0,0,0, 0.5, 1,1,1] -- index 3 is the new knot.
    let removed = inserted
        .remove_knot_u(3, 0.0)
        .expect("should be able to remove freshly inserted knot");
    assert_eq!(removed.num_u, surface.num_u);

    let after = sample_surface(&removed, 5);
    for (p1, p2) in before.iter().zip(after.iter()) {
        let dist = p1.distance_to(p2);
        assert!(
            dist < 1e-9,
            "shape should be preserved after insert+remove: dist={dist}"
        );
    }
}

/// Insert a knot in v then remove it -- shape should be preserved exactly.
#[test]
fn remove_knot_v_roundtrip_biquadratic() {
    let surface = biquadratic_bump();
    let before = sample_surface(&surface, 5);

    let inserted = surface.insert_knot_v(0.5).unwrap();
    assert_eq!(inserted.num_v, surface.num_v + 1);

    // knots_v after insert: [0,0,0, 0.5, 1,1,1] -- index 3.
    let removed = inserted
        .remove_knot_v(3, 0.0)
        .expect("should be able to remove freshly inserted knot");
    assert_eq!(removed.num_v, surface.num_v);

    let after = sample_surface(&removed, 5);
    for (p1, p2) in before.iter().zip(after.iter()) {
        let dist = p1.distance_to(p2);
        assert!(
            dist < 1e-9,
            "shape should be preserved after insert+remove v: dist={dist}"
        );
    }
}

/// Attempting to remove a knot that would change the shape beyond tolerance
/// should return None.
#[test]
fn remove_knot_u_fails_when_exceeds_tolerance() {
    // Start with the biquadratic bump (3x3 control points, no interior knots).
    // Insert a knot at 0.5, then manually perturb a control point so that
    // the reverse removal with zero tolerance fails.
    let surface = biquadratic_bump();
    let mut inserted = surface.insert_knot_u(0.5).unwrap();
    // Perturb one of the new control points significantly.
    inserted.control_points[1 * inserted.num_v + 1].z += 5.0;

    // With zero tolerance, removal should fail because the geometry was altered.
    let result = inserted.remove_knot_u(3, 0.0);
    assert!(
        result.is_none(),
        "removal should fail when control points have been perturbed"
    );
}

/// Attempting to remove a knot that would change the shape beyond tolerance
/// in the v direction should return None.
#[test]
fn remove_knot_v_fails_when_exceeds_tolerance() {
    let surface = biquadratic_bump();
    let mut inserted = surface.insert_knot_v(0.5).unwrap();
    // Perturb a control point.
    inserted.control_points[1 * inserted.num_v + 1].z += 5.0;

    let result = inserted.remove_knot_v(3, 0.0);
    assert!(
        result.is_none(),
        "removal should fail when control points have been perturbed"
    );
}

/// With a large enough tolerance, removal should succeed even on perturbed data.
#[test]
fn remove_knot_u_succeeds_with_large_tolerance() {
    let surface = biquadratic_bump();
    let mut inserted = surface.insert_knot_u(0.5).unwrap();
    inserted.control_points[1 * inserted.num_v + 1].z += 0.5;

    // With a generous tolerance, removal should succeed.
    let result = inserted.remove_knot_u(3, 10.0);
    assert!(
        result.is_some(),
        "removal should succeed with large tolerance"
    );
}

/// remove_knot_u_times reduces multiplicity: insert a knot twice, then
/// remove up to 3 times. Should remove exactly 2.
#[test]
fn remove_knot_u_times_reduces_multiplicity() {
    let surface = biquadratic_bump();
    let inserted = surface.insert_knots_u(&[0.5, 0.5]).unwrap();
    assert_eq!(inserted.num_u, surface.num_u + 2);

    let before = sample_surface(&surface, 5);

    // knots_u = [0,0,0, 0.5,0.5, 1,1,1]. Index 3 is the first 0.5.
    let (result, count) = inserted.remove_knot_u_times(3, 3, 0.0);
    assert_eq!(count, 2, "should remove both inserted knots");
    assert_eq!(result.num_u, surface.num_u);

    let after = sample_surface(&result, 5);
    for (p1, p2) in before.iter().zip(after.iter()) {
        let dist = p1.distance_to(p2);
        assert!(
            dist < 1e-8,
            "shape should be preserved after removing multiplicity"
        );
    }
}

/// remove_knot_v_times reduces multiplicity: insert twice, remove up to 3.
#[test]
fn remove_knot_v_times_reduces_multiplicity() {
    let surface = biquadratic_bump();
    let inserted = surface.insert_knots_v(&[0.5, 0.5]).unwrap();
    assert_eq!(inserted.num_v, surface.num_v + 2);

    let before = sample_surface(&surface, 5);

    let (result, count) = inserted.remove_knot_v_times(3, 3, 0.0);
    assert_eq!(count, 2, "should remove both inserted knots");
    assert_eq!(result.num_v, surface.num_v);

    let after = sample_surface(&result, 5);
    for (p1, p2) in before.iter().zip(after.iter()) {
        let dist = p1.distance_to(p2);
        assert!(
            dist < 1e-8,
            "shape should be preserved after removing v-multiplicity"
        );
    }
}

/// Boundary knots (clamped ends) should not be removable.
#[test]
fn remove_boundary_knots_u_returns_none() {
    let surface = biquadratic_bump();
    // Indices 0, 1, 2 are the left boundary (knot value 0.0).
    // Indices 3, 4, 5 are the right boundary (knot value 1.0).
    for idx in [0, 1, 2, 3, 4, 5] {
        let result = surface.remove_knot_u(idx, 1.0);
        assert!(
            result.is_none(),
            "boundary knot at index {idx} should not be removable"
        );
    }
}

/// Boundary knots (clamped ends) should not be removable in v direction.
#[test]
fn remove_boundary_knots_v_returns_none() {
    let surface = biquadratic_bump();
    for idx in [0, 1, 2, 3, 4, 5] {
        let result = surface.remove_knot_v(idx, 1.0);
        assert!(
            result.is_none(),
            "boundary v-knot at index {idx} should not be removable"
        );
    }
}

/// remove_knot_u_times with zero `times` returns the surface unchanged.
#[test]
fn remove_knot_u_times_zero_is_noop() {
    let surface = biquadratic_bump();
    let inserted = surface.insert_knot_u(0.5).unwrap();
    let (result, count) = inserted.remove_knot_u_times(3, 0, 0.0);
    assert_eq!(count, 0);
    assert_eq!(result.num_u, inserted.num_u);
    assert_eq!(result.knots_u.len(), inserted.knots_u.len());
}

/// remove_knot_v_times with zero `times` returns the surface unchanged.
#[test]
fn remove_knot_v_times_zero_is_noop() {
    let surface = biquadratic_bump();
    let inserted = surface.insert_knot_v(0.5).unwrap();
    let (result, count) = inserted.remove_knot_v_times(3, 0, 0.0);
    assert_eq!(count, 0);
    assert_eq!(result.num_v, inserted.num_v);
    assert_eq!(result.knots_v.len(), inserted.knots_v.len());
}

/// Verify that surface shape is preserved after knot removal by sampling
/// a dense grid of points before and after.
#[test]
fn shape_preserved_after_removal_dense_sampling() {
    let surface = bilinear_patch();
    let inserted = surface.insert_knot_u(0.3).unwrap();
    // knots_u after insert: [0,0, 0.3, 1,1]. Index 2 is the new knot.
    let removed = inserted
        .remove_knot_u(2, 0.0)
        .expect("should remove freshly inserted knot from bilinear patch");

    // Dense 10x10 grid sampling.
    let n = 10;
    for i in 0..n {
        let u = i as f64 / (n - 1) as f64;
        for j in 0..n {
            let v = j as f64 / (n - 1) as f64;
            let p_orig = surface.point_at(u, v);
            let p_removed = removed.point_at(u, v);
            let dist = p_orig.distance_to(&p_removed);
            assert!(dist < 1e-9, "shape deviation at u={u:.2},v={v:.2}: {dist}");
        }
    }
}

/// remove_knot_u_times on a surface with no interior knots should return
/// count=0 and leave the surface unchanged.
#[test]
fn remove_knot_u_times_no_interior_knots() {
    let surface = biquadratic_bump();
    // Index 0 is a boundary knot, removal should fail immediately.
    let (result, count) = surface.remove_knot_u_times(0, 5, 1.0);
    assert_eq!(count, 0);
    assert_eq!(result.num_u, surface.num_u);
}

/// remove_knot_v_times on a surface with no interior knots should return
/// count=0 and leave the surface unchanged.
#[test]
fn remove_knot_v_times_no_interior_knots() {
    let surface = biquadratic_bump();
    let (result, count) = surface.remove_knot_v_times(0, 5, 1.0);
    assert_eq!(count, 0);
    assert_eq!(result.num_v, surface.num_v);
}

/// Insert a single knot and attempt to remove it many times: should succeed
/// once and then stop (index out of range or no more removable knots).
#[test]
fn remove_knot_u_times_stops_after_single_removal() {
    let surface = bilinear_patch();
    let inserted = surface.insert_knot_u(0.5).unwrap();
    // knots_u = [0,0, 0.5, 1,1]. Index 2 is the inserted knot.
    let (result, count) = inserted.remove_knot_u_times(2, 10, 0.0);
    assert_eq!(count, 1, "only one instance to remove");
    assert_eq!(result.num_u, surface.num_u);
}

/// Same test for v direction.
#[test]
fn remove_knot_v_times_stops_after_single_removal() {
    let surface = bilinear_patch();
    let inserted = surface.insert_knot_v(0.5).unwrap();
    let (result, count) = inserted.remove_knot_v_times(2, 10, 0.0);
    assert_eq!(count, 1, "only one instance to remove");
    assert_eq!(result.num_v, surface.num_v);
}

/// Perturbed control point removal with intermediate tolerance: removal
/// succeeds but produces a slightly different surface.
#[test]
fn remove_knot_u_with_tolerance_allows_approximation() {
    let surface = biquadratic_bump();
    let mut inserted = surface.insert_knot_u(0.5).unwrap();
    // Small perturbation -- within a reasonable tolerance.
    inserted.control_points[1 * inserted.num_v + 0].z += 0.01;

    // Zero tolerance should fail.
    assert!(inserted.remove_knot_u(3, 0.0).is_none());
    // Generous tolerance should succeed.
    let result = inserted.remove_knot_u(3, 1.0);
    assert!(result.is_some(), "tolerance=1.0 should allow removal");
}
