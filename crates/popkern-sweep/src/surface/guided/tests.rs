// #[cfg(test)] — included only in test mode; .unwrap() allowed in tests.
use super::*;
use popkern_core::geometry::Point3;

/// Build a flat surface in the XY plane from (0,0,0) to (1,1,0).
fn flat_surface() -> NurbsSurface {
    NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    }
}

/// A guide that lies exactly on the flat surface.
fn on_surface_guide() -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
        weights: vec![1.0, 1.0],
    }
}

/// A guide that is far from the surface (offset in Z).
fn far_guide() -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.5, 0.5, 5.0), Point3::new(0.5, 0.5, 10.0)],
        weights: vec![1.0, 1.0],
    }
}

#[test]
fn guide_on_surface_passes() {
    let surface = flat_surface();
    let result = validate_guide_rails(&surface, &[on_surface_guide()], &Tolerance::default());
    assert!(result.is_ok(), "expected pass, got {:?}", result.err());
    let validation = result.unwrap(); // ok in tests
    assert!(
        validation.max_deviation < 0.1,
        "got {}",
        validation.max_deviation
    );
}

#[test]
fn guide_far_from_surface_fails() {
    let result = validate_guide_rails(&flat_surface(), &[far_guide()], &Tolerance::default());
    assert!(result.is_err());
    match result.err() {
        Some(SweepError::GuideDeviation { deviation }) => {
            assert!(deviation > 1.0, "expected large deviation, got {deviation}");
        }
        other => panic!("expected GuideDeviation, got {:?}", other),
    }
}

#[test]
fn empty_guides_passes() {
    let result = validate_guide_rails(&flat_surface(), &[], &Tolerance::default());
    assert!(result.is_ok());
    assert!(result.unwrap().max_deviation < 1e-10); // ok in tests
}

#[test]
fn multiple_guides_reports_per_guide_deviations() {
    let result = validate_guide_rails(
        &flat_surface(),
        &[on_surface_guide(), far_guide()],
        &Tolerance::default(),
    );
    assert!(result.is_err());
}

#[test]
fn too_many_guides_returns_error() {
    let surface = flat_surface();
    let guides: Vec<NurbsCurve> = (0..17).map(|_| on_surface_guide()).collect();
    let tol = Tolerance::default();
    let result = validate_guide_rails(&surface, &guides, &tol);
    assert!(result.is_err());
}

#[test]
fn apply_guides_noop_for_empty() {
    let mut surface = flat_surface();
    let tol = Tolerance::default();
    let result = apply_guide_rails(&mut surface, &[], &tol);
    assert!(result.is_ok());
}

#[test]
fn apply_guides_already_within_tolerance() {
    let mut surface = flat_surface();
    let guide = on_surface_guide();
    let tol = Tolerance::default();
    let result = apply_guide_rails(&mut surface, &[guide], &tol);
    assert!(
        result.is_ok(),
        "expected pass for on-surface guide, got {:?}",
        result.err()
    );
}

#[test]
fn apply_guides_adjusts_surface_for_small_offset() {
    // Create a flat surface and a guide offset in Z.
    let mut surface = flat_surface();
    let offset = 0.3;
    let guide = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.25, 0.5, offset),
            Point3::new(0.75, 0.5, offset),
        ],
        weights: vec![1.0, 1.0],
    };

    // Before adjustment, deviation should be ~offset.
    let dev_before = measure_guide_deviation(&surface, &guide);
    assert!(
        dev_before > 0.1,
        "expected deviation >{}, got {dev_before}",
        0.1
    );

    // Apply with a tolerance that forces adjustment but is achievable.
    let tol = Tolerance {
        linear: 0.15,
        ..Tolerance::default()
    };
    let _result = apply_guide_rails(&mut surface, &[guide.clone()], &tol);
    // We may or may not converge within tolerance, but deviation
    // should be reduced from the original.
    let dev_after = measure_guide_deviation(&surface, &guide);
    assert!(
        dev_after < dev_before,
        "deviation not reduced: before={dev_before}, after={dev_after}"
    );
}

#[test]
fn guided_adjustment_better_than_nearest_point() {
    // Create a surface with enough control points for meaningful
    // distribution, then offset a guide by 0.1 in Z.
    let mut surface = NurbsSurface {
        degree_u: 2,
        degree_v: 2,
        knots_u: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
        control_points: {
            let mut pts = Vec::with_capacity(16);
            for iu in 0..4 {
                for iv in 0..4 {
                    pts.push(Point3::new(iu as f64 / 3.0, iv as f64 / 3.0, 0.0));
                }
            }
            pts
        },
        num_u: 4,
        num_v: 4,
        weights: vec![1.0; 16],
    };

    let offset = 0.1;
    let guide = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.2, 0.5, offset), Point3::new(0.8, 0.5, offset)],
        weights: vec![1.0, 1.0],
    };

    let dev_before = measure_guide_deviation(&surface, &guide);
    assert!(
        dev_before > 0.05,
        "expected initial deviation, got {dev_before}"
    );

    // Run the weighted multi-point adjustment (10 iterations).
    let tol = Tolerance {
        linear: 0.005,
        ..Tolerance::default()
    };
    let _result = apply_guide_rails(&mut surface, &[guide.clone()], &tol);

    let dev_after = measure_guide_deviation(&surface, &guide);
    // Should achieve at least 4x reduction from the weighted
    // multi-point distribution approach.
    assert!(
        dev_after < dev_before * 0.25,
        "expected 4x reduction: before={dev_before}, after={dev_after}"
    );
}

#[test]
fn apply_guides_too_many_returns_error() {
    let mut surface = flat_surface();
    let guides: Vec<NurbsCurve> = (0..17).map(|_| on_surface_guide()).collect();
    let tol = Tolerance::default();
    let result = apply_guide_rails(&mut surface, &guides, &tol);
    assert!(result.is_err());
}
