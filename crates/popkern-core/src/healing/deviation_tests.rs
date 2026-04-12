use super::*;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};

fn xy_plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

fn x_axis_line() -> CurveGeometry {
    CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    }
}

fn y_axis_line() -> CurveGeometry {
    CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(0.0, 1.0, 0.0),
    }
}

#[test]
fn identical_curves_zero_deviation() {
    let line = x_axis_line();
    let range = (0.0, 10.0);

    let result = curve_curve_deviation(&line, range, &line, range, 50);

    assert!(
        result.max_deviation < 1e-12,
        "identical curves should have ~0 deviation, got {}",
        result.max_deviation
    );
    assert!(
        result.mean_deviation < 1e-12,
        "mean deviation should be ~0 for identical curves, got {}",
        result.mean_deviation
    );
    // 50 uniform + 8 adaptive refinement samples.
    assert!(
        result.sample_count >= 50,
        "sample count should include refinement, got {}",
        result.sample_count
    );
}

#[test]
fn curve_on_plane_zero_deviation() {
    let line = x_axis_line();
    let plane = xy_plane();

    let result = edge_surface_deviation(&line, (0.0, 5.0), &plane, 30);

    assert!(
        result.max_deviation < 1e-6,
        "line on XY plane should have ~0 deviation, got {}",
        result.max_deviation
    );
    assert!(
        result.mean_deviation < 1e-6,
        "mean deviation should be ~0 for line on plane, got {}",
        result.mean_deviation
    );
}

#[test]
fn parallel_curves_known_deviation() {
    let line_a = x_axis_line();
    // Line at z=1, parallel to x_axis_line.
    let line_b = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 1.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    let range = (0.0, 5.0);

    let result = curve_curve_deviation(&line_a, range, &line_b, range, 20);

    assert!(
        (result.max_deviation - 1.0).abs() < 1e-10,
        "parallel lines 1 unit apart should have max deviation ~1.0, got {}",
        result.max_deviation
    );
    assert!(
        (result.mean_deviation - 1.0).abs() < 1e-10,
        "mean deviation should be ~1.0, got {}",
        result.mean_deviation
    );
}

#[test]
fn surface_surface_same_plane() {
    let plane = xy_plane();
    let range = ((-1.0, 1.0), (-1.0, 1.0));

    let result = surface_surface_deviation(&plane, range, &plane, range, 25);

    assert!(
        result.max_deviation < 1e-6,
        "same plane should have ~0 deviation, got {}",
        result.max_deviation
    );
    assert!(result.sample_count > 0);
}

#[test]
fn deviation_result_fields() {
    let line_a = x_axis_line();
    let line_b = y_axis_line();

    let result = curve_curve_deviation(&line_a, (0.0, 1.0), &line_b, (0.0, 1.0), 10);

    // At t=0 both are at origin, so min deviation is 0.
    // At t=1 they are sqrt(2) apart.
    assert!(result.max_deviation > 0.0);
    // 10 uniform + 8 adaptive refinement samples.
    assert!(
        result.sample_count >= 10,
        "sample count should include refinement, got {}",
        result.sample_count
    );
}

#[test]
fn bounded_deviation_line_on_plane_zero() {
    let line = x_axis_line();
    let plane = xy_plane();

    // Use a small range and generous budget so subdivision reaches min_interval.
    let result = edge_surface_deviation_bounded(&line, (0.0, 1.0), &plane, 0.001, 5000);

    // Sampled deviation must be ~0.
    assert!(
        result.result.max_deviation < 1e-6,
        "sampled max deviation should be ~0, got {}",
        result.result.max_deviation
    );
    // The upper bound includes Lipschitz headroom proportional to min_interval.
    // For speed=1.0 and half_interval=0.0005, expect bound ~ 0.0005.
    assert!(
        result.upper_bound < 0.01,
        "bounded deviation should be small for line on plane, got {}",
        result.upper_bound
    );
    assert!(
        result.intervals_processed > 0,
        "should have processed at least one interval"
    );
}

#[test]
fn bounded_deviation_offset_line_known_bound() {
    // Line at z=0.5 over XY plane: deviation should be 0.5 everywhere.
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.5),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    let plane = xy_plane();

    let result = edge_surface_deviation_bounded(&line, (0.0, 5.0), &plane, 0.01, 1000);

    // Upper bound must be >= true max deviation (0.5).
    assert!(
        result.upper_bound >= 0.5 - 1e-6,
        "upper bound should be >= 0.5, got {}",
        result.upper_bound
    );
    // And not wildly above (since the curve is flat, the bound should be tight).
    assert!(
        result.upper_bound < 0.6,
        "upper bound should be tight for a flat curve, got {}",
        result.upper_bound
    );
    assert!(
        (result.result.max_deviation - 0.5).abs() < 1e-6,
        "sampled max should be ~0.5, got {}",
        result.result.max_deviation
    );
}

#[test]
fn bounded_deviation_min_interval_clamp() {
    let line = x_axis_line();
    let plane = xy_plane();

    // Very large min_interval should terminate quickly.
    let result = edge_surface_deviation_bounded(&line, (0.0, 1.0), &plane, 10.0, 100);

    // The full range (1.0) is already below min_interval (10.0), so should process just 1 interval.
    assert!(
        result.intervals_processed <= 2,
        "large min_interval should stop early, got {} intervals",
        result.intervals_processed
    );
}

#[test]
fn bounded_deviation_max_subdivisions_respected() {
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.5),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    let plane = xy_plane();

    // Very small min_interval but low max_subdivisions.
    let result = edge_surface_deviation_bounded(&line, (0.0, 10.0), &plane, 1e-10, 5);

    assert!(
        result.intervals_processed <= 5,
        "should respect max_subdivisions, got {} intervals",
        result.intervals_processed
    );
}
