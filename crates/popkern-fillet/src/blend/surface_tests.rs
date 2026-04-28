//! Tests for blend surface construction.

use super::*;
use crate::blend::boundary::{BlendBoundary, BoundarySegment};
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, VertexHandle};

fn test_tolerance() -> Tolerance {
    Tolerance::default()
        .with_linear(1e-3)
        .with_angular(1e-3)
        .with_parametric(1e-6)
}

fn make_test_curve(start: Point3, end: Point3) -> NurbsCurve {
    crate::center_curve::fit_nurbs_through_points(&[start, end]).expect("fit should succeed")
}

fn make_blend_boundary(vertex: VertexHandle) -> BlendBoundary {
    // Two fillet-end arcs forming a right-angle-ish boundary.
    let seg1 = BoundarySegment::FilletEnd {
        fillet_index: 0,
        curve: make_test_curve(Point3::new(0.1, 0.0, 0.0), Point3::new(0.0, 0.1, 0.0)),
    };
    let seg2 = BoundarySegment::FilletEnd {
        fillet_index: 1,
        curve: make_test_curve(Point3::new(0.0, 0.1, 0.0), Point3::new(0.0, 0.0, 0.1)),
    };

    BlendBoundary {
        vertex,
        vertex_position: Point3::new(0.0, 0.0, 0.0),
        segments: vec![seg1, seg2],
        bridge_fallbacks: Vec::new(),
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
    }
}

#[test]
fn mitred_blend_two_fillets() {
    let tol = test_tolerance();
    let brep = BRep::unit_cube(tol).expect("cube should build");
    let vertex = brep.vertices.keys().next().expect("cube has vertices");

    let boundary = make_blend_boundary(vertex);
    let fillets: Vec<FilletSurfaceResult> = Vec::new();

    let result = build_blend_surface(
        &boundary,
        &fillets,
        &brep,
        VertexMode::Mitred,
        crate::vertex_classifier::VertexBlendClassification::SameSideConvex,
        &tol,
    );
    assert!(result.is_ok(), "mitred blend failed: {:?}", result.err());

    let blend = result.expect("already checked");
    assert!(
        matches!(blend.surface, SurfaceGeometry::Plane { .. }),
        "mitred blend should produce a plane, got {:?}",
        std::mem::discriminant(&blend.surface)
    );
    assert!(
        !blend.boundary_curves.is_empty(),
        "should have boundary curves"
    );
}

#[test]
fn smooth_blend_two_fillets() {
    let tol = test_tolerance();
    let brep = BRep::unit_cube(tol).expect("cube should build");
    let vertex = brep.vertices.keys().next().expect("cube has vertices");

    // Build boundary with points on a sphere of radius 0.1.
    let r = 0.1;
    let seg1 = BoundarySegment::FilletEnd {
        fillet_index: 0,
        curve: make_test_curve(Point3::new(r, 0.0, 0.0), Point3::new(0.0, r, 0.0)),
    };
    let seg2 = BoundarySegment::FilletEnd {
        fillet_index: 1,
        curve: make_test_curve(Point3::new(0.0, r, 0.0), Point3::new(0.0, 0.0, r)),
    };

    let boundary = BlendBoundary {
        vertex,
        vertex_position: Point3::new(0.0, 0.0, 0.0),
        segments: vec![seg1, seg2],
        bridge_fallbacks: Vec::new(),
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
    };

    let fillets: Vec<FilletSurfaceResult> = Vec::new();

    let result = build_blend_surface(
        &boundary,
        &fillets,
        &brep,
        VertexMode::Smooth,
        crate::vertex_classifier::VertexBlendClassification::SameSideConvex,
        &tol,
    );
    assert!(result.is_ok(), "smooth blend failed: {:?}", result.err());

    let blend = result.expect("already checked");
    // Should produce either a sphere or a NURBS patch (or plane fallback).
    // The exact result depends on point distribution; we just verify success.
    assert!(
        !blend.boundary_curves.is_empty(),
        "should have boundary curves"
    );
}

#[test]
fn blend_surface_empty_boundary() {
    let tol = test_tolerance();
    let brep = BRep::unit_cube(tol).expect("cube should build");
    let vertex = brep.vertices.keys().next().expect("cube has vertices");

    let boundary = BlendBoundary {
        vertex,
        vertex_position: Point3::new(0.0, 0.0, 0.0),
        segments: vec![],
        bridge_fallbacks: Vec::new(),
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
    };

    let result = build_blend_surface(
        &boundary,
        &[],
        &brep,
        VertexMode::Mitred,
        crate::vertex_classifier::VertexBlendClassification::SameSideConvex,
        &tol,
    );
    assert!(result.is_err(), "should fail with empty boundary");
}

#[test]
fn fit_plane_coplanar_points() {
    let points = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
    ];

    let (centroid, normal) = fit_plane_least_squares(&points).expect("should fit");

    // Centroid should be at (0.5, 0.5, 0.0).
    assert!((centroid.x - 0.5).abs() < 1e-10);
    assert!((centroid.y - 0.5).abs() < 1e-10);
    assert!(centroid.z.abs() < 1e-10);

    // Normal should be along z (up to sign).
    assert!(
        normal.z.abs() > 0.99,
        "normal should be along z, got {normal:?}"
    );
}
