//! Tests for the blend::gregory builder + twist-vector extraction.
//!
//! Split from `gregory.rs` so the implementation file stays under
//! the 500-line structural cap.

use crate::blend::gregory::{
    extract_triangular_inner_cps, fit_cubic_bezier_through, try_rational_gregory_blend,
};
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3};

use crate::types::{ContinuityLevel, CornerMethod};

/// `fit_cubic_bezier_through` on a straight-line NURBS (4 CPs
/// on a line) reproduces the input CPs exactly.
#[test]
fn fit_straight_line_cubic_recovers_cps() {
    let cps = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0 / 3.0, 0.0, 0.0),
        Point3::new(2.0 / 3.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    ];
    let curve = NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: cps.clone(),
        weights: vec![1.0; 4],
    };
    let fitted = fit_cubic_bezier_through(&curve);
    for (a, b) in fitted.iter().zip(cps.iter()) {
        assert!(
            (a.x - b.x).abs() < 1e-10 && (a.y - b.y).abs() < 1e-10 && (a.z - b.z).abs() < 1e-10,
            "fit mismatch: {:?} vs {:?}",
            a,
            b,
        );
    }
}

/// `fit_cubic_bezier_through` on a non-linear cubic reconstructs
/// the 4 CPs exactly (the 4 Bernstein samples at `t = 0, 1/3, 2/3, 1`
/// uniquely determine a cubic Bezier).
#[test]
fn fit_arbitrary_cubic_reconstructs_cps() {
    let cps = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.2, 0.5, 0.3),
        Point3::new(0.8, 0.5, 0.7),
        Point3::new(1.0, 0.0, 0.0),
    ];
    let curve = NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: cps.clone(),
        weights: vec![1.0; 4],
    };
    let fitted = fit_cubic_bezier_through(&curve);
    for (a, b) in fitted.iter().zip(cps.iter()) {
        let d = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt();
        assert!(d < 1e-8, "fit mismatch: {:?} vs {:?} (d={d:.2e})", a, b);
    }
}

fn make_cube_corner_arc(start: Point3, end: Point3) -> NurbsCurve {
    let mid1 = Point3::new(
        start.x + (end.x - start.x) / 3.0,
        start.y + (end.y - start.y) / 3.0,
        start.z + (end.z - start.z) / 3.0,
    );
    let mid2 = Point3::new(
        start.x + 2.0 * (end.x - start.x) / 3.0,
        start.y + 2.0 * (end.y - start.y) / 3.0,
        start.z + 2.0 * (end.z - start.z) / 3.0,
    );
    NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![start, mid1, mid2, end],
        weights: vec![1.0, 1.0, 1.0, 1.0],
    }
}

/// End-to-end builder acceptance on a synthetic cube-corner-like
/// fixture: 3 arcs around origin at radius r with 3 cylinder
/// strips. The builder returns a triangular Gregory patch with
/// preserved boundary curves and `CornerMethod::Gregory`.
#[test]
fn try_rational_gregory_on_synthetic_cube_corner() {
    use crate::blend::boundary::{BlendBoundary, BoundarySegment};
    use crate::fillet_surface::FilletSurfaceResult;
    use popkern_core::geometry::Vec3;
    use popkern_core::topology::VertexHandle;

    let r = 0.1;
    let vertex = VertexHandle::default();
    let v_pos = Point3::new(0.0, 0.0, 0.0);
    let ax = Point3::new(r, 0.0, 0.0);
    let ay = Point3::new(0.0, r, 0.0);
    let az = Point3::new(0.0, 0.0, r);
    let arc_xy = make_cube_corner_arc(ax, ay);
    let arc_yz = make_cube_corner_arc(ay, az);
    let arc_zx = make_cube_corner_arc(az, ax);
    let boundary_curves = vec![arc_xy.clone(), arc_yz.clone(), arc_zx.clone()];
    let boundary = BlendBoundary {
        vertex,
        vertex_position: v_pos,
        segments: vec![
            BoundarySegment::FilletEnd {
                fillet_index: 0,
                curve: arc_xy,
            },
            BoundarySegment::FilletEnd {
                fillet_index: 1,
                curve: arc_yz,
            },
            BoundarySegment::FilletEnd {
                fillet_index: 2,
                curve: arc_zx,
            },
        ],
        bridge_fallbacks: Vec::new(),
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
        hex_boundary_cps_deg9: Vec::new(),
        fillet_end_arcs: Vec::new(),
    };
    let fillet_surfaces = vec![
        FilletSurfaceResult {
            surface: SurfaceGeometry::Cylinder {
                origin: v_pos,
                axis: Vec3::new(1.0, 0.0, 0.0),
                radius: r,
            },
            trim_curve_1: make_cube_corner_arc(ax, ay),
            trim_curve_2: make_cube_corner_arc(ax, ay),
            continuity_downgrade: None,
        },
        FilletSurfaceResult {
            surface: SurfaceGeometry::Cylinder {
                origin: v_pos,
                axis: Vec3::new(0.0, 1.0, 0.0),
                radius: r,
            },
            trim_curve_1: make_cube_corner_arc(ay, az),
            trim_curve_2: make_cube_corner_arc(ay, az),
            continuity_downgrade: None,
        },
        FilletSurfaceResult {
            surface: SurfaceGeometry::Cylinder {
                origin: v_pos,
                axis: Vec3::new(0.0, 0.0, 1.0),
                radius: r,
            },
            trim_curve_1: make_cube_corner_arc(az, ax),
            trim_curve_2: make_cube_corner_arc(az, ax),
            continuity_downgrade: None,
        },
    ];

    let result = try_rational_gregory_blend(&boundary, &fillet_surfaces, &boundary_curves)
        .expect("Gregory build should succeed on synthetic cube corner");
    match &result.surface {
        SurfaceGeometry::Gregory(GregoryPatch::Triangular(_)) => {}
        other => panic!("expected Triangular Gregory, got {other:?}"),
    }
    assert_eq!(result.boundary_curves.len(), 3);
    assert_eq!(result.method, CornerMethod::Gregory);
    assert_eq!(result.continuity, ContinuityLevel::G2);
}

/// §5.3 — twist-vector extraction displaces each inner
/// candidate from its associated boundary midpoint by ~avg_radius
/// (the defining property vs the earlier `[v_pos; 3]` placeholder
/// which collapsed all 3 onto the vertex position).
#[test]
fn twist_extraction_displaces_inner_cps_from_boundary_midpoint() {
    use crate::fillet_surface::FilletSurfaceResult;
    use popkern_core::geometry::Vec3;

    let r = 0.1;
    let v_pos = Point3::new(0.0, 0.0, 0.0);
    let ax = Point3::new(r, 0.0, 0.0);
    let ay = Point3::new(0.0, r, 0.0);
    let az = Point3::new(0.0, 0.0, r);
    let boundary_curves = vec![
        make_cube_corner_arc(ay, az),
        make_cube_corner_arc(az, ax),
        make_cube_corner_arc(ax, ay),
    ];
    let make_cyl = |axis: Vec3| FilletSurfaceResult {
        surface: SurfaceGeometry::Cylinder {
            origin: v_pos,
            axis,
            radius: r,
        },
        trim_curve_1: make_cube_corner_arc(ax, ay),
        trim_curve_2: make_cube_corner_arc(ax, ay),
        continuity_downgrade: None,
    };
    let fillet_surfaces = vec![
        make_cyl(Vec3::new(1.0, 0.0, 0.0)),
        make_cyl(Vec3::new(0.0, 1.0, 0.0)),
        make_cyl(Vec3::new(0.0, 0.0, 1.0)),
    ];
    let inner = extract_triangular_inner_cps(&boundary_curves, &fillet_surfaces, r, v_pos)
        .expect("extraction should succeed");
    let midpoints = [
        Point3::new(0.0, r / 2.0, r / 2.0),
        Point3::new(r / 2.0, 0.0, r / 2.0),
        Point3::new(r / 2.0, r / 2.0, 0.0),
    ];
    for (k, (p, m)) in inner.iter().zip(midpoints.iter()).enumerate() {
        let d = ((p.x - m.x).powi(2) + (p.y - m.y).powi(2) + (p.z - m.z).powi(2)).sqrt();
        assert!(
            (d - r).abs() < 0.05 * r,
            "inner CP {k} at {p:?} should be ≈r from its boundary midpoint \
             {m:?}; got {d}",
        );
    }
}

/// Builder rejects mixed-radius incidence — guards against
/// firing the Gregory path when the no-twist placement would
/// produce a G0-at-corner result.
#[test]
fn try_rational_gregory_rejects_mixed_radii() {
    use crate::blend::boundary::{BlendBoundary, BoundarySegment};
    use crate::fillet_surface::FilletSurfaceResult;
    use popkern_core::geometry::Vec3;
    use popkern_core::topology::VertexHandle;

    let v_pos = Point3::new(0.0, 0.0, 0.0);
    let vertex = VertexHandle::default();
    let curve = NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.1, 0.0, 0.0),
            Point3::new(0.0, 0.033, 0.0),
            Point3::new(0.0, 0.066, 0.0),
            Point3::new(0.0, 0.1, 0.0),
        ],
        weights: vec![1.0; 4],
    };
    let boundary = BlendBoundary {
        vertex,
        vertex_position: v_pos,
        segments: vec![
            BoundarySegment::FilletEnd {
                fillet_index: 0,
                curve: curve.clone(),
            },
            BoundarySegment::FilletEnd {
                fillet_index: 1,
                curve: curve.clone(),
            },
            BoundarySegment::FilletEnd {
                fillet_index: 2,
                curve: curve.clone(),
            },
        ],
        bridge_fallbacks: Vec::new(),
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
        hex_boundary_cps_deg9: Vec::new(),
        fillet_end_arcs: Vec::new(),
    };
    let fillet_surfaces: Vec<_> = [0.08, 0.10, 0.12]
        .iter()
        .map(|&r| FilletSurfaceResult {
            surface: SurfaceGeometry::Cylinder {
                origin: v_pos,
                axis: Vec3::new(1.0, 0.0, 0.0),
                radius: r,
            },
            trim_curve_1: curve.clone(),
            trim_curve_2: curve.clone(),
            continuity_downgrade: None,
        })
        .collect();
    let boundary_curves = vec![curve.clone(), curve.clone(), curve.clone()];
    let result = try_rational_gregory_blend(&boundary, &fillet_surfaces, &boundary_curves);
    assert!(result.is_none(), "mixed radii should return None");
}
