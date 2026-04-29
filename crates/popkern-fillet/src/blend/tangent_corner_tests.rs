//! Unit tests for the tangent-corner planar cap.

use super::try_tangent_plane_blend;

use crate::blend::boundary::{BlendBoundary, BoundarySegment};
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ContinuityLevel, CornerMethod};

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::VertexHandle;

/// Build a degree-1 NURBS curve from two points. Used to fake
/// `FilletEnd` arcs in these synthetic tests.
fn line_curve(a: Point3, b: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    }
}

#[test]
fn tangent_cap_cylinder_plane_tangent_produces_xz_plane() {
    // Configuration: cylinder axis = +z, radius = 1, touching the
    // plane z = -1 along the line y = 0 (tangent configuration).
    // The fillet between them has a near-π dihedral at the trim
    // edge. The tangent cap should emit a plane whose normal is
    // the fillet surface's normal at the arc midpoint — which for
    // a cylinder at (x=1, y=0) is +x.
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let fillet_arc = line_curve(Point3::new(1.0, 0.0, -0.1), Point3::new(1.0, 0.0, 0.1));
    let fillet = FilletSurfaceResult {
        surface: cyl,
        trim_curve_1: fillet_arc.clone(),
        trim_curve_2: fillet_arc.clone(),
        continuity_downgrade: None,
    };
    let boundary = BlendBoundary {
        vertex: VertexHandle::default(),
        vertex_position: Point3::new(1.0, 0.0, 0.0),
        segments: vec![BoundarySegment::FilletEnd {
            fillet_index: 0,
            curve: fillet_arc.clone(),
        }],
        bridge_fallbacks: Vec::new(),
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
        hex_boundary_cps_deg9: Vec::new(),
    };
    let points = vec![fillet_arc.point_at(0.0), fillet_arc.point_at(1.0)];

    let result = try_tangent_plane_blend(&boundary, &[fillet], &points, &[fillet_arc])
        .expect("tangent cap should build");

    assert_eq!(result.method, CornerMethod::TangentPlane);
    assert_eq!(result.continuity, ContinuityLevel::G1);
    // The cap plane's normal must match the cylinder's normal at
    // the arc midpoint (1, 0, 0), i.e. +x (or −x).
    if let SurfaceGeometry::Plane { normal, .. } = &result.surface {
        let dot_x = normal.x.abs();
        assert!(
            dot_x > 0.999,
            "tangent cap plane normal must align with cylinder radial direction; got {normal:?}"
        );
    } else {
        panic!("tangent cap must be a Plane, got {:?}", result.surface);
    }
}

#[test]
fn tangent_cap_empty_boundary_returns_none() {
    let boundary = BlendBoundary {
        vertex: VertexHandle::default(),
        vertex_position: Point3::new(0.0, 0.0, 0.0),
        segments: vec![],
        bridge_fallbacks: Vec::new(),
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
        hex_boundary_cps_deg9: Vec::new(),
    };
    let result = try_tangent_plane_blend(&boundary, &[], &[], &[]);
    assert!(result.is_none(), "empty boundary should return None");
}

/// §5.3 item #7 — end-to-end dispatch: `build_blend_surface` with
/// `VertexBlendClassification::Tangent` must route to the tangent-
/// plane cap (ahead of the sphere/Gregory/mitred fallback chain).
#[test]
fn build_blend_surface_dispatches_tangent_to_planar_cap() {
    use crate::blend::surface::build_blend_surface;
    use crate::types::VertexMode;
    use crate::vertex_classifier::VertexBlendClassification;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    let tol = Tolerance::default();
    let brep = BRep::unit_cube(tol).expect("unit cube");
    let vertex = brep.vertices.keys().next().expect("has vertex");

    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let arc = line_curve(Point3::new(1.0, 0.0, -0.1), Point3::new(1.0, 0.0, 0.1));
    let fillet = FilletSurfaceResult {
        surface: cyl,
        trim_curve_1: arc.clone(),
        trim_curve_2: arc.clone(),
        continuity_downgrade: None,
    };
    let boundary = BlendBoundary {
        vertex,
        vertex_position: Point3::new(1.0, 0.0, 0.0),
        segments: vec![BoundarySegment::FilletEnd {
            fillet_index: 0,
            curve: arc.clone(),
        }],
        bridge_fallbacks: Vec::new(),
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
        hex_boundary_cps_deg9: Vec::new(),
    };

    let result = build_blend_surface(
        &boundary,
        &[fillet],
        &brep,
        VertexMode::Smooth,
        VertexBlendClassification::Tangent,
        &tol,
    )
    .expect("tangent classification must dispatch successfully");
    assert_eq!(result.method, CornerMethod::TangentPlane);
    assert_eq!(result.continuity, ContinuityLevel::G1);
}

#[test]
fn tangent_cap_averages_multiple_fillet_normals() {
    // Two fillets with nearly-parallel normals (both roughly +z
    // tangent to a flat top). Averaged normal should also be +z.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let arc1 = line_curve(Point3::new(0.0, 0.0, 0.0), Point3::new(0.5, 0.0, 0.0));
    let arc2 = line_curve(Point3::new(0.5, 0.0, 0.0), Point3::new(0.5, 0.5, 0.0));

    let fillet = |curve: NurbsCurve| FilletSurfaceResult {
        surface: plane.clone(),
        trim_curve_1: curve.clone(),
        trim_curve_2: curve,
        continuity_downgrade: None,
    };
    let boundary = BlendBoundary {
        vertex: VertexHandle::default(),
        vertex_position: Point3::new(0.0, 0.0, 0.0),
        segments: vec![
            BoundarySegment::FilletEnd {
                fillet_index: 0,
                curve: arc1.clone(),
            },
            BoundarySegment::FilletEnd {
                fillet_index: 1,
                curve: arc2.clone(),
            },
        ],
        bridge_fallbacks: Vec::new(),
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
        hex_boundary_cps_deg9: Vec::new(),
    };
    let points = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.5, 0.0, 0.0),
        Point3::new(0.5, 0.5, 0.0),
    ];

    let result = try_tangent_plane_blend(
        &boundary,
        &[fillet(arc1.clone()), fillet(arc2.clone())],
        &points,
        &[arc1, arc2],
    )
    .expect("tangent cap should build with multiple fillets");
    if let SurfaceGeometry::Plane { normal, .. } = &result.surface {
        assert!(
            normal.z.abs() > 0.999,
            "averaged tangent-cap normal should align with +z; got {normal:?}"
        );
    } else {
        panic!("expected Plane");
    }
}
