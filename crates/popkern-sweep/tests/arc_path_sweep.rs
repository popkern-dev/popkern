//! PK-440: sweeping along a curved (arc) path must not promote the
//! swept walls to planes. A straight profile segment parallel to the
//! arc's axis passes the extruded-plane recognizer's per-row
//! translation test exactly — every row is a rigid translate of row
//! zero — while the rail curves; before the collinear-rail gate, a 45°
//! arc sweep promoted its inner and outer walls to planes and every
//! corner edge deviated from its "plane" by the arc's sagitta
//! (0.056–0.064 on a 0.8-radius arc).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_sweep::{SweepOptions, SweepProfile, sweep};

/// Exact rational 45° arc about the z axis, starting at angle 0.
fn arc_45(r: f64, z: f64) -> NurbsCurve {
    let half = 22.5_f64.to_radians();
    let w = half.cos();
    let mid_r = r / w;
    NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(r, 0.0, z),
            Point3::new(mid_r * half.cos(), mid_r * half.sin(), z),
            Point3::new(r * (2.0 * half).cos(), r * (2.0 * half).sin(), z),
        ],
        weights: vec![1.0, w, 1.0],
    }
}

fn square_wire(half: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0],
        control_points: vec![
            Point3::new(-half, -half, 0.0),
            Point3::new(half, -half, 0.0),
            Point3::new(half, half, 0.0),
            Point3::new(-half, half, 0.0),
            Point3::new(-half, -half, 0.0),
        ],
        weights: vec![1.0; 5],
    }
}

#[test]
fn arc_path_sweep_walls_stay_curved_and_edges_lie_on_their_faces() {
    let result = sweep(
        &SweepProfile::Wire(square_wire(0.05)),
        &arc_45(0.8, 0.1),
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("arc-path sweep serves");
    let brep = result.brep;
    assert_eq!(brep.faces.len(), 6, "4 swept walls + 2 caps");

    let plane_count = brep
        .faces
        .values()
        .filter(|f| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Plane { .. })
            )
        })
        .count();
    assert_eq!(
        plane_count, 2,
        "only the caps are planar; curved-rail walls must not be \
         promoted to planes"
    );

    let deviations = brep.validate_curve_surface_all();
    assert!(
        deviations.is_empty(),
        "every edge lies on both adjacent faces: {deviations:?}"
    );
}

/// The gate must not regress straight-path sweeps: a square swept
/// along a line still promotes its four walls (and caps) to planes.
#[test]
fn straight_path_sweep_still_promotes_planes() {
    let line = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.7)],
        weights: vec![1.0, 1.0],
    };
    let result = sweep(
        &SweepProfile::Wire(square_wire(0.05)),
        &line,
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("straight sweep serves");
    let brep = result.brep;
    let plane_count = brep
        .faces
        .values()
        .filter(|f| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Plane { .. })
            )
        })
        .count();
    assert_eq!(
        plane_count,
        brep.faces.len(),
        "a straight extrusion is all planes"
    );
    assert!(brep.validate_curve_surface_all().is_empty());
}
