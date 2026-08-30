//! Exactness contract of analytic profiles (scripted-modeling design
//! invariant 6, PK-364): the compiled wire is the form the sweep
//! kernel's canonical recognition maps to true `Cylinder` walls, and
//! the compiled face is the form the extrude operator promotes the
//! same way. Moved here from `popkern-part` with the compilation.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_sketch::{AnalyticProfile, PathSegment, ProfileShape, SketchPlane, Tolerance};
use popkern_sweep::{SweepOptions, SweepProfile, sweep};

fn z_path(z0: f64, z1: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, z0), Point3::new(0.0, 0.0, z1)],
        weights: vec![1.0, 1.0],
    }
}

fn surface_census(brep: &popkern_core::topology::BRep) -> (usize, usize, usize) {
    let (mut planes, mut cylinders, mut other) = (0, 0, 0);
    for f in brep.faces.values() {
        match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Plane { .. }) => planes += 1,
            Some(SurfaceGeometry::Cylinder { .. }) => cylinders += 1,
            _ => other += 1,
        }
    }
    (planes, cylinders, other)
}

/// The rational-quadratic circle wire sweeps to a TRUE `Cylinder` at
/// the exact centre and radius.
#[test]
fn circle_wire_sweeps_to_a_true_cylinder() {
    let wire = AnalyticProfile::new(ProfileShape::Circle { radius: 0.25 })
        .local_wire(&SketchPlane::xy())
        .unwrap();
    let result = sweep(
        &SweepProfile::Wire(wire),
        &z_path(0.0, 0.5),
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("circle sweeps");
    result.brep.validate().expect("cylinder validates");
    let mut found = false;
    for f in result.brep.faces.values() {
        if let Some(SurfaceGeometry::Cylinder { origin, radius, .. }) =
            f.surface.as_ref().map(|s| &s.geometry)
        {
            found = true;
            assert!(origin.x.abs() < 1e-9 && origin.y.abs() < 1e-9);
            assert!((radius - 0.25).abs() < 1e-9);
        }
    }
    assert!(found, "canonical recognition must produce a Cylinder wall");
    assert_eq!(surface_census(&result.brep), (2, 1, 0));
}

/// A mixed line/arc outline (a rectangle whose front side bulges as a
/// circular arc) sweeps to a watertight, χ-clean solid whose arc wall
/// is a TRUE `Cylinder` and whose straight walls are `Plane`s.
#[test]
fn mixed_path_wire_sweeps_to_true_cylinder_wall() {
    let wire = AnalyticProfile::new(ProfileShape::Path {
        start: [-30.0, 0.0],
        segments: vec![
            PathSegment::Line { to: [30.0, 0.0] },
            PathSegment::Arc {
                via: [0.0, 18.0],
                to: [-30.0, 0.0],
            },
        ],
    })
    .local_wire(&SketchPlane::xy())
    .unwrap();
    let result = sweep(
        &SweepProfile::Wire(wire),
        &z_path(0.0, 20.0),
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("mixed sweep");
    let b = &result.brep;
    b.validate().expect("mixed solid validates");
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "mixed-path solid must be genus 0");
    let (planes, cylinders, other) = surface_census(b);
    assert!(
        cylinders >= 1,
        "arc wall must be a true Cylinder (got {cylinders})"
    );
    assert_eq!(other, 0, "no unrecognized walls (got {other})");
    assert!(
        planes >= 3,
        "straight wall + caps must be Planes (got {planes})"
    );
}

/// The face path reaches the same exactness through the extrude
/// operator's `Arc` → `Cylinder` promotion (PK-377).
#[test]
fn circle_face_extrudes_to_the_same_cylinder_as_the_wire() {
    let profile = AnalyticProfile {
        shape: ProfileShape::Circle { radius: 0.25 },
        center: [1.0, 0.0],
    };
    let result = profile
        .to_brep(SketchPlane::xy(), Tolerance::default(), None)
        .unwrap();
    let mut b = result.brep;
    b.extrude(result.regions[0].face, Vec3::new(0.0, 0.0, 1.0), 0.5, None)
        .unwrap();
    assert!(b.validate_all().is_empty(), "{:?}", b.validate_all());
    // Two half-arc edges → two cylinder walls on one axis; the wire
    // path serves one seamed wall. Same surface, different seam census.
    assert_eq!(surface_census(&b), (2, 2, 0));
    for f in b.faces.values() {
        if let Some(SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        }) = f.surface.as_ref().map(|s| &s.geometry)
        {
            assert!((origin.x - 1.0).abs() < 1e-12 && origin.y.abs() < 1e-12);
            assert!((axis.z.abs() - 1.0).abs() < 1e-12);
            assert!((radius - 0.25).abs() < 1e-12);
        }
    }
}
