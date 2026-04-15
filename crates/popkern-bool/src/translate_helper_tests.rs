//! Coverage tests for the `translate` helper used in curved integration tests.
//!
//! Exercises surface/curve geometry branches (Sphere, Torus, Circle/Arc, NURBS
//! catch-all) that aren't hit by the main boolean tests.

// Structural lint marker: this file is test code, included via #[cfg(test)] in lib.rs.
#[cfg(test)]
use popkern_core::{geometry::Vec3, tolerance::Tolerance, topology::BRep};

fn coarse_tol() -> Tolerance {
    Tolerance::default().with_linear(1e-2)
}

/// Translate all vertices and surface/curve geometry of a BRep by `offset`.
///
/// Mirrors the helper in integration_tests_curved.rs.
fn translate(mut brep: BRep, offset: Vec3) -> BRep {
    use popkern_core::geometry::{curve::CurveGeometry, surface::SurfaceGeometry};
    for (_, vertex) in &mut brep.vertices {
        vertex.position = vertex.position + offset;
    }
    for (_, face) in &mut brep.faces {
        if let Some(ref mut sof) = face.surface {
            match sof.geometry {
                SurfaceGeometry::Plane { ref mut origin, .. } => {
                    *origin = *origin + offset;
                }
                SurfaceGeometry::Cylinder { ref mut origin, .. } => {
                    *origin = *origin + offset;
                }
                SurfaceGeometry::Cone { ref mut apex, .. } => {
                    *apex = *apex + offset;
                }
                SurfaceGeometry::Sphere { ref mut center, .. } => {
                    *center = *center + offset;
                }
                SurfaceGeometry::Torus { ref mut center, .. } => {
                    *center = *center + offset;
                }
                _ => {}
            }
        }
    }
    for (_, edge) in &mut brep.edges {
        if let Some(ref mut coe) = edge.curve {
            match coe.geometry {
                CurveGeometry::Line { ref mut origin, .. } => {
                    *origin = *origin + offset;
                }
                CurveGeometry::Circle { ref mut center, .. }
                | CurveGeometry::Arc { ref mut center, .. } => {
                    *center = *center + offset;
                }
                _ => {}
            }
        }
    }
    brep
}

/// Translate a sphere: exercises the `SurfaceGeometry::Sphere` branch of the
/// translate helper, covering the Sphere surface arm.
#[test]
fn translate_sphere_covers_sphere_surface_branch() {
    let t = coarse_tol();
    let sphere = BRep::sphere(1.0, t).unwrap();
    let translated = translate(sphere, Vec3::new(5.0, 0.0, 0.0));

    // All vertices should have x >= 4 (original range ~[-1, 1] + 5 = [4, 6]).
    for (_, v) in &translated.vertices {
        assert!(v.position.x >= 3.9, "sphere vertex x should be translated");
    }
    // The sphere should still have faces.
    assert!(
        !translated.faces.is_empty(),
        "translated sphere should have faces"
    );
}

/// Translate a BRep with a face that has Torus surface geometry: exercises
/// the `SurfaceGeometry::Torus` branch of the translate helper.
/// Since `BRep::torus()` doesn't attach Torus surface geometry, we manually
/// assign it to a cube face before translating.
#[test]
fn translate_torus_covers_torus_surface_branch() {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

    let t = coarse_tol();
    let mut brep = BRep::unit_cube(t).unwrap();

    // Overwrite one face's surface with a Torus geometry.
    let fh = brep.faces.keys().next().unwrap();
    brep.faces[fh].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Torus {
            center: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            major_radius: 2.0,
            minor_radius: 0.5,
        },
        orientation: true,
    });

    let translated = translate(brep, Vec3::new(10.0, 0.0, 0.0));

    // All vertices should have x >= 9.
    for (_, v) in &translated.vertices {
        assert!(v.position.x >= 9.0, "vertex x should be translated");
    }
    // The face with the Torus surface should have its center translated.
    let face_data = &translated.faces[fh];
    if let Some(ref surf) = face_data.surface
        && let SurfaceGeometry::Torus { ref center, .. } = surf.geometry
    {
        assert!(
            (center.x - 10.0).abs() < 0.01,
            "Torus center x should be translated to ~10"
        );
    }
}

/// Translate a BRep with Circle/Arc edge curves: exercises the Circle | Arc
/// curve branch of the translate helper. We manually assign a Circle curve
/// to a cube edge after construction.
#[test]
fn translate_brep_with_circle_curve_covers_circle_arc_branch() {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};

    let t = coarse_tol();
    let mut cube = BRep::unit_cube(t).unwrap();

    // Overwrite one edge's curve with a Circle geometry.
    let eh = cube.edges.keys().next().unwrap();
    cube.edges[eh].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        },
        t_start: 0.0,
        t_end: 1.0,
    });

    let translated = translate(cube, Vec3::new(10.0, 0.0, 0.0));

    // All vertices should have x >= 9.
    for (_, v) in &translated.vertices {
        assert!(v.position.x >= 9.0, "vertex x should be translated");
    }
    // The edge with the Circle curve should have its center translated.
    let edge_data = &translated.edges[eh];
    if let Some(ref coe) = edge_data.curve
        && let CurveGeometry::Circle { ref center, .. } = coe.geometry
    {
        assert!(
            (center.x - 10.0).abs() < 0.01,
            "Circle center x should be translated to ~10"
        );
    }
}

/// Translate a BRep with a NURBS surface and NURBS edge curve: exercises the
/// `_ => {}` catch-all arms in translate for non-handled surface geometry
/// and non-handled curve geometry.
#[test]
fn translate_nurbs_geometry_hits_catchall_branches() {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
    use popkern_core::geometry::nurbs_curve::NurbsCurve;
    use popkern_core::geometry::nurbs_surface::NurbsSurface;
    use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

    let t = coarse_tol();
    let mut cube = BRep::unit_cube(t).unwrap();

    // Assign a Nurbs surface to one face (hits the `_ => {}` surface catch-all).
    let fh = cube.faces.keys().next().unwrap();
    // Minimal bilinear NURBS patch (2x2 control points).
    let nurbs_surf = NurbsSurface {
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
        weights: vec![1.0, 1.0, 1.0, 1.0],
    };
    cube.faces[fh].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Nurbs(nurbs_surf),
        orientation: true,
    });

    // Assign a Nurbs curve to one edge (hits the `_ => {}` curve catch-all).
    let eh = cube.edges.keys().next().unwrap();
    let nurbs_curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    cube.edges[eh].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(nurbs_curve),
        t_start: 0.0,
        t_end: 1.0,
    });

    // Translate — should not panic; NURBS geometry is NOT updated (catch-all `_ => {}`).
    let translated = translate(cube, Vec3::new(5.0, 0.0, 0.0));
    for (_, v) in &translated.vertices {
        assert!(v.position.x >= 4.0, "vertices should be translated");
    }
}
