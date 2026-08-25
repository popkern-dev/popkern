use crate::offset::offset_solid;
use crate::thicken::thicken;
use crate::types::{JoinStyle, ShellOptions};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn arc_options() -> ShellOptions {
    ShellOptions {
        join_style: JoinStyle::Arc,
        ..ShellOptions::default()
    }
}

fn sharp_options() -> ShellOptions {
    ShellOptions {
        join_style: JoinStyle::Intersection,
        ..ShellOptions::default()
    }
}

fn surface_census(brep: &BRep) -> (usize, usize, usize, usize) {
    let mut planes = 0;
    let mut cylinders = 0;
    let mut spheres = 0;
    let mut other = 0;
    for (_, f) in &brep.faces {
        match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Plane { .. }) => planes += 1,
            Some(SurfaceGeometry::Cylinder { .. }) => cylinders += 1,
            Some(SurfaceGeometry::Sphere { .. }) => spheres += 1,
            _ => other += 1,
        }
    }
    (planes, cylinders, spheres, other)
}

/// Row 472 topology half: the outward cube offset under Arc join is the
/// rounded box — 6 planes + 12 cylinder windows + 8 sphere octants.
#[test]
fn cube_outward_arc_offset_is_the_rounded_box() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let d = 0.15;
    let result = offset_solid(&brep, d, &arc_options(), tol()).unwrap();
    let out = &result.brep;

    assert_eq!(out.faces.len(), 26, "6 planes + 12 strips + 8 caps");
    assert_eq!(out.edges.len(), 48, "24 rails + 24 arcs");
    assert_eq!(out.vertices.len(), 24, "3 per corner");
    let (planes, cylinders, spheres, other) = surface_census(out);
    assert_eq!(planes, 6);
    assert_eq!(cylinders, 12);
    assert_eq!(spheres, 8);
    assert_eq!(other, 0);

    let errors = out.validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");
    let _ = d;
}

/// Row 472 sharp half: Intersection join reproduces the extended sharp
/// box (today's behavior) on the same input.
#[test]
fn cube_outward_intersection_offset_stays_sharp() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = offset_solid(&brep, 0.15, &sharp_options(), tol()).unwrap();
    let out = &result.brep;
    assert_eq!(out.faces.len(), 6);
    assert_eq!(out.edges.len(), 12);
    assert_eq!(out.vertices.len(), 8);
}

/// Thicken under Arc: rounded outer skin, sharp (erosion) inner skin.
#[test]
fn cube_thicken_arc_rounds_only_the_outer_skin() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let t = 0.2;
    let result = thicken(&brep, t, &arc_options(), tol()).unwrap();
    let out = &result.brep;
    // Outer skin 26 faces + inner skin 6 faces.
    assert_eq!(out.faces.len(), 32, "rounded outer + sharp inner");
    let errors = out.validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");
}

/// The cylinder thicken under Arc: closed rim edges pre-split and serve
/// torus quarter-round windows; the seam-transparent k=2 corners need
/// no caps.
#[test]
fn cylinder_thicken_arc_serves_torus_rims() {
    let brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let t = 0.2;
    let result = thicken(&brep, t, &arc_options(), tol()).unwrap();
    let out = &result.brep;

    let torus_faces = out
        .faces
        .iter()
        .filter(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Torus { .. })
            )
        })
        .count();
    assert_eq!(torus_faces, 4, "two rims x two half-windows");
    let errors = out.validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");
}
